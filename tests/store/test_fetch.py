"""Tests for genome.store.fetch — the package's one fetch step, and why it lives alone.

``fetch_url`` is exercised for real against a file already sitting at the destination,
which pooch serves without touching the network; nothing here monkeypatches pooch's own
retrieve function. Everything else in the suite replaces this one function instead — see
the ``fake_fetch`` fixture in tests/conftest.py.

The module is a leaf on purpose. Both the downloader and the annotation registration
fetch, and while the fetch lived under the downloader the second of those was an import
back into it — the last edge of the ``download -> chimera -> gtf -> download`` cycle. A
leaf cannot close a cycle, so the guard below is what keeps that true.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from types import ModuleType

import pooch
import pytest

from genome.annotation import registration as annotation_module
from genome.assembly import download as download_module
from genome.store import fetch as fetch_module

from ..assembly.test_source import _module_level_imports


def _sha256(path: Path) -> str:
    """Return the sha256 of ``path`` in the ``algorithm:hexdigest`` form pooch accepts."""
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def test_a_matching_local_file_is_served_without_downloading_and_a_processor_output_returned(
    tmp_path: Path, data_dir: Path
) -> None:
    # A file already at the destination whose hash matches is handed back as-is:
    # no downloader is ever constructed, so this exercises fetch_url offline.
    dest = tmp_path / "tiny.fa"
    dest.write_bytes((data_dir / "tiny.fa").read_bytes())

    result = fetch_module.fetch_url(
        "https://example.org/tiny.fa",
        tmp_path,
        known_hash=_sha256(dest),
        fname="tiny.fa",
        progressbar=False,
    )

    assert result.resolve() == dest.resolve()
    assert result.read_bytes() == (data_dir / "tiny.fa").read_bytes()

    # And when a processor is given, its output — not the raw destination — comes back.
    gz_dest = tmp_path / "sacCer3.fa.gz"
    gz_dest.write_bytes((data_dir / "tiny.fa.gz").read_bytes())

    processed = fetch_module.fetch_url(
        "https://example.org/sacCer3.fa.gz",
        tmp_path,
        known_hash=_sha256(gz_dest),
        fname="sacCer3.fa.gz",
        processor=pooch.Decompress(method="gzip", name="sacCer3.fa"),
        progressbar=False,
    )

    assert processed.resolve() == (tmp_path / "sacCer3.fa").resolve()
    assert processed.read_text() == (data_dir / "tiny.fa").read_text()


def test_the_fetch_step_imports_nothing_from_its_own_package() -> None:
    # What makes this a place two modules can both import: it depends on no other module
    # of the package, so no import of it can ever be half of a cycle. The moment
    # something here needs a name from `genome`, the fetch has stopped being one step and
    # the edge it was split out to close is available again.
    own = {name for name in _module_level_imports(fetch_module) if name.startswith("genome")}

    assert own == set()


@pytest.mark.parametrize(
    "module", [download_module, annotation_module], ids=["download", "annotation"]
)
def test_every_caller_reaches_the_fetch_through_the_module(
    module: ModuleType,
) -> None:
    # The spelling the whole suite rests on. Both callers hold the *module* and look the
    # function up on it at call time, so rebinding this one attribute — which is all
    # `fake_fetch` and the offline guard do — reaches every download in the package. A
    # caller that imported the name would hold a reference no rebinding can follow, and
    # would keep downloading while the suite believed itself offline.
    assert "genome.store.fetch" in _module_level_imports(module)
    assert "genome.store.fetch.fetch_url" not in _module_level_imports(module)
