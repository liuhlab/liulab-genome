"""Tests for genome.io.utils — native tools, checksums, and output-freshness caching.

None of these need the native binaries: ``_run``'s success path is exercised by
the end-to-end tests in test_fasta, and here ``_run`` is either stubbed (for the
caching tests, via the ``run_calls`` fixture) or pointed at the interpreter to
drive its error handling.
"""

from __future__ import annotations

import gzip
import hashlib
import sys
from collections.abc import Callable
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.external import ToolNotFoundError
from genome.io import utils as utils_mod
from genome.io.utils import (
    ChecksumMismatchError,
    _gunzip,
    _is_fresh,
    _run,
    _run_to,
    sha256_file,
)


def test_sha256_file_digests_a_real_fixture(data_dir: Path) -> None:
    # The committed sacCer3 subsample, hashed independently with `shasum -a 256`.
    assert (
        sha256_file(data_dir / "tiny.fa")
        == "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"
    )


def test_sha256_of_a_gzipped_file_differs_from_its_contents(data_dir: Path) -> None:
    # Why the metadata table records the unpacked digest: the archive's is a different
    # number entirely, and changes whenever the file is recompressed.
    assert sha256_file(data_dir / "tiny.fa.gz") != sha256_file(data_dir / "tiny.fa")


def test_sha256_file_spans_many_read_buffers(tmp_path: Path) -> None:
    # Several megabytes, so the digest is assembled from many chunks rather than one.
    payload = b"ACGTN" * 1_000_000
    path = tmp_path / "big.fa"
    path.write_bytes(payload)

    assert sha256_file(path) == hashlib.sha256(payload).hexdigest()


def test_sha256_file_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        sha256_file(tmp_path / "nope.fa")


@pytest.fixture(scope="session")
def digest_scratch(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A reusable file path for the property test (session-scoped: hypothesis-safe)."""
    return tmp_path_factory.mktemp("digest") / "payload.bin"


@given(payload=st.binary(max_size=1 << 16))
def test_chunking_is_invisible(digest_scratch: Path, payload: bytes) -> None:
    # Streaming the file must give the digest of its bytes, whatever they are.
    digest_scratch.write_bytes(payload)
    assert sha256_file(digest_scratch) == hashlib.sha256(payload).hexdigest()


def test_checksum_mismatch_names_the_file_and_both_values() -> None:
    err = ChecksumMismatchError(Path("/data/hg38.fa"), "1a2b3c", "9f8e7d")

    assert isinstance(err, ValueError)
    message = str(err)
    assert "/data/hg38.fa" in message
    assert "1a2b3c" in message  # expected
    assert "9f8e7d" in message  # actual
    assert (err.expected, err.actual) == ("1a2b3c", "9f8e7d")


def test_gunzip_round_trips(tmp_path: Path) -> None:
    payload = b"".join(b">seq%d\nACGTACGTNN\n" % i for i in range(1000))
    src = tmp_path / "data.gz"
    with gzip.open(src, "wb") as fh:
        fh.write(payload)

    dest = tmp_path / "data"
    result = _gunzip(src, dest)

    assert result == dest
    assert dest.read_bytes() == payload


def test_is_fresh_rules(tmp_path: Path, touch_newer_than: Callable[..., None]) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    assert _is_fresh(out, [src]) is False  # missing output
    out.write_text("")
    assert _is_fresh(out, [src]) is False  # empty output
    out.write_text("y")
    touch_newer_than(out, src)
    assert _is_fresh(out, [src]) is True  # non-empty and newer
    touch_newer_than(src, out)
    assert _is_fresh(out, [src]) is False  # input regenerated -> stale


def test_run_to_runs_when_output_missing(
    tmp_path: Path, run_calls: list[tuple[str, list[str]]]
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    result = _run_to("tool", ["build", str(out)], out, [src])

    assert result == out
    assert run_calls == [("tool", ["build", str(out)])]


def test_run_to_skips_when_output_fresh(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    result = _run_to("tool", ["build"], out, [src])

    assert result == out
    assert run_calls == []  # served from the fresh cache


def test_run_to_overwrite_forces_run(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    _run_to("tool", ["build"], out, [src], overwrite=True)

    assert run_calls == [("tool", ["build"])]


def test_run_to_reruns_when_input_is_newer(
    tmp_path: Path,
    run_calls: list[tuple[str, list[str]]],
    touch_newer_than: Callable[..., None],
) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(src, out)  # input regenerated after the output

    _run_to("tool", ["build"], out, [src])

    assert run_calls == [("tool", ["build"])]


def test_run_raises_when_tool_missing() -> None:
    with pytest.raises(ToolNotFoundError):
        _run("definitely-not-a-real-tool-xyz", [])


def test_run_wraps_nonzero_exit_in_runtime_error(monkeypatch: pytest.MonkeyPatch) -> None:
    # Resolve to the interpreter and run a snippet that fails with known stderr;
    # _run should surface that stderr in an actionable RuntimeError.
    monkeypatch.setattr(utils_mod, "_resolve", lambda _name: sys.executable)

    with pytest.raises(RuntimeError, match="boom"):
        _run("python", ["-c", "import sys; sys.stderr.write('boom'); sys.exit(1)"])
