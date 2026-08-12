"""Shared fixtures for the test suite."""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pytest

from genome.io import download as download_mod
from genome.io import utils as utils_mod

#: Committed fixture files — small, subsampled real sacCer3 bytes. See tests/data/README.md.
DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def data_dir() -> Path:
    """Return the directory of committed fixture files (``tests/data``)."""
    return DATA_DIR


@pytest.fixture
def touch_newer_than() -> Callable[..., None]:
    """Return a helper that sets ``path``'s mtime ``delta`` seconds after ``reference``'s."""

    def _touch(path: Path, reference: Path, *, delta: float = 10.0) -> None:
        mtime = reference.stat().st_mtime + delta
        os.utime(path, (mtime, mtime))

    return _touch


@pytest.fixture
def run_calls(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, list[str]]]:
    """Record (and suppress) every ``genome.io.utils._run`` call so caching is observable."""
    calls: list[tuple[str, list[str]]] = []

    def fake_run(name: str, args: Sequence[str]) -> None:
        calls.append((name, list(args)))

    monkeypatch.setattr(utils_mod, "_run", fake_run)
    return calls


@dataclass(frozen=True)
class FetchCall:
    """One recorded call to the replaced fetch step, with every argument it received."""

    url: str
    dest_dir: Path
    known_hash: str | None
    fname: str | None
    processor: object | None
    progressbar: bool


class FakeFetch:
    """Offline stand-in for ``genome.io.download.fetch_url``.

    Copies a file out of ``tests/data`` instead of downloading it, then applies whatever
    pooch processor the caller passed, so a caller sees the same path shape a real fetch
    produces (``fname`` under ``dest_dir``, or the processor's output). Every call is
    recorded in ``calls``, so a test can assert on the URL the package built and the
    hash it expected.

    Serves ``tiny.fa`` unless told otherwise; call ``serve("tiny.fa.gz")`` (or any other
    name under ``tests/data``, or an absolute path) to choose a different file. A caller
    that passes no ``fname`` gets the served file's own name.
    """

    def __init__(self, source: str | Path = "tiny.fa", *, data_dir: Path = DATA_DIR) -> None:
        self.data_dir = data_dir
        self.calls: list[FetchCall] = []
        self.source: Path = self._locate(source)

    def serve(self, source: str | Path) -> None:
        """Serve ``source`` — a name under ``tests/data``, or an absolute path — from now on."""
        self.source = self._locate(source)

    def _locate(self, source: str | Path) -> Path:
        path = Path(source)
        located = path if path.is_absolute() else self.data_dir / path
        if not located.is_file():
            raise FileNotFoundError(f"no such test fixture file: {located}")
        return located

    @property
    def last(self) -> FetchCall:
        """The most recent recorded call."""
        return self.calls[-1]

    def __call__(
        self,
        url: str,
        dest_dir: Path,
        *,
        known_hash: str | None = None,
        fname: str | None = None,
        processor: Callable[..., object] | None = None,
        progressbar: bool = True,
    ) -> Path:
        """Copy the served fixture to ``dest_dir`` as if it had been downloaded."""
        dest_dir = Path(dest_dir)
        self.calls.append(
            FetchCall(
                url=url,
                dest_dir=dest_dir,
                known_hash=known_hash,
                fname=fname,
                processor=processor,
                progressbar=progressbar,
            )
        )
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / (fname or self.source.name)
        shutil.copy2(self.source, dest)
        if processor is None:
            return dest
        return Path(str(processor(str(dest), "download", None)))


@pytest.fixture
def fake_fetch(monkeypatch: pytest.MonkeyPatch) -> FakeFetch:
    """Replace the package's one fetch step with an offline copy from ``tests/data``.

    Every download in the package goes through ``genome.io.download.fetch_url``, so
    patching that one name takes the whole package offline. Use this fixture for any
    test whose code path would otherwise download something::

        def test_something(fake_fetch, tmp_path):
            fake_fetch.serve("tiny.fa.gz")        # what the server "has"
            ...                                   # exercise the caller
            assert fake_fetch.last.url.endswith("sacCer3.fa.gz")
            assert fake_fetch.last.known_hash is None

    See :class:`FakeFetch` for what it records and how to point it at another fixture.
    """
    fake = FakeFetch()
    monkeypatch.setattr(download_mod, "fetch_url", fake)
    return fake
