"""Shared fixtures for the test suite."""

from __future__ import annotations

import os
import shutil
import socket
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NoReturn, Protocol

import pytest
import requests

from genome import Genome
from genome.io import fetch as fetch_mod
from genome.io import utils as utils_mod
from genome.io.fasta import PREPARATION_TOOLS


class NetworkAccessError(RuntimeError):
    """Raised when a test attempts to reach the network."""


#: Address families that leave this machine. ``AF_UNIX`` is local IPC, never the
#: network, so the guard delegates to the real call for it.
_NETWORK_FAMILIES = frozenset({socket.AF_INET, socket.AF_INET6})

#: What to do instead — carried by every blocked call, since an error that only says
#: "no network" leaves the reader to work out how to write the test offline.
_OFFLINE_HELP = (
    "No test may reach the network. Serve a download offline with the `fake_fetch` "
    "fixture, which replaces the package's one fetch step "
    "(genome.io.fetch.fetch_url) with a copy from tests/data. A code path that "
    "also validates an assembly name at UCSC needs `requests.head` stubbed as well — "
    "see `head_recorder` in tests/test_download.py."
)


def _blocked(call: str, target: str) -> NetworkAccessError:
    """Return the error a blocked call raises: what was attempted, where, and the fix."""
    return NetworkAccessError(f"blocked network call: {call} {target}\n\n{_OFFLINE_HELP}")


def _address_text(address: Any) -> str:
    """Render a socket address as ``host:port`` when it has that shape."""
    if isinstance(address, tuple) and len(address) >= 2:
        return f"{address[0]}:{address[1]}"
    return str(address)


@pytest.fixture(autouse=True)
def no_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make reaching the network a failure — in every test, without being asked for.

    "No test ever touches the network" is a guarantee only if nothing can slip past
    it, so this is autouse and cuts in two places:

    - ``requests.Session.request``, which every ``requests`` call funnels through: the
      ``requests.head`` of the UCSC assembly-name check, and the ``requests.get`` inside
      pooch's HTTP downloader — so a download a test forgot to replace trips here,
      naming the URL it was about to pull.
    - ``socket.socket.connect``/``connect_ex`` for the internet address families, as the
      backstop under everything that is not ``requests``: pooch's ftp and sftp
      transports, urllib, a dependency phoning home. ``AF_UNIX`` is local IPC rather
      than the network, and is left alone.

    ``pooch.retrieve`` is deliberately **not** blocked, though it is the package's only
    transport: it serves a file already sitting at the destination without any network
    call, which two tests in test_download exercise for real. Cutting at the transport
    underneath keeps those honest and still catches every download pooch would actually
    perform.

    Nothing opts out. A test that needs a download stands one in with ``fake_fetch``; a
    test that needs the UCSC name check stubs ``requests.head`` itself. Both are
    installed after this one — an autouse fixture is set up first — so they win for the
    test that asked, and the guard is back for the one that did not.
    """
    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex

    def blocked_request(
        self: requests.Session, method: Any, url: Any, *args: Any, **kwargs: Any
    ) -> NoReturn:
        raise _blocked(f"requests {str(method).upper()}", str(url))

    def blocked_connect(sock: socket.socket, address: Any) -> None:
        if sock.family in _NETWORK_FAMILIES:
            raise _blocked("socket connect", _address_text(address))
        real_connect(sock, address)

    def blocked_connect_ex(sock: socket.socket, address: Any) -> int:
        if sock.family in _NETWORK_FAMILIES:
            raise _blocked("socket connect_ex", _address_text(address))
        return real_connect_ex(sock, address)

    monkeypatch.setattr(requests.sessions.Session, "request", blocked_request)
    monkeypatch.setattr(socket.socket, "connect", blocked_connect)
    monkeypatch.setattr(socket.socket, "connect_ex", blocked_connect_ex)


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
    """Offline stand-in for ``genome.io.fetch.fetch_url``.

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

    Every download in the package goes through ``genome.io.fetch.fetch_url``, so
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
    monkeypatch.setattr(fetch_mod, "fetch_url", fake)
    return fake


# ---------------------------------------------------------------------------------------
# Tiny chimera components
# ---------------------------------------------------------------------------------------

#: Where the component fixtures live. ``tests/data/README.md`` carries their table in
#: prose, :data:`CHIMERA_COMPONENTS` carries it as data, and ``test_chimera_fixtures``
#: asserts the committed bytes against the data — so the two cannot drift.
CHIMERA_DATA_DIR = DATA_DIR / "chimera"

#: What a component's own annotation is registered as. Deliberately colourless: what a
#: *merged* annotation is named is the merge's decision, and no fixture should pre-empt it.
COMPONENT_ANNOTATION = "genes"


@dataclass(frozen=True)
class ChimeraComponent:
    """One tiny component assembly and the committed bytes it is cut from.

    What each one is in the set *for* is ``tests/data/README.md``'s table; what is here is
    only what a test can check the bytes against.

    Attributes
    ----------
    name : str
        The component assembly name — alphanumeric, as the model requires.
    slices : dict
        Each chromosome, in file order, mapped to the ``(chromosome, start, end)`` slice
        of ``tests/data/tiny.fa`` it was cut from — 1-based inclusive, as ``samtools
        faidx`` takes them. Slices never overlap, across the whole set.
    line_width : int
        The column the FASTA is wrapped at.
    soft_masked : tuple or None
        The one lower-cased stretch, as ``(chromosome, bases)`` counted from its start, or
        ``None`` for a component that carries no masking.
    has_gtf : bool
        Whether an annotation is shipped beside the FASTA.
    """

    name: str
    slices: dict[str, tuple[str, int, int]]
    line_width: int
    soft_masked: tuple[str, int] | None = None
    has_gtf: bool = True

    @property
    def fasta(self) -> Path:
        """Path to the committed component FASTA."""
        return CHIMERA_DATA_DIR / f"{self.name}.fa"

    @property
    def gtf(self) -> Path | None:
        """Path to the committed component GTF, or ``None`` when it ships none."""
        return CHIMERA_DATA_DIR / f"{self.name}.gtf" if self.has_gtf else None

    @property
    def chromosomes(self) -> list[str]:
        """Chromosome names, in the order the FASTA declares them."""
        return list(self.slices)

    @property
    def lengths(self) -> dict[str, int]:
        """Each chromosome's length, derived from the slice it was cut from."""
        return {chromosome: end - start + 1 for chromosome, (_, start, end) in self.slices.items()}


#: Every tiny component, keyed by name and held in the sorted order a chimera puts its
#: components in — so reading this table is reading the derived name left to right.
CHIMERA_COMPONENTS: dict[str, ChimeraComponent] = {
    "tinyCe": ChimeraComponent(
        name="tinyCe",
        slices={
            "I": ("chrII", 2701, 5200),
            "II": ("chrIII", 6301, 8650),
            "MtDNA": ("chrI", 6901, 9200),
        },
        line_width=60,
        soft_masked=("I", 200),
    ),
    "tinyEc": ChimeraComponent(
        name="tinyEc",
        slices={
            "NZ_TINY01000001.1": ("chrIII", 2601, 6300),
            "NZ_TINY01000002.1": ("chrI", 2401, 4600),
            "chr1_KI270706v1_random": ("chrII", 5701, 7700),
        },
        line_width=80,
    ),
    "tinyEcDub": ChimeraComponent(
        name="tinyEcDub",
        slices={
            "NZ_TINY02000001.1": ("chrII", 7701, 8300),
            "NZ_TINY02__000002.1": ("chrII", 8301, 8800),
        },
        line_width=60,
        has_gtf=False,
    ),
    "tinySc": ChimeraComponent(
        name="tinySc",
        slices={
            "I": ("chrI", 1, 2400),
            "II": ("chrII", 1, 2700),
            "III": ("chrIII", 1, 2600),
        },
        line_width=60,
    ),
}

#: The everyday components: three, sorted, and asking for the ``__`` separator every
#: shipped assembly asks for. A pair from these mirrors the shipped ``ce11_ecHT115``;
#: all three make N > 2 an ordinary case rather than an unexercised one.
CHIMERA_EVERYDAY: tuple[str, ...] = ("tinyCe", "tinyEc", "tinySc")

#: The component that forces the separator to escalate. Kept out of
#: :data:`CHIMERA_EVERYDAY` so the everyday chimera is the ordinary one and escalation is
#: opted into.
CHIMERA_ESCALATION = "tinyEcDub"


class ComponentFactory(Protocol):
    """Registers a tiny component assembly and hands back the opened :class:`Genome`."""

    def __call__(self, name: str, *, with_annotation: bool = False) -> Genome:
        """Register ``name`` — and its annotation when asked — under the test's tmp dir."""
        ...


@pytest.fixture
def chimera_component(tmp_path: Path) -> Iterator[ComponentFactory]:
    """Return a factory that registers a tiny component assembly and opens it.

    Each component lands in its own directory under ``tmp_path`` and is seeded from its
    committed FASTA, so registration runs the real native tools and the real record path
    without reaching the network — which the autouse guard would refuse anyway. Every
    genome opened is closed at teardown::

        def test_something(chimera_component):
            worm = chimera_component("tinyCe", with_annotation=True)
            draft = chimera_component("tinyEc")

    Skips when the preparation tools are not on ``PATH``, so a test using this needs no
    skip marker of its own. See ``tests/data/README.md`` for what each one exercises.
    """
    missing = [tool for tool in PREPARATION_TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")

    opened: list[Genome] = []

    def register(name: str, *, with_annotation: bool = False) -> Genome:
        component = CHIMERA_COMPONENTS[name]
        genome = Genome(
            name,
            path_or_url=component.fasta,
            cache_dir=tmp_path / name,
            progressbar=False,
        )
        opened.append(genome)
        if with_annotation:
            if component.gtf is None:
                raise ValueError(
                    f"{name} ships no annotation; ask for one of "
                    f"{[c.name for c in CHIMERA_COMPONENTS.values() if c.has_gtf]} instead."
                )
            genome.register_gtf(component.gtf, COMPONENT_ANNOTATION)
        return genome

    yield register
    for genome in opened:
        genome.close()
