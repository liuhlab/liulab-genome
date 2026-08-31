"""Shared fixtures for the test suite.

The two autouse guards are **not** here: they are in :mod:`tests._guards`, loaded as a
plugin by the root ``conftest.py``, which says why. Everything in this file is shared by
tests and only by tests.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import pytest
from hypothesis import settings

from genome import Genome
from genome.assembly.fasta import PREPARATION_TOOLS
from genome.assembly.registration import assembly_data_dir
from genome.external import ExternalTool, clear_version_cache
from genome.store import fetch as fetch_mod

from ._guards import install_network_guard

# Hypothesis measures a deadline in wall-clock time, and this lane runs ten workers deep
# with several modules spawning process pools of their own — so what the default 200 ms
# deadline reads here is the scheduler and not the code, which is the same reason a per-test
# wall-clock budget was declined outright:
# docs/research/test-suite-cost-and-parallelism-2026-08-30.md, where the rate is measured.
# Off for the whole suite rather than per test, since which test is unlucky is decided by
# scheduling and not by anything about the test. What a genuinely slow test costs stays
# visible: the CI lane prints `--durations=10` on every run.
settings.register_profile("suite", deadline=None)
settings.load_profile("suite")

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
    """Record (and suppress) every **External tool** invocation, so caching is observable.

    Patches :meth:`genome.external.ExternalTool.run` — the one method every invocation
    converges on, since :meth:`~genome.external.ExternalTool.run_to` reaches its binary
    through ``self.run`` rather than around it. One ``setattr`` on the base class
    therefore catches every tool driven by any adapter anywhere in the package, which is
    the property ``fetch_url`` is spelled for on the download side.

    Freshness is left running: ``run_to`` decides whether to call ``run`` before this
    stub is reached, so an empty list still means *the tool was skipped*.
    """
    calls: list[tuple[str, list[str]]] = []

    def fake_run(
        self: ExternalTool, args: Sequence[str], *, cwd: Path | None = None, capture: bool = True
    ) -> str:
        calls.append((self.name, list(args)))
        return ""

    monkeypatch.setattr(ExternalTool, "run", fake_run)
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
    """Offline stand-in for ``genome.store.fetch.fetch_url``.

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

    Every download in the package goes through ``genome.store.fetch.fetch_url``, so
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
        """Register ``name`` — and its annotation when asked — under the test's data root."""
        ...


class PreparedComponents:
    """Every tiny component registered once, kept to be copied rather than built again.

    Registering one is three native-tool runs over a few kilobytes, and the suite asks
    for a hundred and sixty-odd of them — the same four files, every time, from committed bytes that
    cannot differ. So each ``(name, annotated)`` pair is built at most once and the
    result is handed out as a *copy*: :meth:`directory` returns a path a caller must
    never write to, and :meth:`copy_into` is how a test gets one of its own.

    Nothing here is shared *mutable* state. The suite runs under ``--dist=load``, where
    which test runs beside which changes run to run, so a directory two tests both wrote
    to would fail in whichever order happened to expose it.
    """

    def __init__(self, root: Path) -> None:
        self._root = root
        self._built: dict[tuple[str, bool], Path] = {}

    def directory(self, name: str, *, with_annotation: bool) -> Path:
        """Return the prepared tree for ``name``, registering it on first ask. Read-only."""
        key = (name, with_annotation)
        if key not in self._built:
            self._built[key] = self._prepare(name, with_annotation=with_annotation)
        return self._built[key]

    def copy_into(self, destination: Path, name: str, *, with_annotation: bool) -> None:
        """Copy the prepared tree for ``name`` to ``destination``, mtimes and all."""
        shutil.copytree(self.directory(name, with_annotation=with_annotation), destination)

    def _prepare(self, name: str, *, with_annotation: bool) -> Path:
        component = CHIMERA_COMPONENTS[name]
        path = self._root / ("annotated" if with_annotation else "plain") / name
        # The build runs before any test's setup, so the autouse guard is not up yet:
        # stand one up here rather than leave the one registration nobody watches.
        with pytest.MonkeyPatch.context() as guard:
            install_network_guard(guard)
            genome = Genome(name, path_or_url=component.fasta, cache_dir=path, progressbar=False)
            try:
                if with_annotation:
                    assert component.gtf is not None
                    genome.annotations.register_path(component.gtf, COMPONENT_ANNOTATION)
            finally:
                genome.close()
        return path


@pytest.fixture(scope="session")
def _prepared_components(tmp_path_factory: pytest.TempPathFactory) -> PreparedComponents:
    """Return the session's :class:`PreparedComponents`, in a directory no test can reach."""
    return PreparedComponents(tmp_path_factory.mktemp("prepared-components"))


@pytest.fixture
def chimera_component(_prepared_components: PreparedComponents) -> Iterator[ComponentFactory]:
    """Return a factory that registers a tiny component assembly and opens it.

    Each component lands where the layout puts one — ``<data root>/genome/<name>/``, the
    root being wherever it points when the factory is called, which the :func:`liulab_data`
    fixture makes this test's own directory. So a chimera built from two of them finds
    them beside itself, which is what the staleness comparison does and what a test
    placing each component in a directory of its own would never exercise. Every genome
    opened is closed at teardown::

        def test_something(chimera_component):
            worm = chimera_component("tinyCe", with_annotation=True)
            draft = chimera_component("tinyEc")

    The bytes are a copy of a registration :class:`PreparedComponents` ran once for the
    whole session, not a fresh build: the same four files either way, and the copy is the
    only thing a test is ever handed. Asking twice for the same component opens what is
    there rather than laying the prepared bytes back down over it — so a test that
    re-registers one underneath a chimera sees the new one, as it would on disk.

    Skips when the preparation tools are not on ``PATH``, so a test using this needs no
    skip marker of its own. See ``tests/data/README.md`` for what each one exercises.
    """
    missing = [tool for tool in PREPARATION_TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")

    opened: list[Genome] = []

    def register(name: str, *, with_annotation: bool = False) -> Genome:
        component = CHIMERA_COMPONENTS[name]
        if with_annotation and component.gtf is None:
            raise ValueError(
                f"{name} ships no annotation; ask for one of "
                f"{[c.name for c in CHIMERA_COMPONENTS.values() if c.has_gtf]} instead."
            )
        destination = assembly_data_dir(name)
        if not destination.exists():
            _prepared_components.copy_into(destination, name, with_annotation=with_annotation)
        genome = Genome(name, path_or_url=component.fasta, progressbar=False)
        opened.append(genome)
        if with_annotation and COMPONENT_ANNOTATION not in genome.annotations.registered:
            assert component.gtf is not None
            genome.annotations.register_path(component.gtf, COMPONENT_ANNOTATION)
        return genome

    yield register
    for genome in opened:
        genome.close()


# ---------------------------------------------------------------------------------------
# Stub binaries
# ---------------------------------------------------------------------------------------

#: Every stub binary in the suite is a link to this one script, which sources the ``.sh``
#: written beside the link — so the file that gets *executed* is the same file every time
#: and only the behaviour beside it changes.
#:
#: The indirection buys speed, not style. macOS runs a security check the first time a
#: newly created executable is exec'd, and it is not cheap: measured on this machine, the
#: first exec of a freshly written ``#!/bin/sh`` stub costs 200-900ms and swings four-fold
#: between runs, while a second exec of the same file costs ~6ms. Copying known content to
#: a new path pays it again, so the check is keyed on the file rather than on what is in
#: it. A link to a file that has already been through it is ~6ms. One scanned file per
#: worker, therefore, and the stubs a test installs are free.
#:
#: ``$0`` is the path that was exec'd — the link, not its target, and under a shebang the
#: kernel supplies it whether the caller ran the path directly or the shell found it on
#: ``PATH`` — so a stub's behaviour is looked up under the tool name the test chose.
#: Sourcing keeps the positional parameters, so a body reading ``"$@"`` sees the tool's
#: own arguments.
_STUB_DISPATCHER = '#!/bin/sh\n. "${0}.sh"\n'


class StubBinary(Protocol):
    """Installs one stub binary into a directory, under a name of the caller's choosing."""

    def __call__(self, bin_dir: Path, name: str, body: str, *, executable: bool = True) -> Path:
        """Write ``name`` into ``bin_dir`` running ``body``, and return the path to it."""
        ...


@pytest.fixture(scope="session")
def _stub_dispatcher(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Return the one executable every stub binary links to, exec'd once to warm it."""
    path = tmp_path_factory.mktemp("stub-dispatcher") / "stub"
    path.write_text(_STUB_DISPATCHER)
    path.chmod(0o755)
    Path(f"{path}.sh").write_text("exit 0\n")
    # Pay the platform's first-exec check here, once per worker, rather than in whichever
    # test happened to install a stub first.
    subprocess.run([str(path)], capture_output=True, check=True)
    return path


@pytest.fixture
def stub_binary(_stub_dispatcher: Path) -> StubBinary:
    """Return a helper that installs a stub binary — a shell body under a tool's name.

    A stub is a real binary as far as the package and the lane's guard script are
    concerned: it is found by ``shutil.which``, exec'd by ``PATH`` lookup, and its exit
    status and streams are its own. That is what makes it usable for the cases a developed
    machine cannot produce — a tool that runs and refuses to say what it is, one that
    resolves but cannot be executed::

        def test_something(tmp_path, stub_binary):
            stub_binary(tmp_path / "bin", "samtools", "echo 'samtools 1.21'")

    Pass ``executable=False`` for a file the shell can find but not run, which is a real
    file rather than a link: the shell answers 126 for it and never reaches the body.
    Installing the same name twice replaces it.
    """

    def install(bin_dir: Path, name: str, body: str, *, executable: bool = True) -> Path:
        bin_dir.mkdir(parents=True, exist_ok=True)
        path = bin_dir / name
        path.unlink(missing_ok=True)
        if not executable:
            path.write_text(f"#!/bin/sh\n{body}\n")
            path.chmod(0o644)
            return path
        path.symlink_to(_stub_dispatcher)
        Path(f"{path}.sh").write_text(f"{body}\n")
        return path

    return install


@pytest.fixture(autouse=True)
def _forget_tool_versions() -> Iterator[None]:
    """Empty the process-wide version cache around every test — in every test, unasked.

    ``genome.external`` remembers what each binary answered to ``--version`` for the life
    of the process, keyed on the path it was located at. That is right for a build and
    unsafe for a suite: many tests here put a stub tool on ``PATH``, and several point
    ``PATH`` somewhere else entirely, so an answer that outlived its test could be handed
    to the next one. Under ``--dist=load`` which test that is changes run to run, so the
    failure would be an order-dependent one — hence autouse, and hence no opt-out.

    Cleared on the way in as well as on the way out, so a test inherits nothing from
    whatever ran before it, test or not. The saving this gives up is only the *cross*-test
    one: within a test the four build steps of a chimera still ask each binary once.
    """
    clear_version_cache()
    yield
    clear_version_cache()
