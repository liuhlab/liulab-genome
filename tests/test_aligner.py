"""Tests for genome.aligner — the Aligner abstraction, STAR, chromap, and the mixin.

An aligner is *given* its **External tool**, so the pure-logic tests hand it a
:class:`~genome.external.RecordingTool` and construct nothing else — no resolution to
patch out, no version call to intercept, and no binary anywhere. A handful of integration
tests build a real index — from a toy FASTA + GTF, and from a chimera assembled out of
the tiny components — and those are the only tests in this suite needing a binary the
package does not ship: :func:`_needs` is how one says so.

What a finished index is, is asked of its completion record and of nothing else, so
most of what is worth asserting here — that an unbuilt index raises, that an
interrupted one raises differently, that a damaged one names the file that changed,
and what the record says about the build — needs no aligner installed at all.
"""

from __future__ import annotations

import re
import shutil
import types
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from typing import TypeVar, cast

import pandas as pd
import pytest

import genome.aligner.aligner as aligner_mod
import genome.external as external_mod
from genome import Genome
from genome.aligner.aligner import IndexNotBuiltError
from genome.aligner.chromap import Chromap
from genome.aligner.mixin import AlignerMixin, _resolve_aligner
from genome.aligner.star import STAR
from genome.external import RecordingTool, ToolCall, ToolNotFoundError
from genome.io.completion import (
    RECORD_NAME,
    CompletionRecord,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    write_record,
)
from genome.io.registration import AssemblyDir, assembly_data_dir

from .conftest import CHIMERA_COMPONENTS, CHIMERA_EVERYDAY, ComponentFactory

#: What a mark may be applied to. Bounded, so both marks below resolve to the overload
#: that hands the test function back rather than the one that builds a parametrised mark.
_Test = TypeVar("_Test", bound=Callable[..., None])


def _needs(binary: str) -> Callable[[_Test], _Test]:
    """Mark a test as needing ``binary``, and skip it when that binary is absent.

    One name applies both halves, because they must never come apart. The ``aligner``
    marker is what puts a test in the CI lane that installs STAR and chromap and proves
    they answer before selecting; the skip is what keeps the same test out of the way on
    a machine that has neither. A test carrying only the skip would sit in the other
    lane and report green having run nothing, which is the failure the split ends.
    """
    skip = pytest.mark.skipif(shutil.which(binary) is None, reason=f"{binary} not on PATH")

    def mark(test: _Test) -> _Test:
        return pytest.mark.aligner(skip(test))

    return mark


_needs_star = _needs("STAR")
_needs_chromap = _needs("chromap")

# A 10 kb single-chromosome genome — large enough for STAR to index, small
# enough to build in a second or two.
_SEQ = "ACGTACGTAGGCATCGATCG" * 500

# One gene, two exons (GTF is 1-based inclusive): exon1 [101, 300], exon2
# [601, 800]. STAR should derive the single intron chr1:301-600 (+).
_TOY_GTF = (
    'chr1\ttoy\texon\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chr1\ttoy\texon\t601\t800\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
)

#: Old bookkeeping this build no longer writes: a bare success flag and a parameter
#: sidecar beside it, both absorbed into the one completion record.
_RETIRED_FILES = (".success", "star.index.json", "chromap.index.json")


def _make_genome(
    tmp_path: Path,
    gtfs: dict[str, str] | None = None,
    chrom_sizes: pd.Series | None = None,
) -> Genome:
    """Return a minimal Genome-like stub backed by a real FASTA + GTF(s) on disk.

    ``gtfs`` maps annotation name -> GTF text; each is written to ``tmp_path`` and
    exposed through an ``annotations`` registry stand-in whose ``path`` matches
    :meth:`~genome.io.gtf.AnnotationRegistry.path`. ``default_gtf`` is the first of
    them, as a genome with one registered annotation and no table flag would have it.

    ``chrom_sizes`` replaces the single-chromosome default, which is what lets a
    test dictate the *shape* of a reference — how many sequences and how long —
    without writing one, since the index parameters STAR is given are computed
    from those lengths and from nothing else.
    """
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chr1\n" + _SEQ + "\n")

    gtf_text = gtfs or {"toy": _TOY_GTF}
    gtf_paths: dict[str, Path] = {}
    for name, content in gtf_text.items():
        path = tmp_path / f"{name}.gtf"
        path.write_text(content)
        gtf_paths[name] = path

    sizes = pd.Series({"chr1": len(_SEQ)}) if chrom_sizes is None else chrom_sizes
    stub = types.SimpleNamespace(
        assembly="tiny",
        # Where this genome lives, asked of it rather than re-derived: the stub says so
        # explicitly, exactly as a real Genome does, so a test that moved the data root
        # moves the index dir with it.
        assembly_dir=AssemblyDir.locate("tiny"),
        files=types.SimpleNamespace(fasta=fasta),
        chrom_sizes=sizes,
        annotations=types.SimpleNamespace(path=lambda name: gtf_paths[name]),
        default_gtf=next(iter(gtf_paths)),
    )
    return cast("Genome", stub)


class _Tools:
    """The recording stand-ins this module's stubbed aligners are built with.

    One object owns every call any of them made, in order, because a test asserts on
    *the* command that ran and does not care which tool object carried it. Building an
    aligner is then a constructor call and nothing else — the fifteen ``monkeypatch``
    calls that used to stand in for a binary's resolution and its version are what an
    aligner taking its tool removes.
    """

    def __init__(self) -> None:
        self.calls: list[list[str]] = []
        self.leaves_output = False
        self.exit_code = 0

    def __call__(self, binary: str) -> RecordingTool:
        """Return a recording tool for ``binary`` that reports to this collection."""
        tool = RecordingTool(binary, version="0.0-test", on_run=self._record)
        tool.exit_code = self.exit_code
        return tool

    def fail(self) -> None:
        """Make every tool handed out from now on exit non-zero, as a real one can."""
        self.exit_code = 1

    def _record(self, call: ToolCall) -> None:
        self.calls.append(list(call.args))
        if self.leaves_output:
            _leave_plausible_output(call)


def _leave_plausible_output(call: ToolCall) -> None:
    """Write what a real index build would have left in the directory it ran in.

    A build that writes nothing leaves a record claiming nothing, which cannot disagree
    with anything. Tests about a *damaged* index need a build whose record has real files
    to hold the directory to, so this writes the artifact the command line asks for plus
    the log the tool drops in its working directory — read off the call itself, exactly
    as the real tool reads it.
    """
    assert call.cwd is not None, "an aligner runs in its index dir"
    (call.cwd / "Log.out").write_text("tool log\n")
    args = list(call.args)
    if "--output" in args:  # chromap: one file, named on the command line
        Path(args[args.index("--output") + 1]).write_text("index bytes\n")
    else:  # STAR: the genomeDir it ran in is the artifact
        (call.cwd / "SA").write_text("suffix array bytes\n")


@pytest.fixture
def tools() -> _Tools:
    """The recording tools this test's aligners are built with, and every call they made."""
    return _Tools()


def _star_over(tools: _Tools, tmp_path: Path, chrom_sizes: pd.Series) -> STAR:
    """A stubbed STAR bound to ``toy`` over a genome of exactly ``chrom_sizes``."""
    return STAR(_make_genome(tmp_path, chrom_sizes=chrom_sizes), gtf="toy", tool=tools("STAR"))


def _shape_of(*components: str) -> pd.Series:
    """The ``chrom_sizes`` of a genome shaped like those tiny components together.

    Only the sequence count and the total length reach the parameters under test, so
    the names here are placeholders — what a chimera calls its sequences is the naming
    contract's business, not this one's.
    """
    lengths = [
        length
        for component in components
        for length in CHIMERA_COMPONENTS[component].lengths.values()
    ]
    return pd.Series(lengths, index=[f"seq{i}" for i in range(len(lengths))])


def _flag_value(args: Sequence[str], flag: str) -> str:
    """The single value ``flag`` was given in ``args``, asserting it was given once."""
    assert list(args).count(flag) == 1, f"{flag} appears {list(args).count(flag)}x in {args}"
    return args[list(args).index(flag) + 1]


def _record_of(aligner: aligner_mod.Aligner) -> CompletionRecord:
    """Return the completion record written into ``aligner``'s index directory."""
    record = read_record(aligner.index_dir)
    assert record is not None, f"no {RECORD_NAME} in {aligner.index_dir}"
    return record


@pytest.fixture
def stub_star(tools: _Tools, tmp_path: Path) -> STAR:
    """A STAR (bound to the ``toy`` annotation) driving a recording tool."""
    return STAR(_make_genome(tmp_path), gtf="toy", tool=tools("STAR"))


@pytest.fixture
def stub_chromap(tools: _Tools, tmp_path: Path) -> Chromap:
    """A Chromap driving a recording tool (no annotation — chromap needs none)."""
    return Chromap(_make_genome(tmp_path), tool=tools("chromap"))


@pytest.fixture
def captured_run(tools: _Tools) -> list[list[str]]:
    """Every argument list the stubbed aligners' tools were run with, in order."""
    return tools.calls


@pytest.fixture
def building_run(tools: _Tools) -> list[list[str]]:
    """As ``captured_run``, but each run leaves plausible output files behind."""
    tools.leaves_output = True
    return tools.calls


# -- pure logic -------------------------------------------------------------


def test_one_renderer_spells_each_aligners_long_options_and_expands_list_values(
    stub_star: STAR, stub_chromap: Chromap
) -> None:
    # There is one renderer, and the only thing an aligner varies is the character its
    # long options put between words — declared as a class attribute, not as a body of
    # its own. Handing both the same keywords is what shows the difference is only that.
    kwargs = {"min_frag_length": 30, "read_format": ["r1", "bc"]}

    assert stub_star._flags(kwargs) == ["--min_frag_length", "30", "--read_format", "r1", "bc"]
    assert stub_chromap._flags(kwargs) == [
        "--min-frag-length",
        "30",
        "--read-format",
        "r1",
        "bc",
    ]
    # A list value expands after its own flag, whatever the flag itself is called.
    flags = stub_star._flags({"genomeSAindexNbases": 11, "genomeFastaFiles": ["a.fa", "b.fa"]})
    assert flags == ["--genomeSAindexNbases", "11", "--genomeFastaFiles", "a.fa", "b.fa"]


def test_constructing_an_aligner_runs_nothing_but_indexing_raises_naming_what_installs_it(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Nothing on PATH and no interpreter bin/ either: an aligner that resolved its binary
    # in its constructor could not be built here at all, which is what made every test
    # above patch resolution out. The instructions then travel in the error rather than
    # being printed to stderr on the way past: a library that writes to a console its
    # caller may not have is not an error message, and the caller cannot catch what it
    # cannot see. Driven for both aligners, so neither's message is a copy that quietly
    # stopped matching its own tool.
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external_mod.sys, "executable", str(tmp_path / "bin" / "python"))

    star = STAR(_make_genome(tmp_path), gtf="toy")
    assert star.index_dir.name == "star_toy"  # construction alone runs nothing

    with pytest.raises(ToolNotFoundError) as raised:
        star.index()
    message = str(raised.value)
    assert "STAR is not installed" in message
    assert "pixi add star" in message
    assert "bioconda" in message
    assert message == star.install_instructions()

    chromap = Chromap(_make_genome(tmp_path))
    with pytest.raises(ToolNotFoundError, match="pixi add chromap"):
        chromap.index()
    assert "chromap is not installed" in chromap.install_instructions()


def test_index_dir_naming_is_per_annotation_and_rooted_in_the_genomes_own_assembly_dir(
    stub_star: STAR,
    tools: _Tools,
    tmp_path: Path,
    chimera_component: ComponentFactory,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The annotation key names the genomeDir, so different GTFs never collide.
    assert stub_star.index_dir.parts[-4:] == ("genome", "tiny", "index", "star_toy")

    genome = _make_genome(tmp_path, {"a": _TOY_GTF, "b": _TOY_GTF})
    star_a = STAR(genome, gtf="a", tool=tools("STAR"))
    star_b = STAR(genome, gtf="b", tool=tools("STAR"))
    assert star_a.index_dir != star_b.index_dir
    assert star_a.index_dir.name == "star_a"
    assert star_b.index_dir.name == "star_b"

    # An **Index dir** is `<assembly dir>/index/<name>/`, and the assembly dir is the one
    # this genome was opened in — never one re-derived from the shared data root. Moving
    # the root out from under an already-open genome is the case that tells the two apart:
    # every other test's genome sits under the root, where the two answers agree. The
    # digest guard reads the same assembly dir, and re-deriving it would read a record
    # that is not there — a reference rebuilt underneath the index would stop being
    # noticed, silently, which is what makes this worth pinning down alongside it.
    chimera = chimera_component("tinyCe")
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "elsewhere"))
    assembly = read_record(chimera.fasta_path.parent)
    assert assembly is not None
    assert assembly.sha256 is not None

    chromap = Chromap(chimera, tool=tools("chromap"))
    assert chromap.index_dir == chimera.fasta_path.parent / "index" / "chromap"

    chromap.index()
    assert _record_of(chromap).details[_DIGEST_KEY] == assembly.sha256


# -- asking for an index that is not there ----------------------------------


def test_index_path_raises_for_every_kind_of_unfinished_or_damaged_index(
    stub_star: STAR, stub_chromap: Chromap, building_run: list[list[str]]
) -> None:
    with pytest.raises(IndexNotBuiltError) as raised:
        _ = stub_star.index_path
    message = str(raised.value)
    assert "nothing has been built" in message
    assert str(stub_star.index_dir) in message
    assert "Genome.build_star_index(gtf='toy')" in message

    # Index files with no record: a run that died between writing them and vouching
    # for them. Not the same as never built, and not silently rebuilt either.
    stub_star.index_dir.mkdir(parents=True)
    (stub_star.index_dir / "SA").write_text("half a suffix array\n")
    with pytest.raises(UnfinishedRegistrationError) as raised:
        _ = stub_star.index_path
    message = str(raised.value)
    assert "SA" in message
    assert "Genome.build_star_index(gtf='toy', overwrite=True)" in message

    # A directory prepared by an older version: the flag is not read, so it cannot
    # vouch for anything either, and raises exactly as any other unfinished build.
    (stub_star.index_dir / "SA").unlink()
    (stub_star.index_dir / ".success").touch()
    with pytest.raises(UnfinishedRegistrationError):
        _ = stub_star.index_path
    (stub_star.index_dir / ".success").unlink()

    # A real build, then damaged two different ways: changed in place (STAR) and
    # deleted outright (chromap) — the base class's guard, exercised through both
    # concrete aligners so neither could have special-cased the other's shape.
    stub_star.index()
    (stub_star.index_dir / "SA").write_text("truncated\n")
    with pytest.raises(RegistrationMismatchError) as raised:
        _ = stub_star.index_path
    message = str(raised.value)
    assert "SA" in message
    assert "Genome.build_star_index(gtf='toy', overwrite=True)" in message

    built = stub_chromap.index()
    built.unlink()
    with pytest.raises(RegistrationMismatchError, match=re.escape("chromap.index")):
        _ = stub_chromap.index_path


# -- what a finished build writes -------------------------------------------


def test_index_writes_a_complete_record_and_renders_each_knob_exactly_once(
    stub_star: STAR, stub_chromap: Chromap, captured_run: list[list[str]], tmp_path: Path
) -> None:
    out = stub_star.index(threads=3)

    assert out == stub_star.index_dir == stub_star.index_path
    assert len(captured_run) == 1

    record = _record_of(stub_star)
    assert record.kind == "index"
    assert record.name == "star_toy"
    assert record.package_version
    assert record.completed_at
    assert record.tool_versions == {"STAR": "0.0-test"}
    assert record.details["fasta"] == str(tmp_path / "tiny.fa")
    assert record.details["assembly"] == "tiny"
    assert record.details["aligner"] == "star"

    command = record.details["command"]
    assert command[0] == "STAR"
    assert command[1:] == captured_run[0]
    assert "genomeGenerate" in command
    assert "--genomeFastaFiles" in command

    parameters = record.details["parameters"]
    assert parameters["threads"] == 3
    assert parameters["gtf"] == "toy"
    # Small genome -> a reduced suffix-array size is auto-added.
    assert "genomeSAindexNbases" in parameters

    # `parameters` is every tuning knob that determined the build — the four STAR writes
    # under its own flag names included — and the command line is rendered from the
    # caller's keywords alone. Recording and rendering read the same dict in the old
    # code; a knob appearing under both spellings is what that confusion looked like.
    stub_star.index(threads=3, sjdb_overhang=49, overwrite=True)
    parameters = _record_of(stub_star).details["parameters"]
    assert {"threads", "gtf", "sjdb_gtf_file", "sjdb_overhang"} <= set(parameters)
    args = captured_run[-1]
    assert _flag_value(args, "--runThreadN") == "3"
    assert _flag_value(args, "--sjdbOverhang") == "49"
    assert not {"--threads", "--gtf", "--sjdb_gtf_file", "--sjdb_overhang"} & set(args)

    # The same contract from the other side: chromap now spells its two minimizer knobs
    # on the command line and records them, so rendering the recorded dict as well would
    # emit each flag a second time. `_flag_value` asserts the flag appears exactly once.
    stub_chromap.index(kmer=20, window=10, min_frag_length=25)
    chromap_parameters = _record_of(stub_chromap).details["parameters"]
    assert chromap_parameters == {"kmer": 20, "window": 10, "min_frag_length": 25}
    chromap_args = captured_run[-1]
    assert _flag_value(chromap_args, "--kmer") == "20"
    assert _flag_value(chromap_args, "--window") == "10"
    assert _flag_value(chromap_args, "--min-frag-length") == "25"


def test_each_aligners_record_claims_every_file_and_writes_none_of_the_retired_ones(
    stub_star: STAR, stub_chromap: Chromap, building_run: list[list[str]]
) -> None:
    # An index is claimed whole — the log the tool dropped in its working directory
    # included — so a file going missing later is caught, not just the artifact.
    stub_star.index()
    record = _record_of(stub_star)
    assert set(record.files) == {"SA", "Log.out"}
    assert record.files["SA"] == (stub_star.index_dir / "SA").stat().st_size
    assert RECORD_NAME not in record.files  # the record never claims itself

    built = stub_chromap.index()
    assert _record_of(stub_chromap).files["chromap.index"] == built.stat().st_size

    # Old bookkeeping this build no longer writes: a bare success flag and a parameter
    # sidecar beside it, both absorbed into the one completion record above.
    for retired in _RETIRED_FILES:
        assert not (stub_star.index_dir / retired).exists()
        assert not (stub_chromap.index_dir / retired).exists()


# -- reuse, repair and rebuild ----------------------------------------------


def test_index_is_reused_unless_overwrite_forces_a_rebuild(
    stub_star: STAR, tools: _Tools, tmp_path: Path, building_run: list[list[str]]
) -> None:
    stub_star.index()
    stub_star.index()  # a valid record -> reused, no rebuild
    assert len(building_run) == 1

    stub_star.index(overwrite=True)
    assert len(building_run) == 2

    # A different, half-built index: a plain rebuild refuses it, and only overwrite
    # repairs it — reuse and refusal are the same rule seen from its two sides.
    interrupted = STAR(_make_genome(tmp_path, {"b": _TOY_GTF}), gtf="b", tool=tools("STAR"))
    interrupted.index_dir.mkdir(parents=True)
    (interrupted.index_dir / "SA").write_text("half a suffix array\n")

    with pytest.raises(UnfinishedRegistrationError, match=re.escape("overwrite=True")):
        interrupted.index()
    assert len(building_run) == 2  # the interrupted one was left alone

    out = interrupted.index(overwrite=True)
    assert len(building_run) == 3
    assert out == interrupted.index_path

    # A rebuild that dies leaves nothing vouching for the directory. A second STAR over
    # the already-finished `stub_star` assembly, this one exiting non-zero — the failure
    # comes out of the tool's own error path rather than a raise patched over the call.
    tools.fail()
    dying = STAR(_make_genome(tmp_path), gtf="toy", tool=tools("STAR"))
    with pytest.raises(RuntimeError, match="STAR failed"):
        dying.index(overwrite=True)

    # The earlier record was dropped before the tool ran, so what is left reads as
    # interrupted rather than as a finished index whose sizes happen to still match.
    with pytest.raises(UnfinishedRegistrationError):
        _ = stub_star.index_path


def test_reusing_a_star_index_never_resolves_the_annotation_that_named_it(
    mixin_genome: Genome, tools: _Tools, tmp_path: Path
) -> None:
    # Reuse short-circuits before anything composes a command line, and composing STAR's
    # is what resolves the annotation — so a finished index is handed back without the
    # annotation being looked up at all. The template is what could quietly have moved
    # this: hand `_build` a command line already built and the lookup overtakes the
    # reuse check, turning a returned path into a raise. It composes lazily so that
    # cannot happen, and this pins the order from both sides.
    registered = {"toy": tmp_path / "toy.gtf"}
    mixin_genome.annotations = types.SimpleNamespace(  # type: ignore[misc]
        path=registered.__getitem__
    )

    built = mixin_genome.build_star_index("toy", tool=tools("STAR"))

    registered.clear()

    assert mixin_genome.build_star_index("toy", tool=tools("STAR")) == built
    assert len(tools.calls) == 1  # reused: STAR was not run a second time
    # The read-only way in never composed a command line to begin with.
    assert mixin_genome.get_star_index("toy") == built


def test_index_emits_sjdb_flags_from_bound_gtf_and_an_explicit_bin_size_wins(
    stub_star: STAR, captured_run: list[list[str]], tmp_path: Path
) -> None:
    stub_star.index(sjdb_overhang=49)

    gtf_path = str((tmp_path / "toy.gtf").resolve())
    args = captured_run[0]
    assert "--sjdbGTFfile" in args
    assert gtf_path in args
    assert "--sjdbOverhang" in args
    assert "49" in args

    parameters = _record_of(stub_star).details["parameters"]
    assert parameters["gtf"] == "toy"
    assert parameters["sjdb_gtf_file"] == gtf_path
    assert parameters["sjdb_overhang"] == 49

    # This genome would compute 13 for its bin size; what the caller asked for wins.
    stub_star.index(genomeChrBinNbits=16, overwrite=True)
    assert _flag_value(captured_run[-1], "--genomeChrBinNbits") == "16"
    assert _record_of(stub_star).details["parameters"]["genomeChrBinNbits"] == 16


# -- index parameters computed from the reference's shape -------------------


@pytest.mark.parametrize(
    ("components", "sequences", "total", "expected"),
    [
        (CHIMERA_EVERYDAY, 9, 22_750, 11),
        (("tinyCe",), 3, 7_150, 11),
    ],
    ids=["everyday-chimera", "one-component"],
)
def test_chr_bin_nbits_is_computed_from_the_shape_of_the_reference(
    tools: _Tools,
    tmp_path: Path,
    captured_run: list[list[str]],
    components: tuple[str, ...],
    sequences: int,
    total: int,
    expected: int,
) -> None:
    # The shapes being asked about, spelled out so the expected value is readable.
    sizes = _shape_of(*components)
    assert (len(sizes), int(sizes.sum())) == (sequences, total)

    star = _star_over(tools, tmp_path, sizes)
    star.index()

    assert _flag_value(captured_run[0], "--genomeChrBinNbits") == str(expected)
    assert _record_of(star).details["parameters"]["genomeChrBinNbits"] == expected


def test_chr_bin_nbits_stays_within_its_bounds_for_extreme_shapes(
    tools: _Tools, tmp_path: Path, captured_run: list[list[str]]
) -> None:
    # Two 1 Gb sequences: the recommendation is ~29.9, and the clause is a min, so 18
    # stands. It is still passed, so a record's parameters mean one thing either way.
    star = _star_over(tools, tmp_path, pd.Series({"chr1": 10**9, "chr2": 10**9}))
    star.index()
    assert _flag_value(captured_run[0], "--genomeChrBinNbits") == "18"
    assert _record_of(star).details["parameters"]["genomeChrBinNbits"] == 18

    # 40 scaffolds of 50 bp: the average sequence is far shorter than a read, so the
    # read length is what the bin is sized from — never the 5 that log2(50) would give
    # — the read length being sjdb_overhang + 1. Two overhangs, straddling a log2
    # truncation boundary: log2(100) -> 6 and log2(150) -> 7, so a truncation bug at one
    # magnitude cannot hide behind the other. A distinct annotation key each time, not
    # just a distinct shape: every star here shares one assembly dir, and reusing "toy"
    # would hand back an earlier build instead of a new one.
    scaffolds = pd.Series({f"scaffold{i}": 50 for i in range(40)})
    for overhang, expected in ((99, "6"), (149, "7")):
        genome = _make_genome(
            tmp_path, gtfs={f"scattered{overhang}": _TOY_GTF}, chrom_sizes=scaffolds
        )
        scattered = STAR(genome, gtf=f"scattered{overhang}", tool=tools("STAR"))
        scattered.index(sjdb_overhang=overhang)
        assert _flag_value(captured_run[-1], "--genomeChrBinNbits") == expected
        assert _record_of(scattered).details["parameters"]["genomeChrBinNbits"] == int(expected)


# -- mixin: build/get index entry points ------------------------------------


@pytest.fixture
def mixin_genome(tmp_path: Path) -> Genome:
    """A Genome-like object carrying :class:`AlignerMixin`.

    No patch anywhere: the mixin forwards ``tool=`` to the aligner it builds, so a test
    that wants a recording stand-in hands one to the ``build_*`` call. Constructing an
    aligner without one still runs nothing, which is why the read-only entry points here
    need no tool at all.
    """
    stub = _make_genome(tmp_path)

    class _MixinGenome(AlignerMixin):
        pass

    genome = _MixinGenome()
    for attr in ("assembly", "assembly_dir", "files", "chrom_sizes", "annotations", "default_gtf"):
        setattr(genome, attr, getattr(stub, attr))
    return cast("Genome", genome)


def test_resolve_aligner_and_get_index_shortcuts_all_agree(
    mixin_genome: Genome, tools: _Tools
) -> None:
    assert _resolve_aligner("star") is STAR
    assert _resolve_aligner("STAR") is STAR
    assert _resolve_aligner("chromap") is Chromap
    assert _resolve_aligner("CHROMAP") is Chromap
    with pytest.raises(ValueError, match="Unknown aligner 'bowtie'"):
        _resolve_aligner("bowtie")

    star_built = mixin_genome.build_star_index("toy", tool=tools("STAR"))
    assert mixin_genome.get_index("star", gtf="toy") == star_built
    assert mixin_genome.get_star_index("toy") == star_built

    chromap_built = mixin_genome.build_chromap_index(tool=tools("chromap"))
    assert mixin_genome.get_index("chromap") == chromap_built
    assert mixin_genome.get_chromap_index() == chromap_built
    assert chromap_built.name == "chromap.index"


def test_default_annotation_is_used_when_present_needed_by_neither_chromap_nor_when_absent(
    mixin_genome: Genome, tools: _Tools
) -> None:
    # The everyday call: a chimera, and any assembly the table flags one for, already
    # carries a default, so spelling it out at the call site said nothing the genome did
    # not already know. The index it names is the same one either way.
    named = mixin_genome.build_star_index("toy", tool=tools("STAR"))

    defaulted = mixin_genome.build_star_index(tool=tools("STAR"))

    assert mixin_genome.default_gtf == "toy"
    assert defaulted == named
    assert defaulted.name == "star_toy"
    assert len(tools.calls) == 1  # the second call reused the first's finished index

    # The other half of the same rule: chromap carries no annotation at all, so a genome
    # with no default indexes for it perfectly well.
    mixin_genome.default_gtf = None  # type: ignore[misc]
    assert mixin_genome.build_chromap_index(tool=tools("chromap")).name == "chromap.index"

    # A STAR index is built against one annotation and cannot be built against none, so
    # the refusal is a ValueError naming the per-call fix and the once-and-for-all one.
    calls_before = len(tools.calls)
    with pytest.raises(ValueError, match="no default") as excinfo:
        mixin_genome.build_star_index(tool=tools("STAR"))
    message = str(excinfo.value)
    assert "build_star_index(gtf=<name>)" in message
    assert "default_gtf=<name>" in message
    assert len(tools.calls) == calls_before  # nothing ran


def test_the_tool_is_forwarded_to_the_aligner_and_never_read_as_a_selector(
    mixin_genome: Genome, tools: _Tools
) -> None:
    # The seam `Aligner.__init__` offers is open the whole way from `Genome`. It used to
    # be closed here — the mixin built the aligner and handed it nothing — so the only
    # way to stand a binary in was to patch the fallback out from under it.
    built = mixin_genome.build_star_index("toy", tool=tools("STAR"))
    mixin_genome.build_chromap_index(tool=tools("chromap"))

    assert len(tools.calls) == 2
    assert "genomeGenerate" in tools.calls[0]
    assert "--build-index" in tools.calls[1]

    # `get_index`'s remaining keywords pin down *which* index; `tool` is not one of them
    # and is spelled out so it cannot be mistaken for one.
    assert mixin_genome.get_index("star", gtf="toy", tool=tools("STAR")) == built


def test_get_index_raises_before_a_build_and_for_an_unknown_aligner(mixin_genome: Genome) -> None:
    with pytest.raises(IndexNotBuiltError, match=re.escape("Genome.build_star_index(gtf='toy')")):
        mixin_genome.get_index("star", gtf="toy")

    with pytest.raises(IndexNotBuiltError, match=re.escape("Genome.build_chromap_index()")):
        mixin_genome.get_chromap_index()

    with pytest.raises(ValueError, match="Unknown aligner 'bowtie'"):
        mixin_genome.get_index("bowtie", gtf="toy")


# -- integration (require a real STAR) --------------------------------------


@_needs_star
def test_real_star_index_builds(tmp_path: Path) -> None:
    star = STAR(_make_genome(tmp_path), gtf="toy")

    out = star.index(threads=2)

    assert out == star.index_path  # reopening a fresh index does not raise
    assert out.name == "star_toy"
    names = {p.name for p in out.iterdir()}
    assert {"SA", "SAindex", "Genome", "genomeParameters.txt"} <= names
    assert RECORD_NAME in names
    assert not names & set(_RETIRED_FILES)

    # Every binary file STAR emitted is claimed, so losing one is caught.
    record = _record_of(star)
    assert {"SA", "SAindex", "Genome"} <= set(record.files)
    assert record.tool_versions["STAR"] == star.version

    # The command a real STAR is given carries the bin size computed for this
    # fixture's shape — 10 kb in one sequence asks for 2^13, not STAR's 2^18 — and
    # a real STAR accepts it, which is what the stubbed tests cannot show.
    assert record.details["parameters"]["genomeChrBinNbits"] == 13
    assert _flag_value(record.details["command"], "--genomeChrBinNbits") == "13"
    echoed = [
        line.split("\t")[1]
        for line in (out / "genomeParameters.txt").read_text().splitlines()
        if line.startswith("genomeChrBinNbits\t")
    ]
    assert echoed == ["13"]


@_needs_star
def test_real_star_index_with_gtf(tmp_path: Path) -> None:
    star = STAR(_make_genome(tmp_path), gtf="toy")

    out = star.index(sjdb_overhang=49, threads=2)

    # GTF-derived annotation files are produced.
    assert (out / "geneInfo.tab").is_file()
    sjdb = (out / "sjdbList.out.tab").read_text()
    # The single intron between the two exons: chr1:301-600.
    assert "301" in sjdb
    assert "600" in sjdb

    parameters = _record_of(star).details["parameters"]
    assert parameters["gtf"] == "toy"
    assert parameters["sjdb_gtf_file"] == str((tmp_path / "toy.gtf").resolve())


# ===========================================================================
# Chromap — a second aligner: one index per assembly, no annotation, a file
# (not a directory) as the artifact, and hyphenated long flags.
# ===========================================================================


# -- pure logic -------------------------------------------------------------


def test_chromap_index_dir_naming_default_record_reuse_and_kmer_window_flags(
    stub_chromap: Chromap, building_run: list[list[str]], tmp_path: Path
) -> None:
    # No annotation -> no per-GTF suffix; one chromap index serves the assembly.
    assert stub_chromap.index_dir.parts[-3:] == ("tiny", "index", "chromap")

    with pytest.raises(IndexNotBuiltError) as raised:
        _ = stub_chromap.index_path
    # No annotation selects a chromap index, so the build call takes no arguments.
    assert "Genome.build_chromap_index()" in str(raised.value)

    out = stub_chromap.index()

    # The artifact is a single file inside the index dir, not the dir itself.
    assert out == stub_chromap.index_path
    assert out.name == "chromap.index"
    assert out.parent == stub_chromap.index_dir
    assert len(building_run) == 1

    record = _record_of(stub_chromap)
    assert record.kind == "index"
    assert record.name == "chromap"
    assert record.tool_versions == {"chromap": "0.0-test"}
    assert record.details["assembly"] == "tiny"

    command = record.details["command"]
    assert command[0] == "chromap"
    assert "--build-index" in command
    assert "--ref" in command
    assert "--output" in command

    # With no tuning kwargs, the command is exactly build-index over ref -> output.
    fasta = str(tmp_path / "tiny.fa")
    artifact = str(stub_chromap.index_dir / "chromap.index")
    assert building_run[0] == ["--build-index", "--ref", fasta, "--output", artifact]

    stub_chromap.index()  # a valid record -> reused, no rebuild
    assert len(building_run) == 1

    # Overwrite forces a rebuild, tuned this time with its two minimizer knobs.
    stub_chromap.index(overwrite=True, kmer=20, window=10, min_frag_length=25)
    assert len(building_run) == 2
    args = building_run[-1]
    assert args[args.index("--kmer") + 1] == "20"
    assert args[args.index("--window") + 1] == "10"
    assert args[args.index("--min-frag-length") + 1] == "25"

    parameters = _record_of(stub_chromap).details["parameters"]
    assert parameters["kmer"] == 20
    assert parameters["window"] == 10
    assert parameters["min_frag_length"] == 25


# -- integration (require a real chromap) -----------------------------------


@_needs_chromap
def test_real_chromap_index_builds(tmp_path: Path) -> None:
    chromap = Chromap(_make_genome(tmp_path))

    out = chromap.index()

    assert out == chromap.index_path  # reopening a fresh index does not raise
    assert out.name == "chromap.index"
    assert out.is_file()  # chromap's index is a single file, not a directory
    assert out.stat().st_size > 0

    names = {p.name for p in chromap.index_dir.iterdir()}
    assert RECORD_NAME in names
    assert not names & set(_RETIRED_FILES)

    record = _record_of(chromap)
    assert record.files["chromap.index"] == out.stat().st_size
    assert record.tool_versions["chromap"] == chromap.version


@_needs_chromap
def test_real_chromap_accepts_the_minimizer_knobs_it_now_spells_itself(tmp_path: Path) -> None:
    # kmer and window moved out of the rendered keywords and onto the command line, so
    # a real chromap is what says the two spellings it is given are ones it accepts —
    # the stubbed tests assert the string and would pass on a flag chromap rejects.
    chromap = Chromap(_make_genome(tmp_path))

    out = chromap.index(kmer=17, window=7)

    assert out.is_file()
    assert out.stat().st_size > 0

    record = _record_of(chromap)
    assert record.details["parameters"] == {"kmer": 17, "window": 7}
    assert _flag_value(record.details["command"], "--kmer") == "17"
    assert _flag_value(record.details["command"], "--window") == "7"


# ===========================================================================
# The reference an index was built from. The guard lives in the base class, so
# every case below runs against both aligners: chromap is the one that cannot
# notice a renamed chromosome any other way, its index storing no sequence
# names at all, and STAR is the one that could have been special-cased.
# ===========================================================================

#: What the *assembly* pinned when the index was built, and what it pins after being
#: registered again over different bytes. Shaped like real digests so a message
#: quoting both reads the way it will on disk.
_DIGEST_BUILT_FROM = "a1" * 32
_DIGEST_AFTER_REBUILD = "b2" * 32

#: The key an index's ``details`` pins the assembly's digest under. Asserted as a
#: literal because it is an on-disk contract someone reads by eye months later.
_DIGEST_KEY = "assembly_sha256"


def _pin_assembly_digest(aligner: aligner_mod.Aligner, digest: str | None) -> Path:
    """Write the *assembly*'s completion record, one directory up, pinning ``digest``.

    Registering an assembly again over different bytes is exactly this: a record in
    the assembly dir whose ``sha256`` is a different string. Nothing in the index
    directory is touched, which is the whole failure — for chromap the index file
    stays byte-identical while the names it will be used with change.
    """
    directory = assembly_data_dir(aligner.assembly)
    write_record(
        directory,
        CompletionRecord(
            kind="genome",
            name=aligner.assembly,
            files={},
            source_url=None,
            sha256=digest,
            tool_versions={},
            package_version="0.0-test",
            completed_at="2026-01-01T00:00:00+00:00",
            details={},
        ),
    )
    return directory


@pytest.fixture(params=["stub_star", "stub_chromap"])
def stub_aligner(request: pytest.FixtureRequest) -> tuple[aligner_mod.Aligner, str]:
    """Each concrete aligner in turn, with the forced-rebuild call it must name."""
    repair = {
        "stub_star": "Genome.build_star_index(gtf='toy', overwrite=True)",
        "stub_chromap": "Genome.build_chromap_index(overwrite=True)",
    }[request.param]
    return cast("aligner_mod.Aligner", request.getfixturevalue(request.param)), repair


def test_a_rebuilt_reference_is_caught_raised_and_repaired_only_by_overwrite(
    stub_aligner: tuple[aligner_mod.Aligner, str], building_run: list[list[str]]
) -> None:
    # The lifecycle a digest pin goes through: pinned on a fresh build, still trusted
    # when the assembly is registered again over the same bytes, and a hard stop the
    # moment those bytes actually change — raised on sight, and refused a second time
    # by `index()` itself so a caller cannot rebuild over the mismatch by accident.
    aligner, repair = stub_aligner
    _pin_assembly_digest(aligner, _DIGEST_BUILT_FROM)

    aligner.index()
    details = _record_of(aligner).details
    assert details[_DIGEST_KEY] == _DIGEST_BUILT_FROM
    assert details["fasta"].endswith("tiny.fa")  # the path stays beside the digest

    # Registered again over the same bytes — a repeat of `--force`, not a new reference.
    _pin_assembly_digest(aligner, _DIGEST_BUILT_FROM)
    built = aligner.index_path
    aligner.index()
    assert len(building_run) == 1  # still finished, so still reused

    # The defect this closes: every file the index claims is present at the size it
    # claims, and the reference underneath it is a different genome.
    _pin_assembly_digest(aligner, _DIGEST_AFTER_REBUILD)
    with pytest.raises(RegistrationMismatchError) as raised:
        _ = aligner.index_path
    message = str(raised.value)
    assert _DIGEST_BUILT_FROM in message
    assert _DIGEST_AFTER_REBUILD in message
    assert repair in message

    with pytest.raises(RegistrationMismatchError, match=re.escape("overwrite=True")):
        aligner.index()
    assert len(building_run) == 1  # rebuilding is a deliberate act, as everywhere else

    # The repair the message names must not be blocked by the mismatch it repairs.
    out = aligner.index(overwrite=True)
    assert len(building_run) == 2
    assert out == aligner.index_path == built
    assert _record_of(aligner).details[_DIGEST_KEY] == _DIGEST_AFTER_REBUILD


def test_no_digest_pinned_anywhere_means_no_guard_and_nothing_recorded(
    stub_aligner: tuple[aligner_mod.Aligner, str], building_run: list[list[str]]
) -> None:
    # Absent, not null: a fact that could not be gathered is left out, exactly as
    # tool_versions leaves out a tool that would not identify itself.
    aligner, _ = stub_aligner

    # Built with nothing pinned at all: the record carries no digest...
    built = aligner.index()
    assert _DIGEST_KEY not in _record_of(aligner).details
    # ...and a digest that shows up only afterwards does not retroactively guard it.
    _pin_assembly_digest(aligner, _DIGEST_AFTER_REBUILD)
    assert aligner.index_path == built

    # Symmetric case: a digest that was pinned and then withdrawn — registered again
    # pinning nothing to compare — stops guarding just the same.
    _pin_assembly_digest(aligner, _DIGEST_BUILT_FROM)
    rebuilt = aligner.index(overwrite=True)
    _pin_assembly_digest(aligner, None)
    assert aligner.index_path == rebuilt  # no raise; nothing to compare against


# ===========================================================================
# A chimera is indexed like any other assembly, and that is the claim. There
# is no chimera-specific code here to exercise — the reference is prepared,
# the merged annotation is registered, and an index build sees an assembly —
# so what these tests hold is that the general path survives the reference a
# chimera actually is: nine sequences whose names carry their component, and
# bytes copied from three components without being rewrapped or recased.
# ===========================================================================

#: What the everyday three merge to. Nobody chose this name: each component registers
#: its annotation under the same colourless key, and the merge joins the contributions
#: in component-sorted order — so the derived name is what picks the index directory,
#: `+` included, which no shell has to quote.
_MERGED_ANNOTATION = "genes+genes+genes"

#: The assembly the everyday three build, likewise derived from the component set.
_CHIMERA_ASSEMBLY = "tinyCe_tinyEc_tinySc"

#: The nine chromosomes it carries, in FASTA order: components sorted, then each
#: component's own file order. Spelled out rather than assembled from the fixture table
#: by the same rule the build follows, which would agree with a build that was wrong.
_CHIMERA_CHROMOSOMES = [
    "I__tinyCe",
    "II__tinyCe",
    "MtDNA__tinyCe",
    "NZ_TINY01000001.1__tinyEc",
    "NZ_TINY01000002.1__tinyEc",
    "chr1_KI270706v1_random__tinyEc",
    "I__tinySc",
    "II__tinySc",
    "III__tinySc",
]

#: One gene from each component's own annotation. The fixture slices are disjoint
#: set-wide, so each of these could only have reached the merge from where it says.
_GENE_PER_COMPONENT = {"tinyCe": "YBL111C", "tinyEc": "YCL074W", "tinySc": "YAL069W"}

#: How many genes the three annotations carry between them — 3 + 3 + 7 — which is the
#: count STAR writes on the first line of ``geneInfo.tab``.
_CHIMERA_GENES = 13


@pytest.fixture
def everyday_chimera(chimera_component: ComponentFactory) -> Iterator[Genome]:
    """The everyday chimera, built for real, with its merged annotation registered.

    The components and the chimera land beside each other under the test's own root,
    which is what a chimera's staleness comparison looks for.
    Where the *index* goes no longer depends on it: an index dir and the assembly digest
    it pins both come from the genome's own **Assembly dir**, so a chimera built anywhere
    indexes into itself.

    Skips with :func:`chimera_component` when the preparation tools are absent, so a
    test using this needs no marker of its own beyond the aligner it drives.
    """
    components = [chimera_component(name, with_annotation=True) for name in CHIMERA_EVERYDAY]
    with Genome.chimera(*components) as built:
        yield built


@pytest.fixture
def stub_star_over_chimera(everyday_chimera: Genome, tools: _Tools) -> STAR:
    """A STAR over the everyday chimera, bound to its merged annotation, binary stood in for."""
    return STAR(everyday_chimera, gtf=_MERGED_ANNOTATION, tool=tools("STAR"))


# -- the command a chimera is indexed with ----------------------------------


def test_the_everyday_chimeras_index_is_named_shaped_and_pinned_from_its_own_build(
    everyday_chimera: Genome, stub_star_over_chimera: STAR, captured_run: list[list[str]]
) -> None:
    # The per-annotation layout is the general one; what is new is only that both halves
    # of the path are derived — the assembly from its components, the annotation from
    # what those components contributed.
    assert everyday_chimera.default_gtf == _MERGED_ANNOTATION
    # Nine sequences and 22,750 bases is the whole of what reaches these two knobs:
    # log2(22750) / 2 - 1 -> 6 for the suffix array, and log2(22750 / 9) -> 11 for the
    # bin, the mean sequence being far longer than the read sjdb_overhang implies.
    sizes = everyday_chimera.chrom_sizes
    assert (len(sizes), int(sizes.sum())) == (9, 22_750)

    stub_star_over_chimera.index()

    assert stub_star_over_chimera.index_dir.parts[-4:] == (
        "genome",
        _CHIMERA_ASSEMBLY,
        "index",
        f"star_{_MERGED_ANNOTATION}",
    )
    assert _flag_value(captured_run[0], "--genomeDir") == str(stub_star_over_chimera.index_dir)
    assert _flag_value(captured_run[0], "--genomeSAindexNbases") == "6"
    assert _flag_value(captured_run[0], "--genomeChrBinNbits") == "11"

    record = _record_of(stub_star_over_chimera)
    assert record.name == f"star_{_MERGED_ANNOTATION}"
    # Nothing in the record is chimera-specific: the derived assembly name already
    # carries the components, so a key listing them again could only disagree with it.
    assert record.details["assembly"] == _CHIMERA_ASSEMBLY
    assert "components" not in record.details
    parameters = record.details["parameters"]
    assert (parameters["genomeSAindexNbases"], parameters["genomeChrBinNbits"]) == (6, 11)

    # A chimera is registered like any other assembly, so the digest guard needs no case
    # of its own: it is read from the record the chimera build wrote.
    assembly = read_record(everyday_chimera.fasta_path.parent)
    assert assembly is not None
    assert assembly.sha256 is not None
    assert record.details[_DIGEST_KEY] == assembly.sha256


# -- integration: a real STAR over a real chimera ---------------------------


@_needs_star
def test_real_star_index_over_a_chimera_writes_every_suffixed_name(
    everyday_chimera: Genome,
) -> None:
    # Where the fixtures' prefix traps meet a real tool. The list is asserted whole
    # because membership would pass against a build that dropped a suffix: `I` is a
    # strict prefix of `II` and of `III`, and `tinyEc` of `tinyEcDub`, so a subtly wrong
    # name still reads as present inside a right one.
    star = STAR(everyday_chimera, gtf=_MERGED_ANNOTATION)

    out = star.index(threads=2)

    assert out == star.index_path  # reopening a fresh index does not raise
    assert out.name == f"star_{_MERGED_ANNOTATION}"
    names = (out / "chrName.txt").read_text().splitlines()
    assert names == _CHIMERA_CHROMOSOMES
    # Said out loud, because it is the mistake the equality above exists to refuse: an
    # unsuffixed name is what makes a read's species unattributable.
    assert "I" not in names


@_needs_star
def test_real_star_index_over_a_chimera_carries_the_merged_annotation(
    everyday_chimera: Genome,
) -> None:
    # Every component's genes reach the index, and the one splice junction any of them
    # carries sits on a suffixed name — so STAR read the merged GTF against the merged
    # FASTA and the two agreed, which is the whole of what the merge is for.
    star = STAR(everyday_chimera, gtf=_MERGED_ANNOTATION)

    out = star.index(threads=2)

    genes = (out / "geneInfo.tab").read_text().splitlines()
    assert genes[0] == str(_CHIMERA_GENES)  # a count, then one line per gene
    assert set(_GENE_PER_COMPONENT.values()) <= {line.split("\t")[0] for line in genes[1:]}

    # tinyCe's YBL111C is the set's only two-exon transcript: exons [207, 1416] and
    # [1516, 2309] leave the intron STAR reports here, 1-based inclusive as a GTF is.
    junctions = [line.split("\t") for line in (out / "sjdbList.out.tab").read_text().splitlines()]
    assert junctions == [["I__tinyCe", "1417", "1515", "-"]]


@_needs_star
def test_a_real_star_accepts_the_bin_size_a_chimeras_shape_asks_for(
    everyday_chimera: Genome,
) -> None:
    # 11 rather than STAR's own 18, and a real STAR builds with it and echoes it back.
    # A computed bin size had never reached a real build over a many-sequence reference,
    # which is the shape it was computed for and the one a stub cannot refuse.
    star = STAR(everyday_chimera, gtf=_MERGED_ANNOTATION)

    out = star.index(threads=2)

    assert _record_of(star).details["parameters"]["genomeChrBinNbits"] == 11
    echoed = [
        line.split("\t")[1]
        for line in (out / "genomeParameters.txt").read_text().splitlines()
        if line.startswith("genomeChrBinNbits\t")
    ]
    assert echoed == ["11"]


# -- integration: a real chromap over the same chimera ----------------------


@_needs_chromap
def test_real_chromap_index_over_a_chimera_reads_a_heterogeneous_reference(
    everyday_chimera: Genome, capfd: pytest.CaptureFixture[str]
) -> None:
    # Not symmetry with STAR, and not the naming contract either — a chromap index
    # stores no sequence names, so it can say nothing about those. What only this test
    # reaches is chromap's own bundled kseq.h over the reference a chimera is: nine
    # sequences, wrapped at 60 columns for tinyCe and at 80 for tinyEc, soft-masked in
    # the first and not in the second. The other real chromap test builds over
    # _make_genome — one sequence, on one unwrapped line — so nothing else in this suite
    # has ever handed chromap a second sequence or a wrap at all, and copying component
    # bytes verbatim was verified through faidx, faToTwoBit and twoBitInfo, none of them
    # this parser. Deleting this leaves that property with no witness.
    chromap = Chromap(everyday_chimera)

    out = chromap.index()

    assert out == chromap.index_path
    assert out.is_file()
    assert out.stat().st_size > 0
    record = _record_of(chromap)
    assert record.files["chromap.index"] == out.stat().st_size
    assert record.details["assembly"] == _CHIMERA_ASSEMBLY

    # chromap's own count, read back off its stderr: every sequence and every base
    # arrived, rather than the parser merely not crashing. The wording is chromap
    # 0.3.2's, so a release that changes it fails here loudly instead of quietly
    # proving less.
    reported = capfd.readouterr().err
    assert "number of sequences: 9" in reported
    assert "number of bases: 22750" in reported
