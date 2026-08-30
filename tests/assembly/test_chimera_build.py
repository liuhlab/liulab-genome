"""Tests for genome.assembly.chimera_build and Genome.chimera — building a chimera, and reading one.

The naming contract has its own tests (``test_chimera``) and the components have theirs
(``test_chimera_fixtures``); what is asserted here is the build: that the bytes arrive
verbatim, that the names the native tools read back are the ones the contract predicts,
that the record says what was done, and that the two accessors answer from it.

Every build runs the real ``samtools``/``faToTwoBit``/``twoBitInfo`` over the tiny
components, offline, in a temporary data root — never the lab's own.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

from genome import Genome
from genome.assembly import chimera_build as chimera_mod
from genome.assembly.chimera import ChimeraNamingError, split_suffixed
from genome.assembly.chimera_build import ChimeraBuilder, _extend_header, read_chimera_details
from genome.assembly.download import assembly_data_dir
from genome.assembly.fasta import GenomeFiles
from genome.store import fetch as fetch_mod
from genome.store.checksum import sha256_file
from genome.store.completion import (
    CompletionRecord,
    RegistrationError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
)

from ..conftest import (
    CHIMERA_COMPONENTS,
    CHIMERA_ESCALATION,
    CHIMERA_EVERYDAY,
    ComponentFactory,
)

#: What the everyday three build: nine sequences, one contiguous block per component in
#: component-sorted order, each block in that component's own declared order.
_EVERYDAY_NAME = "tinyCe_tinyEc_tinySc"
_EVERYDAY_CHROMOSOMES = [
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

ChimeraFactory = Callable[..., Genome]


def _fasta_records(path: Path) -> dict[str, list[str]]:
    """Read a FASTA into ``{first header token: [raw sequence lines]}``, in file order.

    The lines are kept as they are written — case and wrapping included — because that is
    exactly what a chimera must not change.
    """
    records: dict[str, list[str]] = {}
    current: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = records.setdefault(line[1:].split()[0], [])
        else:
            current.append(line)
    return records


def _wrap_widths(lines: list[str]) -> set[int]:
    """The widths of every line but the last, which is short whenever the length divides."""
    return {len(line) for line in lines[:-1]}


def _record_of(genome: Genome) -> CompletionRecord:
    """The completion record written in ``genome``'s own assembly directory."""
    record = read_record(genome.fasta_path.parent)
    assert record is not None, f"no record in {genome.fasta_path.parent}"
    return record


@pytest.fixture
def build_chimera(chimera_component: ComponentFactory) -> Iterator[ChimeraFactory]:
    """Return a factory building a chimera of the named tiny components, in a temp root.

    No ``cache_dir``: the default path is the one worth exercising, and the shared
    ``liulab_data`` fixture has already pointed the root at the test's own directory —
    a chimera built into the lab's shared reference data by a test would be a serious
    bug. Every chimera opened is closed at teardown.
    """
    opened: list[Genome] = []

    def build(*names: str, cache_dir: str | Path | None = None, force: bool = False) -> Genome:
        components = [chimera_component(name) for name in names]
        chimera = Genome.chimera(*components, cache_dir=cache_dir, force=force)
        opened.append(chimera)
        return chimera

    yield build
    for chimera in opened:
        chimera.close()


# --------------------------------------------------------------------------------------
# The everyday chimera: it builds, it opens, and it is named as the contract says
# --------------------------------------------------------------------------------------


def test_the_everyday_chimera_builds_in_component_sorted_order_under_a_cache_dir_override(
    build_chimera: ChimeraFactory, tmp_path: Path, liulab_data: Path
) -> None:
    # Built from the reversed order on purpose: identity is the component set, so the
    # published layout contract — component-sorted blocks, each in its own declared order
    # — must hold whatever order the caller handed the components in. Load-bearing
    # off-repo: a consumer filtering one component's sequences back out of an alignment
    # header recovers a single-assembly header only because this holds, and a different
    # concatenation order would hand it a silently wrong header rather than a failure.
    chimera = build_chimera(*reversed(CHIMERA_EVERYDAY))

    assert chimera.assembly == _EVERYDAY_NAME
    assert chimera.chromosomes == _EVERYDAY_CHROMOSOMES
    assert int(chimera.chrom_sizes.sum()) == 22_750
    # Prepared like any other assembly: four files, and the name derived from the set.
    assert chimera.fasta_path.name == f"{_EVERYDAY_NAME}.fa"
    assert all(
        path.is_file()
        for path in (chimera.files.fai, chimera.files.twobit, chimera.files.chrom_sizes)
    )

    # And it lands wherever it is told to rather than always under the shared root.
    elsewhere = tmp_path / "elsewhere"
    pair = build_chimera("tinyCe", "tinySc", cache_dir=elsewhere)
    assert pair.fasta_path == elsewhere / "tinyCe_tinySc.fa"
    assert not (liulab_data / "genome" / "tinyCe_tinySc").exists()


# --------------------------------------------------------------------------------------
# Bytes copied verbatim — the property with no other failing case
# --------------------------------------------------------------------------------------


def test_component_bytes_survive_verbatim_including_wrapping_and_masking(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    built = _fasta_records(chimera.fasta_path)

    # every source line, byte for byte
    for name in CHIMERA_EVERYDAY:
        component = CHIMERA_COMPONENTS[name]
        source = _fasta_records(component.fasta)
        for chromosome, lines in source.items():
            assert built[f"{chromosome}__{name}"] == lines

    # 60 for tinyCe's shape and 80 for tinyEc's, side by side in one file — which is what
    # every tool that could have done this concatenation would have destroyed.
    assert _wrap_widths(built["I__tinyCe"]) == {60}
    assert _wrap_widths(built["NZ_TINY01000001.1__tinyEc"]) == {80}

    # tinyCe carries 200 soft-masked bases and tinyEc carries none, so a build that
    # upper-cased or lower-cased anything loses one of the two.
    worm = "".join(built["I__tinyCe"])
    assert worm[:200].islower()
    assert worm[200:].isupper()
    assert "".join(built["NZ_TINY01000002.1__tinyEc"]).isupper()


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (b">I\n", b">I__tinyCe\n"),
        (b">I some description\n", b">I__tinyCe some description\n"),
        (b">I", b">I__tinyCe"),  # a final header with no line ending
        # Whitespace after '>' is skipped by samtools faidx and faToTwoBit, which name
        # the sequence from the token after it — so `> desc` declares `desc`, and that is
        # what carries the suffix. The skipped bytes are written back where they were.
        (b"> desc\n", b"> desc__tinyCe\n"),
    ],
)
def test_only_the_name_a_header_gives_its_sequence_is_extended(
    line: bytes, expected: bytes
) -> None:
    # The suffix rides on the sequence's name because STAR and chromap both truncate a
    # header at the first whitespace; everything else on the line is description or
    # layout, and is written back untouched.
    assert _extend_header(line, "tinyCe", "__") == expected


@pytest.mark.parametrize("line", [b">", b">   \n"])
def test_a_header_that_names_no_sequence_is_refused(line: bytes) -> None:
    # There is nothing for the suffix to ride on, and writing '>__tinyCe' would file the
    # sequence under a name that is suffix and nothing else — which samtools reads as the
    # empty name and the naming contract refuses outright.
    with pytest.raises(ChimeraNamingError, match="names no sequence"):
        _extend_header(line, "tinyCe", "__")


def test_retrieved_sequence_agrees_with_its_component_however_the_header_was_shaped(
    chimera_component: ComponentFactory, tmp_path: Path, build_chimera: ChimeraFactory
) -> None:
    # End to end, because this is where reading the header differently from samtools shows
    # up: the component registers under `oddChr`, so the chimera must carry
    # `oddChr__tinyOdd` — and a build that suffixed the empty token instead writes
    # `__tinyOdd`, which the cross-check against the components' prediction refuses.
    fasta = tmp_path / "odd.fa"
    fasta.write_text(">  oddChr a description after two spaces\n" + "ACGT" * 15 + "\n")
    odd = Genome("tinyOdd", path_or_url=fasta, cache_dir=tmp_path / "tinyOdd", progressbar=False)

    assert odd.chromosomes == ["oddChr"]  # what samtools and faToTwoBit called it

    with odd, Genome.chimera(odd, chimera_component("tinySc")) as chimera:
        assert "oddChr__tinyOdd" in chimera.chromosomes
        assert chimera["oddChr__tinyOdd"] == odd["oddChr"]

    # And, for ordinary components, the whole 2bit path over concatenated bytes agrees
    # with the component it came from: same bases, same case, same offsets — including a
    # chromosome shaped like hg38's decoy name, which must round-trip like any other.
    everyday = build_chimera(*CHIMERA_EVERYDAY)
    worm = chimera_component("tinyCe")
    draft = chimera_component("tinyEc")
    assert everyday["I__tinyCe:0-100"] == worm["I:0-100"]
    assert everyday["I__tinyCe:0-100"].islower()  # soft-masking survived faToTwoBit
    assert everyday["MtDNA__tinyCe"] == worm["MtDNA"]
    assert everyday["chr1_KI270706v1_random__tinyEc:0-100"] == draft["chr1_KI270706v1_random:0-100"]


# --------------------------------------------------------------------------------------
# What the reference then answers
# --------------------------------------------------------------------------------------


def test_a_bare_chromosome_name_names_every_resolving_spelling_or_falls_back_like_any_assembly(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
    # Resolving a bare name would restore the ambiguity the suffix abolishes, and would
    # make what a name means depend on which components happen to be present (ADR-0009).
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    # `III` is carried by one component: the refusal names the spelling that does
    # resolve — it is the ninth of the nine sequences here, so a message merely listing
    # the first few never mentions it — and that spelling is a next action.
    with pytest.raises(ValueError, match="III__tinySc"):
        chimera["III:0-100"]
    assert len(chimera["III__tinySc:0-100"]) == 100

    # `I` is carried by two components — tinyCe and tinySc — so naming one spelling would
    # be picking for the caller the very thing the suffix exists to make them say; both
    # are named.
    with pytest.raises(ValueError, match="I__tinySc") as refusal:
        chimera["I:0-100"]
    assert "I__tinyCe" in str(refusal.value)

    # `chrZ` is carried by no component: a chimera has one more thing to say only about a
    # bare name it actually carries, so this gets the same general message any assembly
    # would, with no mention of ADR-0009.
    with pytest.raises(ValueError, match="known sequences include") as unknown:
        chimera["chrZ:0-5"]
    assert "ADR-0009" not in str(unknown.value)

    # And a plain assembly — not a chimera at all — has no suffixed spelling to offer for
    # a name it does not carry, so nothing here changes for it either.
    worm = chimera_component("tinyCe")
    with pytest.raises(ValueError, match="known sequences include") as plain:
        worm["III:0-100"]  # a real tinySc chromosome, unknown to tinyCe
    assert "ADR-0009" not in str(plain.value)


def test_two_components_that_collide_get_four_distinct_names_and_escalation_lengthens_the_run(
    build_chimera: ChimeraFactory,
) -> None:
    # tinyCe and tinySc both carry `I` and `II`; the shipped pair never does.
    chimera = build_chimera("tinyCe", "tinySc")
    assert chimera.assembly == "tinyCe_tinySc"
    assert {"I__tinyCe", "II__tinyCe", "I__tinySc", "II__tinySc"} <= set(chimera.chromosomes)
    assert len(set(chimera.chromosomes)) == len(chimera.chromosomes) == 6
    assert chimera["I__tinyCe"] != chimera["I__tinySc"]

    # The separator belongs to one chimera, so it is derived from these components and
    # recorded — a build that wrote a constant would lose the self-announcing property
    # exactly here, and nothing downstream could tell.
    escalated = build_chimera("tinyCe", CHIMERA_ESCALATION)
    details = read_chimera_details(escalated.fasta_path.parent)
    assert escalated.assembly == "tinyCe_tinyEcDub"
    assert details is not None
    assert details.separator == "___"
    # And the same fact through the accessor an off-repo caller has, which is the only one
    # that stops it hardcoding the default and splitting the doubled name in the wrong place.
    assert escalated.separator == "___"
    assert "NZ_TINY02__000002.1___tinyEcDub" in escalated.chromosomes
    assert escalated.chromosomes[0] == "I___tinyCe"
    # And the recorded separator is the one that reads the names back.
    assert [split_suffixed(name, details.separator) for name in escalated.chromosomes] == [
        ("I", "tinyCe"),
        ("II", "tinyCe"),
        ("MtDNA", "tinyCe"),
        ("NZ_TINY02000001.1", "tinyEcDub"),
        ("NZ_TINY02__000002.1", "tinyEcDub"),
    ]


# --------------------------------------------------------------------------------------
# components and chrom_components
# --------------------------------------------------------------------------------------


def test_a_chimeras_components_attribute_every_chromosome(build_chimera: ChimeraFactory) -> None:
    # Reversed on arrival, same as the layout test above: components is the sorted set.
    chimera = build_chimera(*reversed(CHIMERA_EVERYDAY))
    assert chimera.components == list(CHIMERA_EVERYDAY)

    attribution = chimera.chrom_components
    assert list(attribution.index) == chimera.chromosomes
    assert not attribution.isna().to_numpy().any()
    assert attribution["I__tinyCe"] == "tinyCe"
    assert attribution["I__tinySc"] == "tinySc"
    assert attribution["chr1_KI270706v1_random__tinyEc"] == "tinyEc"
    assert Counter(attribution) == {"tinyCe": 3, "tinyEc": 3, "tinySc": 3}


def test_a_plain_assembly_is_not_a_chimera_in_every_attribute_that_shows(
    chimera_component: ComponentFactory,
) -> None:
    worm = chimera_component("tinyCe")

    assert worm.components is None
    assert worm.separator is None
    assert worm.component_annotations is None
    # Totality is what keeps `components` the single is-chimera test: no caller has to
    # read `chrom_components` to find out which kind of assembly it is holding.
    attribution = worm.chrom_components
    assert list(attribution.index) == worm.chromosomes
    assert set(attribution) == {"tinyCe"}


# --------------------------------------------------------------------------------------
# The record, and what it vouches for
# --------------------------------------------------------------------------------------


def test_the_record_says_a_chimera_was_built_and_from_what(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    record = _record_of(chimera)
    digests = {name: _record_of(chimera_component(name)).sha256 for name in CHIMERA_EVERYDAY}

    assert (record.kind, record.name) == ("genome", _EVERYDAY_NAME)
    # Nothing was fetched, so nothing points anywhere — and the digest is of the bytes
    # this build wrote, which is what a later verification has to fall back to.
    assert record.source_url is None
    assert record.sha256 == sha256_file(chimera.fasta_path)
    assert record.details == {
        "separator": "__",
        "components": [{"name": name, "sha256": digests[name]} for name in CHIMERA_EVERYDAY],
    }
    # Record to record: the digests are the components' own, not bytes rehashed here.
    assert all(digest is not None for digest in digests.values())


def test_rebuilding_a_finished_chimera_rewrites_nothing_and_a_broken_one_is_repaired(
    build_chimera: ChimeraFactory,
) -> None:
    first = build_chimera(*CHIMERA_EVERYDAY)
    written = first.fasta_path.stat().st_mtime_ns
    record = read_record(first.fasta_path.parent)

    again = build_chimera(*CHIMERA_EVERYDAY)
    assert again.fasta_path == first.fasta_path
    assert again.fasta_path.stat().st_mtime_ns == written
    assert read_record(again.fasta_path.parent) == record

    # A directory that cannot be trusted is a different matter: the record naming what
    # finished there is gone, so a rebuild raises rather than guess, and names the
    # command that repairs it.
    pair = build_chimera("tinyCe", "tinySc")
    record_path(pair.fasta_path.parent).unlink()
    with pytest.raises(
        UnfinishedRegistrationError, match="genome assembly register tinyCe_tinySc --force"
    ):
        build_chimera("tinyCe", "tinySc")
    assert build_chimera("tinyCe", "tinySc", force=True).chromosomes == pair.chromosomes


def test_a_chrom_sizes_that_disagrees_with_the_components_is_never_recorded_and_nothing_is_fetched(
    chimera_component: ComponentFactory,
    build_chimera: ChimeraFactory,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The cross-check between the tools' answer and the contract's prediction. Standing a
    # wrong chrom.sizes in is the only way to reach it: the build is a verbatim copy under
    # a derived name, so it is meant to be unreachable. Scoped to its own MonkeyPatch
    # context rather than the fixture-provided one: the autouse `liulab_data` fixture
    # patches the data-root env var through that same shared instance, and undoing it
    # here would undo that redirect too, sending the rest of this test at the real one.
    prepare = chimera_mod.prepare_fasta

    def prepare_and_mangle(fasta: Path, *, overwrite: bool = False) -> GenomeFiles:
        files = prepare(fasta, overwrite=overwrite)
        lines = files.chrom_sizes.read_text().splitlines()
        name, length = lines[0].split("\t")
        lines[0] = f"{name}\t{int(length) - 1}"  # one base short, as a short read would be
        files.chrom_sizes.write_text("\n".join(lines) + "\n")
        return files

    components = [chimera_component(name) for name in ("tinyCe", "tinySc")]
    with pytest.MonkeyPatch.context() as scoped:
        scoped.setattr(chimera_mod, "prepare_fasta", prepare_and_mangle)
        with pytest.raises(
            RegistrationError, match="genome assembly register tinyCe_tinySc --force"
        ):
            Genome.chimera(*components)
    assert read_record(assembly_data_dir("tinyCe_tinySc")) is None

    # The autouse guard would raise on a real network call; this says deliberately that
    # the package's one fetch step is not reached at all, which is a stronger claim.
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("a chimera build fetches nothing")

    monkeypatch.setattr(fetch_mod, "fetch_url", refuse)
    assert build_chimera(*CHIMERA_EVERYDAY).components == list(CHIMERA_EVERYDAY)


# --------------------------------------------------------------------------------------
# What a chimera refuses to be built from
# --------------------------------------------------------------------------------------


def test_a_chimera_refuses_an_illegal_component_list(
    chimera_component: ComponentFactory, build_chimera: ChimeraFactory
) -> None:
    with pytest.raises(ChimeraNamingError, match="at least 2 components"):
        Genome.chimera(chimera_component("tinyCe"))

    worm = chimera_component("tinyCe")
    with pytest.raises(ChimeraNamingError, match="must not repeat"):
        Genome.chimera(worm, chimera_component("tinySc"), worm)

    # Nesting is forbidden by the model, and the record on the component's own disk is
    # what says so — the naming contract can only refuse the spelling.
    chimera = build_chimera("tinyCe", "tinySc")
    with pytest.raises(ChimeraNamingError, match="is itself a chimera"):
        Genome.chimera(chimera, chimera_component("tinyEc"))


def test_the_builder_derives_its_name_and_separator_before_anything_is_written(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # Constructing one is pure: the name, the order and the separator are settled from the
    # components' chromosome names, and the directory is not touched until a build.
    builder = ChimeraBuilder(
        [chimera_component("tinySc"), chimera_component("tinyCe")], tmp_path / "unbuilt"
    )

    assert builder.assembly == "tinyCe_tinySc"
    assert builder.separator == "__"
    assert [component.assembly for component in builder.components] == ["tinyCe", "tinySc"]
    assert not (tmp_path / "unbuilt").exists()
