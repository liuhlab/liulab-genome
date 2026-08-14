"""Tests for genome.io.chimera and Genome.chimera — building a chimera, and reading one.

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
from genome.chimera import ChimeraNamingError, split_suffixed
from genome.io import chimera as chimera_mod
from genome.io import fetch as fetch_mod
from genome.io.chimera import ChimeraBuilder, _extend_header, read_chimera_details
from genome.io.completion import (
    CompletionRecord,
    RegistrationError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
)
from genome.io.download import assembly_data_dir
from genome.io.fasta import GenomeFiles
from genome.io.utils import sha256_file

from .conftest import (
    CHIMERA_COMPONENTS,
    CHIMERA_ESCALATION,
    CHIMERA_EVERYDAY,
    ComponentFactory,
)

#: What the everyday three build: nine sequences, in component-sorted then file order.
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
def build_chimera(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Iterator[ChimeraFactory]:
    """Return a factory building a chimera of the named tiny components, in a temp root.

    ``LIULAB_DATA`` is pointed at the test's own directory first: a chimera built into
    the lab's shared reference data by a test would be a serious bug, and the default
    path is the one worth exercising. Every chimera opened is closed at teardown.
    """
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
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


def test_the_everyday_chimera_builds_and_opens_under_its_derived_name(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    assert chimera.assembly == _EVERYDAY_NAME
    assert chimera.chromosomes == _EVERYDAY_CHROMOSOMES
    assert int(chimera.chrom_sizes.sum()) == 22_750
    # Prepared like any other assembly: four files, and the name derived from the set.
    assert chimera.fasta_path.name == f"{_EVERYDAY_NAME}.fa"
    assert all(
        path.is_file()
        for path in (chimera.files.fai, chimera.files.twobit, chimera.files.chrom_sizes)
    )


def test_the_component_order_a_caller_passes_does_not_reach_the_reference(
    build_chimera: ChimeraFactory,
) -> None:
    # Identity is the component set, not the order it arrived in, so the reversed call
    # builds and opens the same directory rather than a second one.
    chimera = build_chimera(*reversed(CHIMERA_EVERYDAY))

    assert chimera.assembly == _EVERYDAY_NAME
    assert chimera.chromosomes == _EVERYDAY_CHROMOSOMES


def test_a_built_chimera_lands_where_it_is_told_to(
    build_chimera: ChimeraFactory, tmp_path: Path
) -> None:
    elsewhere = tmp_path / "elsewhere"
    chimera = build_chimera("tinyCe", "tinySc", cache_dir=elsewhere)

    assert chimera.fasta_path == elsewhere / "tinyCe_tinySc.fa"
    assert not (tmp_path / "data" / "genome" / "tinyCe_tinySc").exists()


# --------------------------------------------------------------------------------------
# Bytes copied verbatim — the property with no other failing case
# --------------------------------------------------------------------------------------


def test_every_component_sequence_line_survives_byte_for_byte(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    built = _fasta_records(chimera.fasta_path)

    for name in CHIMERA_EVERYDAY:
        component = CHIMERA_COMPONENTS[name]
        source = _fasta_records(component.fasta)
        for chromosome, lines in source.items():
            assert built[f"{chromosome}__{name}"] == lines


def test_components_that_disagree_about_wrapping_are_not_rewrapped(
    build_chimera: ChimeraFactory,
) -> None:
    # 60 for ce11's shape and 80 for ecHT115's, side by side in one file — which is what
    # every tool that could have done this concatenation would have destroyed.
    built = _fasta_records(build_chimera(*CHIMERA_EVERYDAY).fasta_path)

    assert _wrap_widths(built["I__tinyCe"]) == {60}
    assert _wrap_widths(built["NZ_TINY01000001.1__tinyEc"]) == {80}


def test_components_that_disagree_about_masking_keep_their_own_case(
    build_chimera: ChimeraFactory,
) -> None:
    # tinyCe carries 200 soft-masked bases and tinyEc carries none, so a build that
    # upper-cased or lower-cased anything loses one of the two.
    built = _fasta_records(build_chimera(*CHIMERA_EVERYDAY).fasta_path)
    worm = "".join(built["I__tinyCe"])

    assert worm[:200].islower()
    assert worm[200:].isupper()
    assert "".join(built["NZ_TINY01000002.1__tinyEc"]).isupper()


@pytest.mark.parametrize(
    ("line", "expected"),
    [
        (b">I\n", b">I__tinyCe\n"),
        (b">I some description\n", b">I__tinyCe some description\n"),
        (b">I\tdesc here\n", b">I__tinyCe\tdesc here\n"),
        (b">I", b">I__tinyCe"),  # a final header with no line ending
        (b">I\r\n", b">I__tinyCe\r\n"),
        # Whitespace after '>' is skipped by samtools faidx and faToTwoBit, which name
        # the sequence from the token after it — so `> desc` declares `desc`, and that is
        # what carries the suffix. The skipped bytes are written back where they were.
        (b"> desc\n", b"> desc__tinyCe\n"),
        (b">\t chrA desc\n", b">\t chrA__tinyCe desc\n"),
        (b">   I\n", b">   I__tinyCe\n"),
    ],
)
def test_only_the_name_a_header_gives_its_sequence_is_extended(
    line: bytes, expected: bytes
) -> None:
    # The suffix rides on the sequence's name because STAR and chromap both truncate a
    # header at the first whitespace; everything else on the line is description or
    # layout, and is written back untouched.
    assert _extend_header(line, "tinyCe", "__") == expected


@pytest.mark.parametrize("line", [b">\n", b">", b">   \n", b">\t\n", b">\r\n"])
def test_a_header_that_names_no_sequence_is_refused(line: bytes) -> None:
    # There is nothing for the suffix to ride on, and writing '>__tinyCe' would file the
    # sequence under a name that is suffix and nothing else — which samtools reads as the
    # empty name and the naming contract refuses outright.
    with pytest.raises(ChimeraNamingError, match="names no sequence"):
        _extend_header(line, "tinyCe", "__")


def test_a_component_header_starting_with_whitespace_is_named_as_the_tools_name_it(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # End to end, because this is where reading the header differently from samtools shows
    # up: the component registers under `oddChr`, so the chimera must carry
    # `oddChr__tinyOdd` — and a build that suffixed the empty token instead writes
    # `__tinyOdd`, which the cross-check against the components' prediction refuses.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    fasta = tmp_path / "odd.fa"
    fasta.write_text(">  oddChr a description after two spaces\n" + "ACGT" * 15 + "\n")
    odd = Genome("tinyOdd", path_or_url=fasta, cache_dir=tmp_path / "tinyOdd", progressbar=False)

    assert odd.chromosomes == ["oddChr"]  # what samtools and faToTwoBit called it

    with odd, Genome.chimera(odd, chimera_component("tinySc")) as chimera:
        assert "oddChr__tinyOdd" in chimera.chromosomes
        assert chimera["oddChr__tinyOdd"] == odd["oddChr"]


# --------------------------------------------------------------------------------------
# What the reference then answers
# --------------------------------------------------------------------------------------


def test_sequence_retrieval_agrees_with_the_component_it_came_from(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
    # The whole 2bit path over concatenated bytes: same bases, same case, same offsets.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    worm = chimera_component("tinyCe")
    draft = chimera_component("tinyEc")

    assert chimera["I__tinyCe:0-100"] == worm["I:0-100"]
    assert chimera["I__tinyCe:0-100"].islower()  # soft-masking survived faToTwoBit
    assert chimera["chr1_KI270706v1_random__tinyEc:0-100"] == draft["chr1_KI270706v1_random:0-100"]
    assert chimera["MtDNA__tinyCe"] == worm["MtDNA"]


def test_a_bare_chromosome_name_does_not_resolve_and_the_refusal_names_the_spelling(
    build_chimera: ChimeraFactory,
) -> None:
    # Resolving `III` too would restore the ambiguity the suffix abolishes, and would make
    # what a name means depend on which components happen to be present (ADR-0009). What
    # makes that bearable is naming the spelling that does resolve — `III` is the ninth of
    # the nine sequences here, so a message merely listing the first few never mentions it.
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    with pytest.raises(ValueError, match="III__tinySc"):
        chimera["III:0-100"]

    # ...and what it named is a name that resolves, which is what makes it a next action.
    assert len(chimera["III__tinySc:0-100"]) == 100


def test_a_bare_name_two_components_carry_is_answered_with_both(
    build_chimera: ChimeraFactory,
) -> None:
    # tinyCe and tinySc both carry `I`. Naming one spelling would be picking for the
    # caller the very thing the suffix exists to make them say, so both are named.
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    with pytest.raises(ValueError, match="I__tinySc") as refusal:
        chimera["I:0-100"]

    assert "I__tinyCe" in str(refusal.value)


def test_an_unknown_name_no_component_carries_is_refused_as_it_always_was(
    build_chimera: ChimeraFactory,
) -> None:
    # A chimera has one more thing to say only about a bare name it actually carries; a
    # name nothing carries gets the general message, on a chimera as anywhere else.
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    with pytest.raises(ValueError, match="known sequences include") as refusal:
        chimera["chrZ:0-5"]

    assert "ADR-0009" not in str(refusal.value)


def test_a_plain_assembly_says_nothing_about_components(
    chimera_component: ComponentFactory,
) -> None:
    # `III` is a real tinySc chromosome and unknown to tinyCe, and an assembly that is not
    # a chimera has no suffixed spelling to offer for it — so nothing here changes.
    worm = chimera_component("tinyCe")

    with pytest.raises(ValueError, match="known sequences include") as refusal:
        worm["III:0-100"]

    assert "ADR-0009" not in str(refusal.value)


def test_two_components_that_collide_get_four_distinct_names(
    build_chimera: ChimeraFactory,
) -> None:
    # tinyCe and tinySc both carry `I` and `II`; the shipped pair never does.
    chimera = build_chimera("tinyCe", "tinySc")

    assert chimera.assembly == "tinyCe_tinySc"
    assert {"I__tinyCe", "II__tinyCe", "I__tinySc", "II__tinySc"} <= set(chimera.chromosomes)
    assert len(set(chimera.chromosomes)) == len(chimera.chromosomes) == 6
    assert chimera["I__tinyCe"] != chimera["I__tinySc"]


def test_a_component_carrying_a_doubled_underscore_escalates_the_separator(
    build_chimera: ChimeraFactory,
) -> None:
    # The separator belongs to one chimera, so it is derived from these components and
    # recorded — a build that wrote a constant would lose the self-announcing property
    # exactly here, and nothing downstream could tell.
    chimera = build_chimera("tinyCe", CHIMERA_ESCALATION)
    details = read_chimera_details(chimera.fasta_path.parent)

    assert chimera.assembly == "tinyCe_tinyEcDub"
    assert details is not None
    assert details.separator == "___"
    assert "NZ_TINY02__000002.1___tinyEcDub" in chimera.chromosomes
    assert chimera.chromosomes[0] == "I___tinyCe"
    # And the recorded separator is the one that reads the names back.
    assert [split_suffixed(name, details.separator) for name in chimera.chromosomes] == [
        ("I", "tinyCe"),
        ("II", "tinyCe"),
        ("MtDNA", "tinyCe"),
        ("NZ_TINY02000001.1", "tinyEcDub"),
        ("NZ_TINY02__000002.1", "tinyEcDub"),
    ]


# --------------------------------------------------------------------------------------
# components and chrom_components
# --------------------------------------------------------------------------------------


def test_components_is_the_sorted_component_names(build_chimera: ChimeraFactory) -> None:
    assert build_chimera(*reversed(CHIMERA_EVERYDAY)).components == list(CHIMERA_EVERYDAY)


def test_a_plain_assembly_is_not_a_chimera_of_nothing(
    chimera_component: ComponentFactory,
) -> None:
    worm = chimera_component("tinyCe")

    assert worm.components is None


def test_chrom_components_attributes_every_chromosome_of_a_chimera(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    attribution = chimera.chrom_components

    assert list(attribution.index) == chimera.chromosomes
    assert not attribution.isna().to_numpy().any()
    assert attribution["I__tinyCe"] == "tinyCe"
    assert attribution["I__tinySc"] == "tinySc"
    assert attribution["chr1_KI270706v1_random__tinyEc"] == "tinyEc"
    assert Counter(attribution) == {"tinyCe": 3, "tinyEc": 3, "tinySc": 3}


def test_chrom_components_is_total_for_a_plain_assembly_too(
    chimera_component: ComponentFactory,
) -> None:
    # Totality is what keeps `components` the single is-chimera test: no caller has to
    # read this one to find out which kind of assembly it is holding.
    worm = chimera_component("tinyCe")
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


def test_rebuilding_a_finished_chimera_rewrites_nothing(build_chimera: ChimeraFactory) -> None:
    first = build_chimera(*CHIMERA_EVERYDAY)
    written = first.fasta_path.stat().st_mtime_ns
    record = read_record(first.fasta_path.parent)

    again = build_chimera(*CHIMERA_EVERYDAY)

    assert again.fasta_path == first.fasta_path
    assert again.fasta_path.stat().st_mtime_ns == written
    assert read_record(again.fasta_path.parent) == record


def test_a_chimera_directory_that_cannot_be_trusted_raises_naming_the_repair(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera("tinyCe", "tinySc")
    record_path(chimera.fasta_path.parent).unlink()

    with pytest.raises(UnfinishedRegistrationError, match="genome register tinyCe_tinySc --force"):
        build_chimera("tinyCe", "tinySc")
    # ...and that command is what repairs it, which is what `force` is.
    assert build_chimera("tinyCe", "tinySc", force=True).chromosomes == chimera.chromosomes


def test_a_chrom_sizes_that_disagrees_with_the_components_is_never_recorded(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # The cross-check between the tools' answer and the contract's prediction. Standing a
    # wrong chrom.sizes in is the only way to reach it: the build is a verbatim copy under
    # a derived name, so it is meant to be unreachable.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    prepare = chimera_mod.prepare_fasta

    def prepare_and_mangle(fasta: Path, *, overwrite: bool = False) -> GenomeFiles:
        files = prepare(fasta, overwrite=overwrite)
        lines = files.chrom_sizes.read_text().splitlines()
        name, length = lines[0].split("\t")
        lines[0] = f"{name}\t{int(length) - 1}"  # one base short, as a short read would be
        files.chrom_sizes.write_text("\n".join(lines) + "\n")
        return files

    monkeypatch.setattr(chimera_mod, "prepare_fasta", prepare_and_mangle)
    components = [chimera_component(name) for name in ("tinyCe", "tinySc")]

    with pytest.raises(RegistrationError, match="genome register tinyCe_tinySc --force"):
        Genome.chimera(*components)
    assert read_record(assembly_data_dir("tinyCe_tinySc")) is None


def test_building_a_chimera_fetches_nothing(
    build_chimera: ChimeraFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The autouse guard would raise on a real network call; this says deliberately that
    # the package's one fetch step is not reached at all, which is a stronger claim.
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("a chimera build fetches nothing")

    monkeypatch.setattr(fetch_mod, "fetch_url", refuse)

    assert build_chimera(*CHIMERA_EVERYDAY).components == list(CHIMERA_EVERYDAY)


# --------------------------------------------------------------------------------------
# What a chimera refuses to be built from
# --------------------------------------------------------------------------------------


def test_one_component_is_not_a_chimera(chimera_component: ComponentFactory) -> None:
    with pytest.raises(ChimeraNamingError, match="at least 2 components"):
        Genome.chimera(chimera_component("tinyCe"))


def test_a_repeated_component_raises(chimera_component: ComponentFactory) -> None:
    worm = chimera_component("tinyCe")

    with pytest.raises(ChimeraNamingError, match="must not repeat"):
        Genome.chimera(worm, chimera_component("tinySc"), worm)


def test_a_chimera_cannot_be_a_component_of_another(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
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
