"""Tests for the merged annotation a chimera build writes and registers.

The FASTA half has its own tests (``test_chimera_build``); what is asserted here is the
annotation half: which annotation each component contributes, what the merge is called,
what its GTF carries, and what the two records say afterwards.

Every build runs the real ``samtools``/``faToTwoBit``/``twoBitInfo`` and a real
``gffutils`` database build over the tiny components, offline, in a temporary data root —
never the lab's own.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable, Iterator
from pathlib import Path

import gffutils
import pytest

from genome import Genome
from genome.chimera import split_suffixed
from genome.io import fetch as fetch_mod
from genome.io import gtf as gtf_mod
from genome.io.chimera import AmbiguousDefaultAnnotationError, read_chimera_details
from genome.io.completion import read_record
from genome.io.gtf import AnnotationNotRegisteredError, ChromosomeMismatchError, annotation_dir
from genome.io.utils import sha256_file

from .conftest import (
    CHIMERA_COMPONENTS,
    CHIMERA_ESCALATION,
    CHIMERA_EVERYDAY,
    COMPONENT_ANNOTATION,
    ComponentFactory,
)

#: What the everyday three, each carrying an annotation named ``genes``, merge to.
_EVERYDAY_MERGED = "genes+genes+genes"

ChimeraFactory = Callable[..., Genome]


def _data_lines(gtf: Path) -> list[str]:
    """Every line of ``gtf`` that carries a tab — the data lines, comments left out."""
    return [line for line in gtf.read_text().splitlines() if "\t" in line]


def _seqnames(gtf: Path) -> list[str]:
    """The first column of every data line, in file order."""
    return [line.split("\t")[0] for line in _data_lines(gtf)]


def _gene_ids(db: Path) -> set[str]:
    """Every ``gene_id`` the built database can be asked for."""
    database = gffutils.FeatureDB(str(db))
    try:
        return {
            value
            for feature in database.all_features()
            for value in feature.attributes.get("gene_id", [])
        }
    finally:
        database.conn.close()


def _seqids(db: Path) -> set[str]:
    """Every sequence name the built database's features sit on."""
    database = gffutils.FeatureDB(str(db))
    try:
        return {feature.seqid for feature in database.all_features()}
    finally:
        database.conn.close()


@pytest.fixture
def build_chimera(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Iterator[ChimeraFactory]:
    """Return a factory building a chimera of the named tiny components, in a temp root.

    Every component that ships a GTF is registered with it, so the everyday call is the
    annotated one; ``annotate=False`` builds the same components with nothing registered.
    ``LIULAB_DATA`` is pointed at the test's own directory first — a chimera built into
    the lab's shared reference data by a test would be a serious bug.
    """
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    opened: list[Genome] = []

    def build(*names: str, annotate: bool = True, force: bool = False) -> Genome:
        components = [
            chimera_component(name, with_annotation=annotate and CHIMERA_COMPONENTS[name].has_gtf)
            for name in names
        ]
        chimera = Genome.chimera(*components, force=force)
        opened.append(chimera)
        return chimera

    yield build
    for chimera in opened:
        chimera.close()


# --------------------------------------------------------------------------------------
# The merge is part of the build, and it is named after what went into it
# --------------------------------------------------------------------------------------


def test_the_everyday_chimera_arrives_with_its_merged_annotation_registered(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    # No second surface: nobody called register_gtf, and the annotation is simply there.
    directory = annotation_dir(chimera.fasta_path.parent, _EVERYDAY_MERGED)
    assert chimera.annotations == [_EVERYDAY_MERGED]
    assert chimera.default_gtf == _EVERYDAY_MERGED
    assert chimera.get_gtf_path(_EVERYDAY_MERGED) == directory / f"{_EVERYDAY_MERGED}.gtf"
    assert (directory / f"{_EVERYDAY_MERGED}.gtf").is_file()
    assert (directory / f"{_EVERYDAY_MERGED}.db").is_file()


def test_the_merged_name_is_the_plus_join_in_sorted_component_order(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Distinct annotation names per component, so the order is visible: sorted by
    # component (tinyCe, tinyEc, tinySc) and not by annotation name, which would give
    # 'ecoli+worm+yeast'.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    named = {"tinyCe": "worm", "tinyEc": "ecoli", "tinySc": "yeast"}
    components = []
    for component, annotation in named.items():
        genome = chimera_component(component)
        gtf = CHIMERA_COMPONENTS[component].gtf
        assert gtf is not None
        genome.register_gtf(gtf, annotation)
        components.append(genome)

    with Genome.chimera(*reversed(components)) as chimera:
        assert chimera.annotations == ["worm+ecoli+yeast"]


def test_the_merged_gtf_carries_exactly_the_chimeras_own_chromosome_names(
    build_chimera: ChimeraFactory,
) -> None:
    # The fixtures' GTF seqnames are set-equal to their assembly's, so a merge that
    # dropped or misspelled one shows up as an inequality here — and a misspelling would
    # not even reach this assertion, check_chromosomes having refused it first.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    merged = chimera.get_gtf_path(_EVERYDAY_MERGED)

    assert set(_seqnames(merged)) == set(chimera.chromosomes)
    assert sorted({split_suffixed(name)[1] for name in _seqnames(merged)}) == list(CHIMERA_EVERYDAY)


def test_the_merged_gtf_is_unsorted_and_in_component_order(
    build_chimera: ChimeraFactory,
) -> None:
    # Unsorted: each component's lines arrive in that component's own order, and the
    # components in the order the chimera's name spells them.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    components = [
        split_suffixed(name)[1] for name in _seqnames(chimera.get_gtf_path(_EVERYDAY_MERGED))
    ]

    first_seen = list(dict.fromkeys(components))
    assert first_seen == list(CHIMERA_EVERYDAY)
    # ...and no component's lines are interleaved with another's.
    assert components == sorted(components, key=first_seen.index)


def test_every_data_line_survives_and_no_pragma_does(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # A #!genome-build pragma names the one assembly its file was built for, so several
    # of them concatenated would each be false about the chimera.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    expected = 0
    components = []
    for name in CHIMERA_EVERYDAY:
        source = CHIMERA_COMPONENTS[name].gtf
        assert source is not None
        body = source.read_text()
        commented = tmp_path / f"{name}-commented.gtf"
        commented.write_text(
            f"#!genome-build {name}\n##provider: fixture\n# a plain comment\n{body}"
        )
        genome = chimera_component(name)
        genome.register_gtf(commented, COMPONENT_ANNOTATION)
        expected += len(_data_lines(source))
        components.append(genome)

    with Genome.chimera(*components) as chimera:
        merged = chimera.get_gtf_path(_EVERYDAY_MERGED)
        lines = merged.read_text().splitlines()

    assert len(lines) == expected
    assert not [line for line in lines if line.startswith("#")]


def test_gene_ids_from_every_component_are_queryable_through_the_database(
    build_chimera: ChimeraFactory,
) -> None:
    # The fixtures' slices are disjoint set-wide, so no gene id is carried by two
    # components and every one of these can only have come from where it says.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    database = (
        annotation_dir(chimera.fasta_path.parent, _EVERYDAY_MERGED) / f"{_EVERYDAY_MERGED}.db"
    )

    assert {"YBL111C", "YCL074W", "YAL069W"} <= _gene_ids(database)
    assert _seqids(database) == set(chimera.chromosomes)


# --------------------------------------------------------------------------------------
# What each component contributes — and what it means to contribute nothing
# --------------------------------------------------------------------------------------


def test_a_component_with_no_annotation_contributes_nothing(
    build_chimera: ChimeraFactory,
) -> None:
    # tinyEcDub ships no GTF at all, which is the shape the one-directional chromosome
    # check already blesses: an assembly may carry sequences an annotation never mentions.
    chimera = build_chimera("tinyCe", CHIMERA_ESCALATION)

    assert chimera.annotations == [COMPONENT_ANNOTATION]
    merged = chimera.get_gtf_path(COMPONENT_ANNOTATION)
    assert {split_suffixed(name, "___")[1] for name in _seqnames(merged)} == {"tinyCe"}
    assert set(_seqnames(merged)) < set(chimera.chromosomes)


def test_no_contributors_at_all_registers_no_annotation(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY, annotate=False)

    assert chimera.annotations == []
    assert chimera.broken_annotations == []
    assert chimera.default_gtf is None
    # Not an empty annotation, and not an empty directory pretending to be one.
    assert not (chimera.fasta_path.parent / "gtf").exists()
    assert chimera.chromosomes  # ...and the chimera itself opened perfectly well


def test_a_component_whose_default_annotation_is_not_registered_raises(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # The cold-machine rule, one level down: a default is a name until something registers
    # it, and a chimera build must not turn that into a silently annotation-less reference.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    prepared = chimera_component("tinyCe")
    yeast = chimera_component("tinySc", with_annotation=True)
    with (
        Genome(
            "tinyCe",
            cache_dir=prepared.fasta_path.parent,
            progressbar=False,
            default_gtf=COMPONENT_ANNOTATION,
        ) as named_only,
        pytest.raises(AnnotationNotRegisteredError, match="tinyCe") as raised,
    ):
        Genome.chimera(named_only, yeast)

    assert "register_gtf" in str(raised.value)
    # Nothing was built: the refusal comes before a byte is written.
    assert not (tmp_path / "data" / "genome" / "tinyCe_tinySc").exists()


def test_several_registered_and_none_flagged_raises_naming_default_gtf(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    worm = chimera_component("tinyCe", with_annotation=True)
    gtf = CHIMERA_COMPONENTS["tinyCe"].gtf
    assert gtf is not None
    worm.register_gtf(gtf, "genes_again")
    yeast = chimera_component("tinySc", with_annotation=True)
    # Reopened, because a component that registered its second annotation in this very
    # session keeps the first as its default; a component opened over two of them, with no
    # table flag to choose between, is the state that has no default at all.
    with Genome("tinyCe", cache_dir=worm.fasta_path.parent, progressbar=False) as undecided:
        assert undecided.default_gtf is None

        with pytest.raises(AmbiguousDefaultAnnotationError, match="default_gtf=") as raised:
            Genome.chimera(undecided, yeast)

    assert "tinyCe" in str(raised.value)
    assert "genes_again" in str(raised.value)
    assert not (tmp_path / "data" / "genome" / "tinyCe_tinySc").exists()


# --------------------------------------------------------------------------------------
# The two records, and repairing what they vouch for
# --------------------------------------------------------------------------------------


def test_the_chimeras_record_carries_each_components_annotation_and_digest(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
    chimera = build_chimera("tinyCe", CHIMERA_ESCALATION)
    details = read_chimera_details(chimera.fasta_path.parent)
    worm = chimera_component("tinyCe", with_annotation=True)
    pinned = read_record(worm.get_gtf_path(COMPONENT_ANNOTATION).parent)

    assert details is not None
    assert details.components == ["tinyCe", CHIMERA_ESCALATION]
    contributor, silent = details.component_details
    assert (contributor.annotation, contributor.annotation_sha256) == (
        COMPONENT_ANNOTATION,
        pinned.sha256 if pinned is not None else None,
    )
    assert contributor.annotation_sha256 is not None
    # Contributed nothing — which is not the same answer as "nobody asked".
    assert (silent.annotation, silent.annotation_sha256) == (None, None)


def test_the_merged_annotations_record_pins_no_source(build_chimera: ChimeraFactory) -> None:
    # Nothing was downloaded, so there is nothing for a table row or a source URL to
    # describe, and the record says so rather than pointing at a file it merged from.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    record = read_record(annotation_dir(chimera.fasta_path.parent, _EVERYDAY_MERGED))

    assert record is not None
    assert (record.kind, record.name) == ("annotation", _EVERYDAY_MERGED)
    assert record.source_url is None
    assert record.sha256 == sha256_file(chimera.get_gtf_path(_EVERYDAY_MERGED))
    assert record.details["merged_from"] == [
        {"component": name, "annotation": COMPONENT_ANNOTATION} for name in CHIMERA_EVERYDAY
    ]
    # check_chromosomes is left on, and this is the record of it having run.
    assert record.details["chromosomes_checked"] is True


def test_rebuilding_a_finished_chimera_rebuilds_no_annotation(
    build_chimera: ChimeraFactory,
) -> None:
    first = build_chimera(*CHIMERA_EVERYDAY)
    database = annotation_dir(first.fasta_path.parent, _EVERYDAY_MERGED) / f"{_EVERYDAY_MERGED}.db"
    built = database.stat().st_mtime_ns

    again = build_chimera(*CHIMERA_EVERYDAY)

    assert again.annotations == [_EVERYDAY_MERGED]
    assert database.stat().st_mtime_ns == built


def test_force_repairs_the_annotation_as_well_as_the_fasta(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    directory = annotation_dir(chimera.fasta_path.parent, _EVERYDAY_MERGED)
    (directory / f"{_EVERYDAY_MERGED}.db").unlink()
    written = chimera.fasta_path.stat().st_mtime_ns

    # One command repairs both halves, because one act built them.
    repaired = build_chimera(*CHIMERA_EVERYDAY, force=True)

    assert repaired.annotations == [_EVERYDAY_MERGED]
    assert (directory / f"{_EVERYDAY_MERGED}.db").is_file()
    assert repaired.fasta_path.stat().st_mtime_ns != written


def test_a_rebuild_whose_contributors_changed_leaves_only_the_annotation_it_wrote(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # The merged name is the +-join of what contributed, so a contributing set that
    # changes across a rebuild changes it — and the build owns both, so the one it no
    # longer owns goes. Left beside the new one it would be a second annotation with no
    # default between them, and a chimera that arrived annotated would come back from a
    # legitimate repair with `default_gtf` of None.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    worm = chimera_component("tinyCe", with_annotation=True)
    yeast = chimera_component("tinySc", with_annotation=True)
    with Genome.chimera(worm, yeast) as first:
        assembly_dir = first.fasta_path.parent
        assert first.default_gtf == "genes+genes"

    gtf = CHIMERA_COMPONENTS["tinyCe"].gtf
    assert gtf is not None
    worm.register_gtf(gtf, "genes_again")
    with (
        Genome(
            "tinyCe",
            cache_dir=worm.fasta_path.parent,
            progressbar=False,
            default_gtf="genes_again",
        ) as renamed,
        Genome.chimera(renamed, yeast, force=True) as rebuilt,
    ):
        assert rebuilt.annotations == ["genes_again+genes"]
        assert rebuilt.default_gtf == "genes_again+genes"

    # Removed rather than left broken: nothing vouches for it and nothing points at it.
    assert not annotation_dir(assembly_dir, "genes+genes").exists()


def test_a_rebuild_with_nothing_left_to_merge_leaves_no_annotation_at_all(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # The same fault at its worst: with nothing to merge, the annotation the previous
    # build wrote would be the *only* one registered here and therefore the default,
    # answering queries with gene models this chimera no longer merges from anything.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    worm = chimera_component("tinyCe", with_annotation=True)
    yeast = chimera_component("tinySc")
    with Genome.chimera(worm, yeast) as first:
        assembly_dir = first.fasta_path.parent
        assert first.default_gtf == COMPONENT_ANNOTATION

    shutil.rmtree(annotation_dir(worm.fasta_path.parent, COMPONENT_ANNOTATION))
    with (
        Genome("tinyCe", cache_dir=worm.fasta_path.parent, progressbar=False) as bare,
        Genome.chimera(bare, yeast, force=True) as rebuilt,
    ):
        assert rebuilt.annotations == []
        assert rebuilt.default_gtf is None

    assert not (assembly_dir / "gtf").exists()


def test_a_broken_merged_annotation_names_the_command_that_rebuilds_the_chimera(
    build_chimera: ChimeraFactory,
) -> None:
    # Neither registering it by name nor handing it a GTF would rebuild a derived
    # annotation, so neither may be the command a broken one prints.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    assembly_dir = chimera.fasta_path.parent
    (annotation_dir(assembly_dir, _EVERYDAY_MERGED) / f"{_EVERYDAY_MERGED}.db").unlink()

    with Genome(chimera.assembly, cache_dir=assembly_dir, progressbar=False) as reopened:
        broken = reopened.broken_annotations

    assert [entry.name for entry in broken] == [_EVERYDAY_MERGED]
    assert broken[0].repair == f"genome register {chimera.assembly} --force"


# --------------------------------------------------------------------------------------
# What the merge is checked against, and what it never touches
# --------------------------------------------------------------------------------------


def test_a_merge_that_misspells_a_name_is_refused_before_it_is_registered(
    chimera_component: ComponentFactory, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # check_chromosomes is the one place the merge's answer and the FASTA build's answer
    # are set against each other, so it is left on and this is what it buys.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "data"))
    components = [chimera_component(name, with_annotation=True) for name in ("tinyCe", "tinySc")]
    monkeypatch.setattr(
        gtf_mod, "suffixed", lambda chromosome, component, separator: f"{chromosome}__wrong"
    )

    with pytest.raises(ChromosomeMismatchError, match="wrong"):
        Genome.chimera(*components)


def test_building_a_chimera_with_annotations_fetches_nothing(
    build_chimera: ChimeraFactory, monkeypatch: pytest.MonkeyPatch
) -> None:
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("a chimera build fetches nothing, annotation included")

    monkeypatch.setattr(fetch_mod, "fetch_url", refuse)

    assert build_chimera(*CHIMERA_EVERYDAY).annotations == [_EVERYDAY_MERGED]
