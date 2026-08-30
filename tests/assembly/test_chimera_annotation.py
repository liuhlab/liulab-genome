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
from genome.annotation import (
    AnnotationNotRegisteredError,
    ChromosomeMismatchError,
    annotation_dir,
)
from genome.annotation import registration as registration_mod
from genome.assembly.chimera import split_suffixed
from genome.assembly.chimera_build import AmbiguousDefaultAnnotationError, read_chimera_details
from genome.store import fetch as fetch_mod
from genome.store.checksum import sha256_file
from genome.store.completion import read_record

from ..conftest import (
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
def build_chimera(chimera_component: ComponentFactory) -> Iterator[ChimeraFactory]:
    """Return a factory building a chimera of the named tiny components, in a temp root.

    Every component that ships a GTF is registered with it, so the everyday call is the
    annotated one; ``annotate=False`` builds the same components with nothing registered.
    The chimera lands where the layout puts it, under the root the shared ``liulab_data``
    fixture pointed at this test's own directory — a chimera built into the lab's shared
    reference data by a test would be a serious bug.
    """
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


def test_the_everyday_chimera_arrives_with_its_merged_annotation_registered_and_queryable(
    build_chimera: ChimeraFactory,
) -> None:
    chimera = build_chimera(*CHIMERA_EVERYDAY)

    # No second surface: nobody registered anything, and the annotation is simply there.
    directory = annotation_dir(chimera.fasta_path.parent, _EVERYDAY_MERGED)
    assert chimera.annotations.registered == [_EVERYDAY_MERGED]
    assert chimera.default_gtf == _EVERYDAY_MERGED
    assert chimera.annotations.path(_EVERYDAY_MERGED) == directory / f"{_EVERYDAY_MERGED}.gtf"
    assert (directory / f"{_EVERYDAY_MERGED}.gtf").is_file()
    assert (directory / f"{_EVERYDAY_MERGED}.db").is_file()

    # And its gene ids are queryable through the built database: the fixtures' slices are
    # disjoint set-wide, so no gene id is carried by two components and every one of
    # these can only have come from where it says.
    assert {"YBL111C", "YCL074W", "YAL069W"} <= _gene_ids(directory / f"{_EVERYDAY_MERGED}.db")
    assert _seqids(directory / f"{_EVERYDAY_MERGED}.db") == set(chimera.chromosomes)


def test_the_merged_name_is_the_plus_join_in_sorted_component_order(
    chimera_component: ComponentFactory,
) -> None:
    # Distinct annotation names per component, so the order is visible: sorted by
    # component (tinyCe, tinyEc, tinySc) and not by annotation name, which would give
    # 'ecoli+worm+yeast'.
    named = {"tinyCe": "worm", "tinyEc": "ecoli", "tinySc": "yeast"}
    components = []
    for component, annotation in named.items():
        genome = chimera_component(component)
        gtf = CHIMERA_COMPONENTS[component].gtf
        assert gtf is not None
        genome.annotations.register_path(gtf, annotation)
        components.append(genome)

    with Genome.chimera(*reversed(components)) as chimera:
        assert chimera.annotations.registered == ["worm+ecoli+yeast"]


def test_the_merged_gtf_carries_exactly_the_chimeras_names_kept_in_component_order(
    build_chimera: ChimeraFactory,
) -> None:
    # The fixtures' GTF seqnames are set-equal to their assembly's, so a merge that
    # dropped or misspelled one shows up as an inequality here — and a misspelling would
    # not even reach this assertion, check_chromosomes having refused it first.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    merged = chimera.annotations.path(_EVERYDAY_MERGED)
    seqnames = _seqnames(merged)

    assert set(seqnames) == set(chimera.chromosomes)
    assert sorted({split_suffixed(name)[1] for name in seqnames}) == list(CHIMERA_EVERYDAY)

    # Unsorted: each component's lines arrive in that component's own order, and the
    # components in the order the chimera's name spells them, with none interleaved.
    components = [split_suffixed(name)[1] for name in seqnames]
    first_seen = list(dict.fromkeys(components))
    assert first_seen == list(CHIMERA_EVERYDAY)
    assert components == sorted(components, key=first_seen.index)


def test_every_data_line_survives_and_no_pragma_does(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # A #!genome-build pragma names the one assembly its file was built for, so several
    # of them concatenated would each be false about the chimera.
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
        genome.annotations.register_path(commented, COMPONENT_ANNOTATION)
        expected += len(_data_lines(source))
        components.append(genome)

    with Genome.chimera(*components) as chimera:
        merged = chimera.annotations.path(_EVERYDAY_MERGED)
        lines = merged.read_text().splitlines()

    assert len(lines) == expected
    assert not [line for line in lines if line.startswith("#")]


# --------------------------------------------------------------------------------------
# What each component contributes — and what it means to contribute nothing
# --------------------------------------------------------------------------------------


def test_a_partial_or_total_absence_of_contributors_is_never_treated_as_an_empty_merge(
    build_chimera: ChimeraFactory,
) -> None:
    # Nobody contributed first, on the everyday three — before either one is touched by
    # the partial case below, since `chimera_component` only ever *adds* a registration,
    # never removes one, so a component annotated first would stay annotated second.
    empty = build_chimera(*CHIMERA_EVERYDAY, annotate=False)
    # `.registered`, and never the registry object: it is always truthy, and reading it
    # as a collection is exactly what would have made this build raise rather than
    # register nothing — see `ChimeraBuilder._contribution`.
    assert empty.annotations.registered == []
    assert empty.annotations.broken == []
    assert empty.default_gtf is None
    # Not an empty annotation, and not an empty directory pretending to be one.
    assert not (empty.fasta_path.parent / "gtf").exists()
    assert empty.chromosomes  # ...and the chimera itself opened perfectly well

    # tinyEcDub ships no GTF at all, which is the shape the one-directional chromosome
    # check already blesses: an assembly may carry sequences an annotation never mentions.
    partial = build_chimera("tinyCe", CHIMERA_ESCALATION)
    assert partial.annotations.registered == [COMPONENT_ANNOTATION]
    # Contributed against contributed-nothing, per component: the distinction that decides
    # which annotation a per-component count is taken against, and the one the merged name
    # cannot carry, since it names only the contributors.
    assert partial.component_annotations == {
        "tinyCe": COMPONENT_ANNOTATION,
        CHIMERA_ESCALATION: None,
    }
    merged = partial.annotations.path(COMPONENT_ANNOTATION)
    assert {split_suffixed(name, "___")[1] for name in _seqnames(merged)} == {"tinyCe"}
    assert set(_seqnames(merged)) < set(partial.chromosomes)


def test_an_undecided_default_annotation_raises_whether_unregistered_or_ambiguous(
    chimera_component: ComponentFactory, liulab_data: Path
) -> None:
    # The cold-machine rule: a default is a name until something registers it, and a
    # chimera build must not turn that into a silently annotation-less reference.
    prepared = chimera_component("tinyCe")
    yeast = chimera_component("tinySc", with_annotation=True)
    with (
        Genome(
            "tinyCe",
            cache_dir=prepared.fasta_path.parent,
            progressbar=False,
            default_gtf=COMPONENT_ANNOTATION,
        ) as named_only,
        pytest.raises(AnnotationNotRegisteredError, match="tinyCe") as unregistered,
    ):
        Genome.chimera(named_only, yeast)
    assert "register_path" in str(unregistered.value)
    # Nothing was built: the refusal comes before a byte is written.
    assert not (liulab_data / "genome" / "tinyCe_tinySc").exists()

    # The state with no default at all is a different failure: reopened, because a
    # component that registered its second annotation in this very session keeps the
    # first as its default; a component opened over two of them, with no table flag to
    # choose between, has none.
    worm = chimera_component("tinyCe", with_annotation=True)
    gtf = CHIMERA_COMPONENTS["tinyCe"].gtf
    assert gtf is not None
    worm.annotations.register_path(gtf, "genes_again")
    with Genome("tinyCe", cache_dir=worm.fasta_path.parent, progressbar=False) as undecided:
        assert undecided.default_gtf is None
        with pytest.raises(AmbiguousDefaultAnnotationError, match="default_gtf=") as ambiguous:
            Genome.chimera(undecided, yeast)
    assert "tinyCe" in str(ambiguous.value)
    assert "genes_again" in str(ambiguous.value)
    assert not (liulab_data / "genome" / "tinyCe_tinySc").exists()


# --------------------------------------------------------------------------------------
# The two records, and repairing what they vouch for
# --------------------------------------------------------------------------------------


def test_the_two_records_carry_each_components_digest_and_the_merged_annotations_provenance(
    build_chimera: ChimeraFactory, chimera_component: ComponentFactory
) -> None:
    partial = build_chimera("tinyCe", CHIMERA_ESCALATION)
    details = read_chimera_details(partial.fasta_path.parent)
    worm = chimera_component("tinyCe", with_annotation=True)
    pinned = read_record(worm.annotations.path(COMPONENT_ANNOTATION).parent)

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

    # And the merged annotation's own record: nothing was downloaded, so there is nothing
    # for a table row or a source URL to describe, and the record says so rather than
    # pointing at a file it merged from.
    everyday = build_chimera(*CHIMERA_EVERYDAY)
    record = read_record(annotation_dir(everyday.fasta_path.parent, _EVERYDAY_MERGED))
    assert record is not None
    assert (record.kind, record.name) == ("annotation", _EVERYDAY_MERGED)
    assert record.source_url is None
    assert record.sha256 == sha256_file(everyday.annotations.path(_EVERYDAY_MERGED))
    assert record.details["merged_from"] == [
        {"component": name, "annotation": COMPONENT_ANNOTATION} for name in CHIMERA_EVERYDAY
    ]
    # check_chromosomes is left on, and this is the record of it having run.
    assert record.details["chromosomes_checked"] is True


def test_rebuilding_a_finished_chimera_rebuilds_nothing_but_force_repairs_both_halves(
    build_chimera: ChimeraFactory,
) -> None:
    first = build_chimera(*CHIMERA_EVERYDAY)
    database = annotation_dir(first.fasta_path.parent, _EVERYDAY_MERGED) / f"{_EVERYDAY_MERGED}.db"
    built = database.stat().st_mtime_ns
    written = first.fasta_path.stat().st_mtime_ns

    again = build_chimera(*CHIMERA_EVERYDAY)
    assert again.annotations.registered == [_EVERYDAY_MERGED]
    assert database.stat().st_mtime_ns == built

    database.unlink()
    # One command repairs both halves, because one act built them.
    repaired = build_chimera(*CHIMERA_EVERYDAY, force=True)
    assert repaired.annotations.registered == [_EVERYDAY_MERGED]
    assert database.is_file()
    assert repaired.fasta_path.stat().st_mtime_ns != written


def test_a_rebuild_never_keeps_a_stale_annotation_whether_renamed_or_gone_entirely(
    chimera_component: ComponentFactory,
) -> None:
    # The merged name is the +-join of what contributed, so a contributing set that
    # changes across a rebuild changes it — and the build owns both, so the one it no
    # longer owns goes. Left beside the new one it would be a second annotation with no
    # default between them, and a chimera that arrived annotated would come back from a
    # legitimate repair with `default_gtf` of None.
    worm = chimera_component("tinyCe", with_annotation=True)
    yeast = chimera_component("tinySc", with_annotation=True)
    with Genome.chimera(worm, yeast) as first:
        assembly_dir = first.fasta_path.parent
        assert first.default_gtf == "genes+genes"

    gtf = CHIMERA_COMPONENTS["tinyCe"].gtf
    assert gtf is not None
    worm.annotations.register_path(gtf, "genes_again")
    with (
        Genome(
            "tinyCe",
            cache_dir=worm.fasta_path.parent,
            progressbar=False,
            default_gtf="genes_again",
        ) as renamed,
        Genome.chimera(renamed, yeast, force=True) as rebuilt,
    ):
        assert rebuilt.annotations.registered == ["genes_again+genes"]
        assert rebuilt.default_gtf == "genes_again+genes"
    # Removed rather than left broken: nothing vouches for it and nothing points at it.
    assert not annotation_dir(assembly_dir, "genes+genes").exists()

    # The same fault at its worst: with nothing left to merge, the annotation the
    # previous build wrote would otherwise be the *only* one registered and therefore the
    # default, answering queries with gene models this chimera no longer merges from
    # anything — so a rebuild with no contributors leaves no annotation at all. A fresh
    # pair, independent of the one above: `renamed` was closed when its `with` exited.
    annotated = chimera_component("tinySc", with_annotation=True)
    plain = chimera_component("tinyEc")
    with Genome.chimera(annotated, plain) as second:
        second_dir = second.fasta_path.parent
        assert second.default_gtf == COMPONENT_ANNOTATION

    shutil.rmtree(annotation_dir(annotated.fasta_path.parent, COMPONENT_ANNOTATION))
    with (
        Genome("tinySc", cache_dir=annotated.fasta_path.parent, progressbar=False) as bare,
        Genome.chimera(bare, plain, force=True) as bare_rebuilt,
    ):
        assert bare_rebuilt.annotations.registered == []
        assert bare_rebuilt.default_gtf is None
    assert not (second_dir / "gtf").exists()


def test_a_broken_merged_annotation_names_the_command_that_rebuilds_the_chimera(
    build_chimera: ChimeraFactory,
) -> None:
    # Neither registering it by name nor handing it a GTF would rebuild a derived
    # annotation, so neither may be the command a broken one prints.
    chimera = build_chimera(*CHIMERA_EVERYDAY)
    assembly_dir = chimera.fasta_path.parent
    (annotation_dir(assembly_dir, _EVERYDAY_MERGED) / f"{_EVERYDAY_MERGED}.db").unlink()

    with Genome(chimera.assembly, cache_dir=assembly_dir, progressbar=False) as reopened:
        broken = reopened.annotations.broken

    assert [entry.name for entry in broken] == [_EVERYDAY_MERGED]
    assert broken[0].repair == f"genome register {chimera.assembly} --force"


# --------------------------------------------------------------------------------------
# What the merge is checked against, and what it never touches
# --------------------------------------------------------------------------------------


def test_a_misspelled_merge_is_refused_before_registration_and_nothing_is_ever_fetched(
    chimera_component: ComponentFactory,
    build_chimera: ChimeraFactory,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # check_chromosomes is the one place the merge's answer and the FASTA build's answer
    # are set against each other, so it is left on and this is what it buys. Scoped to its
    # own MonkeyPatch context rather than the fixture-provided one: the autouse
    # `liulab_data` fixture patches the data-root env var through that same shared
    # instance, and undoing it here would undo that redirect too, sending the rest of
    # this test at the real data root.
    components = [chimera_component(name, with_annotation=True) for name in ("tinyCe", "tinySc")]
    with pytest.MonkeyPatch.context() as scoped:
        scoped.setattr(
            registration_mod,
            "suffixed",
            lambda chromosome, component, separator: f"{chromosome}__wrong",
        )
        with pytest.raises(ChromosomeMismatchError, match="wrong"):
            Genome.chimera(*components)

    # And the autouse network guard would raise on a real call; this says deliberately
    # that the package's one fetch step is not reached at all, annotation included.
    def refuse(*args: object, **kwargs: object) -> None:
        raise AssertionError("a chimera build fetches nothing, annotation included")

    monkeypatch.setattr(fetch_mod, "fetch_url", refuse)
    assert build_chimera(*CHIMERA_EVERYDAY).annotations.registered == [_EVERYDAY_MERGED]
