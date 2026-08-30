"""Tests for genome.annotation.stems — a **Gene id stem** against an annotation's ids.

The seam the Xref, Orthology and TF contexts all cross, and the first surface here to open
the **Annotation database**, so every test registers a fixture GTF for real and asks the
registry. What is asserted is what a caller holding a published table keyed by unversioned
ids gets back: both of the ids a stem names rather than one of them, the stems this
annotation has nothing for said out loud, and the walk those answers ride on staying one
row at a time.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from genome.annotation.database import gene_ids
from genome.annotation.registration import MergeSource, register_merged_gtf
from genome.annotation.registry import AnnotationNotRegisteredError, AnnotationRegistry
from genome.annotation.stems import NoGeneFeaturesError

from .conftest import (
    _BARE_GTF,
    _CHIMERA,
    _CURATED,
    _CURATED_ASSEMBLY,
    _FOOD,
    _FOOD_COMPONENT,
    _GTF,
    _WORM,
    _WORM_COMPONENT,
    _register_by_path,
    _write_chrom_sizes,
)

# ---------------------------------------------------------------------------------------
# Which of an annotation's own gene ids a Gene id stem names
# ---------------------------------------------------------------------------------------


#: A pseudoautosomal gene as GENCODE spells it on a lift: one stem, two gene ids, the
#: second the Y copy. Real ids — PLCXD1 and TP53 — because the shape of the collision is
#: the whole point and an invented id would not have it.
_PAR_X = "ENSG00000182378.14"


_PAR_Y = "ENSG00000182378.14_PAR_Y"


_ALONE = "ENSG00000141510.18"


_PAR_STEM = "ENSG00000182378"


_ALONE_STEM = "ENSG00000141510"


#: A stem no fixture below carries a gene for — what the census holds and the annotation
#: does not.
_ABSENT_STEM = "ENSG00000288541"


def _gtf_declaring(*gene_ids: str) -> str:
    """A GTF declaring one gene, transcript and exon per id, in the shape GENCODE has.

    The ids are the only thing that varies: nothing reading it looks at a coordinate, so
    every gene sits on the same interval rather than inviting anyone to.
    """
    lines: list[str] = []
    for gene_id in gene_ids:
        attributes = f'gene_id "{gene_id}"; transcript_id "{gene_id}_t";'
        lines.extend(
            f"chrI\ttest\t{feature}\t1\t100\t.\t+\t.\t{attributes}"
            for feature in ("gene", "transcript", "exon")
        )
    return "\n".join(lines) + "\n"


class TestResolveGeneIds:
    """``AnnotationRegistry.resolve_gene_ids`` — a **Gene id stem** against real gene ids.

    The first surface here to open the **Annotation database**, so every test registers a
    fixture GTF for real and asks the registry, exactly as the gene-category tests do. What
    is asserted is what a caller holding a published table keyed by unversioned ids gets
    back: both of the ids a stem names rather than one of them, and the stems this
    annotation has nothing for said out loud.
    """

    def _registry(self, tmp_path: Path, gtf: str, name: str = "mine") -> AnnotationRegistry:
        """Register ``gtf`` under ``name`` for the ``tiny`` assembly and open the registry."""
        tmp_path.mkdir(parents=True, exist_ok=True)
        source = tmp_path / f"{name}.gtf"
        source.write_text(gtf)
        _register_by_path(tmp_path, source, name)
        return AnnotationRegistry.locate("tiny", tmp_path)

    def test_a_stem_resolves_to_the_versioned_id_or_to_itself_when_unversioned(
        self, tmp_path: Path
    ) -> None:
        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM], "mine")

        assert answer.resolved == {_ALONE_STEM: (_ALONE,)}
        assert answer.unresolved == ()
        assert (answer.assembly, answer.annotation) == ("tiny", "mine")

        # WormBase and SGD never versioned a gene id, and an Ensembl-shaped assumption
        # would leave both unresolvable. An id with no version is its own stem.
        unversioned_dir = tmp_path / "unversioned"
        unversioned_dir.mkdir()
        unversioned = self._registry(unversioned_dir, _GTF)
        unversioned_answer = unversioned.resolve_gene_ids(["g1", "g2"], "mine")
        assert unversioned_answer.resolved == {"g1": ("g1",)}
        assert unversioned_answer.unresolved == ("g2",)

    def test_an_unresolvable_stem_rides_back_and_order_and_repeats_are_handled(
        self, tmp_path: Path
    ) -> None:
        # The `gencode_v50lift37` case: nine stems name two genes each, eight of them a
        # pseudoautosomal pair. Answering with the first would hand back the X copy of a Y
        # gene without ever saying a choice had been made.
        par_stem_registry = self._registry(
            tmp_path / "par-stem", _gtf_declaring(_PAR_Y, _PAR_X, _ALONE)
        )
        par_stem_answer = par_stem_registry.resolve_gene_ids([_PAR_STEM], "mine")
        assert par_stem_answer.resolved[_PAR_STEM] == (_PAR_X, _PAR_Y)

        registry = self._registry(tmp_path, _gtf_declaring(_ALONE))

        answer = registry.resolve_gene_ids([_ALONE_STEM, _ABSENT_STEM], "mine")
        assert answer.unresolved == (_ABSENT_STEM,)
        assert _ABSENT_STEM not in answer.resolved
        # …and what did resolve is still there, so an absence costs the caller nothing else.
        assert answer.gene_ids == [_ALONE]

        empty_answer = registry.resolve_gene_ids([], "mine")
        assert (dict(empty_answer.resolved), empty_answer.unresolved) == ({}, ())

        # A caller passing a few thousand at once reads its own list against the answer.
        par_dir = tmp_path / "par"
        par_dir.mkdir()
        par_registry = self._registry(par_dir, _gtf_declaring(_ALONE, _PAR_X))
        par_answer = par_registry.resolve_gene_ids(
            [_PAR_STEM, _ABSENT_STEM, _ALONE_STEM, _PAR_STEM], "mine"
        )
        assert list(par_answer.resolved) == [_PAR_STEM, _ALONE_STEM]
        assert par_answer.gene_ids == [_PAR_X, _ALONE]
        assert par_answer.unresolved == (_ABSENT_STEM,)

    def test_an_unregistered_name_or_no_name_is_resolved_by_the_same_lookup(
        self, tmp_path: Path
    ) -> None:
        # An exon-level GTF registers as exons alone, since reconstructing the features
        # above them is off by default. Every stem would come back unresolved, which reads
        # as *this annotation has none of your genes* and is a different fact entirely.
        bare_registry = self._registry(tmp_path / "bare", _BARE_GTF)
        with pytest.raises(NoGeneFeaturesError) as no_genes:
            bare_registry.resolve_gene_ids(["g1"], "mine")
        no_genes_message = str(no_genes.value)
        assert "mine" in no_genes_message
        assert "--infer-genes" in no_genes_message  # the argument that rebuilds it with genes

        # Resolved through the same lookup every other question goes through, so the
        # message and its repair are the ones that surface already has.
        registry = AnnotationRegistry.locate(_CURATED_ASSEMBLY, tmp_path)
        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            registry.resolve_gene_ids([_ALONE_STEM], _CURATED)
        assert f"genome register-annotation {_CURATED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

        default_dir = tmp_path / "default"
        default_dir.mkdir()
        default_registry = self._registry(default_dir, _gtf_declaring(_ALONE))
        assert default_registry.default == "mine"
        assert default_registry.resolve_gene_ids([_ALONE_STEM]).annotation == "mine"

    def test_a_merged_annotation_answers_from_its_own_database_for_either_component(
        self, tmp_path: Path
    ) -> None:
        # A merge rewrites seqnames and never a `gene_id`, so its one database holds both
        # components' genes under the ids their own annotations gave them — and resolution
        # reads that database, exactly as it reads any other annotation's.
        worm = tmp_path / "worm.gtf"
        worm.write_text(_gtf_declaring("WBGene00004512"))
        food = tmp_path / "food.gtf"
        food.write_text(_gtf_declaring(_ALONE))
        chrom_sizes = _write_chrom_sizes(
            tmp_path,
            *(f"chrI__{component}" for component in (_WORM_COMPONENT, _FOOD_COMPONENT)),
            assembly=_CHIMERA,
        )
        register_merged_gtf(
            tmp_path,
            "merged",
            [
                MergeSource(_WORM_COMPONENT, _WORM, worm),
                MergeSource(_FOOD_COMPONENT, _FOOD, food),
            ],
            separator="__",
            chrom_sizes=chrom_sizes,
        )
        registry = AnnotationRegistry.locate(_CHIMERA, tmp_path)

        answer = registry.resolve_gene_ids(["WBGene00004512", _ALONE_STEM], "merged")

        assert answer.gene_ids == ["WBGene00004512", _ALONE]
        assert answer.unresolved == ()


def test_the_stem_pass_walks_the_databases_gene_rows_one_at_a_time(tmp_path: Path) -> None:
    # What keeps a GENCODE-sized annotation out of memory: the stem match rides on a walk
    # that yields, so the first gene row is answered before the rest of the table has been
    # read and only what matched is ever held. A read that returned a list would pass every
    # assertion about the answer above and quietly hold 78,000 ids to do it.
    source = tmp_path / "mine.gtf"
    source.write_text(_gtf_declaring(_ALONE, _PAR_X, _PAR_Y))
    annotation = _register_by_path(tmp_path, source, "mine")

    assert inspect.isgeneratorfunction(gene_ids)
    walk = gene_ids(annotation.db)
    assert next(walk) == _ALONE  # answered without the rest of the rows being read
    walk.close()  # and closing it early releases the SQLite connection behind it

    assert list(gene_ids(annotation.db)) == [_ALONE, _PAR_X, _PAR_Y]  # ascending, whole
