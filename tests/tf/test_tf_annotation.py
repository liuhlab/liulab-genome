"""Tests for the TF halves' annotation crossing — a shipped table met with real gene ids.

``genome.tf.gene.annotation`` and ``genome.tf.cofactor.annotation`` together, because the
two are one crossing asked twice: a census against an annotation and a **Cofactor table**
against one differ in which shipped file is read and in what a row of it says, and in
nothing about the annotation, so they register their fixtures the same way and one test
below needs both surfaces at once — the worm assembly one half answers and the other
refuses.

The shipped files answer throughout, as the curated lists do for gene categories: what a
census holds and what a publisher lists are their publishers' business and no fixture
pretends otherwise. What is asserted is the crossing — the **Gene id stem**s arriving as
this annotation's own gene ids, what rode back unresolved, and that an assembly nothing
published can answer for raises rather than answering with nothing.

Nothing here prepares an assembly, fetches anything or reads the **Data dir**: every test
registers one small GTF into ``tmp_path`` and asks. Registration itself is tested under
``tests/annotation/``, whose fixture helpers are reused rather than spelled again.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genome.annotation import AnnotationNotRegisteredError, AnnotationRegistry
from genome.assembly.metadata import assembly_metadata
from genome.tf.cofactor import UNIFORM_COLUMNS as COFACTOR_UNIFORM_COLUMNS
from genome.tf.cofactor import (
    CofactorTable,
    NoCofactorTableError,
    cofactor_table,
    resolve_tf_cofactors,
    tf_cofactor_list,
)
from genome.tf.gene import (
    UNIFORM_COLUMNS,
    NoTFCensusError,
    TFGeneTable,
    resolve_tf_genes,
    tf_gene_list,
    tf_gene_table,
)
from genome.tf.species import UnknownSpeciesError

from ..annotation.conftest import _CHIMERA, _CURATED, _CURATED_ASSEMBLY, _register_by_path
from ..annotation.test_stems import _gtf_declaring


def _registry_declaring(
    tmp_path: Path, *gene_ids: str, assembly: str, name: str = "mine"
) -> AnnotationRegistry:
    """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry.

    The setup both shipped-table surfaces share — a census against an annotation and a
    cofactor table against one differ in which file is read and in nothing about the
    annotation, so they register the same way. ``assembly`` is required: which species an
    assembly is is exactly what those tests are about, and a default would hide it.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    source = tmp_path / f"{name}.gtf"
    source.write_text(_gtf_declaring(*gene_ids))
    _register_by_path(tmp_path, source, name, assembly=assembly)
    return AnnotationRegistry.locate(assembly, tmp_path)


# ---------------------------------------------------------------------------------------
# Which of an annotation's genes a published census judges transcription factors
# ---------------------------------------------------------------------------------------

#: An assembly whose species has a census, one whose species has none, and one nothing
#: names a species for at all — the chimera, which is more than one species by
#: construction. Named rather than derived: which species a census ships for is exactly
#: what these tests are about, so the pairing is written down where the test controls it.
_CENSUSED_ASSEMBLY = _CURATED_ASSEMBLY
_UNCENSUSED_ASSEMBLY = "ce11"

#: ZBED1 as GENCODE spells it on a pseudoautosomal region: one **Gene id stem**, two gene
#: ids, and a gene Lambert judges a transcription factor — so the collision a resolver must
#: not silently pick from is a **TF gene**'s own, rather than an invented one.
_TF_PAR_X = "ENSG00000214717.13"
_TF_PAR_Y = "ENSG00000214717.13_PAR_Y"
_TF_PAR_STEM = "ENSG00000214717"


def _census() -> TFGeneTable:
    """The census shipped for the species of the assembly the tests below register against.

    Read off the shipped file rather than named, exactly as the curated gene lists are:
    which genes are transcription factors is the census's judgement, and no fixture may
    stand in for it.
    """
    species = assembly_metadata(_CENSUSED_ASSEMBLY).species
    assert species is not None, f"{_CENSUSED_ASSEMBLY} has no species in the assembly table"
    census = tf_gene_table(species)
    assert census is not None, f"no census ships for {species}"
    return census


def _census_row(stem: str) -> dict[str, str | None]:
    """The census's own row for one **Gene id stem**, keyed by its own column names."""
    census = _census()
    row = census.rows[census.gene_id_stems.index(stem)]
    return dict(zip(census.columns, row, strict=True))


def _rejected_stem() -> str:
    """A stem the census assessed and judged *not* to be a transcription factor."""
    census = _census()
    positive = set(census.assessed_positive)
    return next(stem for stem in census.gene_id_stems if stem not in positive)


def _stem_assessed(assessment: str) -> str:
    """A stem the census judged a TF under one of its own **TF assessment** grades."""
    return next(
        stem
        for stem in _census().assessed_positive
        if _census_row(stem)["tf_assessment"] == assessment
    )


def _versioned(stem: str) -> str:
    """One gene id an annotation might spell that stem with — a version it never carries."""
    return f"{stem}.1"


class TestTFGeneList:
    """``resolve_tf_genes`` — where a published census meets a registered annotation.

    The shipped census answers here, as the shipped curated lists do for gene categories:
    what a census holds is its publisher's business and no fixture pretends otherwise. What
    is asserted is the crossing — the census's **Gene id stem**s arriving as this
    annotation's own gene ids, what rode back unresolved, and that an assembly no census
    can answer for raises rather than answering with nothing.
    """

    def _registry(
        self,
        tmp_path: Path,
        *gene_ids: str,
        assembly: str = _CENSUSED_ASSEMBLY,
        name: str = "mine",
    ) -> AnnotationRegistry:
        """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry."""
        return _registry_declaring(tmp_path, *gene_ids, assembly=assembly, name=name)

    def test_the_census_arrives_filtered_widened_or_tightened_with_judgement_and_provenance(
        self, tmp_path: Path
    ) -> None:
        # The whole point of the surface: the answer joins to a counts matrix keyed by this
        # annotation's ids, with no normalisation left for the caller.
        stems = _census().assessed_positive[:2]
        gene_ids = [_versioned(stem) for stem in stems]
        registry = self._registry(tmp_path / "two", *gene_ids)

        answer = resolve_tf_genes(registry, "mine")

        assert (answer.assembly, answer.annotation) == (_CENSUSED_ASSEMBLY, "mine")
        assert [gene.gene_id_stem for gene in answer.genes] == list(stems)
        assert answer.gene_ids == gene_ids

        # The common case is not 2,765 rows to filter down to 1,639: only the genes the
        # census judged transcription factors are carried.
        positive, rejected = _census().assessed_positive[0], _rejected_stem()
        filtered = self._registry(tmp_path / "filtered", _versioned(positive), _versioned(rejected))
        filtered_answer = resolve_tf_genes(filtered, "mine")
        assert [gene.gene_id_stem for gene in filtered_answer.genes] == [positive]
        assert [gene.is_tf for gene in filtered_answer.genes] == [True]

        census = _census()
        stem = census.assessed_positive[0]
        single = self._registry(tmp_path / "single", _versioned(stem))
        single_answer = resolve_tf_genes(single, "mine")
        gene = single_answer.genes[0]
        cells = _census_row(stem)
        assert gene.symbol == cells["symbol"]
        assert gene.dbd_family == cells["dbd_family"]
        # Everything the publisher recorded beyond the uniform four, under its own names:
        # the assessment, the binding mode, the motif status, the KRAB flag and the votes.
        assert dict(gene.judgements) == {
            name: cells[name] for name in census.columns[len(UNIFORM_COLUMNS) :]
        }
        # The verdict travels with the census that reached it.
        assert single_answer.provenance == census.provenance
        assert single_answer.provenance.publisher
        assert single_answer.provenance.version
        assert single_answer.provenance.pubmed_id
        assert single_answer.species == assembly_metadata(_CENSUSED_ASSEMBLY).species

        # Widening carries the verdict rather than dropping it: a gene the census
        # assessed and turned down arrives saying so, which is not the same fact as a gene
        # it never looked at, and that one is absent from both answers.
        registry = self._registry(tmp_path / "widen", _versioned(positive), _versioned(rejected))

        widened = resolve_tf_genes(registry, "mine", include_rejected=True)

        assert {gene.gene_id_stem: gene.is_tf for gene in widened.genes} == {
            positive: True,
            rejected: False,
        }

        # The **TF assessment** is graded, and tightening to `Known motif` or loosening to
        # include `Inferred motif` is a re-filter on what the answer already carries rather
        # than a second flag this package invents.
        known, inferred = _stem_assessed("Known motif"), _stem_assessed("Inferred motif")
        graded = self._registry(tmp_path / "graded", _versioned(known), _versioned(inferred))

        answer = resolve_tf_genes(graded, "mine")

        assert {gene.gene_id_stem for gene in answer.genes} == {known, inferred}
        assert [
            gene.gene_id_stem
            for gene in answer.genes
            if gene.judgements["tf_assessment"] == "Known motif"
        ] == [known]

    def test_unresolved_stems_ride_back_and_a_par_stem_answers_with_both_ids(
        self, tmp_path: Path
    ) -> None:
        census = _census()
        carried, absent = census.assessed_positive[0], census.assessed_positive[1]
        registry = self._registry(tmp_path / "unresolved", _versioned(carried))

        answer = resolve_tf_genes(registry, "mine")

        assert [gene.gene_id_stem for gene in answer.genes] == [carried]
        assert absent in answer.unresolved
        # Every stem the census judged a transcription factor is accounted for one way or
        # it holds and this annotation does not is visible rather than dropped.
        assert len(answer.genes) + len(answer.unresolved) == len(census.assessed_positive)

        assert _TF_PAR_STEM in census.assessed_positive
        par_registry = self._registry(tmp_path / "par", _TF_PAR_Y, _TF_PAR_X)
        par_answer = resolve_tf_genes(par_registry, "mine")
        assert [gene.gene_ids for gene in par_answer.genes] == [(_TF_PAR_X, _TF_PAR_Y)]
        assert par_answer.gene_ids == [_TF_PAR_X, _TF_PAR_Y]

    def test_the_species_follows_the_assembly_never_the_gtf_and_the_two_absences_differ(
        self, tmp_path: Path
    ) -> None:
        # Human gene ids registered for a worm assembly. Asking for one species'
        # transcription factors while holding another's assembly is not expressible, so this
        # is answered about the assembly's own species and never about what is in the GTF.
        uncensused = self._registry(tmp_path / "worm", _TF_PAR_X, assembly=_UNCENSUSED_ASSEMBLY)
        with pytest.raises(NoTFCensusError) as no_census:
            resolve_tf_genes(uncensused, "mine")
        message = str(no_census.value)
        assert str(assembly_metadata(_UNCENSUSED_ASSEMBLY).species) in message
        assert str(assembly_metadata(_CENSUSED_ASSEMBLY).species) in message

        # And an assembly nothing names a species for — the chimera, or a plain local
        # key — says so rather than guessing, one representative of that whole class.
        unnamed = self._registry(tmp_path / "chimera", _TF_PAR_X, assembly=_CHIMERA)
        with pytest.raises(UnknownSpeciesError) as no_species:
            resolve_tf_genes(unnamed, "mine")
        unnamed_message = str(no_species.value)
        assert _CHIMERA in unnamed_message
        assert str(assembly_metadata(_CENSUSED_ASSEMBLY).species) in unnamed_message

        # As the curated gene lists' pair already are: *no census ships for this species*
        # and *nothing says what species this is* are different answers, both lookups, and
        # neither is an empty collection.
        assert isinstance(no_census.value, NoTFCensusError)
        assert isinstance(no_species.value, UnknownSpeciesError)
        assert not isinstance(no_census.value, UnknownSpeciesError)
        assert not isinstance(no_species.value, NoTFCensusError)

    def test_an_unregistered_name_earns_its_error_and_naming_none_asks_the_default(
        self, tmp_path: Path
    ) -> None:
        registry = AnnotationRegistry.locate(_CENSUSED_ASSEMBLY, tmp_path / "unregistered")

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            resolve_tf_genes(registry, _CURATED)

        assert f"genome register-annotation {_CENSUSED_ASSEMBLY} {_CURATED}" in str(excinfo.value)

        default_registry = self._registry(tmp_path / "default", _TF_PAR_X, name=_CURATED)
        assert default_registry.default == _CURATED
        assert resolve_tf_genes(default_registry).annotation == _CURATED

    def test_the_json_record_carries_the_genes_and_provenance_and_it_answers_by_assembly_name(
        self, tmp_path: Path
    ) -> None:
        # What ``--json`` has to be able to emit: the genes with their **TF assessment** and
        # **DBD family**, the census's provenance, and the stems that resolved to nothing.
        stem = _census().assessed_positive[0]
        registry = self._registry(tmp_path, _versioned(stem))

        payload = resolve_tf_genes(registry, "mine").as_json()

        assert payload["assembly"] == _CENSUSED_ASSEMBLY
        assert payload["gene_ids"] == [_versioned(stem)]
        assert payload["provenance"]["pubmed_id"] == _census().provenance.pubmed_id
        assert payload["unresolved"]
        gene = payload["genes"][0]
        assert (gene["gene_id_stem"], gene["gene_ids"]) == (stem, [_versioned(stem)])
        assert gene["dbd_family"] == _census_row(stem)["dbd_family"]
        assert gene["judgements"]["tf_assessment"] == _census_row(stem)["tf_assessment"]
        assert json.loads(json.dumps(payload)) == payload  # serializes as it stands

        answer = tf_gene_list(_CENSUSED_ASSEMBLY, annotation="mine", cache_dir=tmp_path)
        assert [gene.gene_id_stem for gene in answer.genes] == [stem]
        assert answer.provenance == _census().provenance


# ---------------------------------------------------------------------------------------
# Which of an annotation's genes a publisher lists as transcription cofactors
# ---------------------------------------------------------------------------------------

#: An assembly whose species has a cofactor table, one whose species has none, and the
#: worm — which has a table although no publisher has released a worm TF census, so it is
#: the one assembly the two halves answer differently for. Named rather than derived, for
#: the reason the census pairing above is: which species ships what is what these tests
#: are about. Yeast rather than human is the untabled one on purpose: human has no table
#: only until one is built for it, and a test pinned to that would go green by accident.
_TABLED_ASSEMBLY = "mm39"
_UNTABLED_ASSEMBLY = "sacCer3"
_WORM_ASSEMBLY = "ce11"


def _cofactors(assembly: str = _TABLED_ASSEMBLY) -> CofactorTable:
    """The **Cofactor table** shipped for that assembly's species.

    Read off the shipped file rather than named, exactly as the census is: which genes a
    publisher lists as cofactors is the publisher's judgement, and no fixture may stand
    in for it.
    """
    species = assembly_metadata(assembly).species
    assert species is not None, f"{assembly} has no species in the assembly table"
    table = cofactor_table(species)
    assert table is not None, f"no cofactor table ships for {species}"
    return table


def _cofactor_row(stem: str, assembly: str = _TABLED_ASSEMBLY) -> dict[str, str | None]:
    """The table's own row for one **Gene id stem**, keyed by its own column names."""
    table = _cofactors(assembly)
    row = table.rows[table.gene_id_stems.index(stem)]
    return dict(zip(table.columns, row, strict=True))


class TestTFCofactorList:
    """``resolve_tf_cofactors`` — a shipped cofactor table meets a registered annotation.

    The counterpart of :class:`TestTFGeneList`, registering its fixture annotations the
    same way and asserting the same crossing: the table's **Gene id stem**s arriving as
    this annotation's own gene ids, what rode back unresolved, and that an assembly no
    published table can answer for raises rather than answering with nothing. The shipped
    table answers throughout; what it holds is its publisher's business.
    """

    def _registry(
        self,
        tmp_path: Path,
        *gene_ids: str,
        assembly: str = _TABLED_ASSEMBLY,
        name: str = "mine",
    ) -> AnnotationRegistry:
        """Register a GTF declaring ``gene_ids`` under ``name`` and open the registry."""
        return _registry_declaring(tmp_path, *gene_ids, assembly=assembly, name=name)

    def test_the_table_arrives_carrying_uniform_and_publisher_columns_and_provenance(
        self, tmp_path: Path
    ) -> None:
        # The whole point of the surface: the answer joins to a counts matrix keyed by this
        # annotation's ids, with no normalisation left for the caller.
        stems = _cofactors().cofactor_stems[:2]
        gene_ids = [_versioned(stem) for stem in stems]
        registry = self._registry(tmp_path / "two", *gene_ids)

        answer = resolve_tf_cofactors(registry, "mine")

        assert (answer.assembly, answer.annotation) == (_TABLED_ASSEMBLY, "mine")
        assert [entry.gene_id_stem for entry in answer.cofactors] == list(stems)
        assert answer.gene_ids == gene_ids

        table = _cofactors()
        stem = table.cofactor_stems[0]
        single = self._registry(tmp_path / "single", _versioned(stem))
        single_answer = resolve_tf_cofactors(single, "mine")
        entry = single_answer.cofactors[0]
        cells = _cofactor_row(stem)
        assert (entry.symbol, entry.source) == (cells["symbol"], cells["source"])
        assert entry.is_cofactor is True
        # Everything the publisher classified it with, under that publisher's own
        # namespaced name — the AnimalTFDB family and the category joined onto it.
        assert dict(entry.classifications) == {
            name: cells[name] for name in table.columns[len(COFACTOR_UNIFORM_COLUMNS) :]
        }
        assert "animaltfdb_family" in entry.classifications

        # Membership travels with the publishers that listed the gene.
        assert single_answer.provenance == table.provenance
        assert single_answer.provenance.sources
        assert all(
            source.publisher and source.pubmed_id for source in single_answer.provenance.sources
        )
        assert single_answer.species == assembly_metadata(_TABLED_ASSEMBLY).species

    def test_unresolved_stems_ride_back_and_a_par_stem_answers_with_both_ids(
        self, tmp_path: Path
    ) -> None:
        table = _cofactors()
        carried, absent = table.cofactor_stems[0], table.cofactor_stems[1]
        registry = self._registry(tmp_path / "unresolved", _versioned(carried))

        answer = resolve_tf_cofactors(registry, "mine")

        assert [entry.gene_id_stem for entry in answer.cofactors] == [carried]
        assert absent in answer.unresolved
        # Every stem the table lists is accounted for one way or the other, so what the
        # publisher holds and this annotation does not is visible rather than dropped.
        assert len(answer.cofactors) + len(answer.unresolved) == len(table.cofactor_stems)

        # The collision is what is under test rather than the biology: two gene ids that
        # reduce to one stem, in the shape GENCODE's pseudoautosomal copies have, so a
        # resolver taking the first would hand back one of them without saying it chose.
        stem = table.cofactor_stems[0]
        first, second = _versioned(stem), f"{_versioned(stem)}_PAR_Y"
        par_registry = self._registry(tmp_path / "par", second, first)
        par_answer = resolve_tf_cofactors(par_registry, "mine")
        assert [entry.gene_ids for entry in par_answer.cofactors] == [(first, second)]
        assert par_answer.gene_ids == [first, second]

    def test_the_species_follows_the_assembly_never_the_gtf_and_the_two_absences_differ(
        self, tmp_path: Path
    ) -> None:
        # Mouse gene ids registered for a yeast assembly. Asking for one species' cofactors
        # while holding another's assembly is not expressible, so this is answered about the
        # assembly's own species and never about what is in the GTF.
        stem = _cofactors().cofactor_stems[0]
        untabled = self._registry(
            tmp_path / "untabled", _versioned(stem), assembly=_UNTABLED_ASSEMBLY
        )

        with pytest.raises(NoCofactorTableError) as no_table:
            resolve_tf_cofactors(untabled, "mine")

        message = str(no_table.value)
        assert str(assembly_metadata(_UNTABLED_ASSEMBLY).species) in message
        assert str(assembly_metadata(_TABLED_ASSEMBLY).species) in message

        # And an assembly nothing names a species for — the chimera, or a plain local
        # key — says so rather than guessing, one representative of that whole class.
        unnamed = self._registry(tmp_path / "chimera", _versioned(stem), assembly=_CHIMERA)
        with pytest.raises(UnknownSpeciesError) as no_species:
            resolve_tf_cofactors(unnamed, "mine")
        unnamed_message = str(no_species.value)
        assert _CHIMERA in unnamed_message
        assert "cofactor table" in unnamed_message
        assert str(assembly_metadata(_TABLED_ASSEMBLY).species) in unnamed_message

        # As the census half's pair is: *nobody has published a table for this species* and
        # *nothing says what species this is* are different answers, both lookups, and
        # neither is an empty collection.
        assert isinstance(no_table.value, NoCofactorTableError)
        assert isinstance(no_species.value, UnknownSpeciesError)
        assert not isinstance(no_table.value, UnknownSpeciesError)
        assert not isinstance(no_species.value, NoCofactorTableError)

        # The asymmetry, pinned: AnimalTFDB assessed worm cofactors and nobody has released
        # a worm TF census, so one assembly gets two different answers. That is the
        # publishers' shape and not a defect, and a test says so where it would otherwise
        # be filed as one.
        worm_stem = _cofactors(_WORM_ASSEMBLY).cofactor_stems[0]
        worm = self._registry(tmp_path / "worm", _versioned(worm_stem), assembly=_WORM_ASSEMBLY)
        worm_answer = resolve_tf_cofactors(worm, "mine")
        assert [entry.gene_id_stem for entry in worm_answer.cofactors] == [worm_stem]
        with pytest.raises(NoTFCensusError):
            resolve_tf_genes(worm, "mine")

    def test_an_unregistered_name_earns_its_error_and_naming_none_asks_the_default(
        self, tmp_path: Path
    ) -> None:
        registry = AnnotationRegistry.locate(_TABLED_ASSEMBLY, tmp_path / "unregistered")

        with pytest.raises(AnnotationNotRegisteredError) as excinfo:
            resolve_tf_cofactors(registry, "gencode_vM39")

        assert f"genome register-annotation {_TABLED_ASSEMBLY} gencode_vM39" in str(excinfo.value)

        stem = _cofactors().cofactor_stems[0]
        default_registry = self._registry(
            tmp_path / "default", _versioned(stem), name="gencode_vM39"
        )
        assert default_registry.default == "gencode_vM39"
        assert resolve_tf_cofactors(default_registry).annotation == "gencode_vM39"

    def test_the_json_record_carries_provenance_and_it_answers_by_assembly_name_too(
        self, tmp_path: Path
    ) -> None:
        # What ``--json`` has to be able to emit, in the shape the TF gene list's answer
        # already uses: the entries with the publisher's own classification of each, the
        # provenance to cite, and the stems that resolved to nothing.
        stem = _cofactors().cofactor_stems[0]
        registry = self._registry(tmp_path, _versioned(stem))

        payload = resolve_tf_cofactors(registry, "mine").as_json()

        assert payload["assembly"] == _TABLED_ASSEMBLY
        assert payload["gene_ids"] == [_versioned(stem)]
        assert payload["provenance"]["sources"][0]["pubmed_id"] == (
            _cofactors().provenance.sources[0].pubmed_id
        )
        assert payload["unresolved"]
        entry = payload["cofactors"][0]
        assert (entry["gene_id_stem"], entry["gene_ids"]) == (stem, [_versioned(stem)])
        assert entry["source"] == _cofactor_row(stem)["source"]
        assert (
            entry["classifications"]["animaltfdb_family"]
            == (_cofactor_row(stem)["animaltfdb_family"])
        )
        assert json.loads(json.dumps(payload)) == payload  # serializes as it stands

        # One code path, asserted whole rather than sampled: the module-level function
        # and the method reach the same answer, so a shell surface over either says one
        # thing.
        answer = tf_cofactor_list(_TABLED_ASSEMBLY, annotation="mine", cache_dir=tmp_path)
        assert [entry.gene_id_stem for entry in answer.cofactors] == [stem]
        assert answer == resolve_tf_cofactors(registry, "mine")
