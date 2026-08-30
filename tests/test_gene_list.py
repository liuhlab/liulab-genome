"""Tests for genome.gene_list — the curated gene lists shipped inside the package.

The shipped JSON files answer here; no fixture stands in for them, since guarding what
ships is most of what this module exists to do. Nothing touches the network, and nothing
here writes into the package's own data directory: the malformed cases are handed to the
reader as text rather than laid down as files.

The category names are **data**. They come from whatever a curated list declares and
differ per annotation, so what is asserted below is structure and never a closed set of
names — a curator adding one must not break a test.

One test below walks every shipped file's categories and gene ids — the whole-table
structural invariant that catches a newly dropped-in file with no code change here,
alongside the coverage check that every annotation the table lists ships one at all.
"""

from __future__ import annotations

import json
from typing import Any

import pytest

from genome.gene_list import (
    GENE_LIST_SUFFIX,
    CuratedCategory,
    CuratedGeneList,
    CuratedGeneListError,
    GeneCategoryNotDeclaredError,
    GeneListAssemblyMismatchError,
    NoGeneCategoriesError,
    _read_gene_list,
    curated_annotations,
    curated_gene_list,
)
from genome.metadata import annotation_table

#: Annotations the shipped table lists that deliberately ship no curated gene list.
#: Empty, and that is the point: adding a table row without a list is then a visible
#: decision rather than an annotation that silently cannot answer.
_WITHOUT_A_CURATED_LIST: tuple[str, ...] = ()

#: The one category every shipped list declares. Not a vocabulary — the rest differ per
#: annotation and are read off the files — but rRNA is what #111 was opened for, so an
#: annotation that cannot answer for it is one this feature does not serve.
_UNIVERSAL_CATEGORY = "rRNA"

#: Which assembly each shipped annotation belongs to, read from the annotation table —
#: the same fact a curated list carries, and the one it is checked against.
_TABLE_ASSEMBLIES: dict[str, str] = {record.name: record.assembly for record in annotation_table()}


def _payload(**overrides: Any) -> str:
    """Return one well-formed curated gene list as text, with ``overrides`` applied."""
    body: dict[str, Any] = {
        "annotation": "mine",
        "assembly": "tiny",
        "categories": {
            "rRNA": {
                "description": "the mature ribosomal RNA genes",
                "source": "curated by hand from the GTF",
                "gene_ids": ["g1", "g2"],
            }
        },
    }
    body.update(overrides)
    return json.dumps(body)


def _read(text: str) -> CuratedGeneList:
    """Read ``text`` as the curated gene list shipped for the annotation ``mine``."""
    return _read_gene_list(text, annotation="mine", origin=f"mine{GENE_LIST_SUFFIX}")


def _shipped(annotation: str) -> CuratedGeneList:
    """Return the shipped list for ``annotation``, failing by name when none ships."""
    listed = curated_gene_list(annotation)
    assert listed is not None, f"no {annotation}{GENE_LIST_SUFFIX} ships in the package"
    return listed


# ---------------------------------------------------------------------------------------
# What ships: every file valid, and every listed annotation covered
# ---------------------------------------------------------------------------------------


def test_the_package_ships_sorted_deduplicated_annotations_covering_every_table_row() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is exactly the silent zero #111 is about.
    annotations = curated_annotations()

    assert annotations
    assert list(annotations) == sorted(annotations)
    assert len(set(annotations)) == len(annotations)

    # A row added without a list is then a decision someone made out loud, rather than an
    # annotation that answers *absent* for every category and nobody noticed.
    missing = sorted(set(_TABLE_ASSEMBLIES) - set(annotations))
    assert missing == sorted(_WITHOUT_A_CURATED_LIST)


@pytest.mark.parametrize("annotation", curated_annotations())
def test_every_shipped_file_is_valid_matches_its_assembly_and_answers_the_universal_category(
    annotation: str,
) -> None:
    # Loading is what validates, so a file that cannot be trusted raises here rather than
    # answering. The assertions below are the same invariants said again, so that a
    # loosened loader is caught by a failure naming the offending file. The fact the
    # assembly guard rests on: a curated list is pinned to one assembly's annotation, and
    # the table is where that pairing is written down.
    where = f"{annotation}{GENE_LIST_SUFFIX}"
    listed = _shipped(annotation)

    assert listed.annotation == annotation, f"{where}: names another annotation"
    assert listed.assembly == _TABLE_ASSEMBLIES.get(annotation), (
        f"{where}: names an assembly the annotation table does not file it under"
    )
    assert _UNIVERSAL_CATEGORY in listed.categories

    assert listed.categories, f"{where}: declares no categories"
    seen: dict[str, str] = {}
    for name, category in listed.categories.items():
        assert category.category == name, f"{where}: {name} is filed under another name"
        assert category.description.strip(), f"{where}: {name} says nothing about membership"
        assert category.source.strip(), f"{where}: {name} says nothing about where it came from"
        assert category.gene_ids, f"{where}: {name} declares no gene ids"
        assert len(set(category.gene_ids)) == len(category.gene_ids), (
            f"{where}: {name} lists a gene id twice"
        )
        for gene_id in category.gene_ids:
            assert gene_id not in seen, (
                f"{where}: {gene_id} is in both {seen.get(gene_id)} and {name}"
            )
            seen[gene_id] = name


# ---------------------------------------------------------------------------------------
# The raw absence: one function's ``None``, and nobody else's
# ---------------------------------------------------------------------------------------


def test_an_annotation_no_list_ships_for_and_a_path_shaped_name_both_reach_nothing() -> None:
    assert curated_gene_list("no_such_annotation") is None

    # Names are looked up among what ships rather than joined onto the resource
    # directory, so a name shaped like a path finds nothing instead of walking out of it.
    assert curated_gene_list("../annotation_metadata") is None
    assert curated_gene_list("") is None


# ---------------------------------------------------------------------------------------
# A file that is there and cannot be read is a packaging defect
# ---------------------------------------------------------------------------------------


def test_a_well_formed_file_reads_back_as_its_categories_in_file_order() -> None:
    listed = _read(_payload())
    assert listed.annotation == "mine"
    assert listed.assembly == "tiny"
    assert list(listed.categories) == ["rRNA"]
    assert listed.categories["rRNA"].gene_ids == ("g1", "g2")

    spelled = {
        name: {"description": "d", "source": "s", "gene_ids": [f"{name}-1"]}
        for name in ("zeta", "alpha", "mu")
    }
    assert list(_read(_payload(categories=spelled)).categories) == ["zeta", "alpha", "mu"]


def test_invalid_json_and_a_mismatched_annotation_field_both_raise_and_name_the_file() -> None:
    with pytest.raises(CuratedGeneListError, match=f"mine{GENE_LIST_SUFFIX}"):
        _read("{not json at all")

    with pytest.raises(CuratedGeneListError, match="theirs"):
        _read(_payload(annotation="theirs"))


def test_a_missing_top_level_key_names_the_key() -> None:
    # Three keys, three independent failure modes: each is checked in its own block so a
    # future refusal that only catches two of the three still fails here.
    for missing in ("annotation", "assembly", "categories"):
        body = json.loads(_payload())
        del body[missing]

        with pytest.raises(CuratedGeneListError, match=missing):
            _read(json.dumps(body))


def test_various_malformed_category_payloads_each_raise_and_name_the_offender() -> None:
    # Absence is spelled by shipping no file. A file that declares nothing would be a
    # second spelling of it, and the one that reads as emptiness.
    with pytest.raises(CuratedGeneListError, match="categories"):
        _read(_payload(categories={}))

    empty = {"rRNA": {"description": "d", "source": "s", "gene_ids": []}}
    with pytest.raises(CuratedGeneListError, match="rRNA"):
        _read(_payload(categories=empty))

    doubled = {"rRNA": {"description": "d", "source": "s", "gene_ids": ["g1", "g1"]}}
    with pytest.raises(CuratedGeneListError, match="g1"):
        _read(_payload(categories=doubled))

    # A caller summing two categories would count that gene twice, which is the swing
    # #111 measured rather than a tidiness complaint.
    both = {
        "rRNA": {"description": "d", "source": "s", "gene_ids": ["g1"]},
        "Mt_rRNA": {"description": "d", "source": "s", "gene_ids": ["g1"]},
    }
    with pytest.raises(CuratedGeneListError, match="g1"):
        _read(_payload(categories=both))

    numeric = {"rRNA": {"description": "d", "source": "s", "gene_ids": [1, 2]}}
    with pytest.raises(CuratedGeneListError, match="rRNA"):
        _read(_payload(categories=numeric))


def test_a_category_that_explains_nothing_or_is_shaped_wrong_raises() -> None:
    for blank in ("description", "source"):
        entry: dict[str, Any] = {"description": "d", "source": "s", "gene_ids": ["g1"]}
        entry[blank] = "   "
        with pytest.raises(CuratedGeneListError, match=blank):
            _read(_payload(categories={"rRNA": entry}))

    for not_a_mapping in ([], 3):
        with pytest.raises(CuratedGeneListError, match="categories"):
            _read(_payload(categories=not_a_mapping))

    for not_an_object in ({"rRNA": ["g1"]}, {"rRNA": "g1"}):
        with pytest.raises(CuratedGeneListError, match="rRNA"):
            _read(_payload(categories=not_an_object))


# ---------------------------------------------------------------------------------------
# The assembly a list was curated against
# ---------------------------------------------------------------------------------------


def test_a_list_answers_for_its_own_assembly_and_refuses_a_mismatch_as_a_non_lookup_defect() -> (
    None
):
    assert _read(_payload()).check_assembly("tiny") is None

    # A same-named annotation registered against a different assembly: answering would
    # hand back another species' gene ids, so it raises rather than answering.
    with pytest.raises(GeneListAssemblyMismatchError) as excinfo:
        _read(_payload()).check_assembly("hg38")
    message = str(excinfo.value)
    assert "tiny" in message
    assert "hg38" in message

    # It is not one of the two absences below: nothing is missing, and a caller catching
    # LookupError for absence must not swallow this.
    assert issubclass(GeneListAssemblyMismatchError, CuratedGeneListError)
    assert not issubclass(GeneListAssemblyMismatchError, LookupError)


# ---------------------------------------------------------------------------------------
# Absence is not emptiness: two errors, told apart
# ---------------------------------------------------------------------------------------


def test_the_two_absences_are_distinguishable_lookups_with_helpful_messages() -> None:
    # A caller may catch the pair and still act differently on each, which is exactly what
    # #111 asks for: *no categories at all* and *not this category* are different facts.
    assert issubclass(NoGeneCategoriesError, LookupError)
    assert issubclass(GeneCategoryNotDeclaredError, LookupError)
    assert not issubclass(NoGeneCategoriesError, GeneCategoryNotDeclaredError)
    assert not issubclass(GeneCategoryNotDeclaredError, NoGeneCategoriesError)

    error = NoGeneCategoriesError("mine", "tiny", ("wormbase_ws298", "gencode_v50"))
    message = str(error)
    assert "mine" in message
    assert "wormbase_ws298" in message
    assert "gencode_v50" in message
    assert (error.annotation, error.assembly) == ("mine", "tiny")

    # For a merged annotation nothing would be fixed by shipping a list under the merged
    # name: it is the contributing annotations that need one.
    merged = NoGeneCategoriesError("a+b", "x_y", ("gencode_v50",), contributors=("a", "b"))
    assert "a+b" in str(merged)
    assert "a, b" in str(merged)

    declared_error = GeneCategoryNotDeclaredError("mine", "tiny", "tRNA", ("rRNA", "Mt_rRNA"))
    declared_message = str(declared_error)
    assert "tRNA" in declared_message
    assert "rRNA" in declared_message
    assert "Mt_rRNA" in declared_message
    assert declared_error.category == "tRNA"
    assert declared_error.declared == ("rRNA", "Mt_rRNA")


# ---------------------------------------------------------------------------------------
# The values themselves
# ---------------------------------------------------------------------------------------


def test_a_category_carries_its_fields_and_a_curated_list_is_frozen() -> None:
    category = CuratedCategory(
        category="rRNA",
        description="the mature ribosomal RNA genes",
        source="WormBase WS298",
        gene_ids=("b", "a"),
    )
    assert category.category == "rRNA"
    assert category.gene_ids == ("b", "a")

    listed = _read(_payload())
    with pytest.raises(AttributeError):
        listed.assembly = "hg38"  # type: ignore[misc]
