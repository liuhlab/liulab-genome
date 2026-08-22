"""Tests for genome.gene_list — the curated gene lists shipped inside the package.

The shipped JSON files answer here; no fixture stands in for them, since guarding what
ships is most of what this module exists to do. Nothing touches the network, and nothing
here writes into the package's own data directory: the malformed cases are handed to the
reader as text rather than laid down as files.

The category names are **data**. They come from whatever a curated list declares and
differ per annotation, so what is asserted below is structure and never a closed set of
names — a curator adding one must not break a test.
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


def test_the_package_ships_curated_gene_lists_at_all() -> None:
    # The guard under every parametrized test below: with no files at all each of those
    # would collect zero cases and pass, which is exactly the silent zero #111 is about.
    assert curated_annotations()


def test_the_shipped_annotations_are_sorted_and_each_named_once() -> None:
    annotations = curated_annotations()

    assert list(annotations) == sorted(annotations)
    assert len(set(annotations)) == len(annotations)


@pytest.mark.parametrize("annotation", curated_annotations())
def test_every_shipped_file_is_valid(annotation: str) -> None:
    # Loading is what validates, so a file that cannot be trusted raises here rather than
    # answering. The assertions below are the same invariants said again, so that a
    # loosened loader is caught by a failure naming the offending file.
    where = f"{annotation}{GENE_LIST_SUFFIX}"
    listed = _shipped(annotation)

    assert listed.annotation == annotation, f"{where}: names another annotation"
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


@pytest.mark.parametrize("annotation", curated_annotations())
def test_every_shipped_file_names_the_assembly_the_table_files_it_under(annotation: str) -> None:
    # The fact the assembly guard rests on: a curated list is pinned to one assembly's
    # annotation, and the table is where that pairing is written down.
    where = f"{annotation}{GENE_LIST_SUFFIX}"

    assert _shipped(annotation).assembly == _TABLE_ASSEMBLIES.get(annotation), (
        f"{where}: names an assembly the annotation table does not file it under"
    )


@pytest.mark.parametrize("annotation", curated_annotations())
def test_every_shipped_file_can_answer_the_category_this_was_opened_for(annotation: str) -> None:
    assert _UNIVERSAL_CATEGORY in _shipped(annotation).categories


def test_every_annotation_the_table_lists_ships_a_curated_list() -> None:
    # A row added without a list is then a decision someone made out loud, rather than an
    # annotation that answers *absent* for every category and nobody noticed.
    missing = sorted(set(_TABLE_ASSEMBLIES) - set(curated_annotations()))

    assert missing == sorted(_WITHOUT_A_CURATED_LIST)


# ---------------------------------------------------------------------------------------
# The raw absence: one function's ``None``, and nobody else's
# ---------------------------------------------------------------------------------------


def test_an_annotation_no_list_ships_for_is_none() -> None:
    assert curated_gene_list("no_such_annotation") is None


def test_a_name_that_is_a_path_reaches_nothing() -> None:
    # Names are looked up among what ships rather than joined onto the resource directory,
    # so a name shaped like a path finds nothing instead of walking out of it.
    assert curated_gene_list("../annotation_metadata") is None
    assert curated_gene_list("") is None


# ---------------------------------------------------------------------------------------
# A file that is there and cannot be read is a packaging defect
# ---------------------------------------------------------------------------------------


def test_a_well_formed_file_reads_back_as_its_categories() -> None:
    listed = _read(_payload())

    assert listed.annotation == "mine"
    assert listed.assembly == "tiny"
    assert list(listed.categories) == ["rRNA"]
    assert listed.categories["rRNA"].gene_ids == ("g1", "g2")


def test_categories_keep_the_order_the_file_spells_them_in() -> None:
    spelled = {
        name: {"description": "d", "source": "s", "gene_ids": [f"{name}-1"]}
        for name in ("zeta", "alpha", "mu")
    }

    assert list(_read(_payload(categories=spelled)).categories) == ["zeta", "alpha", "mu"]


def test_text_that_is_not_json_names_the_file() -> None:
    with pytest.raises(CuratedGeneListError, match=f"mine{GENE_LIST_SUFFIX}"):
        _read("{not json at all")


@pytest.mark.parametrize("missing", ["annotation", "assembly", "categories"])
def test_a_missing_key_names_the_key(missing: str) -> None:
    body = json.loads(_payload())
    del body[missing]

    with pytest.raises(CuratedGeneListError, match=missing):
        _read(json.dumps(body))


def test_an_annotation_field_disagreeing_with_the_file_name_raises() -> None:
    with pytest.raises(CuratedGeneListError, match="theirs"):
        _read(_payload(annotation="theirs"))


def test_a_file_declaring_no_categories_at_all_raises() -> None:
    # Absence is spelled by shipping no file. A file that declares nothing would be a
    # second spelling of it, and the one that reads as emptiness.
    with pytest.raises(CuratedGeneListError, match="categories"):
        _read(_payload(categories={}))


def test_a_declared_category_with_no_gene_ids_raises() -> None:
    empty = {"rRNA": {"description": "d", "source": "s", "gene_ids": []}}

    with pytest.raises(CuratedGeneListError, match="rRNA"):
        _read(_payload(categories=empty))


@pytest.mark.parametrize("blank", ["description", "source"])
def test_a_category_that_explains_nothing_raises(blank: str) -> None:
    entry: dict[str, Any] = {"description": "d", "source": "s", "gene_ids": ["g1"]}
    entry[blank] = "   "

    with pytest.raises(CuratedGeneListError, match=blank):
        _read(_payload(categories={"rRNA": entry}))


def test_the_same_gene_id_twice_in_one_category_raises() -> None:
    doubled = {"rRNA": {"description": "d", "source": "s", "gene_ids": ["g1", "g1"]}}

    with pytest.raises(CuratedGeneListError, match="g1"):
        _read(_payload(categories=doubled))


def test_the_same_gene_id_in_two_categories_raises() -> None:
    # A caller summing two categories would count that gene twice, which is the swing
    # #111 measured rather than a tidiness complaint.
    both = {
        "rRNA": {"description": "d", "source": "s", "gene_ids": ["g1"]},
        "Mt_rRNA": {"description": "d", "source": "s", "gene_ids": ["g1"]},
    }

    with pytest.raises(CuratedGeneListError, match="g1"):
        _read(_payload(categories=both))


@pytest.mark.parametrize("malformed", [[], "rRNA", 3])
def test_categories_that_are_not_a_mapping_raise(malformed: Any) -> None:
    with pytest.raises(CuratedGeneListError, match="categories"):
        _read(_payload(categories=malformed))


@pytest.mark.parametrize("malformed", [{"rRNA": ["g1"]}, {"rRNA": "g1"}])
def test_a_category_that_is_not_an_object_raises(malformed: Any) -> None:
    with pytest.raises(CuratedGeneListError, match="rRNA"):
        _read(_payload(categories=malformed))


def test_a_gene_id_that_is_not_text_raises() -> None:
    numeric = {"rRNA": {"description": "d", "source": "s", "gene_ids": [1, 2]}}

    with pytest.raises(CuratedGeneListError, match="rRNA"):
        _read(_payload(categories=numeric))


# ---------------------------------------------------------------------------------------
# The assembly a list was curated against
# ---------------------------------------------------------------------------------------


def test_a_list_answers_for_the_assembly_it_was_curated_against() -> None:
    assert _read(_payload()).check_assembly("tiny") is None


def test_a_list_refuses_to_answer_for_another_assembly() -> None:
    # A same-named annotation registered against a different assembly: answering would
    # hand back another species' gene ids, so it raises rather than answering.
    with pytest.raises(GeneListAssemblyMismatchError) as excinfo:
        _read(_payload()).check_assembly("hg38")

    message = str(excinfo.value)
    assert "tiny" in message
    assert "hg38" in message


def test_the_assembly_mismatch_is_a_defect_in_what_ships_rather_than_a_lookup() -> None:
    # It is not one of the two absences: nothing is missing, and a caller catching
    # LookupError for absence must not swallow this.
    assert issubclass(GeneListAssemblyMismatchError, CuratedGeneListError)
    assert not issubclass(GeneListAssemblyMismatchError, LookupError)


# ---------------------------------------------------------------------------------------
# Absence is not emptiness: two errors, told apart
# ---------------------------------------------------------------------------------------


def test_the_two_absences_are_distinguishable_and_both_are_lookups() -> None:
    # A caller may catch the pair and still act differently on each, which is exactly what
    # #111 asks for: *no categories at all* and *not this category* are different facts.
    assert issubclass(NoGeneCategoriesError, LookupError)
    assert issubclass(GeneCategoryNotDeclaredError, LookupError)
    assert not issubclass(NoGeneCategoriesError, GeneCategoryNotDeclaredError)
    assert not issubclass(GeneCategoryNotDeclaredError, NoGeneCategoriesError)


def test_declaring_no_categories_says_who_does_and_what_would_fix_it() -> None:
    error = NoGeneCategoriesError("mine", "tiny", ("wormbase_ws298", "gencode_v50"))

    message = str(error)
    assert "mine" in message
    assert "wormbase_ws298" in message
    assert "gencode_v50" in message
    assert (error.annotation, error.assembly) == ("mine", "tiny")


def test_declaring_no_categories_names_the_contributors_of_a_merge() -> None:
    # For a merged annotation nothing would be fixed by shipping a list under the merged
    # name: it is the contributing annotations that need one.
    error = NoGeneCategoriesError("a+b", "x_y", ("gencode_v50",), contributors=("a", "b"))

    assert "a+b" in str(error)
    assert "a, b" in str(error)


def test_a_category_nobody_declared_lists_the_ones_that_are() -> None:
    error = GeneCategoryNotDeclaredError("mine", "tiny", "tRNA", ("rRNA", "Mt_rRNA"))

    message = str(error)
    assert "tRNA" in message
    assert "rRNA" in message
    assert "Mt_rRNA" in message
    assert error.category == "tRNA"
    assert error.declared == ("rRNA", "Mt_rRNA")


# ---------------------------------------------------------------------------------------
# The values themselves
# ---------------------------------------------------------------------------------------


def test_a_category_carries_its_own_name_and_its_ids_in_file_order() -> None:
    category = CuratedCategory(
        category="rRNA",
        description="the mature ribosomal RNA genes",
        source="WormBase WS298",
        gene_ids=("b", "a"),
    )

    assert category.category == "rRNA"
    assert category.gene_ids == ("b", "a")


def test_a_curated_list_is_frozen() -> None:
    listed = _read(_payload())

    with pytest.raises(AttributeError):
        listed.assembly = "hg38"  # type: ignore[misc]
