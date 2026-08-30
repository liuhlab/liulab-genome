"""The curated gene lists — which genes one annotation puts in a named **Gene category**.

One small hand-maintained JSON file ships inside the package per annotation, under
``data/gene_list/<annotation>.curated_gene_list.json``. Each declares the categories that
annotation can answer for — ``rRNA`` and whatever else was curated for it — and, per
category, what membership means, where membership came from, and the gene ids in it.

A **Curated gene list** rather than the GTF's own biotype attribute, and that is the whole
point of the file existing: the attribute is spelled ``gene_type`` by GENCODE and
``gene_biotype`` by WormBase and RefSeq, the taxonomies behind those spellings do not
agree, and ``sacCer3``'s carries no biotype attribute at all — so a caller deriving
categories from the GTF reports no ribosomal RNA for yeast and never finds out (ADR-0011).
Which categories exist is therefore **data**: it differs per annotation, and nothing here
knows a category vocabulary.

This module sits beside :mod:`genome.metadata` and for the same reason: it reads a shipped
package resource, never the **Data dir** and never the network. It is pure — it holds no
opinion about what is registered on this machine, which is
:class:`~genome.io.annotation.registry.AnnotationRegistry`'s question, and it is where that registry gets
the answer it attributes and returns.

**Absence is not emptiness**, and every accessor here is built around it.
:func:`curated_gene_list` answers ``None`` for an annotation no list ships for — the raw
absence, and the one place ``None`` is how it is said, because this is the layer below the
one a caller touches. Everything above turns that into
:class:`NoGeneCategoriesError`, whose sibling :class:`GeneCategoryNotDeclaredError` says
the other thing: this annotation declares categories, and not that one. Both subclass
:class:`LookupError`, so a caller may catch the pair and still tell them apart, and
neither is ever an empty collection. A file that ships is held to the same rule: a
declared category always carries at least one gene id, and a file declaring nothing at all
is a packaging defect — absence is spelled by shipping no file, and a second spelling of
it would be the one that reads as emptiness.

A shipped file that cannot be trusted never answers. It is validated as it is loaded and
raises :class:`CuratedGeneListError` naming the file and what is wrong with it, since a
hand-curated file that ships broken is a defect in this package rather than anything the
caller did.

Examples
--------
>>> from genome.gene_list import curated_annotations, curated_gene_list
>>> "gencode_v50" in curated_annotations()
True
>>> curated_gene_list("gencode_v50").assembly
'hg38'
>>> curated_gene_list("no_such_annotation") is None
True
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from types import MappingProxyType
from typing import Any

#: Directory inside the package holding one curated gene list per annotation.
GENE_LIST_SUBDIR = "data/gene_list"

#: What one of those files is called: the **Registered name**, then this. The name is what
#: a list is looked up by, so the two halves are never spelled apart.
GENE_LIST_SUFFIX = ".curated_gene_list.json"

#: Top-level keys every curated gene list carries. Anything else in a file is ignored, so
#: a curator may keep notes beside these without this module having to know about them.
_REQUIRED_KEYS = ("annotation", "assembly", "categories")


class CuratedGeneListError(ValueError):
    """A shipped curated gene list cannot be read, so it is not allowed to answer.

    A **packaging defect** and not a caller error: these files ship inside the package,
    so bad JSON, a missing key, a category with no genes, a gene id in two categories or
    an ``annotation`` field disagreeing with the file name are all faults in what was
    committed here. A :class:`ValueError`, because a hand-curated file that says something
    the format does not is a bad value rather than a broken program.

    The message names the file and what is wrong with it, since fixing the file is the
    only thing anyone can do about it.

    Examples
    --------
    >>> try:
    ...     _read_gene_list("{", annotation="mine", origin="mine.curated_gene_list.json")
    ... except CuratedGeneListError as error:
    ...     print("mine.curated_gene_list.json" in str(error))
    True
    """


class GeneListAssemblyMismatchError(CuratedGeneListError):
    """A curated list was asked for by an annotation registered against another assembly.

    The list names the **Assembly** it was curated against, and that is checked rather
    than assumed: an annotation registered under a name whose curated list belongs to a
    different reference would otherwise answer with another species' gene ids, and mixing
    builds is an error and not a warning (ADR-0003, ADR-0005).

    Not one of the two absences: nothing is missing, so a caller catching
    :class:`LookupError` for absence does not swallow this.

    Examples
    --------
    >>> listed = CuratedGeneList("mine", "ce11", {})
    >>> try:
    ...     listed.check_assembly("hg38")
    ... except GeneListAssemblyMismatchError as error:
    ...     print("ce11" in str(error) and "hg38" in str(error))
    True
    """


class NoGeneCategoriesError(LookupError):
    """Nothing is known about this annotation's categories, so none can be asked of it.

    The first of the two absences, and the one a caller must never read as *this
    annotation has no rRNA genes*: no curated list ships for it, so the question was not
    answered rather than answered with nothing. The message says which annotations do
    declare categories and that answering for this one means shipping a curated list.

    A :class:`LookupError`, so a caller may catch it together with
    :class:`GeneCategoryNotDeclaredError` and still act differently on each.

    Parameters
    ----------
    annotation : str
        The **Registered name** that was asked about — the merged name for a **Merged
        annotation**.
    assembly : str
        The **Assembly** it is registered for.
    curated : iterable of str
        The annotations a curated list does ship for.
    contributors : iterable of str, optional
        For a **Merged annotation**, the annotations it merges. Shipping a list under the
        merged name would fix nothing; it is these that need one.

    Attributes
    ----------
    annotation : str
        The name asked about.
    assembly : str
        The assembly it is registered for.
    curated : tuple of str
        The annotations that do ship a curated list.
    contributors : tuple of str
        The contributing annotations, empty for anything but a merge.

    Examples
    --------
    >>> try:
    ...     raise NoGeneCategoriesError("mine", "tiny", ["gencode_v50"])
    ... except LookupError as error:
    ...     print("gencode_v50" in str(error))
    True
    """

    def __init__(
        self,
        annotation: str,
        assembly: str,
        curated: Iterable[str],
        *,
        contributors: Iterable[str] = (),
    ) -> None:
        self.annotation = annotation
        self.assembly = assembly
        self.curated: tuple[str, ...] = tuple(curated)
        self.contributors: tuple[str, ...] = tuple(contributors)
        listed = ", ".join(self.curated) or "(none)"
        if self.contributors:
            merged = ", ".join(self.contributors)
            next_step = (
                f"It merges {merged}, and no curated list ships for any of them either — a "
                f"list under the merged name would fix nothing, so ship one for whichever "
                f"of those should answer."
            )
        else:
            next_step = (
                f"Answering for it means shipping a curated list for it, as "
                f"{GENE_LIST_SUBDIR}/{annotation}{GENE_LIST_SUFFIX}."
            )
        super().__init__(
            f"no curated gene list ships for the annotation {annotation!r} of {assembly!r}, "
            f"so no gene category can be answered for it — which is not the same answer as "
            f"its declaring a category with no genes in it. Curated lists ship for: "
            f"{listed}. {next_step}"
        )


class GeneCategoryNotDeclaredError(LookupError):
    """The annotation declares categories, and not the one asked for.

    The second of the two absences. It is a real answer about a real curated list — the
    category was not curated for this annotation, which is not the same as its being
    curated and empty, and no shipped category is ever empty. The message lists the ones
    it does declare, since those are what a caller may ask for instead.

    Parameters
    ----------
    annotation : str
        The **Registered name** that was asked about.
    assembly : str
        The **Assembly** it is registered for.
    category : str
        The category that is not declared.
    declared : iterable of str
        The categories that are, in the order the curated lists spell them.

    Attributes
    ----------
    annotation : str
        The name asked about.
    assembly : str
        The assembly it is registered for.
    category : str
        The category nobody declared.
    declared : tuple of str
        The categories that are declared.

    Examples
    --------
    >>> try:
    ...     raise GeneCategoryNotDeclaredError("mine", "tiny", "tRNA", ["rRNA"])
    ... except LookupError as error:
    ...     print("rRNA" in str(error))
    True
    """

    def __init__(
        self, annotation: str, assembly: str, category: str, declared: Iterable[str]
    ) -> None:
        self.annotation = annotation
        self.assembly = assembly
        self.category = category
        self.declared: tuple[str, ...] = tuple(declared)
        listed = ", ".join(self.declared) or "(none)"
        super().__init__(
            f"the curated gene list for {annotation!r} of {assembly!r} declares no category "
            f"{category!r}. It declares: {listed}. Ask for one of those, or add {category!r} "
            f"to the curated list — a category nobody curated is unanswered here rather than "
            f"empty."
        )


@dataclass(frozen=True)
class CuratedCategory:
    """One **Gene category** a curated list declares, and the genes curated into it.

    Attributes
    ----------
    category : str
        The category's name, as the file spells it — ``"rRNA"``, ``"Mt_rRNA"``. Which
        names exist is a property of the annotation and not of this package.
    description : str
        What membership means for this annotation, since two annotations spelling a
        category the same way need not have curated it the same way.
    source : str
        Where membership came from, and the caveats on using it.
    gene_ids : tuple of str
        The gene ids, in the order the file lists them. Never empty: a declared category
        carries at least one gene, so a caller never has to tell an empty one from an
        absent one.

    Examples
    --------
    >>> category = CuratedCategory(
    ...     category="rRNA",
    ...     description="the mature ribosomal RNA genes",
    ...     source="WormBase WS298 gene_biotype",
    ...     gene_ids=("WBGene00004512",),
    ... )
    >>> category.category, len(category.gene_ids)
    ('rRNA', 1)
    """

    category: str
    description: str
    source: str
    gene_ids: tuple[str, ...]


@dataclass(frozen=True)
class CuratedGeneList:
    """One annotation's **Curated gene list**: every category it declares, and for what.

    What one shipped file says, read back. It carries the **Assembly** it was curated
    against as well as the annotation it is for, because a curated list is pinned to
    both — see :meth:`check_assembly`.

    Attributes
    ----------
    annotation : str
        The **Registered name** this list is for, which is also its file's name.
    assembly : str
        The **Assembly** it was curated against.
    categories : mapping of str to CuratedCategory
        Every category declared, keyed by name and in the order the file spells them.
        Never empty in one read from a shipped file: a file declaring nothing at all is
        refused as a defect, since absence is spelled by shipping no file.

    Examples
    --------
    >>> listed = CuratedGeneList("gencode_v50", "hg38", {})
    >>> listed.annotation, listed.assembly
    ('gencode_v50', 'hg38')
    >>> listed.check_assembly("hg38") is None
    True
    """

    annotation: str
    assembly: str
    categories: Mapping[str, CuratedCategory]

    def check_assembly(self, assembly: str) -> None:
        """Prove this list was curated against ``assembly``, or refuse to answer for it.

        The one guard between a curated list and the wrong reference. An annotation is
        addressed by name, and a name is unique only within its assembly, so a list found
        by name is not yet known to be about the reference asking — and answering anyway
        would hand back another species' gene ids under this one's name.

        Parameters
        ----------
        assembly : str
            The **Assembly** the annotation is registered for. For a contributor to a
            **Merged annotation** this is that contributor's **Component**, since the
            genes it contributed are its own.

        Raises
        ------
        GeneListAssemblyMismatchError
            If this list was curated against a different assembly.

        Examples
        --------
        >>> CuratedGeneList("wormbase_ws298", "ce11", {}).check_assembly("ce11") is None
        True
        """
        if assembly != self.assembly:
            raise GeneListAssemblyMismatchError(
                f"the curated gene list for {self.annotation!r} was curated against "
                f"{self.assembly!r}, and it was asked for {assembly!r}. A curated list names "
                f"the gene ids of one assembly's annotation, so answering from it here would "
                f"hand back another reference's genes. Either {assembly!r} has an annotation "
                f"registered under a name that belongs to {self.assembly!r}, in which case "
                f"register it under a name of its own, or a curated list has to be shipped "
                f"for {assembly!r}."
            )


@cache
def curated_annotations() -> tuple[str, ...]:
    """Return every annotation a **Curated gene list** ships for, sorted.

    What can be asked about at all, and the answer an error names when an annotation
    cannot be. Read once from the package's own resources and cached; the names are the
    **Registered name**s the lists are filed under, which is how they are looked up.

    Returns
    -------
    tuple of str
        The annotation names, sorted. Empty only if the package ships no lists.

    Examples
    --------
    >>> "ensgene_v101" in curated_annotations()
    True
    """
    directory = files("genome").joinpath(GENE_LIST_SUBDIR)
    if not directory.is_dir():
        return ()
    return tuple(
        sorted(
            entry.name[: -len(GENE_LIST_SUFFIX)]
            for entry in directory.iterdir()
            if entry.name.endswith(GENE_LIST_SUFFIX)
        )
    )


@cache
def curated_gene_list(annotation: str) -> CuratedGeneList | None:
    """Return the **Curated gene list** shipped for ``annotation``, or ``None``.

    The raw absence, and the only place ``None`` is an acceptable way to say it: this is
    the layer below the one a caller touches, and everything above turns it into
    :class:`NoGeneCategoriesError` so that *unanswerable* can never be read as *empty*.

    The name is looked up among what :func:`curated_annotations` found rather than joined
    onto the resource directory, so a name shaped like a path finds nothing instead of
    walking out of it. The file is validated as it is read, so a list that cannot be
    trusted raises rather than answering. Read once per name and cached; everything it
    returns is frozen.

    Parameters
    ----------
    annotation : str
        The **Registered name** to look up, e.g. ``"gencode_v50"``.

    Returns
    -------
    CuratedGeneList or None
        The list, or ``None`` when none ships for that name. ``None`` is legal and
        ordinary — an annotation nobody has curated categories for.

    Raises
    ------
    CuratedGeneListError
        If a list ships under that name and cannot be read; the message names the file.

    Examples
    --------
    >>> list(curated_gene_list("refseq_rs_2025_06_26").categories)
    ['rRNA']
    >>> curated_gene_list("no_such_annotation") is None
    True
    """
    if annotation not in curated_annotations():
        return None
    resource = files("genome").joinpath(GENE_LIST_SUBDIR, f"{annotation}{GENE_LIST_SUFFIX}")
    return _read_gene_list(
        resource.read_text(encoding="utf-8"), annotation=annotation, origin=str(resource)
    )


def _read_gene_list(text: str, *, annotation: str, origin: str) -> CuratedGeneList:
    """Read one curated gene list from ``text``, validating it as it goes.

    Separate from the resource it came out of, so every way a file can be wrong is
    reachable without writing a broken one into the package. ``origin`` is where the text
    came from and is named in every message, since fixing that file is the only repair.
    """
    payload = _payload(text, origin=origin)
    for key in _REQUIRED_KEYS:
        if key not in payload:
            raise CuratedGeneListError(
                f"{origin} declares no {key!r}. A curated gene list names the annotation it "
                f"is for, the assembly it was curated against, and the categories it "
                f"declares. Add {key!r} to the file."
            )
    declared = payload["annotation"]
    if declared != annotation:
        raise CuratedGeneListError(
            f"{origin} says it is the curated gene list for {declared!r}, but its file name "
            f"says {annotation!r}. A list is looked up by the name its annotation is "
            f"registered under, so the two must agree: rename the file, or fix the "
            f"'annotation' field."
        )
    categories = payload["categories"]
    if not isinstance(categories, dict) or not categories:
        raise CuratedGeneListError(
            f"{origin} declares no categories, or declares them as something other than an "
            f"object keyed by category name. A curated list that declares nothing says no "
            f"more than an absent file does, and absence is how *no categories at all* is "
            f"spelled — declare at least one category, or remove the file."
        )
    built: dict[str, CuratedCategory] = {}
    seen: dict[str, str] = {}
    for name, entry in categories.items():
        built[str(name)] = _read_category(str(name), entry, origin=origin, seen=seen)
    return CuratedGeneList(
        annotation=annotation,
        assembly=str(payload["assembly"]),
        categories=MappingProxyType(built),
    )


def _payload(text: str, *, origin: str) -> dict[str, Any]:
    """Parse one curated gene list's JSON, or say which file refused and why."""
    try:
        payload = json.loads(text)
    except ValueError as error:
        raise CuratedGeneListError(
            f"{origin} is not readable as JSON: {error}. A curated gene list ships inside "
            f"this package, so this is a defect in what was committed here rather than "
            f"anything a caller did — fix the file."
        ) from error
    if not isinstance(payload, dict):
        raise CuratedGeneListError(
            f"{origin} holds {type(payload).__name__} where a curated gene list is an object "
            f"with 'annotation', 'assembly' and 'categories'. Fix the file."
        )
    return payload


def _read_category(name: str, entry: Any, *, origin: str, seen: dict[str, str]) -> CuratedCategory:
    """Read one declared category, holding it to what a category must be.

    ``seen`` maps each gene id already read to the category that claimed it, and is
    updated here: a gene in two categories is what would make a caller summing both count
    it twice, which is the judgment this package refuses to get wrong for them.
    """
    if not isinstance(entry, dict):
        raise CuratedGeneListError(
            f"{origin} declares the category {name!r} as {type(entry).__name__} where a "
            f"category is an object with 'description', 'source' and 'gene_ids'. Fix the file."
        )
    description = _prose(entry, "description", name=name, origin=origin)
    source = _prose(entry, "source", name=name, origin=origin)
    gene_ids = entry.get("gene_ids")
    if not isinstance(gene_ids, list) or not gene_ids:
        raise CuratedGeneListError(
            f"{origin} declares the category {name!r} with no gene ids. A declared category "
            f"always carries at least one gene, since an empty one reads as *no genes in it* "
            f"where this package says *nobody curated it* — fill it in, or drop the category."
        )
    ids: list[str] = []
    for gene_id in gene_ids:
        if not isinstance(gene_id, str) or not gene_id.strip():
            raise CuratedGeneListError(
                f"{origin} lists {gene_id!r} among the gene ids of {name!r}, which is not a "
                f"gene id. Every entry is the annotation's own gene_id, as text."
            )
        if gene_id in ids:
            raise CuratedGeneListError(
                f"{origin} lists the gene id {gene_id!r} twice in the category {name!r}. The "
                f"ids are returned as they are listed and never de-duplicated, so a caller "
                f"summing them would count it twice — list it once."
            )
        claimed = seen.get(gene_id)
        if claimed is not None:
            raise CuratedGeneListError(
                f"{origin} lists the gene id {gene_id!r} in both {claimed!r} and {name!r}. A "
                f"gene belongs to one category, so a caller summing two of them would count "
                f"it twice — decide which category it is in."
            )
        seen[gene_id] = name
        ids.append(gene_id)
    return CuratedCategory(
        category=name, description=description, source=source, gene_ids=tuple(ids)
    )


def _prose(entry: Mapping[str, Any], key: str, *, name: str, origin: str) -> str:
    """Return one of a category's two sentences, or say which file left it blank."""
    value = entry.get(key)
    if not isinstance(value, str) or not value.strip():
        raise CuratedGeneListError(
            f"{origin} gives the category {name!r} no {key!r}. A category says what "
            f"membership means for this annotation and where membership came from, because "
            f"two annotations spelling a category the same way need not have curated it the "
            f"same way. Fill {key!r} in."
        )
    return value
