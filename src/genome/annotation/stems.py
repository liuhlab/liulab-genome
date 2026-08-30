"""Resolving a **Gene id stem** against an **Annotation**'s own gene ids.

The seam :mod:`genome.xref`, :mod:`genome.homology` and :mod:`genome.tf` all cross, and
the only one they do. A stem is a gene id with its version dropped —
``ENSG00000123456`` for ``ENSG00000123456.7`` — which is how every published table keyed by
gene arrives, and never how a GENCODE **Annotation** spells the same gene.
:meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids` is the crossing
and this is what it is made of: every gene id in the **Annotation database** is reduced to
its own stem, and a stem answers with every gene id that reduced to it.

**Every id, and never a chosen one.** One stem naming two gene ids is not a malformed
annotation: ``gencode_v50lift37`` has nine such stems, eight of them pseudoautosomal genes
carrying a ``_PAR_Y`` copy, and a resolver taking the first would hand back the X copy of a
Y gene without saying it had chosen. So :class:`ResolvedGeneIds` maps a stem to *all* of
them, and the stems that named nothing ride back on it rather than being dropped.

**The pass stays streaming.** :func:`gene_ids_by_stem` walks
:func:`~genome.annotation.database.gene_ids` a row at a time and keeps only what
matched, so a GENCODE-sized annotation costs an index walk of its gene rows and is never
held in memory.

Nothing here knows what the stems it is handed are a list *of*. Which species selects a
shipped file, and what a row of one says, are facts this module has no stake in, so it
holds no import of any of them and gains nothing when a fourth topic arrives.

Examples
--------
>>> from genome.annotation.stems import ResolvedGeneIds
>>> ResolvedGeneIds("hg38", "gencode_v50", {"A": ("A.1", "A.2")}, ("B",)).gene_ids
['A.1', 'A.2']
"""

from __future__ import annotations

from collections.abc import Container, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genome.annotation.database import gene_ids

#: What separates a gene id from its version — ``ENSG00000123456.7`` — and therefore what
#: a **Gene id stem** is everything before. An id carrying none is its own stem.
_VERSION_SEPARATOR = "."


class NoGeneFeaturesError(LookupError):
    """An annotation's database holds no gene at all, so no gene id can be resolved.

    The absence a caller must never read as *this annotation carries none of my genes*.
    A GTF that declares only exons registers as exons alone —
    :meth:`~genome.annotation.registry.AnnotationRegistry.register` leaves **Feature
    inference** off, and rightly, since it is the library's slow path and the publishers
    who matter declare their genes — so an
    annotation like that would answer every stem with *not found* while looking perfectly
    healthy. It says so instead, and names the argument that rebuilds it with the genes in.

    A :class:`LookupError`, as the other absences on this surface are.

    Parameters
    ----------
    annotation : str
        The **Registered name** that was asked about.
    assembly : str
        The **Assembly** it is registered for.

    Attributes
    ----------
    annotation : str
        The name asked about.
    assembly : str
        The assembly it is registered for.

    Examples
    --------
    >>> try:
    ...     raise NoGeneFeaturesError("mine", "tiny")
    ... except LookupError as error:
    ...     print("--infer-genes" in str(error))
    True
    """

    def __init__(self, annotation: str, assembly: str) -> None:
        self.annotation = annotation
        self.assembly = assembly
        super().__init__(
            f"the annotation {annotation!r} registered for {assembly!r} has no gene features "
            f"in its database, so there is nothing to resolve gene ids against. This is not "
            f"a gene it happens to lack: its GTF declares no gene lines at all, and "
            f"reconstructing them from the exons is off by default. Register it again with "
            f"that turned on — --infer-genes from a shell, disable_infer_genes=False from "
            f"Python — or register an annotation whose GTF declares its genes."
        )


@dataclass(frozen=True)
class ResolvedGeneIds:
    """The gene ids one **Annotation** carries for the **Gene id stem**s it was asked about.

    :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`'s answer,
    and the one result type two
    contexts share: it is defined beside the call that returns it, and
    :mod:`genome.tf.gene`, :mod:`genome.tf.cofactor` and :mod:`genome.homology.annotation`
    import it from here (ADR-0022). A stem is a gene id with its version dropped, and
    inside one annotation it may name more than one gene id — nine do in
    ``gencode_v50lift37``, eight of them pseudoautosomal-Y — so every stem answers with
    **all** of them and nothing here picks one. Two stems never name the same gene id,
    since an id has exactly one stem.

    **What was asked about and is not there rides back on the answer.** A caller holding a
    few thousand stems gets the ones this annotation carries no gene for in
    :attr:`unresolved` rather than a shorter list than it passed, so what the thing it was
    holding contains and this annotation does not is visible instead of dropped.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    resolved : mapping of str to tuple of str
        Every stem that named at least one gene id, in the order the stems were asked
        about, to the ids it names, in ascending order. No value is ever an empty tuple —
        a stem that named nothing is in :attr:`unresolved` instead.
    unresolved : tuple of str
        The stems no gene id in the annotation is of, in the order they were asked about.

    Examples
    --------
    >>> answer = ResolvedGeneIds(
    ...     assembly="hg19",
    ...     annotation="gencode_v50lift37",
    ...     resolved={
    ...         "ENSG00000182378": ("ENSG00000182378.14", "ENSG00000182378.14_PAR_Y"),
    ...         "ENSG00000141510": ("ENSG00000141510.18",),
    ...     },
    ...     unresolved=("ENSG00000288541",),
    ... )
    >>> answer.gene_ids
    ['ENSG00000182378.14', 'ENSG00000182378.14_PAR_Y', 'ENSG00000141510.18']
    >>> answer.as_json()["unresolved"]
    ['ENSG00000288541']
    """

    assembly: str
    annotation: str
    resolved: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id resolved, stem order and then id order — a fresh list each call.

        **Every** id, not one per stem. Flattening is exactly where a reader would take
        the first id of each stem and lose the second, which is the pseudoautosomal gene
        this answer's shape exists to keep; :attr:`resolved` is what says which stem an id
        came from.
        """
        return [gene_id for ids in self.resolved.values() for gene_id in ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``resolved`` as a mapping of stem to a list of
            gene ids, ``unresolved`` as a list, and the flattened ``gene_ids``. The last
            is written out beside the mapping it is read from for the reason
            :attr:`gene_ids` gives: a reader assembling it is a reader who might take one
            id per stem.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "resolved": {stem: list(ids) for stem, ids in self.resolved.items()},
            "unresolved": list(self.unresolved),
            "gene_ids": self.gene_ids,
        }


def _gene_id_stem(gene_id: str) -> str:
    """Return ``gene_id`` with its version dropped — its **Gene id stem**.

    Everything before the first separator, and the whole id when it carries none, which is
    what makes an unversioned annotation's gene its own stem. GENCODE's pseudoautosomal
    ``ENSG00000182378.14_PAR_Y`` stems to ``ENSG00000182378`` alongside the copy it is of,
    which is the collision the caller is handed both halves of.
    """
    stem, separator, _version = gene_id.partition(_VERSION_SEPARATOR)
    return stem if separator else gene_id


def gene_ids_by_stem(database: Path, wanted: Container[str]) -> tuple[dict[str, list[str]], bool]:
    """Return the gene ids in ``database`` under each wanted stem, and whether it has any.

    One pass, matching stems as it walks :func:`~genome.annotation.database.gene_ids` a
    row at a time: the database is queried rather than read, so a GENCODE-sized annotation
    costs an index walk of its gene rows and only what matched is ever held. The ids under
    a stem arrive ascending because the walk does, so two machines answer in one order.

    The second half of the answer is *whether the database holds a gene at all*, because
    **no genes** and **no matching genes** are different facts that the mapping alone
    cannot tell apart, and the caller says something different about each.

    Parameters
    ----------
    database : pathlib.Path
        The **Annotation database** to walk.
    wanted : container of str
        The **Gene id stem**s to keep. Membership is all that is asked of it, so a caller
        with a few thousand passes a set and pays nothing per gene row for it.

    Returns
    -------
    dict of str to list of str
        Each wanted stem that named at least one gene id, to those ids, ascending.
    bool
        Whether the database holds any gene feature at all.

    Examples
    --------
    >>> from pathlib import Path
    >>> gene_ids_by_stem(Path("gencode_v50.db"), {"ENSG00000141510"})   # doctest: +SKIP
    ({'ENSG00000141510': ['ENSG00000141510.18']}, True)
    """
    found: dict[str, list[str]] = {}
    any_genes = False
    for gene_id in gene_ids(database):
        any_genes = True
        stem = _gene_id_stem(gene_id)
        if stem in wanted:
            found.setdefault(stem, []).append(gene_id)
    return found, any_genes
