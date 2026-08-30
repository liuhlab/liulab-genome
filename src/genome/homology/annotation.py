"""Putting a homology answer into one registered **Annotation**'s own gene ids.

The crossing a caller makes once they have homologs and want to join them to their own
counts matrix, and the whole of it is one existing call:
:meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`, used **unchanged**. Nothing is
added to :class:`~genome.assembly.genome.Genome` for it and no mixin is introduced — a **Homology
set** is usable with no genome open, and a ``Genome``-level convenience would quietly
re-introduce the assembly this design removed.

This module is a function, not an object, because it owns no state: it takes an answer and
a registry, asks the registry the one question it already answers, and reports what the
crossing cost. It reaches no file itself — the registry does. What it answers with,
:class:`ResolvedHomologs`, is defined here for the same reason (ADR-0022), and the answer
it consumes is defined beside the set that builds it, in :mod:`genome.homology.compara`.

**The Homology type is the publisher's and stands** (ADR-0020). An annotation that spells
one gene of an ``ortholog_one2many`` link and not the other leaves a view that looks
one-to-one, and the label still reads ``ortholog_one2many``; what fell away is counted as a
**Dropped partner** rather than folded into the label.

Examples
--------
>>> from genome.homology import resolve_homologs
>>> resolve_homologs                                     # doctest: +ELLIPSIS
<function resolve_homologs at ...>
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

from genome.homology.compara import HomologyAnswer, HomologyLink

if TYPE_CHECKING:  # pragma: no cover - imported for typing only, and the cycle it avoids
    from genome.annotation import AnnotationRegistry


@dataclass(frozen=True)
class ResolvedHomologs:
    """A :class:`~genome.homology.compara.HomologyAnswer` put into one **Annotation**'s gene ids.

    :func:`resolve_homologs`'s answer, and the crossing a caller makes once they have
    homologs and want to join them to their own counts matrix. The hop itself is
    :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`, used unchanged.

    **The Homology type is the publisher's and stands** (ADR-0020). An annotation that
    spells one gene of an ``ortholog_one2many`` link and not the other leaves a view that
    looks one-to-one, and the label still reads ``ortholog_one2many``; what the crossing
    removed is in :attr:`dropped_partners` rather than folded into the label.

    **Both qualifications the answer carried ride through.** A **Dropped partner** counts
    partners lost to either step — a **Homology type** filter before the crossing, or an
    annotation missing the gene during it — so the count a caller reads is what the whole
    path cost rather than what its last step did. And
    :attr:`~genome.homology.compara.HomologyAnswer.null_quality_scores` is a fact about the
    set the answer came from, not about the crossing, so it is repeated here for a caller
    who filters on ``goc_score`` after resolving.

    Attributes
    ----------
    species : str
        The species the stems asked about belong to.
    other_species : str
        The species the homologous genes belong to, and the one the annotation annotates.
    release : str
        The Ensembl Compara **Release** that asserted these links.
    assembly : str
        The **Assembly** whose annotation these gene ids belong to.
    annotation : str
        The **Registered name** whose own gene ids these are.
    resolved : mapping of str to tuple of HomologyLink
        Every asked stem that still names at least one link, in ask order, to the links
        whose partner this annotation carries a gene for. No value is ever empty.
    gene_ids : mapping of str to tuple of str
        Every partner **Gene id stem** that survived, to the gene ids this annotation
        spells it with, ascending. Keyed by partner and not by asked stem, because two
        asked stems may name one partner and its ids are the same ids.
    unresolved : tuple of str
        The asked stems left naming nothing here: first those the crossing emptied — every
        partner missing from this annotation — in ask order, then those the set already
        named no homolog for, in ask order. Two groups rather than one interleaved list,
        because *this annotation is missing every partner* and *this release knows no
        homolog* are different facts about a gene.
    dropped_partners : tuple of str
        The **Dropped partner**s: every partner **Gene id stem** this answer no longer
        names, ascending — those a **Homology type** filter removed before the crossing
        and those this annotation carries no gene for, in one count, since the definition
        covers both and a caller wants what the whole path cost.
    null_quality_scores : tuple of str
        The names of the confidence fields the **Homology set** behind this holds no value
        in, carried through unchanged. The crossing neither adds a score nor removes one.

    Examples
    --------
    >>> from genome.homology.compara import HomologyLink
    >>> link = HomologyLink(
    ...     "ENSG00000141510", "ENSMUSG00000059552", "ortholog_one2many", True, 100, 96.79
    ... )
    >>> crossed = ResolvedHomologs(
    ...     species="Homo sapiens",
    ...     other_species="Mus musculus",
    ...     release="116",
    ...     assembly="mm39",
    ...     annotation="gencode_vM39",
    ...     resolved={"ENSG00000141510": (link,)},
    ...     gene_ids={"ENSMUSG00000059552": ("ENSMUSG00000059552.5",)},
    ...     unresolved=(),
    ...     dropped_partners=("ENSMUSG00000000001",),
    ...     null_quality_scores=(),
    ... )
    >>> crossed.homolog_gene_ids
    ['ENSMUSG00000059552.5']
    >>> crossed.resolved["ENSG00000141510"][0].homology_type
    'ortholog_one2many'
    """

    species: str
    other_species: str
    release: str
    assembly: str
    annotation: str
    resolved: Mapping[str, tuple[HomologyLink, ...]]
    gene_ids: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]
    dropped_partners: tuple[str, ...]
    null_quality_scores: tuple[str, ...]

    @property
    def homolog_gene_ids(self) -> list[str]:
        """Every homologous gene id named, partner order then id order — a fresh list.

        **Flattening loses what the mapping carries**: which asked stem reached the gene,
        and under what **Homology type**. It keeps every id rather than one per partner,
        since one stem may be spelled by two gene ids — the pseudoautosomal case
        :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids` answers with both of.
        """
        return [gene_id for ids in self.gene_ids.values() for gene_id in ids]

    def as_json(self) -> dict[str, Any]:
        """Return this crossing as ``--json`` serializes it.

        Returns
        -------
        dict
            The species pair, the ``release``, the ``assembly`` and ``annotation``,
            ``resolved`` as a mapping of stem to a list of
            :meth:`~genome.homology.compara.HomologyLink.as_json` links, ``gene_ids`` as a
            plain mapping, ``unresolved``, ``dropped_partners`` and ``null_quality_scores``
            as lists, and the flattened ``homolog_gene_ids``.
        """
        return {
            "species": self.species,
            "other_species": self.other_species,
            "release": self.release,
            "assembly": self.assembly,
            "annotation": self.annotation,
            "resolved": {
                stem: [link.as_json() for link in links] for stem, links in self.resolved.items()
            },
            "gene_ids": {stem: list(ids) for stem, ids in self.gene_ids.items()},
            "unresolved": list(self.unresolved),
            "dropped_partners": list(self.dropped_partners),
            "null_quality_scores": list(self.null_quality_scores),
            "homolog_gene_ids": self.homolog_gene_ids,
        }


def resolve_homologs(
    answer: HomologyAnswer, registry: AnnotationRegistry, name: str | None = None
) -> ResolvedHomologs:
    """Return ``answer``'s homologous genes in one registered annotation's own gene ids.

    The other species' **Gene id stem**s are resolved in one call against ``registry``,
    which is the annotation hop this package already had; every link whose partner that
    annotation carries a gene for is kept, with its **Homology type** untouched, and every
    partner it carries none for is reported in :attr:`~ResolvedHomologs.dropped_partners`
    rather than dropped in silence — *added to* whatever the answer had already dropped,
    since a **Dropped partner** is one the answer no longer names whichever step removed
    it. Which quality columns the set holds nothing in is a fact about the set rather than
    about the crossing, and rides through unchanged.

    **The registry must annotate the answer's other species.** Nothing here checks that —
    an assembly's species is the assembly's own metadata and a registry does not carry
    one — so passing a mouse answer a worm registry is a question about the wrong genome
    and will simply resolve nothing. Reach the registry from the assembly whose species is
    :attr:`~HomologyAnswer.other_species`.

    Parameters
    ----------
    answer : HomologyAnswer
        What :meth:`~genome.homology.compara.HomologySet.homologs` returned. Whatever it
        was filtered to is what is crossed: this adds nothing back and removes no more.
    registry : genome.annotation.registry.AnnotationRegistry
        The registry of the **Assembly** whose annotation the ids should be in.
    name : str, optional
        The **Registered name** to resolve against. Omitted, that assembly's **Default
        annotation** answers.

    Returns
    -------
    ResolvedHomologs
        The links whose partners this annotation spells, the gene ids it spells them with,
        the asked stems left naming nothing, every partner dropped along the way, and the
        quality columns the set behind it holds nothing in.

    Raises
    ------
    ValueError
        If ``name`` is omitted and no **Default annotation** is decided.
    genome.annotation.registry.AnnotationNotRegisteredError
        If nothing of that name is registered there.
    genome.annotation.stems.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.homology import HomologySet, resolve_homologs
    >>> from genome.annotation import AnnotationRegistry
    >>> answer = HomologySet("Homo sapiens", "Mus musculus").homologs(   # doctest: +SKIP
    ...     ["ENSG00000141510"]
    ... )
    >>> crossed = resolve_homologs(                                      # doctest: +SKIP
    ...     answer, AnnotationRegistry.locate("mm39"), "gencode_vM39"
    ... )
    >>> crossed.homolog_gene_ids                                         # doctest: +SKIP
    ['ENSMUSG00000059552.13']
    """
    partners = tuple(
        dict.fromkeys(
            link.homolog_gene_id_stem for links in answer.resolved.values() for link in links
        )
    )
    crossed = registry.resolve_gene_ids(partners, name)
    kept: dict[str, tuple[HomologyLink, ...]] = {}
    for stem, links in answer.resolved.items():
        surviving = tuple(link for link in links if link.homolog_gene_id_stem in crossed.resolved)
        if surviving:
            kept[stem] = surviving
    named = {link.homolog_gene_id_stem for links in kept.values() for link in links}
    return ResolvedHomologs(
        species=answer.species,
        other_species=answer.other_species,
        release=answer.release,
        assembly=crossed.assembly,
        annotation=crossed.annotation,
        resolved=kept,
        gene_ids={stem: ids for stem, ids in crossed.resolved.items() if stem in named},
        # The stems the crossing emptied first, in the order they were asked about, then
        # the ones the set already named no homolog for, in theirs. Two groups rather than
        # one interleaved list because an answer keeps them apart, exactly as
        # `ResolvedGeneIds` does — and because *this annotation is missing every partner*
        # and *this release knows no homolog* are two different facts about a gene.
        unresolved=tuple(stem for stem in answer.resolved if stem not in kept) + answer.unresolved,
        # Both causes in one count. A **Dropped partner** is a partner the answer no
        # longer names, and the definition does not care which step removed it — so what
        # a Homology type filter dropped before the crossing is added to what the
        # annotation dropped during it, rather than being replaced by it. ADR-0020 turns
        # on the count being reported, not on the label being corrected.
        dropped_partners=tuple(sorted(set(answer.dropped_partners) | set(crossed.unresolved))),
        # A fact about the set the answer came from, not about the crossing, and the
        # module documentation says it rides on every answer. A crossing is an answer.
        null_quality_scores=answer.null_quality_scores,
    )
