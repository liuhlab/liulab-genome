"""Putting a homology answer into one registered **Annotation**'s own gene ids.

The crossing a caller makes once they have homologs and want to join them to their own
counts matrix, and the whole of it is one existing call:
:meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`, used **unchanged**. Nothing is
added to :class:`~genome.genome.Genome` for it and no mixin is introduced — a **Homology
set** is usable with no genome open, and a ``Genome``-level convenience would quietly
re-introduce the assembly this design removed.

This module is a function, not an object, because it owns no state: it takes an answer and
a registry, asks the registry the one question it already answers, and reports what the
crossing cost. It reaches no file itself — the registry does.

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

from typing import TYPE_CHECKING

from genome.io.results import HomologyAnswer, HomologyLink, ResolvedHomologs

if TYPE_CHECKING:  # pragma: no cover - imported for typing only, and the cycle it avoids
    from genome.io.gtf import AnnotationRegistry


def resolve_homologs(
    answer: HomologyAnswer, registry: AnnotationRegistry, name: str | None = None
) -> ResolvedHomologs:
    """Return ``answer``'s homologous genes in one registered annotation's own gene ids.

    The other species' **Gene id stem**s are resolved in one call against ``registry``,
    which is the annotation hop this package already had; every link whose partner that
    annotation carries a gene for is kept, with its **Homology type** untouched, and every
    partner it carries none for is reported in
    :attr:`~genome.io.results.ResolvedHomologs.dropped_partners` rather than dropped in
    silence.

    **The registry must annotate the answer's other species.** Nothing here checks that —
    an assembly's species is the assembly's own metadata and a registry does not carry
    one — so passing a mouse answer a worm registry is a question about the wrong genome
    and will simply resolve nothing. Reach the registry from the assembly whose species is
    :attr:`~genome.io.results.HomologyAnswer.other_species`.

    Parameters
    ----------
    answer : genome.io.results.HomologyAnswer
        What :meth:`~genome.homology.compara.HomologySet.homologs` returned. Whatever it
        was filtered to is what is crossed: this adds nothing back and removes no more.
    registry : genome.io.gtf.AnnotationRegistry
        The registry of the **Assembly** whose annotation the ids should be in.
    name : str, optional
        The **Registered name** to resolve against. Omitted, that assembly's **Default
        annotation** answers.

    Returns
    -------
    genome.io.results.ResolvedHomologs
        The links whose partners this annotation spells, the gene ids it spells them with,
        the asked stems left naming nothing, and the partners it carries no gene for.

    Raises
    ------
    ValueError
        If ``name`` is omitted and no **Default annotation** is decided.
    genome.io.gtf.AnnotationNotRegisteredError
        If nothing of that name is registered there.
    genome.io.gtf.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.homology import HomologySet, resolve_homologs
    >>> from genome.io.gtf import AnnotationRegistry
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
        dropped_partners=tuple(sorted(crossed.unresolved)),
    )
