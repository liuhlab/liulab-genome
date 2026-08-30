"""Where a shipped **Cofactor table** meets one registered **Annotation**'s gene ids.

:mod:`genome.tf.gene.annotation` for the cofactor half, in the same shape and against the
same one call — :meth:`~genome.io.annotation.registry.AnnotationRegistry.resolve_gene_ids`, used
**unchanged**. It is a second caller of one resolver and not a second crossing: the two
halves differ in which shipped file is read and in what a row of it says, and in nothing
else, so a caller who has read one answer has read both.

**Membership is this package's and classification is each publisher's.** A table built
from two publishers is a union nobody else published (ADR-0016), which is why each entry
says who listed the gene and keeps every publisher's own vocabulary under that publisher's
namespaced column name.

**The absences are the census half's pair, and they do not fire for the same assemblies.**
The species is the assembly's own here too, and worm is answered here while
:func:`~genome.tf.gene.annotation.resolve_tf_genes` raises for it — a publisher assessed
worm cofactors and none has released a worm TF census.

Examples
--------
>>> from genome.tf.cofactor import tf_cofactor_list
>>> tf_cofactor_list("mm39").species                     # doctest: +SKIP
'Mus musculus'
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from genome.io.annotation import AnnotationRegistry, ResolvedGeneIds
from genome.metadata import assembly_metadata, species_slug
from genome.tf.cofactor.table import (
    TRUE_CELL,
    CofactorProvenance,
    CofactorTable,
    cofactor_metadata,
    cofactor_species,
    cofactor_table,
)
from genome.tf.cofactor.table import UNIFORM_COLUMNS as COFACTOR_UNIFORM_COLUMNS
from genome.tf.species import NoCofactorTableError, UnknownSpeciesError


@dataclass(frozen=True)
class TFCofactor:
    """One gene a publisher lists as a cofactor, named in one **Annotation**'s gene ids.

    An entry of a :class:`TFCofactorList`, and the counterpart of :class:`TFGene`. The
    table's four uniform columns are fields of their own — the **Gene id stem** it is
    keyed by, the symbol, the cofactor flag and which publisher listed the gene — and
    everything each publisher classified it with stays under that publisher's own
    namespaced name in :attr:`classifications`, because no two publishers carry the same
    columns and nothing here compares one's vocabulary with another's (ADR-0014).

    ``gene_ids`` is a tuple for the reason :class:`TFGene`'s is: one stem may name more
    than one gene id in one annotation and this never picks one. It is never empty — a
    stem the annotation carries no gene for is in :attr:`TFCofactorList.unresolved`
    instead of here.

    Attributes
    ----------
    gene_id_stem : str
        The **Gene id stem** the **Cofactor table** is keyed by.
    gene_ids : tuple of str
        Every gene id this annotation spells that stem with, ascending.
    symbol : str or None
        The gene symbol the table records, or ``None`` where it records none.
    is_cofactor : bool
        The table's own cofactor flag. ``True`` for every entry of a **TF cofactor
        list**, since no publisher shipping here releases a rejected set — a source that
        did would ship rejected rows, and they would be excluded rather than arriving
        here saying ``False``.
    source : str
        Which publisher listed the gene, in the table's own closed vocabulary —
        ``animaltfdb``, ``epifactors``, or ``both`` for a gene two of them listed. ``both``
        is agreement on **membership only** and never on how either classified it.
    classifications : mapping of str to (str or None)
        Every other column the table records for this gene, under the publisher's own
        namespaced snake_case name: ``animaltfdb_family`` and ``animaltfdb_category``
        for a gene AnimalTFDB listed, and the EpiFactors function, target, modification
        and complex name for one EpiFactors did. ``None`` is a cell that publisher left
        blank — which, for a gene the other publisher listed, is that publisher saying
        nothing rather than a value being lost.

    Examples
    --------
    >>> cofactor = TFCofactor(
    ...     gene_id_stem="ENSMUSG00000000085",
    ...     gene_ids=("ENSMUSG00000000085.16",),
    ...     symbol="Scmh1",
    ...     is_cofactor=True,
    ...     source="animaltfdb",
    ...     classifications={"animaltfdb_category": "Other Cofactors"},
    ... )
    >>> cofactor.classifications["animaltfdb_category"]
    'Other Cofactors'
    >>> cofactor.as_json()["gene_ids"]
    ['ENSMUSG00000000085.16']
    """

    gene_id_stem: str
    gene_ids: tuple[str, ...]
    symbol: str | None
    is_cofactor: bool
    source: str
    classifications: Mapping[str, str | None]

    def as_json(self) -> dict[str, Any]:
        """Return this cofactor as ``--json`` serializes it.

        Returns
        -------
        dict
            The fields above under their own names, with ``gene_ids`` as a list and
            ``classifications`` as a plain mapping under the publishers' own column names.
        """
        return {
            "gene_id_stem": self.gene_id_stem,
            "gene_ids": list(self.gene_ids),
            "symbol": self.symbol,
            "is_cofactor": self.is_cofactor,
            "source": self.source,
            "classifications": dict(self.classifications),
        }


@dataclass(frozen=True)
class TFCofactorList:
    """One **Assembly**'s **Transcription cofactor**s, in its annotation's own gene ids.

    :func:`resolve_tf_cofactors`'s answer, and the counterpart of
    :class:`~genome.tf.gene.annotation.TFGeneList` in the same shape: the **Cofactor
    table**'s **Gene id stem**s resolved against one registered annotation, so the ids
    join to a counts matrix with nothing left to normalise.

    **Membership is this package's and classification is each publisher's.** A table built
    from two publishers is a union nobody else published, which is why :attr:`provenance`
    carries a record per publisher rather than one, and why a ``source`` of ``both`` on an
    entry says the two agreed the gene is a cofactor and nothing about how either
    classified it (ADR-0016).

    **What the table holds and this annotation does not is visible.** A stem no gene id
    here is of comes back in :attr:`unresolved` rather than being dropped.

    There is no empty one, for the reasons an absent table would give: an assembly whose
    species has no cofactor table, and one nothing names a species for, each raise a
    :class:`LookupError` of their own.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    species : str
        The species the assembly's own metadata row names, which is what selected the
        table. Never passed in by a caller, so asking for one species' cofactors while
        holding another species' assembly is not expressible (ADR-0003).
    provenance : genome.tf.cofactor.table.CofactorProvenance
        Where the table came from: one record per publisher that contributed to it, plus
        the digest of the shipped bytes.
        :meth:`~genome.tf.cofactor.table.CofactorProvenance.attribution` renders the line
        to print beside anything it answered.
    cofactors : tuple of TFCofactor
        One entry per **Gene id stem** that named at least one gene id here, in the
        table's own row order.
    unresolved : tuple of str
        The stems this annotation carries no gene for, in table row order.

    Examples
    --------
    >>> from genome.tf.cofactor import cofactor_table
    >>> answer = TFCofactorList(
    ...     assembly="mm39",
    ...     annotation="gencode_vM39",
    ...     species="Mus musculus",
    ...     provenance=cofactor_table("Mus musculus").provenance,
    ...     cofactors=(
    ...         TFCofactor(
    ...             "ENSMUSG00000000085",
    ...             ("ENSMUSG00000000085.16",),
    ...             "Scmh1",
    ...             True,
    ...             "animaltfdb",
    ...             {},
    ...         ),
    ...     ),
    ...     unresolved=("ENSMUSG00000000275",),
    ... )
    >>> answer.gene_ids
    ['ENSMUSG00000000085.16']
    >>> answer.provenance.sources[0].publisher
    'AnimalTFDB'
    >>> answer.as_json()["unresolved"]
    ['ENSMUSG00000000275']
    """

    assembly: str
    annotation: str
    species: str
    provenance: CofactorProvenance
    cofactors: tuple[TFCofactor, ...]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id, cofactor order then id order — a fresh list each call.

        **Every** id, not one per gene, for the reason
        :attr:`~genome.tf.gene.annotation.TFGeneList.gene_ids` gives: flattening is
        where a reader would take the first id of a stem that names two and lose the
        other. :attr:`cofactors` is what says which gene an id came from, and what the
        publisher said about it.
        """
        return [gene_id for cofactor in self.cofactors for gene_id in cofactor.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``species``, the table's ``provenance`` under
            its own field names with one entry per publisher under ``sources``,
            ``cofactors`` as a list of :meth:`TFCofactor.as_json` entries, the flattened
            ``gene_ids``, and ``unresolved`` as a list — the keys
            :class:`~genome.tf.gene.annotation.TFGeneList` uses, with the entries named
            for what they are.
        """
        provenance = asdict(self.provenance)
        # ``asdict`` leaves a tuple field a tuple and JSON has no tuple, so the ragged
        # per-publisher records — the one such field here — are written out as the list
        # they serialize to. A payload that did not survive its own round trip would be
        # one whose shape depended on whether anybody had serialized it yet.
        provenance["sources"] = list(provenance["sources"])
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "species": self.species,
            "provenance": provenance,
            "cofactors": [cofactor.as_json() for cofactor in self.cofactors],
            "gene_ids": self.gene_ids,
            "unresolved": list(self.unresolved),
        }


def resolve_tf_cofactors(registry: AnnotationRegistry, name: str | None = None) -> TFCofactorList:
    """Return the genes a publisher lists as cofactors, in one annotation's own gene ids.

    The **Cofactor table** :mod:`genome.tf.cofactor` ships for ``registry``'s assembly's
    species, met with one registered annotation, exactly as
    :func:`~genome.tf.gene.annotation.resolve_tf_genes` meets a census: every **Gene id
    stem** the table is keyed by is resolved through
    :meth:`~genome.io.annotation.registry.AnnotationRegistry.resolve_gene_ids` into the gene ids that
    annotation actually spells, so the answer joins to a counts matrix with nothing left
    for the caller to normalise. A stem naming two gene ids answers with both, and the
    stems the annotation carries no gene for ride back on
    :attr:`TFCofactorList.unresolved` rather than being dropped.

    **The species is the assembly's own** — the curated metadata table's, read here and
    never passed in, so asking for human cofactors while holding a mouse assembly is not
    expressible (ADR-0003).

    **Membership is this package's and classification is each publisher's.** A table built
    from two publishers is a union nobody else published (ADR-0016), which is why each
    entry says which publisher listed the gene in :attr:`TFCofactor.source` and keeps
    every publisher's own vocabulary under that publisher's namespaced column name in
    :attr:`TFCofactor.classifications`. A ``source`` of ``both`` asserts agreement on
    membership only, never on classification.

    There is no widening flag here, as there is on
    :func:`~genome.tf.gene.annotation.resolve_tf_genes`: no publisher shipping today
    releases a rejected set, so the table's listed genes are the whole table. A source
    that did record rejections would ship them, and they would be left out here rather
    than needing a second argument.

    Parameters
    ----------
    registry : genome.io.annotation.registry.AnnotationRegistry
        The registry of the **Assembly** whose species selects the table and whose
        annotation the ids should be in. :class:`~genome.genome.Genome` holds one already.
    name : str, optional
        The **Registered name** to answer in the gene ids of. Omitted, that assembly's
        **Default annotation** answers.

    Returns
    -------
    TFCofactorList
        The cofactors, the publishers' provenance, and the stems that resolved to nothing.

    Raises
    ------
    genome.tf.species.UnknownSpeciesError
        If nothing names that assembly's species — a chimera, or an assembly no curated
        row lists.
    genome.tf.species.NoCofactorTableError
        If no cofactor table ships for that species; the message names the ones that do.
        Worm has a table although no TF census covers it, so this and
        :class:`~genome.tf.species.NoTFCensusError` do not raise for the same set of
        assemblies.
    ValueError
        If ``name`` is omitted and no **Default annotation** is decided.
    genome.io.annotation.registry.AnnotationNotRegisteredError
        If nothing of that name is registered there.
    genome.io.annotation.stems.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.io.annotation import AnnotationRegistry
    >>> from genome.tf.cofactor import resolve_tf_cofactors
    >>> registry = AnnotationRegistry.locate("mm39")                    # doctest: +SKIP
    >>> answer = resolve_tf_cofactors(registry, "gencode_vM39")         # doctest: +SKIP
    >>> first = answer.cofactors[0]                                     # doctest: +SKIP
    >>> first.symbol, first.classifications["animaltfdb_category"]      # doctest: +SKIP
    ('Scmh1', 'Other Cofactors')
    >>> answer.provenance.sources[0].publisher                          # doctest: +SKIP
    'AnimalTFDB'
    """
    species = assembly_metadata(registry.assembly).species
    if species is None:
        raise UnknownSpeciesError(
            registry.assembly, _cofactor_species(), shipped_table="cofactor table"
        )
    table = cofactor_table(species)
    if table is None:
        raise NoCofactorTableError(registry.assembly, species, _cofactor_species())
    resolved = registry.resolve_gene_ids(table.cofactor_stems, name)
    return TFCofactorList(
        assembly=registry.assembly,
        annotation=resolved.annotation,
        species=species,
        provenance=table.provenance,
        cofactors=_tf_cofactors(table, resolved),
        unresolved=resolved.unresolved,
    )


def tf_cofactor_list(
    assembly: str,
    *,
    annotation: str | None = None,
    cache_dir: str | Path | None = None,
) -> TFCofactorList:
    """Return the genes a publisher lists as transcription cofactors in ``assembly``.

    :func:`resolve_tf_cofactors` for an assembly named rather than opened, built the way
    :func:`~genome.tf.gene.annotation.tf_gene_list` is: a registry for the length of the
    call, so a shell surface over it adds no second code path. Nothing is prepared,
    fetched or built to answer it, and the table is read from inside the package.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"mm39"``. Its own metadata row names the
        species, which is what selects the table.
    annotation : str, optional
        The **Registered name** to answer in the gene ids of; the **Default annotation**
        when omitted.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected, as
        :func:`~genome.tf.gene.annotation.tf_gene_list` takes it.

    Returns
    -------
    TFCofactorList
        The answer :func:`resolve_tf_cofactors` describes.

    Raises
    ------
    genome.tf.species.UnknownSpeciesError
        If nothing names the assembly's species.
    genome.tf.species.NoCofactorTableError
        If no cofactor table ships for that species.
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    genome.io.annotation.registry.AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.io.annotation.stems.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.tf.cofactor import tf_cofactor_list
    >>> tf_cofactor_list("mm39").species                     # doctest: +SKIP
    'Mus musculus'
    """
    return resolve_tf_cofactors(AnnotationRegistry.locate(assembly, cache_dir), annotation)


def _cofactor_species() -> tuple[str, ...]:
    """Return the species a **Cofactor table** ships for, in an assembly row's spelling.

    :func:`_censused_species` for the cofactor half, and answering the same question for
    the same reason: it is what both absences name as the thing a caller can ask about
    instead. The shipped files are what is enumerated, since one is what makes a species
    answerable; the provenance table beside them is only read for the publisher's own
    spelling, so a table shipping without a row is still named — badly, as its slug,
    which is the state :func:`~genome.tf.cofactor.table.cofactor_table` raises over
    anyway.
    """
    named = {species_slug(record.species): record.species for record in cofactor_metadata()}
    return tuple(named.get(slug, slug) for slug in cofactor_species())


def _tf_cofactors(table: CofactorTable, resolved: ResolvedGeneIds) -> tuple[TFCofactor, ...]:
    """Return one entry per resolved stem, carrying the table's own row for that gene.

    :func:`~genome.tf.gene.annotation._tf_genes` for the cofactor half. The table's four
    uniform columns become fields and everything after them stays under the publisher's
    own namespaced name, because no two publishers carry the same columns and nothing here
    compares one's vocabulary with another's (ADR-0014). The order is the table's own row
    order, which is the order the stems were asked about.
    """
    publisher_columns = table.columns[len(COFACTOR_UNIFORM_COLUMNS) :]
    rows = {row[0]: row for row in table.rows}
    cofactors: list[TFCofactor] = []
    for stem, gene_ids in resolved.resolved.items():
        cells = dict(zip(table.columns, rows[stem], strict=True))
        cofactors.append(
            TFCofactor(
                gene_id_stem=stem,
                gene_ids=gene_ids,
                symbol=cells["symbol"],
                is_cofactor=cells["is_cofactor"] == TRUE_CELL,
                # ``source`` is one of a closed vocabulary, checked as the file is read,
                # so ``or ""`` only narrows the type of a cell that is always text.
                source=cells["source"] or "",
                classifications={name: cells[name] for name in publisher_columns},
            )
        )
    return tuple(cofactors)
