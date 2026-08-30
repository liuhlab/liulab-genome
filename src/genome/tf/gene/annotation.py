"""Where a published census meets one registered **Annotation**'s own gene ids.

The crossing a caller makes once they have a census and want to join it to their own
counts matrix, and the whole of the crossing is one existing call:
:meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`, used **unchanged**. This
module holds the half the registry has no stake in — which species selects the file, which
file was read, what a row of it says, and what the two absences are — so the Annotation
context answers *which gene ids does this stem name here* and nothing about transcription
factors.

**The species is the assembly's own**, read from its curated metadata row and never passed
in, so asking for human transcription factors while holding a mouse assembly is not
expressible (ADR-0003). It is the only thing here that asks an assembly a question, which
is why :mod:`genome.tf.species` owns the refusals rather than :mod:`genome.annotation`.

**Nothing here decides what a transcription factor is.** The verdict is the census's and
travels with it: a **TF gene** carries the publisher's own **DBD family** and every
**TF assessment** it recorded, and the answer carries the provenance that says whose
judgement it all is.

Two ways in and one code path. :func:`resolve_tf_genes` takes a registry, which is what
:meth:`~genome.assembly.genome.Genome.tf_gene_list` already holds; :func:`tf_gene_list` takes an
assembly name and builds one for the length of the call, which is what ``genome
tf gene-list`` runs. Nothing is prepared, fetched or built to answer either, and the
census is read from inside the package.

Examples
--------
>>> from genome.tf.gene import tf_gene_list
>>> tf_gene_list("hg38").provenance.publisher            # doctest: +SKIP
'Lambert et al. 2018'
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from genome.annotation import AnnotationRegistry, ResolvedGeneIds
from genome.assembly.metadata import assembly_metadata, species_slug
from genome.tf.gene.census import (
    TRUE_CELL,
    UNIFORM_COLUMNS,
    CensusProvenance,
    TFGeneTable,
    census_metadata,
    census_species,
    tf_gene_table,
)
from genome.tf.species import NoTFCensusError, UnknownSpeciesError


@dataclass(frozen=True)
class TFGene:
    """One gene a census assessed, named in one **Annotation**'s own gene ids.

    An entry of a :class:`TFGeneList`. The census's four uniform columns are fields of
    their own — the **Gene id stem** it is keyed by, the symbol, the TF flag and the **DBD
    family** — and everything the publisher recorded beyond them stays under the
    publisher's own name in :attr:`judgements`, because beyond those four no two censuses
    carry the same columns and nothing here compares one publisher's with another's
    (ADR-0014).

    ``gene_ids`` is a tuple because one stem may name more than one gene id in one
    annotation and this never picks one — see
    :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids`. It is never empty: a stem
    the annotation carries no gene for is in :attr:`TFGeneList.unresolved` instead of here.

    Attributes
    ----------
    gene_id_stem : str
        The **Gene id stem** the census is keyed by.
    gene_ids : tuple of str
        Every gene id this annotation spells that stem with, ascending.
    symbol : str or None
        The gene symbol the census records, or ``None`` where it records none.
    is_tf : bool
        The census's own TF flag. ``True`` for every gene of a **TF gene list** left at
        its default, and the field that tells a rejected gene from an accepted one when a
        caller widened to carry both.
    dbd_family : str or None
        The **DBD family** this census classifies the gene under, in the publisher's own
        vocabulary — group by it within a species and never across two.
    judgements : mapping of str to (str or None)
        Every other column the census records for this gene, under the census's own
        snake_case spelling of its published name: the **TF assessment** a caller tightens
        or loosens on, the binding mode, the motif status, the KRAB flag and the
        third-party votes, for a census that records them. ``None`` is a cell the
        publisher left blank.

    Examples
    --------
    >>> gene = TFGene(
    ...     gene_id_stem="ENSG00000214717",
    ...     gene_ids=("ENSG00000214717.13", "ENSG00000214717.13_PAR_Y"),
    ...     symbol="ZBED1",
    ...     is_tf=True,
    ...     dbd_family="BED ZF",
    ...     judgements={"tf_assessment": "Known motif"},
    ... )
    >>> gene.judgements["tf_assessment"]
    'Known motif'
    >>> gene.as_json()["gene_ids"]
    ['ENSG00000214717.13', 'ENSG00000214717.13_PAR_Y']
    """

    gene_id_stem: str
    gene_ids: tuple[str, ...]
    symbol: str | None
    is_tf: bool
    dbd_family: str | None
    judgements: Mapping[str, str | None]

    def as_json(self) -> dict[str, Any]:
        """Return this gene as ``--json`` serializes it.

        Returns
        -------
        dict
            The fields above under their own names, with ``gene_ids`` as a list and
            ``judgements`` as a plain mapping under the census's own column names.
        """
        return {
            "gene_id_stem": self.gene_id_stem,
            "gene_ids": list(self.gene_ids),
            "symbol": self.symbol,
            "is_tf": self.is_tf,
            "dbd_family": self.dbd_family,
            "judgements": dict(self.judgements),
        }


@dataclass(frozen=True)
class TFGeneList:
    """One **Assembly**'s **TF gene**s, in its registered annotation's own gene ids.

    :func:`resolve_tf_genes`'s answer, and what a ``--json`` surface over it serializes.
    The census's **Gene id stem**s resolved against one annotation, so the ids join to a
    counts matrix with nothing left to normalise, and assessed-positive by default: the
    common case is not 2,765 rows to filter down to 1,639.

    **Nothing here decides what a transcription factor is.** Every verdict is the census's
    and travels with :attr:`provenance`, which names the publisher to cite. Two censuses
    that classify one factor differently are two answers rather than a contradiction, and
    this says which one is speaking.

    **What the census holds and this annotation does not is visible.** A stem no gene id
    here is of comes back in :attr:`unresolved` rather than being dropped, so a caller can
    count what the crossing cost instead of wondering.

    There is no empty one for the reasons an absent census would give: an assembly whose
    species has no census, and one nothing names a species for, each raise a
    :class:`LookupError` of their own.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    species : str
        The species the assembly's own metadata row names, which is what selected the
        census. Never passed in by a caller, so asking for one species' transcription
        factors while holding another species' assembly is not expressible (ADR-0003).
    provenance : genome.tf.gene.census.CensusProvenance
        Where the census came from: publisher, version, PubMed id, source URL and digest.
        :meth:`~genome.tf.gene.census.CensusProvenance.attribution` renders the line to
        print beside anything it answered.
    genes : tuple of TFGene
        One entry per **Gene id stem** that named at least one gene id here, in the
        census's own row order.
    unresolved : tuple of str
        The stems this annotation carries no gene for, in census row order.

    Examples
    --------
    >>> from genome.tf.gene import tf_gene_table
    >>> answer = TFGeneList(
    ...     assembly="hg38",
    ...     annotation="gencode_v50",
    ...     species="Homo sapiens",
    ...     provenance=tf_gene_table("Homo sapiens").provenance,
    ...     genes=(
    ...         TFGene("ENSG00000137203", ("ENSG00000137203.12",), "TFAP2A", True, "AP-2", {}),
    ...     ),
    ...     unresolved=("ENSG00000214717",),
    ... )
    >>> answer.gene_ids
    ['ENSG00000137203.12']
    >>> answer.provenance.publisher
    'Lambert et al. 2018'
    >>> answer.as_json()["unresolved"]
    ['ENSG00000214717']
    """

    assembly: str
    annotation: str
    species: str
    provenance: CensusProvenance
    genes: tuple[TFGene, ...]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id, gene order then id order — a fresh list each call.

        **Every** id, not one per gene, for the reason
        :attr:`~genome.annotation.stems.ResolvedGeneIds.gene_ids` gives: flattening is where a
        reader would take the first id of a stem that names two and lose the other.
        :attr:`genes` is what says which gene an id came from, and what the census said
        about it.
        """
        return [gene_id for gene in self.genes for gene_id in gene.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``species``, the census's ``provenance`` under
            its own field names, ``genes`` as a list of :meth:`TFGene.as_json` entries,
            the flattened ``gene_ids``, and ``unresolved`` as a list. The ids are written
            out beside the genes they are read from for the reason :attr:`gene_ids` gives.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "species": self.species,
            "provenance": asdict(self.provenance),
            "genes": [gene.as_json() for gene in self.genes],
            "gene_ids": self.gene_ids,
            "unresolved": list(self.unresolved),
        }


def resolve_tf_genes(
    registry: AnnotationRegistry, name: str | None = None, *, include_rejected: bool = False
) -> TFGeneList:
    """Return the genes a census judges transcription factors, in one annotation's ids.

    The **TF gene table** :mod:`genome.tf.gene` ships for ``registry``'s assembly's
    species, met with one registered annotation: every **Gene id stem** the census is
    keyed by is resolved through
    :meth:`~genome.annotation.registry.AnnotationRegistry.resolve_gene_ids` into the gene ids that
    annotation actually spells, so the answer joins to a counts matrix with nothing left
    for the caller to normalise. A stem naming two gene ids answers with both, and the
    stems the annotation carries no gene for ride back on :attr:`TFGeneList.unresolved`
    rather than being dropped.

    **The species is the assembly's own** — the curated metadata table's, read here and
    never passed in, so asking for human transcription factors while holding a mouse
    assembly is not expressible (ADR-0003).

    **Nothing here decides what a transcription factor is.** The verdict is the census's,
    and the answer carries the publisher, version and PubMed id that reached it.
    Assessed-positive by default, because the common case is not 2,765 rows to filter down
    to 1,639; ``include_rejected`` widens to the genes the census assessed and turned
    down, which is the whole census and still not every gene there is — one it never
    assessed is absent from both answers, and that is a third fact.

    Wanting only ``Known motif``, or wanting ``Inferred motif`` included, is a re-filter
    on the **TF assessment** each gene already carries in :attr:`TFGene.judgements` rather
    than a second flag here.

    Parameters
    ----------
    registry : genome.annotation.registry.AnnotationRegistry
        The registry of the **Assembly** whose species selects the census and whose
        annotation the ids should be in. :class:`~genome.assembly.genome.Genome` holds one already.
    name : str, optional
        The **Registered name** to answer in the gene ids of. Omitted, that assembly's
        **Default annotation** answers.
    include_rejected : bool, default False
        Carry the genes the census assessed and judged *not* to be transcription factors
        as well, each saying so in :attr:`TFGene.is_tf`. A census that records no
        rejections answers the same either way.

    Returns
    -------
    TFGeneList
        The genes, the census's provenance, and the stems that resolved to nothing.

    Raises
    ------
    genome.tf.species.UnknownSpeciesError
        If nothing names that assembly's species — a chimera, or an assembly no curated
        row lists.
    genome.tf.species.NoTFCensusError
        If no census ships for that species; the message names the ones that do.
    ValueError
        If ``name`` is omitted and no **Default annotation** is decided.
    genome.annotation.registry.AnnotationNotRegisteredError
        If nothing of that name is registered there.
    genome.annotation.stems.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.annotation import AnnotationRegistry
    >>> from genome.tf.gene import resolve_tf_genes
    >>> registry = AnnotationRegistry.locate("hg38")               # doctest: +SKIP
    >>> answer = resolve_tf_genes(registry, "gencode_v50")         # doctest: +SKIP
    >>> answer.genes[0].symbol, answer.genes[0].dbd_family         # doctest: +SKIP
    ('TFAP2A', 'AP-2')
    >>> answer.provenance.publisher                                # doctest: +SKIP
    'Lambert et al. 2018'
    """
    species = assembly_metadata(registry.assembly).species
    if species is None:
        raise UnknownSpeciesError(registry.assembly, _censused_species(), shipped_table="TF census")
    census = tf_gene_table(species)
    if census is None:
        raise NoTFCensusError(registry.assembly, species, _censused_species())
    resolved = registry.resolve_gene_ids(
        census.gene_id_stems if include_rejected else census.assessed_positive, name
    )
    return TFGeneList(
        assembly=registry.assembly,
        annotation=resolved.annotation,
        species=species,
        provenance=census.provenance,
        genes=_tf_genes(census, resolved),
        unresolved=resolved.unresolved,
    )


def tf_gene_list(
    assembly: str,
    *,
    annotation: str | None = None,
    include_rejected: bool = False,
    cache_dir: str | Path | None = None,
) -> TFGeneList:
    """Return the genes a published census judges transcription factors in ``assembly``.

    :func:`resolve_tf_genes` for an assembly named rather than opened, which is what
    ``genome tf gene-list`` runs. A registry is built for the length of the call, so a
    shell surface over it adds no second code path. Nothing is prepared, fetched or built
    to answer it, and the census is read from inside the package.

    Parameters
    ----------
    assembly : str
        The assembly to ask about, e.g. ``"hg38"``. Its own metadata row names the
        species, which is what selects the census.
    annotation : str, optional
        The **Registered name** to answer in the gene ids of; the **Default annotation**
        when omitted.
    include_rejected : bool, default False
        Carry the genes the census assessed and turned down as well.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory is inspected, as
        :meth:`~genome.annotation.registry.AnnotationRegistry.locate` takes it.

    Returns
    -------
    TFGeneList
        The answer :func:`resolve_tf_genes` describes.

    Raises
    ------
    genome.tf.species.UnknownSpeciesError
        If nothing names the assembly's species.
    genome.tf.species.NoTFCensusError
        If no census ships for that species.
    ValueError
        If ``annotation`` is omitted and no **Default annotation** is decided.
    genome.annotation.registry.AnnotationNotRegisteredError
        If that annotation is not registered here.
    genome.annotation.stems.NoGeneFeaturesError
        If its database holds no gene at all.

    Examples
    --------
    >>> from genome.tf.gene import tf_gene_list
    >>> tf_gene_list("hg38").provenance.publisher            # doctest: +SKIP
    'Lambert et al. 2018'
    """
    return resolve_tf_genes(
        AnnotationRegistry.locate(assembly, cache_dir),
        annotation,
        include_rejected=include_rejected,
    )


def _censused_species() -> tuple[str, ...]:
    """Return the species a census ships for, in the spelling an assembly's row uses.

    What both absences name as the thing a caller can ask about instead. The shipped files
    are what is enumerated, since one is what makes a species answerable; the provenance
    table beside them is only read for the publisher's own spelling, so a census shipping
    without a row is still named — badly, as its slug, which is the state
    :func:`~genome.tf.gene.census.tf_gene_table` raises over anyway.
    """
    named = {species_slug(record.species): record.species for record in census_metadata()}
    return tuple(named.get(slug, slug) for slug in census_species())


def _tf_genes(census: TFGeneTable, resolved: ResolvedGeneIds) -> tuple[TFGene, ...]:
    """Return one entry per resolved stem, carrying the census's own row for that gene.

    The census's four uniform columns become fields and everything after them stays under
    the publisher's own name, because beyond those four no two censuses carry the same
    columns and nothing here compares one publisher's with another's (ADR-0014). The order
    is the census's own row order, which is the order the stems were asked about.
    """
    publisher_columns = census.columns[len(UNIFORM_COLUMNS) :]
    rows = {row[0]: row for row in census.rows}
    genes: list[TFGene] = []
    for stem, gene_ids in resolved.resolved.items():
        cells = dict(zip(census.columns, rows[stem], strict=True))
        genes.append(
            TFGene(
                gene_id_stem=stem,
                gene_ids=gene_ids,
                symbol=cells["symbol"],
                is_tf=cells["is_tf"] == TRUE_CELL,
                dbd_family=cells["dbd_family"],
                judgements={name: cells[name] for name in publisher_columns},
            )
        )
    return tuple(genes)
