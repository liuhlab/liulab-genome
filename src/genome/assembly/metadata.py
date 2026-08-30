"""The curated assembly table — which references the lab supports, and how each is named.

One small hand-maintained TSV ships inside the package. ``data/assembly_metadata.tsv``
records, for each known reference assembly, its canonical name and the cross-references
used to talk about it across databases: species, UCSC name, NCBI name, NCBI assembly
accession, and NCBI taxonomy id. A row may additionally pin where that assembly's FASTA
is fetched from and the sha256 of the **unpacked** FASTA it yields, which is what makes
preparing the assembly reproducible, and it may carry the longest gap a spliced aligner
should take for an intron on that assembly together with the reason that number was
chosen — a bound for a consumer to apply, hand-set and never derived here (ADR-0010).

Two accessors, because there are two questions. :func:`assembly_metadata` answers *what
is known about this assembly* and always answers with a record: the table's row, or one
carrying the name with every identifier unknown, so nothing downstream guards a missing
record before reading a field. :func:`lookup_assembly` answers *does the curated table
list this name*, and only that question needs ``None`` — it is what tells a chimera's
derived name from a free-form local key on a machine holding neither (ADR-0003).

The annotation table is the same idea one level down and belongs to the other context that
reads it: it is :mod:`genome.annotation.metadata`, keyed by assembly plus **Registered
name**. The two tables share their shape and nothing else — no lookup here reads that one,
and none there reads this one — so each sits with the context whose callers ask it
questions rather than both sitting between them.

The table is not an allow-list: a row means *officially supported* — pinned source, pinned
checksum — and an assembly with no row is perfectly legal (ADR-0003).

:class:`AssemblyMetadata` declares its field list once: the table is read through those
fields column by column, parsed by each field's own declared type, and a whole record is
what :class:`~genome.assembly.genome.Genome` and the registration functions take to
override the table. It is a **Shipped table** and is read by the one loader every such
table goes through, :mod:`genome.shipped`: this module declares where the table lives,
what its columns are and what repairs it, and the header validation, the cell parsing and
the shape of the error a broken file raises all happen there.

:func:`parse_cell` and :func:`species_slug` live with that reader — the first because five
other tables parse their cells the same way, the second because naming a file after a
species is what every shipped-data directory in the package does and none of them owns.
Both are re-exported here, so a caller that reads a curated row and then names a file after
its species imports the pair from one place.

A row is a record's other spelling, and both directions are public:
:func:`format_table_row` writes one, :meth:`AssemblyMetadata.from_row` reads one. A table
is the records it lists, so every lookup here takes ``table=`` and a caller with rows of
their own hands those over instead of the shipped ones — which is theirs to do, the table
being a cross-reference and not an allow-list (ADR-0003). The shipped TSV is read in
:func:`assembly_table`, and nowhere else.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from typing import Any, get_type_hints

from genome.shipped import MetadataRowError, ShippedTable, parse_cell, species_slug

__all__ = [
    "METADATA_FIELDS",
    "AssemblyMetadata",
    "MetadataRowError",
    "assembly_metadata",
    "assembly_table",
    "format_table_row",
    "lookup_assembly",
    "parse_cell",
    "species_slug",
]

#: Location of the curated assembly metadata table within the package.
_METADATA_RESOURCE = "data/assembly_metadata.tsv"

#: What repairs it: the table is hand-maintained rather than generated, so the repair is an
#: edit and the message says which file to make it in.
_ASSEMBLY_REPAIR = f"edit {_METADATA_RESOURCE}, which is maintained by hand"


@dataclass(frozen=True)
class AssemblyMetadata:
    """Identifiers for one reference assembly (one row of the metadata table).

    The single declaration of what an assembly's metadata consists of: the table
    is parsed through these fields, and a complete record is what
    :class:`~genome.assembly.genome.Genome` accepts in place of the table's own row.

    Only ``assembly_name`` is required. Every other column may be left blank, and a
    blank cell reads back as ``None`` rather than as text: the table fills in over
    time, and a freshly prepared assembly pins its source and digest well before
    anyone supplies its species, its UCSC and NCBI names or its taxonomy id.

    ``source_url`` and ``sha256`` are what makes preparing an assembly reproducible.
    ``source_url`` pins where its FASTA is fetched from, so nothing has to be derived
    or guessed; ``sha256`` pins the digest of the **unpacked** FASTA that source
    yields — not of the compressed archive it arrives in, so a copy taken from a
    mirror or recompressed elsewhere still matches (ADR-0006). A row with no digest is
    unverified rather than wrong.

    ``intron_length_cap`` is the longest gap a spliced aligner should take for an
    intron on this assembly, and ``intron_length_cap_rationale`` says why that number
    and not another. It is a deliberately loose round number set by hand, never
    computed from an annotation — an annotation catalogues the transcripts someone
    observed, so its longest intron is a floor on what the organism does rather than a
    ceiling on it (ADR-0010). Nothing in this package reads either field: they are
    curated here so that a consumer choosing aligner parameters reads a fact about the
    assembly from the same row as its identifiers. A blank cap is an assembly nobody
    has characterised, which is legal and says *no bound has been chosen* — the reading
    that leaves such an assembly aligning exactly as it did before.

    A **Chimera**'s row pins neither, and that is deliberate rather than pending: its
    bytes are not fetched from anywhere, and they are derived by a pure function from
    components whose own rows are pinned, so it is proven transitively. The table pins
    what was downloaded; what was derived is pinned by a test, where a change in this
    package's own concatenation code belongs (ADR-0008). Such a row carries its name and
    nothing else — its identifiers are its components', reachable through
    :attr:`~genome.assembly.genome.Genome.chrom_components` — and it exists so that a machine
    holding none of them can still tell a chimera's name from a local key someone chose.

    ``ucsc_name`` is blank permanently rather than pending in some rows, for the same
    reason the assembly id is a local key rather than a UCSC one (ADR-0003): the lab
    supports references UCSC has never carried, and such a row simply has no name in
    that namespace to give.

    Examples
    --------
    >>> record = AssemblyMetadata(
    ...     "hg38", "Homo sapiens", "hg38", "GRCh38", "GCF_000001405.40", 9606
    ... )
    >>> record.species
    'Homo sapiens'
    >>> record.sha256 is None            # nothing pinned unless the row says so
    True
    """

    assembly_name: str
    species: str | None
    ucsc_name: str | None
    ncbi_name: str | None
    ncbi_assembly_id: str | None
    ncbi_taxid: int | None
    source_url: str | None = None
    sha256: str | None = None
    intron_length_cap: int | None = None
    intron_length_cap_rationale: str | None = None

    @classmethod
    def unknown(cls, assembly_name: str) -> AssemblyMetadata:
        """Return the record for an assembly the table does not list: the name, and nothing else.

        What *unlisted* looks like as a record rather than as a missing one. Every
        identifier is genuinely unknown, which a blank cell already means everywhere
        else, so a caller reads a field and gets ``None`` instead of first asking
        whether there is a record to read it off. It is exactly the line
        :func:`format_table_row` emits for an assembly nobody has curated yet.

        The name is the local key the caller asked for, since that is the only
        identifier an unlisted assembly has — the assembly id is a local key and the
        table is a cross-reference rather than an allow-list (ADR-0003). This answers
        *what is known about this assembly*; whether the table lists it at all is
        :func:`lookup_assembly`'s question and stays a separate one.

        Parameters
        ----------
        assembly_name : str
            The assembly the record is for.

        Returns
        -------
        AssemblyMetadata
            A record carrying ``assembly_name`` with every other field ``None``.

        Examples
        --------
        >>> record = AssemblyMetadata.unknown("my_ref")
        >>> record.assembly_name
        'my_ref'
        >>> record.species is None and record.sha256 is None
        True
        """
        return cls(
            assembly_name=assembly_name,
            species=None,
            ucsc_name=None,
            ncbi_name=None,
            ncbi_assembly_id=None,
            ncbi_taxid=None,
        )

    @classmethod
    def from_row(cls, row: Mapping[str, object]) -> AssemblyMetadata:
        """Build a record from one row of a metadata table.

        The reader half of the register-then-paste flow :func:`format_table_row` writes:
        a row is a mapping of column name to cell, which is how the shipped TSV is read
        and how :func:`dataclasses.asdict` of a record spells one. Each column is parsed
        by its own declared type, and a blank cell — empty, absent, or the NaN pandas
        reads a blank as — means unknown, which only a column that has an unknown takes.

        Parameters
        ----------
        row : mapping of str to object
            Column name to cell. A cell is the table's own text, but a value already of
            the column's type is taken as it stands, so a record's fields are a row.
            Keys outside :data:`METADATA_FIELDS` are ignored.

        Returns
        -------
        AssemblyMetadata
            The record the row spells.

        Raises
        ------
        MetadataRowError
            If a cell cannot be read as its column's type, or a column that has no
            unknown is blank. The record is built from parsed cells or not at all, so a
            caller is handed a whole record or an error naming the column — never a
            record carrying the columns that happened to come before the bad one.

        Examples
        --------
        >>> row = {"assembly_name": "sacCer3", "ncbi_taxid": "559292"}
        >>> AssemblyMetadata.from_row(row).ncbi_taxid
        559292
        >>> AssemblyMetadata.from_row(row).species is None   # blank is unknown
        True
        """
        return cls(**{name: parse_cell(name, row, _FIELD_TYPES) for name in METADATA_FIELDS})


#: Each metadata field's declared type, which parses that field's column of the table.
_FIELD_TYPES: dict[str, Any] = get_type_hints(AssemblyMetadata)

#: The metadata field names, in table-column order — the columns every row carries.
METADATA_FIELDS: tuple[str, ...] = tuple(_FIELD_TYPES)

#: The assembly table as a **Shipped table**: where it lives, what its header is, what it
#: is called and what repairs it. Every check the header and the file are held to lives in
#: :mod:`genome.shipped`; a blank cell is the record's own business, since which columns
#: have an unknown is declared by the field types and not here.
_ASSEMBLY_TABLE = ShippedTable(
    resource=_METADATA_RESOURCE,
    columns=METADATA_FIELDS,
    noun="assembly metadata table",
    repair=_ASSEMBLY_REPAIR,
    error=MetadataRowError,
    identify=("assembly_name",),
)


def format_table_row(row: Mapping[str, object]) -> str:
    r"""Render one metadata row as a tab-separated line, in table-column order.

    Blank means unknown: a field that is ``None``, or missing from ``row`` altogether,
    renders as an empty cell — exactly how the shipped table spells a value nobody has
    pinned yet. The result is the line to paste into ``data/assembly_metadata.tsv``,
    with no trailing newline.

    Parameters
    ----------
    row : mapping of str to object
        Field name to value, such as :func:`dataclasses.asdict` of an
        :class:`AssemblyMetadata`. Keys outside :data:`METADATA_FIELDS` are ignored.

    Returns
    -------
    str
        The row's values joined by tabs, in :data:`METADATA_FIELDS` order.

    Examples
    --------
    >>> format_table_row({"assembly_name": "sacCer3", "ncbi_taxid": 559292})
    'sacCer3\t\t\t\t\t559292\t\t\t\t'
    """
    return "\t".join("" if row.get(name) is None else str(row[name]) for name in METADATA_FIELDS)


@cache
def assembly_table() -> tuple[AssemblyMetadata, ...]:
    """Return every assembly the shipped table lists, in table order.

    What the lab officially supports — a pinned source and a pinned checksum each — and
    what every lookup here reads when it is handed no ``table=`` of its own. Read once
    and cached; the records are frozen, so the tuple is safe to hold on to.

    Returns
    -------
    tuple of AssemblyMetadata
        One record per row of ``data/assembly_metadata.tsv``.

    Raises
    ------
    MetadataRowError
        If the shipped file is empty, its header is not :data:`METADATA_FIELDS`, a row
        holds the wrong number of cells, or a cell cannot be read as its column's type.

    Examples
    --------
    >>> "hg38" in {record.assembly_name for record in assembly_table()}
    True
    """
    return tuple(AssemblyMetadata.from_row(row) for row in _ASSEMBLY_TABLE.read().mappings())


def lookup_assembly(
    assembly: str, *, table: Sequence[AssemblyMetadata] | None = None
) -> AssemblyMetadata | None:
    """Return the :class:`AssemblyMetadata` for a UCSC (or canonical) assembly name, or ``None``.

    Parameters
    ----------
    assembly : str
        The name to look up, matched against each record's ``ucsc_name`` and
        ``assembly_name``.
    table : sequence of AssemblyMetadata, optional
        The rows to look in; the shipped table (:func:`assembly_table`) when omitted.
        Curating rows of your own is ordinary rather than an override — the table is a
        cross-reference and never an allow-list (ADR-0003) — and nothing is installed
        by passing them: the sequence is read for this call and no other.

    Returns
    -------
    AssemblyMetadata or None
        The row for ``assembly``, or ``None`` when the table does not list it. An
        unlisted assembly is legal and its identifiers are simply unknown. A blank cell
        in an optional column reads back as ``None``.

    Examples
    --------
    >>> lookup_assembly("hg38").ncbi_name
    'GRCh38'
    >>> lookup_assembly("hg38").source_url
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    >>> lookup_assembly("no_such_assembly") is None
    True
    >>> mine = AssemblyMetadata.unknown("my_ref")
    >>> lookup_assembly("my_ref", table=[mine]) == mine
    True
    """
    rows = assembly_table() if table is None else table
    # A blank ``ucsc_name`` is ``None`` and equals no name asked for, so a row for a
    # reference UCSC never carried is found by its own name and never by an empty one.
    return next(
        (row for row in rows if assembly in (row.ucsc_name, row.assembly_name)),
        None,
    )


def assembly_metadata(
    assembly: str, *, table: Sequence[AssemblyMetadata] | None = None
) -> AssemblyMetadata:
    """Return what is known about ``assembly`` — the table's row, or an unknown record.

    The **total** accessor, and the one to reach for when the question is *what are this
    assembly's identifiers*: there is always a record, so a caller reads a field rather
    than a record and then a field. An assembly the table does not list has every
    identifier ``None`` and its own name, which is what an unlisted assembly knows about
    itself.

    Deliberately not the same function as :func:`lookup_assembly`, which answers the
    other question — *does the curated table list this name?* — and keeps its ``None``
    for it. That answer is what separates a chimera's derived name from a free-form local
    key on a machine holding neither (ADR-0003, ADR-0008), so it cannot be made total
    without reading every name as listed.

    Parameters
    ----------
    assembly : str
        The name to look up, matched as :func:`lookup_assembly` matches it.
    table : sequence of AssemblyMetadata, optional
        The rows to look in, as :func:`lookup_assembly` takes them.

    Returns
    -------
    AssemblyMetadata
        The table's row for ``assembly``, else
        :meth:`AssemblyMetadata.unknown(assembly) <AssemblyMetadata.unknown>`.

    Examples
    --------
    >>> assembly_metadata("hg38").ncbi_name
    'GRCh38'
    >>> assembly_metadata("no_such_assembly").ncbi_name is None
    True
    >>> assembly_metadata("no_such_assembly").assembly_name
    'no_such_assembly'
    """
    listed = lookup_assembly(assembly, table=table)
    return listed if listed is not None else AssemblyMetadata.unknown(assembly)
