"""The curated tables — which assemblies the lab supports, and what annotates them.

Two small hand-maintained TSVs ship inside the package. ``data/assembly_metadata.tsv``
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

``data/annotation_metadata.tsv`` is the same idea one level down: keyed by assembly plus
the name an annotation is registered under, it says who publishes that annotation, which
release it is, where to fetch it and what the unpacked **GTF** hashes to, so naming one
is enough to register it. :func:`lookup_annotation` answers for one name and
:func:`list_annotation_metadata` for everything an assembly offers.

Neither table is an allow-list: a row means *officially supported* — pinned source,
pinned checksum — and an assembly or an annotation with no row is perfectly legal
(ADR-0003).

Each dataclass declares its field list once: its table is read through those fields
column by column, parsed by each field's own declared type, and a whole record is what
:class:`~genome.genome.Genome` and the registration functions take to override the table.

A row is a record's other spelling, and both directions are public:
:func:`format_table_row` writes one, :meth:`AssemblyMetadata.from_row` and
:meth:`AnnotationMetadata.from_row` read one. A table is the records it lists, so every
lookup here takes ``table=`` and a caller with rows of their own hands those over instead
of the shipped ones — which is theirs to do, the table being a cross-reference and not an
allow-list (ADR-0003). The shipped TSVs are read in :func:`assembly_table` and
:func:`annotation_table`, and nowhere else.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any, get_args, get_type_hints

import pandas as pd

#: Location of the curated assembly metadata table within the package.
_METADATA_RESOURCE = "data/assembly_metadata.tsv"

#: Location of the curated annotation metadata table within the package.
_ANNOTATION_RESOURCE = "data/annotation_metadata.tsv"

#: Cell spellings a boolean column accepts, lower-cased. Anything else in one is a
#: typo in a hand-maintained table and says so rather than reading as ``False``.
_TRUE_CELLS = frozenset({"yes", "true", "1"})
_FALSE_CELLS = frozenset({"no", "false", "0"})


class MetadataRowError(ValueError):
    """A row cannot be read as a record, and the message names the column that refused.

    Raised by :meth:`AssemblyMetadata.from_row` and :meth:`AnnotationMetadata.from_row`,
    for a cell no column's type can read and for a blank cell in a column that has no
    unknown. A :class:`ValueError`, because a hand-maintained row that says something
    the columns do not is a bad value rather than a broken program.

    Examples
    --------
    >>> try:
    ...     AssemblyMetadata.from_row({"assembly_name": "tiny", "ncbi_taxid": "many"})
    ... except MetadataRowError as error:
    ...     print("ncbi_taxid" in str(error))
    True
    """


@dataclass(frozen=True)
class AssemblyMetadata:
    """Identifiers for one reference assembly (one row of the metadata table).

    The single declaration of what an assembly's metadata consists of: the table
    is parsed through these fields, and a complete record is what
    :class:`~genome.genome.Genome` accepts in place of the table's own row.

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
    :attr:`~genome.genome.Genome.chrom_components` — and it exists so that a machine
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
        return cls(**{name: _parse_cell(name, row, _FIELD_TYPES) for name in METADATA_FIELDS})


#: Each metadata field's declared type, which parses that field's column of the table.
_FIELD_TYPES: dict[str, Any] = get_type_hints(AssemblyMetadata)

#: The metadata field names, in table-column order — the columns every row carries.
METADATA_FIELDS: tuple[str, ...] = tuple(_FIELD_TYPES)


@dataclass(frozen=True)
class AnnotationMetadata:
    """One annotation the lab supports for one assembly (one row of the annotation table).

    Keyed by ``assembly`` plus ``name`` — the **Registered name** the annotation is
    addressed by everywhere — and carrying enough to register it from that name alone:
    who publishes it, which release, where to fetch it, and what the **unpacked** GTF
    that source yields hashes to (ADR-0006). A complete record is also what the
    registration functions accept in place of the table's own row.

    Attributes
    ----------
    assembly : str
        The assembly this annotation belongs to — an annotation belongs to exactly one.
    name : str
        The **Registered name**, unique within the assembly, e.g. ``"gencode_v50"``.
    provider : str
        Who publishes it: ``"GENCODE"``, ``"RefSeq"``, ``"WormBase"``, ``"UCSC"``.
    version : str
        The provider's own release identifier, e.g. ``"v50"`` or ``"WS298"``.
    url : str
        Where the GTF is fetched from.
    sha256 : str or None
        Digest of the **unpacked** GTF, or ``None`` when the row pins none — in which
        case whatever is fetched is recorded but nothing is compared.
    default : bool
        Whether this is the assembly's **Default annotation**.

    Examples
    --------
    >>> record = AnnotationMetadata(
    ...     assembly="sacCer3",
    ...     name="ensgene_v101",
    ...     provider="UCSC",
    ...     version="ensGene.v101",
    ...     url="https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz",
    ... )
    >>> record.provider
    'UCSC'
    >>> record.default                   # not the default unless the row says so
    False
    """

    assembly: str
    name: str
    provider: str
    version: str
    url: str
    sha256: str | None = None
    default: bool = False

    @classmethod
    def from_row(cls, row: Mapping[str, object]) -> AnnotationMetadata:
        """Build a record from one row of an annotation table.

        :meth:`AssemblyMetadata.from_row` for the annotation table's own columns, and
        the same rules — with one more: ``default`` is a flag column, where a blank cell
        is the real answer *no* rather than an unknown.

        Parameters
        ----------
        row : mapping of str to object
            Column name to cell, as :meth:`AssemblyMetadata.from_row` takes one. Keys
            outside :data:`ANNOTATION_FIELDS` are ignored.

        Returns
        -------
        AnnotationMetadata
            The record the row spells.

        Raises
        ------
        MetadataRowError
            As :meth:`AssemblyMetadata.from_row` raises it, and additionally for a flag
            cell spelled a way no row spells one.

        Examples
        --------
        >>> AnnotationMetadata.from_row(
        ...     {
        ...         "assembly": "sacCer3",
        ...         "name": "ensgene_v101",
        ...         "provider": "UCSC",
        ...         "version": "ensGene.v101",
        ...         "url": "https://example.org/sacCer3.ensGene.gtf.gz",
        ...         "default": "yes",
        ...     }
        ... ).default
        True
        """
        return cls(
            **{name: _parse_cell(name, row, _ANNOTATION_FIELD_TYPES) for name in ANNOTATION_FIELDS}
        )


#: Each annotation field's declared type, which parses that field's column of the table.
_ANNOTATION_FIELD_TYPES: dict[str, Any] = get_type_hints(AnnotationMetadata)

#: The annotation field names, in table-column order.
ANNOTATION_FIELDS: tuple[str, ...] = tuple(_ANNOTATION_FIELD_TYPES)


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

    Examples
    --------
    >>> "hg38" in {record.assembly_name for record in assembly_table()}
    True
    """
    return tuple(AssemblyMetadata.from_row(row) for row in _rows(_METADATA_RESOURCE))


@cache
def annotation_table() -> tuple[AnnotationMetadata, ...]:
    """Return every annotation the shipped table lists, for all assemblies, in table order.

    The annotation half of :func:`assembly_table`, and the default the annotation
    lookups read. :func:`list_annotation_metadata` is the usual way in, since an
    annotation is asked for by the assembly it belongs to.

    Returns
    -------
    tuple of AnnotationMetadata
        One record per row of ``data/annotation_metadata.tsv``.

    Examples
    --------
    >>> "ensgene_v101" in {record.name for record in annotation_table()}
    True
    """
    return tuple(AnnotationMetadata.from_row(row) for row in _rows(_ANNOTATION_RESOURCE))


def _rows(resource_name: str) -> list[dict[str, Any]]:
    """Read one shipped TSV as rows of text — every column, every cell, unparsed."""
    resource = files("genome").joinpath(resource_name)
    with resource.open("r", encoding="utf-8") as handle:
        frame = pd.read_csv(handle, sep="\t", dtype=str)
    return [dict(row) for _, row in frame.iterrows()]


def _cell_text(name: str, row: Mapping[str, object]) -> str:
    """Return the ``name`` cell of ``row`` as stripped text, empty when it is blank.

    Blank is a missing column, ``None``, the NaN pandas reads a blank cell as, or
    whitespace. Anything else is that cell's own text, so a row that spells a value in
    its column's own type — an ``int`` taxid, a ``bool`` flag — reads back the same as
    the table's text for it.
    """
    cell = row.get(name)
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return ""
    return str(cell).strip()


def _parse_flag(name: str, row: Mapping[str, object]) -> bool:
    """Parse a boolean column — blank is ``False``, and a spelling nobody uses raises."""
    text = _cell_text(name, row).lower()
    if not text or text in _FALSE_CELLS:
        return False
    if text in _TRUE_CELLS:
        return True
    accepted = ", ".join(sorted(_TRUE_CELLS | _FALSE_CELLS))
    raise MetadataRowError(
        f"the {name!r} column holds {text!r}, which is not a flag. Fix that cell to one "
        f"of: {accepted} — or leave it blank, which reads as false."
    )


def _parse_cell(name: str, row: Mapping[str, object], types: Mapping[str, Any]) -> Any:
    """Parse the ``name`` column of ``row`` with that field's own declared type.

    A field declared optional (``T | None``) is parsed by ``T`` when its cell carries
    text, and is ``None`` when the cell is blank or its column is absent — a union is
    not callable, so the type inside it does the parsing. A field declared ``bool`` is a
    flag column, where an empty cell is a real answer (see :func:`_parse_flag`) rather
    than an unknown. Any other field is parsed by its declared type and has no unknown,
    so a blank cell there is a malformed row rather than the text ``'nan'``.
    """
    declared = types[name]
    if declared is bool:
        return _parse_flag(name, row)
    inside = [arg for arg in get_args(declared) if arg is not type(None)]
    text = _cell_text(name, row)
    if inside:
        return _parse_text(name, text, inside[0]) if text else None
    if not text:
        raise MetadataRowError(
            f"the {name!r} column is blank, and it is one no row may leave blank. Fill "
            f"that cell in: a blank cell reads back as unknown, and {name!r} has none."
        )
    return _parse_text(name, text, declared)


def _parse_text(name: str, text: str, declared: Any) -> Any:
    """Parse one cell's text with its column's type, or say which cell refused."""
    try:
        return declared(text)
    except (TypeError, ValueError) as error:
        raise MetadataRowError(
            f"the {name!r} column holds {text!r}, which {declared.__name__} cannot read. "
            f"Fix that cell to a value {declared.__name__} accepts."
        ) from error


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


def lookup_annotation(
    assembly: str, name: str, *, table: Sequence[AnnotationMetadata] | None = None
) -> AnnotationMetadata | None:
    """Return the annotation ``name`` the table lists for ``assembly``, or ``None``.

    The lookup that makes registering an annotation by name possible: the row it
    returns says where to fetch the GTF and what it must hash to.

    Parameters
    ----------
    assembly : str
        The assembly the annotation belongs to, matched against the table's
        ``assembly`` column.
    name : str
        The **Registered name** to look up, e.g. ``"gencode_v50"``.
    table : sequence of AnnotationMetadata, optional
        The rows to look in; the shipped table (:func:`annotation_table`) when omitted,
        and read for this call only, as :func:`lookup_assembly` reads its own.

    Returns
    -------
    AnnotationMetadata or None
        The row, or ``None`` when the table lists no such annotation for that
        assembly. The table is a cross-reference, not an allow-list, so an unlisted
        annotation is legal — it is registered by path instead of by name.

    Examples
    --------
    >>> lookup_annotation("sacCer3", "ensgene_v101").provider
    'UCSC'
    >>> lookup_annotation("sacCer3", "no_such_annotation") is None
    True
    """
    rows = annotation_table() if table is None else table
    return next(
        (row for row in rows if row.assembly == assembly and row.name == name),
        None,
    )


def list_annotation_metadata(
    assembly: str, *, table: Sequence[AnnotationMetadata] | None = None
) -> list[AnnotationMetadata]:
    """Return every annotation the table offers for ``assembly``, in table order.

    What the lab supports for an assembly — a different question from what is
    registered on this machine, which is :func:`~genome.io.gtf.list_annotations`.

    Parameters
    ----------
    assembly : str
        The assembly to list, matched against each record's ``assembly``.
    table : sequence of AnnotationMetadata, optional
        The rows to list from, as :func:`lookup_annotation` takes them.

    Returns
    -------
    list of AnnotationMetadata
        One record per listed annotation; empty for an assembly the table offers
        nothing for. A fresh list each call, so a caller may sort or filter it.

    Examples
    --------
    >>> [record.name for record in list_annotation_metadata("sacCer3")]
    ['ensgene_v101']
    >>> list_annotation_metadata("no_such_assembly")
    []
    """
    rows = annotation_table() if table is None else table
    return [row for row in rows if row.assembly == assembly]
