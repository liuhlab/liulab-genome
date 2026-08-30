"""The curated annotation table — what the lab supports for each **Assembly**.

One small hand-maintained TSV ships inside the package, and it is the assembly table's
idea one level down. ``data/annotation_metadata.tsv`` is keyed by assembly plus the
**Registered name** an annotation is addressed by, and says who publishes that annotation,
which release it is, where to fetch it and what the unpacked **GTF** hashes to — so naming
one is enough to register it. :func:`lookup_annotation` answers for one name and
:func:`list_annotation_metadata` for everything an assembly offers.

The assembly half of the pair is :mod:`genome.assembly.metadata`, and the two share their
shape and nothing else: no lookup here reads that table and none there reads this one, so
each sits with the context whose callers ask it questions.

The table is not an allow-list: a row means *officially supported* — pinned source, pinned
checksum — and an annotation with no row is perfectly legal, registered by path instead of
by name (ADR-0003).

:class:`AnnotationMetadata` declares its field list once: the table is read through those
fields column by column, parsed by each field's own declared type, and a whole record is
what the registration functions take to override the table. It is a **Shipped table** read
by the one loader every such table goes through, :mod:`genome.shipped`, which owns the
header validation, the cell parsing and the shape of the error a broken file raises.

Examples
--------
>>> from genome.annotation.metadata import lookup_annotation
>>> lookup_annotation("sacCer3", "ensgene_v101").provider
'UCSC'
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from typing import Any, get_type_hints

from genome.shipped import MetadataRowError, ShippedTable, parse_cell

__all__ = [
    "ANNOTATION_FIELDS",
    "AnnotationMetadata",
    "MetadataRowError",
    "annotation_table",
    "list_annotation_metadata",
    "lookup_annotation",
]

#: Location of the curated annotation metadata table within the package.
_ANNOTATION_RESOURCE = "data/annotation_metadata.tsv"

#: What repairs it: hand-maintained rather than generated, so the repair is an edit and the
#: message says which file to make it in.
_ANNOTATION_REPAIR = f"edit {_ANNOTATION_RESOURCE}, which is maintained by hand"


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
            **{name: parse_cell(name, row, _ANNOTATION_FIELD_TYPES) for name in ANNOTATION_FIELDS}
        )


#: Each annotation field's declared type, which parses that field's column of the table.
_ANNOTATION_FIELD_TYPES: dict[str, Any] = get_type_hints(AnnotationMetadata)

#: The annotation field names, in table-column order.
ANNOTATION_FIELDS: tuple[str, ...] = tuple(_ANNOTATION_FIELD_TYPES)

#: The annotation table as a **Shipped table**, declared exactly as the assembly one is.
_ANNOTATION_TABLE = ShippedTable(
    resource=_ANNOTATION_RESOURCE,
    columns=ANNOTATION_FIELDS,
    noun="annotation metadata table",
    repair=_ANNOTATION_REPAIR,
    error=MetadataRowError,
    identify=("assembly", "name"),
)


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

    Raises
    ------
    MetadataRowError
        As :func:`assembly_table` raises it, for that table's own columns.

    Examples
    --------
    >>> "ensgene_v101" in {record.name for record in annotation_table()}
    True
    """
    return tuple(AnnotationMetadata.from_row(row) for row in _ANNOTATION_TABLE.read().mappings())


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
    registered on this machine, which is :func:`~genome.annotation.registry.list_annotations`.

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
