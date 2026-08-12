"""Assembly metadata — the curated table mapping assemblies to their identifiers.

A small hand-maintained TSV (``data/assembly_metadata.tsv``) records, for each
known reference assembly, its canonical name and the cross-references used to
talk about it across databases: species, UCSC name, NCBI name, NCBI assembly
accession, and NCBI taxonomy id. A row may additionally pin where that assembly's
FASTA is fetched from and the sha256 of the **unpacked** FASTA it yields, which is
what makes preparing the assembly reproducible. :func:`lookup_assembly` resolves a
UCSC assembly name (the identifier :class:`~genome.genome.Genome` is built from) to
its :class:`AssemblyMetadata` record, or ``None`` when the assembly is not in the
table.

:class:`AssemblyMetadata` declares the field list once: the table is read
through it column by column, and a whole record is what
:class:`~genome.genome.Genome` takes to override the table.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any, get_args, get_type_hints

import pandas as pd

#: Location of the curated metadata table within the package.
_METADATA_RESOURCE = "data/assembly_metadata.tsv"


@dataclass(frozen=True)
class AssemblyMetadata:
    """Identifiers for one reference assembly (one row of the metadata table).

    The single declaration of what an assembly's metadata consists of: the table
    is parsed through these fields, and a complete record is what
    :class:`~genome.genome.Genome` accepts in place of the table's own row.

    The last two fields are what makes preparing an assembly reproducible.
    ``source_url`` pins where its FASTA is fetched from, so nothing has to be derived
    or guessed; ``sha256`` pins the digest of the **unpacked** FASTA that source
    yields — not of the compressed archive it arrives in, so a copy taken from a
    mirror or recompressed elsewhere still matches (ADR-0006). Either may be ``None``:
    the table fills in over time, and a row with no digest is unverified rather than
    wrong.

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
    species: str
    ucsc_name: str
    ncbi_name: str
    ncbi_assembly_id: str
    ncbi_taxid: int
    source_url: str | None = None
    sha256: str | None = None


#: Each metadata field's declared type, which parses that field's column of the table.
_FIELD_TYPES: dict[str, Any] = get_type_hints(AssemblyMetadata)

#: The metadata field names, in table-column order — the columns every row carries.
METADATA_FIELDS: tuple[str, ...] = tuple(_FIELD_TYPES)


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
    'sacCer3\t\t\t\t\t559292\t\t'
    """
    return "\t".join("" if row.get(name) is None else str(row[name]) for name in METADATA_FIELDS)


@cache
def _metadata_table() -> pd.DataFrame:
    """Load and cache the curated assembly metadata table as a DataFrame of text."""
    resource = files("genome").joinpath(_METADATA_RESOURCE)
    with resource.open("r", encoding="utf-8") as handle:
        return pd.read_csv(handle, sep="\t", dtype=str)


def _parse_cell(name: str, row: pd.Series) -> Any:
    """Parse the ``name`` column of ``row`` with that field's own declared type.

    A field declared optional (``T | None``) is parsed by ``T`` when its cell carries
    text, and is ``None`` when the cell is blank or its column is absent — a union is
    not callable, so the type inside it does the parsing. A required field is parsed by
    its declared type, as every column was before any of them became optional.
    """
    declared = _FIELD_TYPES[name]
    inside = [arg for arg in get_args(declared) if arg is not type(None)]
    if not inside:
        return declared(row[name])
    # The whole table is read as text, so anything that is not a non-empty string
    # (a missing column, or the NaN pandas reads a blank cell as) is an empty cell.
    cell = row.get(name)
    text = cell.strip() if isinstance(cell, str) else ""
    return inside[0](text) if text else None


@cache
def lookup_assembly(assembly: str) -> AssemblyMetadata | None:
    """Return the :class:`AssemblyMetadata` for a UCSC (or canonical) assembly name, or ``None``.

    Parameters
    ----------
    assembly : str
        The name to look up, matched against the table's ``ucsc_name`` and
        ``assembly_name`` columns.

    Returns
    -------
    AssemblyMetadata or None
        The row for ``assembly``, or ``None`` when the table does not list it.
        The table is a cross-reference, not an allow-list, so an unlisted
        assembly is legal and its identifiers are simply unknown. A blank cell in
        an optional column reads back as ``None``.

    Examples
    --------
    >>> lookup_assembly("hg38").ncbi_name
    'GRCh38'
    >>> lookup_assembly("hg38").source_url
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    >>> lookup_assembly("no_such_assembly") is None
    True
    """
    table = _metadata_table()
    match = table[(table["ucsc_name"] == assembly) | (table["assembly_name"] == assembly)]
    if match.empty:
        return None
    row = match.iloc[0]
    # Every column is read as text; each field's declared type parses its own column.
    return AssemblyMetadata(**{name: _parse_cell(name, row) for name in METADATA_FIELDS})
