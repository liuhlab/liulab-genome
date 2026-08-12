"""Assembly metadata — the curated table mapping assemblies to their identifiers.

A small hand-maintained TSV (``data/assembly_metadata.tsv``) records, for each
known reference assembly, its canonical name and the cross-references used to
talk about it across databases: species, UCSC name, NCBI name, NCBI assembly
accession, and NCBI taxonomy id. :func:`lookup_assembly` resolves a UCSC
assembly name (the identifier :class:`~genome.genome.Genome` is built from) to
its :class:`AssemblyMetadata` record, or ``None`` when the assembly is not in
the table.

:class:`AssemblyMetadata` declares the field list once: the table is read
through it column by column, and a whole record is what
:class:`~genome.genome.Genome` takes to override the table.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any, get_type_hints

import pandas as pd

#: Location of the curated metadata table within the package.
_METADATA_RESOURCE = "data/assembly_metadata.tsv"


@dataclass(frozen=True)
class AssemblyMetadata:
    """Identifiers for one reference assembly (one row of the metadata table).

    The single declaration of what an assembly's metadata consists of: the table
    is parsed through these fields, and a complete record is what
    :class:`~genome.genome.Genome` accepts in place of the table's own row.

    Examples
    --------
    >>> record = AssemblyMetadata(
    ...     "hg38", "Homo sapiens", "hg38", "GRCh38", "GCF_000001405.40", 9606
    ... )
    >>> record.species
    'Homo sapiens'
    """

    assembly_name: str
    species: str
    ucsc_name: str
    ncbi_name: str
    ncbi_assembly_id: str
    ncbi_taxid: int


#: Each metadata field's declared type, which parses that field's column of the table.
_FIELD_TYPES: dict[str, Any] = get_type_hints(AssemblyMetadata)

#: The metadata field names, in table-column order — the columns every row carries.
METADATA_FIELDS: tuple[str, ...] = tuple(_FIELD_TYPES)


@cache
def _metadata_table() -> pd.DataFrame:
    """Load and cache the curated assembly metadata table as a DataFrame of text."""
    resource = files("genome").joinpath(_METADATA_RESOURCE)
    with resource.open("r", encoding="utf-8") as handle:
        return pd.read_csv(handle, sep="\t", dtype=str)


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
        assembly is legal and its identifiers are simply unknown.

    Examples
    --------
    >>> lookup_assembly("hg38").ncbi_name
    'GRCh38'
    >>> lookup_assembly("no_such_assembly") is None
    True
    """
    table = _metadata_table()
    match = table[(table["ucsc_name"] == assembly) | (table["assembly_name"] == assembly)]
    if match.empty:
        return None
    row = match.iloc[0]
    # Every column is read as text; each field's declared type parses its own column.
    return AssemblyMetadata(**{name: _FIELD_TYPES[name](row[name]) for name in METADATA_FIELDS})
