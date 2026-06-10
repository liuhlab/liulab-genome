"""Assembly metadata — the curated table mapping assemblies to their identifiers.

A small hand-maintained TSV (``data/assembly_metadata.tsv``) records, for each
known reference assembly, its canonical name and the cross-references used to
talk about it across databases: species, UCSC name, NCBI name, NCBI assembly
accession, and NCBI taxonomy id. :func:`lookup_assembly` resolves a UCSC
assembly name (the identifier :class:`~genome.genome.Genome` is built from) to
its :class:`AssemblyMetadata` record, or ``None`` when the assembly is not in
the table.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from functools import cache
from importlib.resources import files

import pandas as pd

#: Location of the curated metadata table within the package.
_METADATA_RESOURCE = "data/assembly_metadata.tsv"


@dataclass(frozen=True)
class AssemblyMetadata:
    """Identifiers for one reference assembly (one row of the metadata table)."""

    assembly_name: str
    species: str
    ucsc_name: str
    ncbi_name: str
    ncbi_assembly_id: str
    ncbi_taxid: int


@cache
def _metadata_table() -> pd.DataFrame:
    """Load and cache the curated assembly metadata table as a DataFrame."""
    resource = files("genome").joinpath(_METADATA_RESOURCE)
    with resource.open("r", encoding="utf-8") as handle:
        return pd.read_csv(handle, sep="\t", dtype={"ncbi_taxid": "int64"})


@cache
def lookup_assembly(assembly: str) -> AssemblyMetadata | None:
    """Return the :class:`AssemblyMetadata` for a UCSC (or canonical) assembly name, or ``None``."""
    table = _metadata_table()
    match = table[(table["ucsc_name"] == assembly) | (table["assembly_name"] == assembly)]
    if match.empty:
        return None
    row = match.iloc[0]
    return AssemblyMetadata(
        assembly_name=str(row["assembly_name"]),
        species=str(row["species"]),
        ucsc_name=str(row["ucsc_name"]),
        ncbi_name=str(row["ncbi_name"]),
        ncbi_assembly_id=str(row["ncbi_assembly_id"]),
        ncbi_taxid=int(row["ncbi_taxid"]),
    )


#: The metadata field names, in table-column order — the kwargs ``Genome`` mirrors.
METADATA_FIELDS: tuple[str, ...] = tuple(f.name for f in fields(AssemblyMetadata))
