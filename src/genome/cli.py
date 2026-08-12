"""Command-line interface — a thin Typer wrapper over the genome API.

Logic lives in :mod:`genome.seq`, :mod:`genome.external`, etc.; this module
only translates arguments, dispatches, and chooses an output format.
"""

from __future__ import annotations

import json as _json

import typer

from genome import __version__ as _package_version
from genome.external import ToolNotFoundError
from genome.external import doctor as _doctor
from genome.io.download import assembly_table_row as _assembly_table_row
from genome.metadata import format_table_row as _format_table_row
from genome.seq import DNA

app = typer.Typer(help="Tools for handling genomic files.", no_args_is_help=True)


@app.command()
def version() -> None:
    """Print the installed package version."""
    typer.echo(_package_version)


@app.command()
def revcomp(
    sequence: str = typer.Argument(..., help="A DNA sequence over A/C/G/T (case is preserved)."),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Reverse-complement a DNA sequence.

    Exits with code 2 on invalid input.
    """
    # The DNA constructor no longer validates (too costly on large sequences),
    # so reject non-A/C/G/T characters here, at the I/O boundary.
    invalid = sorted({c for c in sequence if c.upper() not in "ACGT"})
    if invalid:
        typer.echo(
            f"error: sequence contains characters outside alphabet {{ACGT}}: {invalid!r}",
            err=True,
        )
        raise typer.Exit(code=2)

    result = DNA(sequence).reverse_complement()

    if json:
        typer.echo(_json.dumps({"input": sequence, "reverse_complement": str(result)}))
    else:
        typer.echo(str(result))


@app.command()
def doctor(
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Report availability and versions of required native tools.

    Exits with code 1 if any required tool is missing from PATH.
    """
    try:
        versions = _doctor()
    except ToolNotFoundError as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(versions))
    else:
        for name, ver in versions.items():
            typer.echo(f"{name}: {ver}")


# --- assembly commands -------------------------------------------------------


@app.command("table-row")
def table_row(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'sacCer3'."),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of the TSV line."),
) -> None:
    """Download an assembly and print its finished metadata table row.

    Fetches the FASTA — from the pinned source when the table has one for this
    assembly, otherwise from the UCSC golden path — unpacks it, and computes the sha256
    of the unpacked file. Prints one tab-separated line to paste into the shipped
    metadata table, or the same row as a JSON object.

    Exits with code 1 if the download fails or the FASTA does not match a checksum the
    table already pins.
    """
    try:
        row = _assembly_table_row(assembly, progressbar=not json)
    except (ValueError, OSError, RuntimeError) as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    typer.echo(_json.dumps(row) if json else _format_table_row(row))
