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
from genome.io.download import register_assembly as _register_assembly
from genome.io.download import verify_assembly as _verify_assembly
from genome.io.gtf import annotation_status as _annotation_status
from genome.io.gtf import register_annotation as _register_annotation
from genome.metadata import format_table_row as _format_table_row
from genome.seq import DNA

#: What a failed assembly command raises, in one place. Every one of them is already
#: actionable — a checksum mismatch, a registration that cannot be trusted, a missing
#: native tool — so the CLI prints the message and exits non-zero rather than adding to
#: it. ``RegistrationError`` and ``ToolNotFoundError`` are ``RuntimeError``s;
#: ``ChecksumMismatchError`` is a ``ValueError``; a failed download is an ``OSError``.
_ASSEMBLY_ERRORS = (ValueError, OSError, RuntimeError)

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


@app.command()
def register(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'hg38'."),
    source: str | None = typer.Option(
        None,
        "--source",
        help="Seed from this FASTA (local path or http(s)/ftp/sftp URL) instead of the "
        "source the metadata table pins.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Register again from scratch — the repair for a directory that raises.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Prepare an assembly on disk: fetch, verify, index, and record it.

    Downloads the FASTA from the source the metadata table pins (or the UCSC golden
    path for an assembly no row lists), checks it against the pinned sha256, derives
    the `.fai`, `.2bit` and `chrom.sizes`, and writes the registration record that says
    all of it finished. An assembly that is already registered is reported from its
    record without fetching anything.

    Exits with code 1 when the directory holds a registration that cannot be trusted —
    files with no record, or a record that disagrees with what is on disk. Re-run with
    `--force` to repair it: an unpacked FASTA that still matches the pinned checksum is
    kept and only the derived files are rebuilt.
    """
    try:
        payload = _register_assembly(assembly, source=source, force=force, progressbar=not json)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(payload))
        return
    claimed = payload["files"]
    names = sorted(claimed) if isinstance(claimed, dict) else []
    typer.echo(f"registered {payload['assembly']} in {payload['directory']}")
    typer.echo(f"  source  {payload['source_url']}")
    typer.echo(f"  sha256  {payload['sha256']}")
    typer.echo(f"  files   {', '.join(names)}")


@app.command("register-annotation")
def register_annotation(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'hg38'."),
    name: str = typer.Argument(..., help="Registered annotation name, e.g. 'gencode_v50'."),
    force: bool = typer.Option(
        False,
        "--force",
        help="Register again from scratch — the repair for a directory that raises.",
    ),
    check_chromosomes: bool = typer.Option(
        True,
        "--check-chromosomes/--no-check-chromosomes",
        help="Refuse a GTF naming sequences this assembly does not carry, before paying "
        "for the database build. Pass --no-check-chromosomes to register one whose "
        "mismatch you have inspected and accept.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Register one of an assembly's annotations by name: fetch, verify, build, record.

    Downloads the GTF from the URL the annotation table pins for this assembly, checks
    the unpacked file against the pinned sha256 and its chromosome names against the
    assembly's, builds the gffutils database, and writes the registration record that
    says all of it finished. An annotation that is already registered is reported from
    its record without fetching anything.

    The chromosome check needs the assembly's ``chrom.sizes``, so it can only run once
    the assembly itself is registered. Registering an annotation first is allowed and
    reports ``chromosomes not checked``; register the assembly first to have the names
    verified.

    Exits with code 1 when the table lists no such annotation, when the GTF is not the
    digest pinned for it, when it names chromosomes the assembly does not carry, or when
    the directory holds a registration that cannot be trusted — files with no record, or
    a record that disagrees with what is on disk. Re-run with `--force` to repair it.
    """
    try:
        payload = _register_annotation(
            assembly,
            name,
            force=force,
            progressbar=not json,
            check_chromosomes=check_chromosomes,
        )
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(payload))
        return
    claimed = payload["files"]
    names = sorted(claimed) if isinstance(claimed, dict) else []
    typer.echo(f"registered {payload['name']} for {payload['assembly']} in {payload['directory']}")
    typer.echo(f"  source  {payload['source_url']}")
    typer.echo(f"  sha256  {payload['sha256']}")
    typer.echo(f"  files   {', '.join(names)}")
    # Whether the names were actually verified is not something to leave implicit: the
    # check cannot run before the assembly is registered, and silence would read as a pass.
    details = payload.get("details")
    checked = details.get("chromosomes_checked") if isinstance(details, dict) else None
    if checked is not True:
        typer.echo("  chromosomes not checked — register the assembly first to verify them")


# Named for what it does rather than for the command it serves: ``from __future__ import
# annotations`` already binds that name in this module.
@app.command("annotations")
def list_annotations(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'hg38'."),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """List what the annotation table offers for an assembly against what is registered here.

    Two questions, two answers, side by side: which annotations the lab supports for
    this assembly, and which are actually registered on this machine. The default
    annotation is named last, with the command that registers it when it is one of the
    ones this machine does not have — which is the ordinary state of a fresh install.

    Nothing is downloaded, prepared or built to answer this, so it works for an assembly
    that has never been registered here.
    """
    try:
        payload = _annotation_status(assembly)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(payload))
        return

    rows = payload["annotations"] if isinstance(payload["annotations"], list) else []
    typer.echo(f"annotations for {payload['assembly']} in {payload['directory']}")
    if not rows:
        typer.echo("  (the table offers none, and none is registered here)")
    name_width = max((len(str(row["name"])) for row in rows), default=0)
    state_width = max((len(_state(row)) for row in rows), default=0)
    for row in rows:
        provider = f"  {row['provider']} {row['version']}" if row["offered"] else ""
        line = f"  {row['name']!s:<{name_width}}  {_state(row):<{state_width}}{provider}"
        typer.echo(line.rstrip())
    typer.echo(_default_line(payload, rows))


def _state(row: dict[str, object]) -> str:
    """Return which of the two questions a row answers: offered, registered, or both."""
    if not row["offered"]:
        return "registered, not offered"
    return "registered" if row["registered"] else "offered, not registered"


def _default_line(payload: dict[str, object], rows: list[dict[str, object]]) -> str:
    """Return the closing line naming the default annotation, and how to get it if absent."""
    default = payload["default_annotation"]
    if default is None:
        return "default: (none)"
    registered = any(row["name"] == default and row["registered"] for row in rows)
    if registered:
        return f"default: {default}"
    return (
        f"default: {default} — not registered here; register it with "
        f"`genome register-annotation {payload['assembly']} {default}`"
    )


@app.command()
def verify(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'sacCer3'."),
    fasta: str | None = typer.Option(
        None,
        "--fasta",
        help="Check this FASTA instead of the assembly's registered one — a copy from a "
        "mirror or handed over by hand, checkable before anything is built on it.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Re-read an assembly's FASTA and check its sha256 against the table's row.

    Registering an assembly and reopening it go by presence and size, which is what
    makes them instant. This is the deliberate re-verification for when integrity is
    actually in doubt: it reads the whole file and computes its digest.

    Exits with code 1 when the digest is not the one the row pins, when there is
    nothing registered to verify, or when the assembly's directory cannot be trusted.
    """
    try:
        payload = _verify_assembly(assembly, fasta=fasta)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(payload))
        return
    tail = (
        "matches the digest pinned for it"
        if payload["verified"]
        else f"({assembly} pins no digest, so there was nothing to check it against)"
    )
    typer.echo(f"{payload['fasta']}: sha256 {payload['sha256']} {tail}")


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

    A checksum the table already pins is **reported, never enforced** — this is the
    command to reach for when an upstream file has legitimately changed and the pin has
    to be regenerated, so refusing on a mismatch would refuse exactly when it is needed.
    Use ``genome verify`` to check a FASTA you already hold against the official row.

    Exits with code 1 if the download fails.
    """
    try:
        row = _assembly_table_row(assembly, progressbar=not json)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    typer.echo(_json.dumps(row) if json else _format_table_row(row))
