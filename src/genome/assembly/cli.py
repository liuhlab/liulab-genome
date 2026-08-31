"""The ``genome assembly`` sub-app — prepare a reference on disk, see it, and check it.

A thin Typer wrapper over :mod:`genome.assembly`: it translates arguments, dispatches to
the by-name module functions and chooses an output format. It ships from this package so
that what a command prints and what the result type holds change in one place.

Examples
--------
>>> from genome.assembly.cli import app
>>> [command.name for command in app.registered_commands]
['register', 'list', 'verify', 'table-row']
"""

from __future__ import annotations

import json as _json
from dataclasses import asdict as _asdict

import typer

from genome.assembly import COMPONENTS_UNCHANGED as _COMPONENTS_UNCHANGED
from genome.assembly import COMPONENTS_UNKNOWN as _COMPONENTS_UNKNOWN
from genome.assembly import EXPECTED_FROM_RECORD as _EXPECTED_FROM_RECORD
from genome.assembly import EXPECTED_FROM_TABLE as _EXPECTED_FROM_TABLE
from genome.assembly import AssemblyStatusRow as _AssemblyStatusRow
from genome.assembly import ChimeraDetails as _ChimeraDetails
from genome.assembly import VerifiedAssembly as _VerifiedAssembly
from genome.assembly import assembly_status as _assembly_status
from genome.assembly import assembly_table_row as _assembly_table_row
from genome.assembly import format_table_row as _format_table_row
from genome.assembly import register_assembly as _register_assembly
from genome.assembly import verify_assembly as _verify_assembly

#: What a failed assembly command raises, in one place. Every one of them is already
#: actionable — a checksum mismatch, a registration that cannot be trusted, a missing
#: native tool — so the CLI prints the message and exits non-zero rather than adding to
#: it. ``RegistrationError`` and ``ToolNotFoundError`` are ``RuntimeError``s;
#: ``ChecksumMismatchError`` is a ``ValueError``; a failed download is an ``OSError``.
_ASSEMBLY_ERRORS = (ValueError, OSError, RuntimeError)

#: What answered *which digest should this FASTA have?*, one sentence per answer. Three
#: states and three sentences: the curated row pinned it, this machine's own registration
#: recorded it, or nothing pinned it at all. One wording covering two of them would report
#: the weakest result in the words of the strongest. ``{assembly}`` is the name asked about.
#: Keyed on the constants the API answers with, never on the strings spelled again here:
#: the render falls back to the raw status, so a re-spelling would print it silently.
_EXPECTED_SENTENCES: dict[str | None, str] = {
    _EXPECTED_FROM_TABLE: "matches the digest the metadata table pins for it",
    _EXPECTED_FROM_RECORD: (
        "matches the digest {assembly}'s own registration recorded — the table pins none "
        "for it, so this is what this machine last produced and not an independent pin"
    ),
    None: (
        "({assembly} pins no digest and no record here holds one, so there was nothing to "
        "check it against)"
    ),
}

#: What the component check proved, one sentence per answer. Both print: a chimera whose
#: components could not be compared is unproven rather than passed, and a surface that
#: says nothing in that case says exactly what it says when everything agreed. Keyed on
#: the API's own constants, for the reason :data:`_EXPECTED_SENTENCES` is.
_COMPONENT_SENTENCES: dict[str, str] = {
    _COMPONENTS_UNCHANGED: (
        "unchanged — every component is still the one this chimera was built from"
    ),
    _COMPONENTS_UNKNOWN: (
        "unknown — a digest was missing on one side or the other, so nothing was actually "
        "compared; this is unproven rather than a pass"
    ),
}

#: The escape hatch the listing closes with, in the one place the words live. Fixed text
#: rather than a property of the report, since it reads nothing off it: the metadata table
#: is a cross-reference and never an allow-list, so any UCSC assembly name registers from
#: the golden path. What such a registration does not get is a pinned digest, which is the
#: fact a reader choosing between a listed assembly and an unlisted one needs.
_UNLISTED_SENTENCE = (
    "an assembly the table does not list registers too, from the UCSC golden path — with "
    "no pinned checksum behind it"
)

app = typer.Typer(
    help="Prepare an assembly on disk, check one, and compute its metadata row.",
    no_args_is_help=True,
)


@app.command("register")
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
    """Prepare an assembly on disk: fetch or build, verify, index, and record it.

    Downloads the FASTA from the source the metadata table pins (or the UCSC golden
    path for an assembly no row lists), checks it against the pinned sha256, derives
    the `.fai`, `.2bit` and `chrom.sizes`, and writes the registration record that says
    all of it finished. An assembly that is already registered is reported from its
    record without fetching anything.

    **Naming a chimera builds it**, and there is no flag for listing its parts: an
    assembly whose name is its component assemblies sorted and joined by `_` —
    `ce11_ecHT115` — is concatenated from components already prepared here, merged
    annotation included, and nothing is downloaded. The name carries the fact, so typing
    the components in the wrong order is refused by naming the canonical spelling, and a
    component this machine has not prepared is refused by naming the command that
    prepares it. What an assembly already registered here *is* comes from its record and
    never from its name, so a plain `hg38_mm10` seeded years ago stays what it was.

    Exits with code 1 when the directory holds a registration that cannot be trusted —
    files with no record, or a record that disagrees with what is on disk — and when a
    chimera cannot be built from this machine's parts. Re-run with `--force` to repair a
    directory: an unpacked FASTA that still matches the pinned checksum is kept and only
    the derived files are rebuilt, a chimera is concatenated again from its components,
    and neither is a way past the checks above.
    """
    try:
        registered = _register_assembly(assembly, source=source, force=force, progressbar=not json)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(registered.as_json()))
        return
    # The record answers whether this is a chimera, as it does everywhere else — and what
    # the registration answered with carries that record, so the fact is read from it
    # rather than from disk a second time. One read, one fact: the two surfaces cannot
    # disagree about what just happened.
    chimera = registered.chimera
    typer.echo(f"registered {registered.assembly} in {registered.directory}")
    if chimera is None:
        typer.echo(f"  source  {registered.source_url}")
    else:
        # A chimera was fetched from nowhere, so the source line would read `None`. What
        # it was made of is the fact that line was reaching for.
        typer.echo(f"  components  {', '.join(chimera.components)}")
    typer.echo(f"  sha256  {registered.sha256}")
    typer.echo(f"  files   {', '.join(registered.file_names)}")
    if chimera is not None:
        typer.echo(f"  {_merged_annotation_summary(chimera)}")


def _merged_annotation_summary(details: _ChimeraDetails) -> str:
    """Return the closing line naming the merged annotation this same command registered.

    A chimera's annotation is written by its build rather than by a second command, so the
    registration that just happened has one to name — and a build with nothing to merge
    says so, since silence there would read as an annotation nobody mentioned.
    """
    merged = details.merged_annotation
    if merged is None:
        return (
            "annotation  none — no component carries a default annotation, so none was "
            "merged rather than an empty one registered"
        )
    return f"annotation  {merged} — the components' own, merged and registered by this build"


# Named for what it does rather than for the command it serves, which would shadow the
# builtin ``list`` this module's own annotations are written in terms of.
@app.command("list")
def list_assemblies(
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """List the assemblies the metadata table offers against the ones prepared here.

    Two questions, two answers, side by side: which assemblies the lab pins a source and
    a checksum for, and which are actually prepared on this machine — including ones no
    row lists, which is what a registration from the UCSC golden path is. It takes no
    assembly name, since neither half has one to scope it.

    A directory in the assembly tree with no registration record beside it reads as
    `here, not registered`: nothing vouches for what is in it, and reading it as an
    assembly would be a claim this command has not checked. Registering it again from
    scratch with `--force` is what makes it trustworthy.

    **Whether a registration is still intact is not asked here.** That is `genome assembly
    verify`, which re-reads the FASTA and recomputes its digest; doing it cheaply enough
    to belong on this command would report *unchecked* in the words of *checked*. This one
    reads the table, the directory names and one small record each, so it stays instant.

    Nothing is downloaded, prepared or built to answer, and exit is `0` on a machine where
    nothing is registered — which is a fresh install's ordinary state, and is printed as
    such rather than as an error.
    """
    status = _assembly_status()

    if json:
        typer.echo(_json.dumps(status.as_json()))
        return

    rows = status.assemblies
    typer.echo(f"assemblies in {status.directory}")
    # Every fact below the heading comes off the report: which state a row is in, what the
    # table knows about it, and the closing lines. What is left here is the drawing — the
    # column widths, and dropping a column the table has nothing in.
    name_width = max((len(row.assembly_name) for row in rows), default=0)
    state_width = max((len(row.state) for row in rows), default=0)
    for row in rows:
        line = f"  {row.assembly_name:<{name_width}}  {row.state:<{state_width}}  {_names_of(row)}"
        typer.echo(line.rstrip())
    typer.echo(status.summary)
    if status.unregistered_note is not None:
        typer.echo(status.unregistered_note)
    typer.echo(_UNLISTED_SENTENCE)


def _names_of(row: _AssemblyStatusRow) -> str:
    """Return the table's names for one row — what a reader chooses between assemblies by.

    Whatever of the species and the NCBI name that row carries, since a row the table does
    not list carries neither and a **Chimera**'s row carries only its name.
    """
    return " ".join(part for part in (row.species, row.ncbi_name) if part)


@app.command("verify")
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
    """Re-read an assembly's FASTA and check its sha256 against the digest expected of it.

    Registering an assembly and reopening it go by presence and size, which is what
    makes them instant. This is the deliberate re-verification for when integrity is
    actually in doubt: it reads the whole file and computes its digest.

    What that digest is held to is the metadata table's pin, and failing that the one
    this assembly's own registration recorded — three states with three sentences, since
    being held to the lab's pin and being held to what this machine last produced are
    different results.

    A **chimera** is also checked against its components, which is the one failure no
    digest of its own bytes can show: those bytes do not change when a component is
    registered again underneath them. That line always prints for a chimera, whether it
    proved the components unchanged or found nothing to compare, so silence is never read
    as a pass.

    Exits with code 1 when the digest is not the one expected, when a component is no
    longer the one this was built from, when there is nothing registered to verify, or
    when the assembly's directory cannot be trusted.
    """
    try:
        checked = _verify_assembly(assembly, fasta=fasta)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(checked.as_json()))
        return
    typer.echo(f"{checked.fasta}: sha256 {checked.sha256} {_digest_summary(checked)}")
    components = checked.components
    if components is not None:
        typer.echo(f"  components  {_COMPONENT_SENTENCES.get(components, components)}")


def _digest_summary(checked: _VerifiedAssembly) -> str:
    """Return the sentence saying what a verified FASTA's digest was held to."""
    sentence = _EXPECTED_SENTENCES.get(checked.expected_from, _EXPECTED_SENTENCES[None])
    return sentence.format(assembly=checked.assembly)


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
    Use ``genome assembly verify`` to check a FASTA you already hold against the official
    row.

    A **chimera** is refused, and refused before anything is downloaded: it pins nothing,
    so this command has no job to do for one. The refusal describes the row — the name,
    and every other column blank — rather than printing it, because a row this command
    printed would look like one it had computed something for.

    Exits with code 1 if the download fails, or if the assembly named is a chimera.
    """
    try:
        computed = _assembly_table_row(assembly, progressbar=not json)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    # One row, rendered two ways: the metadata module owns the column order for both,
    # since the JSON object and the line to paste carry the same fields.
    row = _asdict(computed)
    typer.echo(_json.dumps(row) if json else _format_table_row(row))
