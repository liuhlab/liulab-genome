"""Command-line interface — a thin Typer wrapper over the genome API.

Logic lives in :mod:`genome.seq`, :mod:`genome.external`, etc.; this module
only translates arguments, dispatches, and chooses an output format.
"""

from __future__ import annotations

import json as _json
from dataclasses import asdict as _asdict

import typer

from genome import __version__ as _package_version
from genome.external import ToolNotFoundError
from genome.external import doctor as _doctor
from genome.io.chimera import COMPONENTS_UNCHANGED as _COMPONENTS_UNCHANGED
from genome.io.chimera import COMPONENTS_UNKNOWN as _COMPONENTS_UNKNOWN
from genome.io.chimera import ChimeraDetails as _ChimeraDetails
from genome.io.download import EXPECTED_FROM_RECORD as _EXPECTED_FROM_RECORD
from genome.io.download import EXPECTED_FROM_TABLE as _EXPECTED_FROM_TABLE
from genome.io.download import VerifiedAssembly as _VerifiedAssembly
from genome.io.download import assembly_table_row as _assembly_table_row
from genome.io.download import register_assembly as _register_assembly
from genome.io.download import verify_assembly as _verify_assembly
from genome.io.gtf import AnnotationStatus as _AnnotationStatus
from genome.io.gtf import AnnotationStatusRow as _AnnotationStatusRow
from genome.io.gtf import RegisteredAnnotation as _RegisteredAnnotation
from genome.io.gtf import annotation_status as _annotation_status
from genome.io.gtf import register_annotation as _register_annotation
from genome.io.gtf import register_annotation_by_path as _register_annotation_by_path
from genome.metadata import format_table_row as _format_table_row
from genome.seq import DNA

#: What a failed assembly command raises, in one place. Every one of them is already
#: actionable — a checksum mismatch, a registration that cannot be trusted, a missing
#: native tool — so the CLI prints the message and exits non-zero rather than adding to
#: it. ``RegistrationError`` and ``ToolNotFoundError`` are ``RuntimeError``s;
#: ``ChecksumMismatchError`` is a ``ValueError``; a failed download is an ``OSError``.
_ASSEMBLY_ERRORS = (ValueError, OSError, RuntimeError)

#: Help for the two feature-inference switches, in one place: both registration
#: commands offer them and a reader of either help text deserves the same explanation.
#: ``{feature}`` is the GTF feature being reconstructed.
_INFER_HELP = (
    "Reconstruct {feature} features from exon lines. Off by default: GENCODE, Ensembl "
    "and RefSeq GTFs declare them already and inferring is the slow path. Turn it on for "
    "a bare exon-level GTF — one whose only lines are exons — which otherwise registers "
    "as a database of exons and nothing else."
)

#: Help for the chromosome-name check, likewise shared by both registration commands.
_CHECK_CHROMOSOMES_HELP = (
    "Refuse a GTF naming sequences this assembly does not carry, before paying for the "
    "database build. Pass --no-check-chromosomes to register one whose mismatch you have "
    "inspected and accept — the annotation is built as it stands and the record says the "
    "check was stood down, so nothing later mistakes it for a verified one, and nothing "
    "mistakes it for an annotation registered before its assembly either."
)

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
        True, "--check-chromosomes/--no-check-chromosomes", help=_CHECK_CHROMOSOMES_HELP
    ),
    infer_genes: bool = typer.Option(
        False, "--infer-genes/--no-infer-genes", help=_INFER_HELP.format(feature="gene")
    ),
    infer_transcripts: bool = typer.Option(
        False,
        "--infer-transcripts/--no-infer-transcripts",
        help=_INFER_HELP.format(feature="transcript"),
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Register one of an assembly's annotations by name: fetch, verify, build, record.

    Downloads the GTF from the URL the annotation table pins for this assembly, checks
    the unpacked file against the pinned sha256 and its chromosome names against the
    assembly's, builds the gffutils database, and writes the registration record that
    says all of it finished. An annotation that is already registered is reported from
    its record without fetching anything. A GTF the table does not list is registered
    from a path instead, with `genome register-gtf`.

    The chromosome check needs the assembly's ``chrom.sizes``, so it can only run once
    the assembly itself is registered. Registering an annotation first is allowed and
    reports that there was nothing to check the names against, which registering the
    assembly first is what fixes. Standing the check down yourself reports that instead,
    and advises nothing: you meant it.

    Exits with code 1 when the table lists no such annotation, when the GTF is not the
    digest pinned for it, when it names chromosomes the assembly does not carry, or when
    the directory holds a registration that cannot be trusted — files with no record, or
    a record that disagrees with what is on disk. Re-run with `--force` to repair it.
    """
    try:
        registered = _register_annotation(
            assembly,
            name,
            force=force,
            progressbar=not json,
            check_chromosomes=check_chromosomes,
            disable_infer_genes=not infer_genes,
            disable_infer_transcripts=not infer_transcripts,
        )
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    _report_annotation(registered, json=json)


@app.command("register-gtf")
def register_gtf(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'hg38'."),
    gtf: str = typer.Argument(..., help="Path to the GTF to register, plain or .gz."),
    name: str = typer.Argument(
        ..., help="Registered name to address it by afterwards, e.g. 'wormbase_ws298'."
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Register again from scratch — the repair for a directory that raises.",
    ),
    check_chromosomes: bool = typer.Option(
        True, "--check-chromosomes/--no-check-chromosomes", help=_CHECK_CHROMOSOMES_HELP
    ),
    infer_genes: bool = typer.Option(
        False, "--infer-genes/--no-infer-genes", help=_INFER_HELP.format(feature="gene")
    ),
    infer_transcripts: bool = typer.Option(
        False,
        "--infer-transcripts/--no-infer-transcripts",
        help=_INFER_HELP.format(feature="transcript"),
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Register an annotation from a local GTF: place, check, build, record.

    The escape hatch for a GTF the annotation table does not list — `genome
    register-annotation` is the way in for one it does. Nothing is downloaded and no
    checksum is compared against, because an unlisted GTF has none pinned for it: the
    file you name is placed under `<assembly dir>/gtf/<name>/` (a `.gz` is decompressed
    on the way), its chromosome names are checked against the assembly's, the gffutils
    database is built, and the record that says all of it finished is written last. It
    is addressed by `<name>` from then on, exactly as a listed annotation is, and shows
    up in `genome annotations` as registered but not offered.

    The assembly is named rather than inferred — it is what says which reference these
    gene models describe — and naming it is also what lets the chromosome check find the
    assembly's `chrom.sizes` without being told where it is.

    Exits with code 1 when the GTF is not there, when it names chromosomes the assembly
    does not carry, or when the directory holds a registration that cannot be trusted.
    Re-run with `--force` to repair it.
    """
    try:
        registered = _register_annotation_by_path(
            assembly,
            gtf,
            name,
            force=force,
            check_chromosomes=check_chromosomes,
            disable_infer_genes=not infer_genes,
            disable_infer_transcripts=not infer_transcripts,
        )
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    _report_annotation(registered, json=json)


def _report_annotation(registered: _RegisteredAnnotation, *, json: bool) -> None:
    """Print a finished annotation registration, as JSON or as the human summary."""
    if json:
        typer.echo(_json.dumps(registered.as_json()))
        return
    typer.echo(f"registered {registered.name} for {registered.assembly} in {registered.directory}")
    typer.echo(f"  source  {registered.source_url}")
    typer.echo(f"  sha256  {registered.sha256}")
    typer.echo(f"  files   {', '.join(registered.file_names)}")
    # Whether the names were actually verified is not something to leave implicit, and it
    # is printed whichever way it went: silence would read as a pass. Which sentence that
    # is belongs to the record and to the API that reads it, not to this surface.
    typer.echo(f"  {registered.chromosome_check}")


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

    An annotation whose directory is here but cannot be trusted — files with no record,
    or a record that disagrees with what is on disk — reads as `broken` rather than as
    one nobody has fetched, and the line under it says what is wrong and names the
    command that repairs it. This is where such a thing is discovered, so it is reported
    and not raised over: exit is still `0`, and one broken annotation never hides the
    ones beside it.

    Nothing is downloaded, prepared or built to answer this, so it works for an assembly
    that has never been registered here.
    """
    try:
        payload = _annotation_status(assembly)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(payload.as_json()))
        return

    rows = payload.annotations
    typer.echo(f"annotations for {payload.assembly} in {payload.directory}")
    if not rows:
        typer.echo("  (the table offers none, and none is registered here)")
    name_width = max((len(row.name) for row in rows), default=0)
    state_width = max((len(_state(row)) for row in rows), default=0)
    for row in rows:
        provider = f"  {row.provider} {row.version}" if row.offered else ""
        line = f"  {row.name:<{name_width}}  {_state(row):<{state_width}}{provider}"
        typer.echo(line.rstrip())
        # Verbatim from the API, which is the same text re-registering it would print:
        # what is wrong, and the command that fixes it. One wording, two surfaces.
        if row.broken:
            typer.echo(f"      {row.problem}")
    typer.echo(_default_line(payload))


def _state(row: _AnnotationStatusRow) -> str:
    """Return the state a row is in: registered, broken, or merely offered.

    ``broken`` comes first because it is the one that needs acting on, and because a
    broken annotation is not registered — no record vouches for it — so reporting it as
    the absence of one would be true and useless.
    """
    if row.broken:
        return "broken"
    if not row.offered:
        return "registered, not offered"
    return "registered" if row.registered else "offered, not registered"


def _default_line(payload: _AnnotationStatus) -> str:
    """Return the closing line naming the default annotation, and how to get it if absent."""
    default = payload.default_annotation
    if default is None:
        return "default: (none)"
    row = payload.default_row
    if row is not None and row.broken:
        return f"default: {default} — broken here; repair it with `{row.repair}`"
    if row is not None and row.registered:
        return f"default: {default}"
    return (
        f"default: {default} — not registered here; register it with "
        f"`genome register-annotation {payload.assembly} {default}`"
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
    Use ``genome verify`` to check a FASTA you already hold against the official row.

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
