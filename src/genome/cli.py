"""Command-line interface — a thin Typer wrapper over the genome API.

Logic lives in :mod:`genome.seq`, :mod:`genome.external`, etc.; this module
only translates arguments, dispatches, and chooses an output format.
"""

from __future__ import annotations

import json as _json
from collections.abc import Iterator as _Iterator
from dataclasses import asdict as _asdict
from pathlib import Path as _Path
from typing import Any as _Any

import typer

from genome import __version__ as _package_version
from genome.external import ToolNotFoundError
from genome.external import doctor as _doctor
from genome.io.components import COMPONENTS_UNCHANGED as _COMPONENTS_UNCHANGED
from genome.io.components import COMPONENTS_UNKNOWN as _COMPONENTS_UNKNOWN
from genome.io.components import ChimeraDetails as _ChimeraDetails
from genome.io.download import assembly_table_row as _assembly_table_row
from genome.io.download import register_assembly as _register_assembly
from genome.io.download import verify_assembly as _verify_assembly
from genome.io.gtf import annotation_status as _annotation_status
from genome.io.gtf import gene_list as _gene_list
from genome.io.gtf import gene_lists as _gene_lists
from genome.io.gtf import register_annotation as _register_annotation
from genome.io.gtf import register_gtf as _register_gtf
from genome.io.gtf import tf_cofactor_list as _tf_cofactor_list
from genome.io.gtf import tf_gene_list as _tf_gene_list
from genome.io.results import EXPECTED_FROM_RECORD as _EXPECTED_FROM_RECORD
from genome.io.results import EXPECTED_FROM_TABLE as _EXPECTED_FROM_TABLE
from genome.io.results import GeneList as _GeneList
from genome.io.results import GeneListSource as _GeneListSource
from genome.io.results import RegisteredAnnotation as _RegisteredAnnotation
from genome.io.results import ResolvedStems as _ResolvedStems
from genome.io.results import ResolvedXrefIds as _ResolvedXrefIds
from genome.io.results import VerifiedAssembly as _VerifiedAssembly
from genome.metadata import format_table_row as _format_table_row
from genome.seq import DNA
from genome.tf.motif.background import BackgroundMode as _BackgroundMode
from genome.tf.motif.jaspar import DEFAULT_RELEASE as _DEFAULT_RELEASE
from genome.tf.motif.jaspar import DEFAULT_TAX_GROUP as _DEFAULT_TAX_GROUP
from genome.tf.motif.jaspar import JASPAR_RELEASES as _JASPAR_RELEASES
from genome.tf.motif.jaspar import JASPAR_TAX_GROUPS as _JASPAR_TAX_GROUPS
from genome.tf.motif.jaspar import JasparDatabase as _JasparDatabase
from genome.tf.motif.motif import DEFAULT_THRESHOLD as _DEFAULT_THRESHOLD
from genome.tf.motif.motif import MIN_MOTIF_LENGTH as _MIN_MOTIF_LENGTH
from genome.tf.motif.parquet import hit_count as _hit_count
from genome.tf.motif.parquet import provenance_of as _provenance_of
from genome.tf.motif.scan import read_fasta as _read_fasta
from genome.tf.motif.scan import scan_stream as _scan_stream
from genome.tf.motif.workers import resolve_workers as _resolve_workers
from genome.xref import NAMESPACES as _NAMESPACES
from genome.xref import XrefSet as _XrefSet

#: What a failed assembly command raises, in one place. Every one of them is already
#: actionable — a checksum mismatch, a registration that cannot be trusted, a missing
#: native tool — so the CLI prints the message and exits non-zero rather than adding to
#: it. ``RegistrationError`` and ``ToolNotFoundError`` are ``RuntimeError``s;
#: ``ChecksumMismatchError`` is a ``ValueError``; a failed download is an ``OSError``.
_ASSEMBLY_ERRORS = (ValueError, OSError, RuntimeError)

#: What a failed gene-list command raises, on top of the above — the gene-category pair
#: and both TF lists alike. An annotation that ships no curated gene list, one whose
#: list does not declare the category asked for, a species no census or cofactor table
#: ships for and an assembly nothing names a species for are lookups that found nothing
#: rather than bad values — so they are ``LookupError``s, which :data:`_ASSEMBLY_ERRORS`
#: does not cover, and so is the ``KeyError`` an unregistered annotation has always
#: raised. One list for all four commands because they are alike by construction: each
#: asks one registered annotation a question a shipped file answers, and each absence on
#: the way is a lookup.
_GENE_LIST_ERRORS = (*_ASSEMBLY_ERRORS, LookupError)

#: What a failed motif scan raises. The same three, and each already says what to do: a
#: release this package does not prepare, a threshold that is not a p-value and a file that
#: is not FASTA are ``ValueError``s; a FASTA that is not there and a release that cannot be
#: fetched are ``OSError``s; a process pool that died under the scan is a ``RuntimeError``.
#: Its own name rather than :data:`_ASSEMBLY_ERRORS` reused, because a motif belongs to no
#: assembly and the two lists are alike by coincidence rather than by construction.
_MOTIF_SCAN_ERRORS = (ValueError, OSError, RuntimeError)

#: What a failed xref hop raises, and every one of them already names its next action: a
#: species, a source or a **Namespace** that resolves to nothing is a ``LookupError``; a set
#: that could not be fetched, and a directory an interrupted download left unfinished, are
#: ``RuntimeError``s; a publisher's file that does not match its pin and a stored slice this
#: package did not write are ``ValueError``s; and a file that went away under the read is an
#: ``OSError``. Its own name rather than :data:`_GENE_LIST_ERRORS` reused, because an **Xref
#: set** belongs to no assembly and no annotation — the two lists are alike by coincidence
#: rather than by construction.
_XREF_ERRORS = (ValueError, OSError, RuntimeError, LookupError)

#: What to do when neither direction is named, or both are. Spelled once because the two
#: mistakes are one mistake, and it names both flags because which one is wanted is exactly
#: what the caller has not said.
_XREF_DIRECTION_HELP = (
    "error: name exactly one direction. `--to-stems NAMESPACE` reads the ids as that "
    "namespace and answers in gene id stems; `--from-stems NAMESPACE` reads them as stems "
    "and answers in that namespace. An identifier does not say which system it belongs to, "
    "so nothing here infers the direction from the strings you passed."
)

#: Help for the ``--annotation`` option every command asking one annotation a question takes.
_ANNOTATION_HELP = (
    "Ask about this registered annotation instead of the assembly's default one. An "
    "assembly with no default and none named has nothing to answer about, and says so."
)

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
def version(
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the installed package version."""
    if json:
        typer.echo(_json.dumps({"version": _package_version}))
        return
    typer.echo(_package_version)


@app.command()
def revcomp(
    sequence: str = typer.Argument(
        ...,
        # The third place this alphabet used to be spelled by hand, and the one a reader
        # meets first. Rendered from the type like the check and the error below it.
        help=f"A DNA sequence over {'/'.join(sorted(DNA.ALPHABET))} (case is preserved).",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Reverse-complement a DNA sequence.

    Exits with code 2 on invalid input.
    """
    # The DNA constructor no longer validates (too costly on large sequences), so reject
    # non-alphabet characters here, at the I/O boundary. Both halves of that ask the type —
    # which characters offend, and what to call the alphabet they offended against — because
    # an edge that spells `ACGT` itself is a second copy of `DNA.ALPHABET` that drifts from it
    # silently.
    invalid = DNA.outside_alphabet(sequence)
    if invalid:
        alphabet = "".join(sorted(DNA.ALPHABET))
        typer.echo(
            f"error: sequence contains characters outside alphabet {{{alphabet}}}: {invalid!r}",
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
        registered = _register_gtf(
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
        status = _annotation_status(assembly)
    except _ASSEMBLY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(status.as_json()))
        return

    rows = status.annotations
    typer.echo(f"annotations for {status.assembly} in {status.directory}")
    if not rows:
        typer.echo("  (the table offers none, and none is registered here)")
    # Every word below the heading comes off the report: which state a row is in, what is
    # wrong with a broken one, and the closing line about the default. This chooses the
    # column widths and nothing else.
    name_width = max((len(row.name) for row in rows), default=0)
    state_width = max((len(row.state) for row in rows), default=0)
    for row in rows:
        provider = f"  {row.provider} {row.version}" if row.offered else ""
        line = f"  {row.name:<{name_width}}  {row.state:<{state_width}}{provider}"
        typer.echo(line.rstrip())
        if row.broken:
            typer.echo(f"      {row.problem}")
    typer.echo(status.default_summary)


@app.command("gene-list")
def gene_list(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'ce11'."),
    category: str = typer.Argument(..., help="Gene category, e.g. 'rRNA'."),
    annotation: str | None = typer.Option(None, "--annotation", help=_ANNOTATION_HELP),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the gene ids an annotation puts in one category, one per line.

    The genes come from a curated gene list shipped inside this package, not from the
    GTF's own biotype attribute — which is spelled two ways across four publishers, sorts
    genes by three taxonomies that do not agree, and is missing altogether from some
    annotations, so a caller deriving categories itself reports none for those and never
    finds out.

    Only the ids go to stdout, so the output pipes: the heading and the per-source
    attribution go to stderr. For a merged annotation there is one source per contributing
    component, and the ids are concatenated in that order and never de-duplicated — two
    components carrying one gene id is a real ambiguity, not something to collapse here.
    `--json` carries the same answer with the sources kept apart.

    `genome gene-categories <assembly>` says which categories may be asked for; they are
    the curated list's to declare and differ between annotations.

    Exits with code 1 when the annotation is not registered here, when no curated gene
    list ships for it, and when it declares categories but not this one — three different
    facts, each with its own message, and none of them an empty list of genes.
    """
    try:
        listed = _gene_list(assembly, category, annotation=annotation)
    except _GENE_LIST_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(listed.as_json()))
        return
    # Attribution to stderr, ids to stdout: a bare id list is what a shell pipeline wants,
    # and who contributed what is what a reader wants. Neither has to cost the other.
    typer.echo(f"{listed.category} for {listed.assembly} / {listed.annotation}", err=True)
    for source in listed.sources:
        typer.echo(f"  {_source_label(source)}  {len(source.gene_ids)}", err=True)
    for gene_id in listed.gene_ids:
        typer.echo(gene_id)


@app.command("gene-categories")
def gene_categories(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'ce11'."),
    annotation: str | None = typer.Option(None, "--annotation", help=_ANNOTATION_HELP),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """List the gene categories an annotation declares, with a gene count for each.

    What may be asked of `genome gene-list`, since the categories are a property of the
    annotation rather than of this package: human GENCODE splits rRNA pseudogenes out and
    mouse does not, yeast carries a precursor category nothing else has, and a bacterium
    has no mitochondrial category at all.

    A merged annotation shows the per-component split beside each count, so a category one
    component declares and another does not is visible as such rather than as a smaller
    number. `--json` carries every category with its gene ids and its sources — the same
    answer `genome gene-list` gives for one of them, for all of them at once.

    Exits with code 1 when the annotation is not registered here, and when no curated gene
    list ships for it — which is not the same answer as its declaring no genes, and is why
    an empty list is never printed.
    """
    try:
        answers = _gene_lists(assembly, annotation=annotation)
    except _GENE_LIST_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps([answer.as_json() for answer in answers]))
        return
    first = answers[0]
    typer.echo(f"categories for {first.assembly} / {first.annotation}")
    name_width = max(len(answer.category) for answer in answers)
    count_width = max(len(str(len(answer.gene_ids))) for answer in answers)
    for answer in answers:
        typer.echo(_category_row(answer, name_width=name_width, count_width=count_width))


def _source_label(source: _GeneListSource) -> str:
    """Return who contributed a category's genes: the component and its annotation, or it alone."""
    if source.component is None:
        return source.annotation
    return f"{source.component}: {source.annotation}"


def _category_row(answer: _GeneList, *, name_width: int, count_width: int) -> str:
    """Return one category's line: its name, its gene count, and the per-component split.

    The split is printed only where there is one to print — a merged annotation's — since
    for any other the single component would repeat the annotation already in the heading.
    """
    line = f"  {answer.category:<{name_width}}  {len(answer.gene_ids):>{count_width}}"
    split = ", ".join(
        f"{source.component}: {len(source.gene_ids)}"
        for source in answer.sources
        if source.component is not None
    )
    return f"{line}  ({split})" if split else line


@app.command("tf-gene-list")
def tf_gene_list(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'hg38'."),
    annotation: str | None = typer.Option(None, "--annotation", help=_ANNOTATION_HELP),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the gene ids a published census judges transcription factors, one per line.

    Nothing here decides what a transcription factor is. The verdict is the census's —
    Lambert et al. 2018 for human, AnimalTFDB 4.0 for mouse — and which one spoke is
    printed beside the answer, since citing it is the condition on shipping it. The
    species comes from the assembly's own metadata row and is never passed in, so asking
    for human transcription factors while holding a mouse assembly is not expressible.

    A census is keyed by gene id stems — gene ids with the version suffix dropped — and a
    registered annotation is not, so every stem is resolved into the ids that annotation
    actually spells and the output joins to a counts matrix with nothing left to
    normalise. A stem naming two genes prints both rather than one of them.

    Only the ids go to stdout, so the output pipes: the heading, the census's attribution
    and the counts go to stderr, the last of them saying how many stems this annotation
    carries no gene for. `--json` carries the whole record — every gene with the census's
    own assessment and DBD family, the provenance to cite, and those unresolved stems.

    Assessed-positive genes only, and there is no flag to widen it: a gene the census
    assessed and turned down is a verdict too, but a bare id list has nowhere to say which
    of the two an id is, and a pipeline would read the rejected ones as transcription
    factors. `Genome(<assembly>).tf_gene_list(include_rejected=True)` is where that answer
    is expressible, because there each id travels with the verdict reached on it.

    Exits with code 1 when the annotation is not registered here, when no census ships for
    the assembly's species, and when nothing says what species the assembly is — three
    different facts, each with its own message, and none of them an empty list of genes.
    """
    try:
        listed = _tf_gene_list(assembly, annotation=annotation)
    except _GENE_LIST_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(listed.as_json()))
        return
    # Attribution to stderr, ids to stdout, for the reason `gene-list` splits them: a bare
    # id list is what a shell pipeline wants, and a reader needs to know whose judgement
    # it is. The counts join the attribution because what the census holds and this
    # annotation does not is the one thing a plain id list cannot show.
    gene_ids = listed.gene_ids
    typer.echo(f"TF genes for {listed.assembly} / {listed.annotation} ({listed.species})", err=True)
    typer.echo(f"  {listed.provenance.attribution()}", err=True)
    typer.echo(
        f"  {len(listed.genes)} genes, {len(gene_ids)} gene ids, "
        f"{len(listed.unresolved)} stems this annotation carries no gene for",
        err=True,
    )
    for gene_id in gene_ids:
        typer.echo(gene_id)


@app.command("tf-cofactor-list")
def tf_cofactor_list(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'mm39'."),
    annotation: str | None = typer.Option(None, "--annotation", help=_ANNOTATION_HELP),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the gene ids a publisher lists as transcription cofactors, one per line.

    `genome tf-gene-list` for the other half of the machinery, and the same shape: a
    cofactor — a chromatin remodeller, a histone-modifying enzyme, a Mediator subunit —
    recognises no sequence of its own and so has no motif to scan for, but which genes are
    cofactors is published, and this is that list met with one annotation. Nothing here
    decides what a cofactor is: membership and classification both travel with the
    publisher, and who to cite is printed beside the answer.

    The species comes from the assembly's own metadata row and is never passed in, so
    asking for mouse cofactors while holding a worm assembly is not expressible. A table
    is keyed by gene id stems — gene ids with the version suffix dropped — and a
    registered annotation is not, so every stem is resolved into the ids that annotation
    actually spells and the output joins to a counts matrix with nothing left to
    normalise. A stem naming two genes prints both rather than one of them.

    Only the ids go to stdout, so the output pipes: the heading, the publishers'
    attribution and the counts go to stderr, the last of them saying how many stems this
    annotation carries no gene for. `--json` carries the whole record — every gene with
    the publisher that listed it and that publisher's own classification, one provenance
    entry per publisher to cite, and those unresolved stems.

    A worm assembly is answered here and refused by `genome tf-gene-list`: a publisher
    assessed worm cofactors and none has released a worm TF census. That is what the
    publishers have done rather than a defect here.

    Exits with code 1 when the annotation is not registered here, when no cofactor table
    ships for the assembly's species, and when nothing says what species the assembly is —
    three different facts, each with its own message, and none of them an empty list of
    genes.
    """
    try:
        listed = _tf_cofactor_list(assembly, annotation=annotation)
    except _GENE_LIST_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(listed.as_json()))
        return
    # Attribution to stderr, ids to stdout, for the reason `tf-gene-list` splits them: a
    # bare id list is what a shell pipeline wants, and a reader needs to know whose list it
    # is. The counts join the attribution because what the table holds and this annotation
    # does not is the one thing a plain id list cannot show.
    gene_ids = listed.gene_ids
    typer.echo(
        f"TF cofactors for {listed.assembly} / {listed.annotation} ({listed.species})", err=True
    )
    typer.echo(f"  {listed.provenance.attribution()}", err=True)
    typer.echo(
        f"  {len(listed.cofactors)} cofactors, {len(gene_ids)} gene ids, "
        f"{len(listed.unresolved)} stems this annotation carries no gene for",
        err=True,
    )
    for gene_id in gene_ids:
        typer.echo(gene_id)


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


# --- xref commands -----------------------------------------------------------


@app.command("xref")
def xref(
    species: str = typer.Argument(
        ...,
        help="Species an xref set exists for, e.g. 'Homo sapiens' — the slug "
        "'homo_sapiens' names the same one. A species none exists for names the ones "
        "that do rather than answering nothing.",
    ),
    ids: list[str] = typer.Argument(
        ...,
        help="The identifiers to convert. Each comes back in the order you passed it, "
        "with its version suffix and its namespace's CURIE prefix accepted either way.",
    ),
    to_stems: str | None = typer.Option(
        None,
        "--to-stems",
        metavar="NAMESPACE",
        help=f"Read the ids as this namespace and answer in gene id stems: "
        f"{', '.join(_NAMESPACES)}, whichever of them this set carries. Exactly one of "
        f"this and --from-stems is named.",
    ),
    from_stems: str | None = typer.Option(
        None,
        "--from-stems",
        metavar="NAMESPACE",
        help="Read the ids as gene id stems and answer in this namespace. A versioned "
        "gene id is accepted and reduced to its stem, so an annotation's own ids go "
        "straight in.",
    ),
    source: str | None = typer.Option(
        None,
        "--source",
        help="Answer from this xref source rather than the species' default one. Which "
        "publisher answers is a scientific choice and not a detail: NCBI and Ensembl "
        "agree on 57.5% of human gene-level (GeneID, ENSG) pairs.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Convert identifiers to and from gene id stems against one published xref set.

    The way a column of Entrez GeneIDs from a GEO series, UniProt accessions from a
    mass-spec run or HGNC ids from a curated resource reaches this package's answers,
    without writing Python and without the hand-built join everyone in the lab writes
    slightly differently. No assembly is named and no genome is opened: an identifier is a
    name and not a place.

    **The direction is named, never inferred.** `--to-stems NAMESPACE` reads the ids as
    that namespace and answers in gene id stems; `--from-stems NAMESPACE` reads them as
    stems and answers in that namespace. A string does not say which system it belongs to,
    so `HGNC:11998` asked the wrong way answers *nothing found* rather than quietly turning
    around. There is no third direction: Entrez to HGNC is two calls and the join is yours,
    which keeps the hop visible in your pipeline rather than invisible in ours.

    **The pairs go to stdout, tab-separated, so the output pipes** — `cut -f2` is the
    answer, `cut -f1` says what asked for it — and the heading, the publisher's URL and the
    counts go to stderr. An id naming two genes prints two rows rather than whichever came
    first, and **an id that resolved to nothing gets a row too, with an empty second
    column**: what your list holds and this release does not is the one thing a hand-rolled
    join drops silently. `--json` carries the same answer, keyed by what was asked about,
    with those ids under `unresolved`.

    Omitting `--source` answers from the species' default xref source, so everyone in the
    lab reaches for the same one without discussing it. It is a default and not a
    recommendation: naming a source is how the scientific choice gets made deliberately,
    and every answer names the source and the release that produced it either way.

    Naming a species prepares its set, which the first time is a download. **The lab's CPU
    cluster compute nodes have no internet**, so a set must be constructed once from a login
    node — by running this there, or from Python — before a job that needs it is submitted;
    after that it is read from the **Data dir** and shared by every project on the machine.

    Exits with code 2 when no direction is named or both are, and with code 1 when no set
    exists for the species — the message names the ones that do — when the source is not one
    this package prepares, when the set is not here and cannot be fetched, when the
    namespace is not one the set carries, and when a directory holds a set left unfinished.
    """
    # Which way the hop goes and which namespace it is in arrive on one flag, so neither can
    # be given without the other. Naming both, or neither, is the one thing left to check.
    if to_stems is not None and from_stems is None:
        to_hub, namespace = True, to_stems
    elif from_stems is not None and to_stems is None:
        to_hub, namespace = False, from_stems
    else:
        typer.echo(_XREF_DIRECTION_HELP, err=True)
        raise typer.Exit(code=2)

    try:
        xrefs = _XrefSet(species, source, progressbar=not json)
        answer = xrefs.to_stems(ids, namespace) if to_hub else xrefs.from_stems(ids, namespace)
    except _XREF_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(answer.as_json()))
        return
    _report_xref(answer, source_url=xrefs.source_url, to_hub=to_hub)


def _report_xref(
    answer: _ResolvedStems | _ResolvedXrefIds, *, source_url: str, to_hub: bool
) -> None:
    """Print one hop as the pairs it found, with what answered and what it missed beside them.

    The pairs to stdout and everything else to stderr, for the reason `gene-list` splits
    them: two tab-separated columns are what a shell pipeline wants, and which publisher
    asserted them is what a reader wants. **Every id asked about gets at least one row**,
    the ones that resolved to nothing getting one with an empty second column, so nothing
    leaves this command shorter than it arrived.

    Which noun each column holds is the surface's to say and not the answer's: the answer
    carries the namespace and the direction is what the caller named, so it is passed in
    rather than guessed back out of which type came through.
    """
    asked_as = f"{answer.namespace} ids" if to_hub else "gene id stems"
    answered_as = "gene id stems" if to_hub else f"{answer.namespace} ids"
    rows = [(asked, found) for asked, values in answer.resolved.items() for found in values]
    typer.echo(
        f"{asked_as} -> {answered_as} for {answer.species} ({answer.source} {answer.release})",
        err=True,
    )
    typer.echo(f"  source  {source_url}", err=True)
    typer.echo(
        f"  {len(answer.resolved)} resolved, {len(rows)} {answered_as}, "
        f"{len(answer.unresolved)} this release names none for",
        err=True,
    )
    for asked, found in rows:
        typer.echo(f"{asked}\t{found}")
    for asked in answer.unresolved:
        typer.echo(f"{asked}\t")


# --- motif commands ----------------------------------------------------------


@app.command("motif-scan")
def motif_scan(
    fasta: str = typer.Argument(
        ...,
        help="FASTA to scan, plain or .gz. A record is named by its header up to the first "
        "whitespace — what the aligners write into an alignment made from the same file.",
    ),
    output: str = typer.Argument(
        ...,
        help="Where to write the hits, as Parquet. Read it back with "
        "genome.tf.motif.read_hits, which restores the dtypes and the provenance both.",
    ),
    release: str = typer.Option(
        _DEFAULT_RELEASE,
        "--release",
        help=f"JASPAR release to scan with: {', '.join(_JASPAR_RELEASES)}. Recorded on the "
        f"hits, so a table opened months later still says what made it.",
    ),
    tax_group: str = typer.Option(
        _DEFAULT_TAX_GROUP,
        "--tax-group",
        help=f"JASPAR taxonomic group: {', '.join(_JASPAR_TAX_GROUPS)}. It chooses which file "
        f"is fetched rather than filtering one afterwards, so a worm scan never pays for a "
        f"thousand plant matrices.",
    ),
    threshold: float = typer.Option(
        _DEFAULT_THRESHOLD,
        "--threshold",
        help="The per-position p-value each motif's cutoff is converted from — one number "
        "meaning the same stringency for a short matrix and a long one. A motif that cannot "
        "reach it is left out and named among the skipped, never called at something looser.",
    ),
    background: _BackgroundMode = typer.Option(
        "auto",
        "--background",
        help="The base composition scores are taken against, and the parameter that decides "
        "the answer most: auto derives it from the input above 10 000 unambiguous bases and "
        "stays uniform below that, uniform pins it, derive derives whatever the input holds. "
        "Four frequencies of your own are a Python call rather than a flag.",
    ),
    workers: int | None = typer.Option(
        None,
        "--workers",
        help="How many processes to shard the scan across. Every core the allocation granted "
        "by default — the Slurm allocation first, then this process's CPU affinity, then the "
        "machine — where the library defaults to one. More than one produces the identical "
        "table.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Scan a FASTA with a JASPAR release and write the hits to Parquet.

    The batch case, and the one motif operation that belongs in a shell script and a
    scheduler job: a FASTA in, a Parquet file out, a summary of the run on standard
    output. Listing, plotting and comparing motifs are notebook work and get no command.

    **The hits go to the named file and the summary to standard output**, so `--json` is
    never corrupted by table data. The summary says which release was scanned with, how
    many motifs it scanned and which it left out, the background actually used, how many
    sequences were read, and how many hits were written where — the same facts the table
    itself carries, so a pipeline consuming the summary and a reader opening the file
    months later agree about what happened.

    **It defaults to every core the allocation granted**, where the library defaults to
    one: a console script is a proper entry point, so the process-pool hazard that
    justifies the serial default does not apply here. `--workers 1` scans serially and
    answers with the identical table.

    Naming a release prepares it, which the first time is a download. **The lab's CPU
    cluster compute nodes have no internet**, so a release must be constructed once from a
    login node — by running this there, or from Python — before a job that needs it is
    submitted; after that it is read from the **Data dir** and shared by every project on
    the machine.

    Exits with code 1 when the FASTA is not there or is not FASTA, when the release or tax
    group is not one this package prepares, when the threshold is not a p-value in `(0, 1)`,
    when the worker count is below 1, and when the release is not cached here and cannot be
    fetched — which is what a compute node with no internet looks like.
    """
    # Built before the try and read nowhere yet: `read_fasta` is a generator, so the file
    # is opened when the scan draws its first record and not here.
    records = _CountedRecords(_read_fasta(fasta))
    try:
        motifs = _JasparDatabase(release, tax_group, progressbar=not json)
        # Resolved once, here, rather than left to the scan: the count the summary reports
        # has to be the count the scan used, and this is also what turns `--workers 0` into
        # a refusal before a base is read.
        used = _resolve_workers(workers)
        written = _scan_stream(
            motifs,
            records,
            threshold=threshold,
            background=background,
            output=output,
            workers=used,
        )
    except _MOTIF_SCAN_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    summary = _scan_summary(written, sequences=records.count, workers=used)
    if json:
        typer.echo(_json.dumps(summary))
        return
    _report_scan(summary)


class _CountedRecords:
    """A FASTA's records, counting them as the scan draws them.

    How many sequences were scanned is a property of the stream the scan consumed, and
    nothing else knows it: the **Hit table** names only the sequences that had a hit, and
    asking the file afterwards would mean reading a genome FASTA a second time to answer
    one line of a summary.
    """

    def __init__(self, records: _Iterator[tuple[str, str]]) -> None:
        self._records = records
        self.count = 0

    def __iter__(self) -> _Iterator[tuple[str, str]]:
        """Yield every record, counting it on the way past."""
        for record in self._records:
            self.count += 1
            yield record


def _scan_summary(written: _Path, *, sequences: int, workers: int) -> dict[str, _Any]:
    """Describe a finished scan, reading its provenance and its row count off the file.

    **Never by reading the hits back.** A genome-scale scan is exactly when that is fatal —
    550 million rows is what hg38 against a full vertebrate release comes to — and both
    facts are in the Parquet footer, which is the same size for an empty scan and that one.
    """
    scanned = _provenance_of(written)
    return {
        "release": scanned["release"],
        "tax_group": scanned["tax_group"],
        "motifs_scanned": len(scanned["motifs_scanned"]),
        # The ids and not a count, because *which* motifs were left out is what explains an
        # absent factor; how many were scanned is all a log needs of the other side.
        "motifs_skipped": list(scanned["motifs_skipped"]),
        "background": list(scanned["background"]),
        "threshold": scanned["threshold"],
        "sequences_scanned": sequences,
        "hits_written": _hit_count(written),
        "workers": workers,
        "output": str(written),
    }


def _report_scan(summary: dict[str, _Any]) -> None:
    """Print a finished scan the way a person reads it, off the same summary `--json` emits."""
    typer.echo(
        f"scanned {summary['sequences_scanned']} sequences with {summary['motifs_scanned']} "
        f"motifs from JASPAR {summary['release']} {summary['tax_group']}"
    )
    typer.echo(f"  background  {', '.join(str(value) for value in summary['background'])}")
    typer.echo(f"  threshold   {summary['threshold']}")
    typer.echo(f"  skipped     {_skipped_summary(summary['motifs_skipped'])}")
    typer.echo(f"  workers     {summary['workers']}")
    typer.echo(f"  hits        {summary['hits_written']} -> {summary['output']}")


def _skipped_summary(skipped: list[str]) -> str:
    """Return which motifs were left out of a scan, `none` being an answer and not a silence."""
    if not skipped:
        return f"none — every motif here is at least {_MIN_MOTIF_LENGTH} positions long"
    return (
        f"{len(skipped)} under {_MIN_MOTIF_LENGTH} positions, so not scanned: {', '.join(skipped)}"
    )
