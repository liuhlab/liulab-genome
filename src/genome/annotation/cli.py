"""The ``genome annotation`` sub-app — register a GTF, list what is here, read a category.

A thin Typer wrapper over :mod:`genome.annotation`: it translates arguments, dispatches to
the by-assembly-name module functions and chooses an output format. It ships from this
package so that what ``gene-list`` prints and what a **Gene list** holds change in one
place, and no **Genome** is constructed to answer any of it.

Examples
--------
>>> from genome.annotation.cli import app
>>> [command.name for command in app.registered_commands]
['register', 'register-gtf', 'list', 'gene-list', 'gene-categories']
"""

from __future__ import annotations

import json as _json

import typer

from genome.annotation import GeneList as _GeneList
from genome.annotation import GeneListSource as _GeneListSource
from genome.annotation import RegisteredAnnotation as _RegisteredAnnotation
from genome.annotation import annotation_status as _annotation_status
from genome.annotation import gene_list as _gene_list
from genome.annotation import gene_lists as _gene_lists
from genome.annotation import register_annotation as _register_annotation
from genome.annotation import register_gtf as _register_gtf

#: What a failed annotation registration raises, and each of them is already actionable —
#: a checksum mismatch, a chromosome name the assembly does not carry, a registration that
#: cannot be trusted — so the command prints the message and exits non-zero rather than
#: adding to it. ``RegistrationError`` is a ``RuntimeError``, ``ChecksumMismatchError`` a
#: ``ValueError``, a failed download an ``OSError``. Its own name rather than the assembly
#: sub-app's list reused, because an annotation is its own context and the two are alike by
#: coincidence rather than by construction.
_ANNOTATION_ERRORS = (ValueError, OSError, RuntimeError)

#: What a failed gene-list command raises, on top of the above — the gene-category pair
#: here and both TF lists in :mod:`genome.tf.cli`, which imports this rather than spelling
#: it again. An annotation that ships no curated gene list, one whose list does not declare
#: the category asked for, a species no census or cofactor table ships for and an assembly
#: nothing names a species for are lookups that found nothing rather than bad values — so
#: they are ``LookupError``s, which :data:`_ANNOTATION_ERRORS` does not cover, and so is the
#: ``KeyError`` an unregistered annotation has always raised. One list for all four commands
#: because they are alike by construction: each asks one registered annotation a question a
#: shipped file answers, and each absence on the way is a lookup.
_GENE_LIST_ERRORS = (*_ANNOTATION_ERRORS, LookupError)

#: Help for the ``--annotation`` option every command asking one annotation a question
#: takes — the two below, and the two in :mod:`genome.tf.cli`, which imports this so that a
#: reader of either help text gets the same explanation.
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

app = typer.Typer(
    help="Register an assembly's annotations, list them, and read a gene category.",
    no_args_is_help=True,
)


@app.command("register")
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
    from a path instead, with `genome annotation register-gtf`.

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
    except _ANNOTATION_ERRORS as err:
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

    The escape hatch for a GTF the annotation table does not list — `genome annotation
    register` is the way in for one it does. Nothing is downloaded and no checksum is
    compared against, because an unlisted GTF has none pinned for it: the file you name is
    placed under `<assembly dir>/gtf/<name>/` (a `.gz` is decompressed on the way), its
    chromosome names are checked against the assembly's, the gffutils database is built,
    and the record that says all of it finished is written last. It is addressed by
    `<name>` from then on, exactly as a listed annotation is, and shows up in `genome
    annotation list` as registered but not offered.

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
    except _ANNOTATION_ERRORS as err:
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
@app.command("list")
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
    except _ANNOTATION_ERRORS as err:
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

    `genome annotation gene-categories <assembly>` says which categories may be asked for;
    they are the curated list's to declare and differ between annotations.

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

    What may be asked of `genome annotation gene-list`, since the categories are a property
    of the annotation rather than of this package: human GENCODE splits rRNA pseudogenes
    out and mouse does not, yeast carries a precursor category nothing else has, and a
    bacterium has no mitochondrial category at all.

    A merged annotation shows the per-component split beside each count, so a category one
    component declares and another does not is visible as such rather than as a smaller
    number. `--json` carries every category with its gene ids and its sources — the same
    answer `genome annotation gene-list` gives for one of them, for all of them at once.

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
