"""The ``genome tf`` sub-app — which genes a publisher judges factors, and which cofactors.

A thin Typer wrapper over :mod:`genome.tf.gene` and :mod:`genome.tf.cofactor`: it
translates arguments, dispatches to the by-assembly-name module functions and chooses an
output format. It ships from this package so that what ``gene-list`` prints and what a
**TF gene table** holds change in one place, and no **Genome** is constructed to answer
either command.

The two option-help strings and the error list come from :mod:`genome.annotation.cli`
rather than being spelled again: both commands ask one registered annotation a question a
**Shipped table** answers, which is the same shape the gene-category pair has.

Examples
--------
>>> from genome.tf.cli import app
>>> [command.name for command in app.registered_commands]
['gene-list', 'cofactor-list']
"""

from __future__ import annotations

import json as _json

import typer

from genome.annotation.cli import _ANNOTATION_HELP, _GENE_LIST_ERRORS
from genome.tf.cofactor import tf_cofactor_list as _tf_cofactor_list
from genome.tf.gene import tf_gene_list as _tf_gene_list

app = typer.Typer(
    help="Read a species' transcription factors and cofactors in one annotation's gene ids.",
    no_args_is_help=True,
)


@app.command("gene-list")
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
    # Attribution to stderr, ids to stdout, for the reason `annotation gene-list` splits
    # them: a bare id list is what a shell pipeline wants, and a reader needs to know whose
    # judgement it is. The counts join the attribution because what the census holds and
    # this annotation does not is the one thing a plain id list cannot show.
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


@app.command("cofactor-list")
def tf_cofactor_list(
    assembly: str = typer.Argument(..., help="Assembly name, e.g. 'mm39'."),
    annotation: str | None = typer.Option(None, "--annotation", help=_ANNOTATION_HELP),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the gene ids a publisher lists as transcription cofactors, one per line.

    `genome tf gene-list` for the other half of the machinery, and the same shape: a
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

    A worm assembly is answered here and refused by `genome tf gene-list`: a publisher
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
    # Attribution to stderr, ids to stdout, for the reason `tf gene-list` splits them: a
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
