"""The ``genome homology`` sub-app — which genes of another species a gene is homologous to.

A thin Typer wrapper over :mod:`genome.homology`: it translates arguments, makes one call
on a :class:`~genome.homology.compara.HomologySet` and renders. It ships from this package
so that what a link prints and what a **Homology link** holds change in one place. No
**Assembly** is named and no **Genome** is opened — a set is anchored to a species pair and
a **Release**, not to a build.

Examples
--------
>>> from genome.homology.cli import app
>>> [command.name for command in app.registered_commands]
['links']
"""

from __future__ import annotations

import json as _json

import typer

from genome.homology import DEFAULT_RELEASE as _HOMOLOGY_RELEASE
from genome.homology import NULL_CELL as _NULL_CELL
from genome.homology import HomologyMetadata as _HomologyMetadata
from genome.homology import HomologySet as _HomologySet
from genome.homology.compara import HomologyAnswer as _HomologyAnswer

#: What a failed homology lookup raises, and every one of them already names its next
#: action: a species and a species pair nothing is pinned for are ``LookupError``s; a release
#: that is not pinned, both species being the same one, a fetched file that is not Compara's
#: and a stored slice that changed after it was prepared are ``ValueError``s; the partition
#: guard, an unfinished directory and a set that could not be fetched are ``RuntimeError``s;
#: and a file that went away under the read is an ``OSError``. Its own name rather than the
#: xref sub-app's list reused, because a **Homology set** and an **Xref set** are two
#: contexts and the two lists are alike by coincidence rather than by construction.
_HOMOLOGY_ERRORS = (ValueError, OSError, RuntimeError, LookupError)

#: The columns a **Homology link** is printed as, which are the keys
#: :meth:`~genome.homology.compara.HomologyLink.as_json` writes and in its order — so the text
#: rendering and ``--json`` cannot drift apart, and every cell printed is a value the API
#: put in the answer rather than one assembled here. ``homology_type`` and ``is_ortholog``
#: are what **mark** a **Paralogy link**: a duplication label prints where a speciation one
#: would, which is what keeps *not an ortholog* distinguishable from *absent*.
_HOMOLOG_COLUMNS: tuple[str, ...] = (
    "gene_id_stem",
    "homolog_gene_id_stem",
    "homology_type",
    "is_ortholog",
    "is_high_confidence",
    "goc_score",
    "wga_coverage",
)

app = typer.Typer(
    help="Read Ensembl Compara's links between the genes of two species.",
    no_args_is_help=True,
)


@app.command("links")
def homologs(
    species: str = typer.Argument(
        ...,
        help="Species the gene id stems belong to, e.g. 'Homo sapiens' — the slug "
        "'homo_sapiens' names the same one. A species no set is pinned for names the ones "
        "that are rather than answering nothing.",
    ),
    other_species: str = typer.Argument(
        ...,
        help="Species the homologous genes come back in. Any pairing among human, mouse "
        "and worm; the same species twice is refused, since a gene's paralogs within one "
        "species is a different question this does not answer.",
    ),
    stems: list[str] = typer.Argument(
        ...,
        help="The gene id stems to ask about, answered in the order you passed them. "
        "Compara writes its gene ids bare, so a versioned id is refused by name rather "
        "than answered emptily.",
    ),
    release: str = typer.Option(
        _HOMOLOGY_RELEASE,
        "--release",
        help="Ensembl Compara release to answer from. Recorded on the answer, so a result "
        "is reproducible a year later.",
    ),
    paralogs: bool = typer.Option(
        False,
        "--paralogs",
        help="Return every link the publisher wrote for these genes rather than only the "
        "ones its own label calls a speciation event. A paralogy link is marked by that "
        "label in the `homology_type` column, never excluded.",
    ),
    json: bool = typer.Option(False, "--json", help="Emit JSON instead of plain text."),
) -> None:
    """Print the genes of another species a gene id stem's gene is homologous to.

    The way a hit carries across species without leaving the package and without the
    Ensembl BioMart web API, whose intermittent failures make a pipeline built on it fail
    irreproducibly. Everything here is a bulk file fetched once and read locally. No
    assembly is named and no genome is opened: a **Homology set** is anchored to a species
    pair and a release, not to a build.

    **Every cell is Ensembl Compara's.** The `homology_type` — `ortholog_one2one`,
    `ortholog_one2many`, `ortholog_many2many` — is the publisher's own tree-derived label
    printed verbatim, and it is never recomputed from what came back: an answer can show
    one partner and still read `ortholog_one2many`, which is the point of carrying the
    label rather than counting rows. The high-confidence flag and both quality scores come
    through the same way, and this package publishes no score, ranking or "best ortholog"
    of its own.

    **The links go to stdout, tab-separated, so the output pipes** — `cut -f2` is the
    answer, `cut -f1` says what asked for it — and the heading, the attribution, the counts
    and the two qualifications below go to stderr. A gene with three homologs prints three
    rows rather than whichever came first, and **a stem this release names no homolog for
    gets a row too, with every other column empty**: what your list holds and this release
    does not is visible rather than dropped. An empty cell there is not `NULL`, which is
    the publisher's own word for a cell it recorded nothing in on a link that does exist.

    **Two qualifications ride on every answer.** The **Dropped partner**s — the homologous
    genes a filter removed — are counted and named, so a link that merely *looks*
    one-to-one in your view stays distinguishable from one the publisher called one-to-one.
    And whichever quality columns this set holds no value in *anywhere* are named up front:
    Compara records neither `goc_score` nor `wga_coverage` on any link of either worm
    pairing, so a filter written against one would empty itself in silence.

    Orthologs are the answer by default and `--paralogs` returns every link the publisher
    wrote; a **Paralogy link** is marked by its own `homology_type` rather than excluded,
    so *not an ortholog* stays distinguishable from *absent*. Release 116 publishes no
    cross-species paralogy for these three species, so on it the flag changes nothing.

    Naming a pair prepares its set, which the first time is a download. **The lab's CPU
    cluster compute nodes have no internet**, so a set must be constructed once from a
    login node — by running this there, or from Python — before a job that needs it is
    submitted; after that it is read from the **Data dir** and shared by every project on
    the machine.

    Exits with code 1 when no set is pinned for the species — the message names the ones
    that are — when the release is not pinned, when both species are the same one, when a
    stem carries a version, when the set is not here and cannot be fetched, when a
    directory holds a set left unfinished, and when the file that was recorded as holding
    this pair holds none of its rows, which means Compara re-partitioned and the message
    names the other file.
    """
    try:
        homology = _HomologySet(species, other_species, release, progressbar=not json)
        answer = homology.homologs(stems, paralogs=paralogs)
    except _HOMOLOGY_ERRORS as err:
        typer.echo(f"error: {err}", err=True)
        raise typer.Exit(code=1) from err

    if json:
        typer.echo(_json.dumps(answer.as_json()))
        return
    _report_homologs(answer, provenance=homology.provenance, paralogs=paralogs)


def _report_homologs(
    answer: _HomologyAnswer, *, provenance: _HomologyMetadata, paralogs: bool
) -> None:
    """Print one pair's links, with what asserted them and what qualifies them beside.

    The links to stdout and everything else to stderr, for the reason `xref ids` splits
    them: tab-separated columns are what a shell pipeline wants, and which publisher
    asserted them is what a reader wants. **Every stem asked about gets at least one row**,
    the ones naming no homolog getting one with every other column empty, so nothing leaves
    this command shorter than it arrived.

    Which question was asked is the surface's to say and not the answer's — the answer
    holds what came back, not which switch the caller set — so it is passed in rather than
    guessed back out of the labels that survived.
    """
    kind = "homologs, paralogy included" if paralogs else "orthologs"
    typer.echo(
        f"{answer.species} -> {answer.other_species} {kind} "
        f"({provenance.publisher} {answer.release})",
        err=True,
    )
    typer.echo(f"  source   {provenance.attribution()}", err=True)
    typer.echo(f"  columns  {', '.join(_HOMOLOG_COLUMNS)}", err=True)
    typer.echo(
        f"  {len(answer.resolved)} resolved, {len(answer.links)} links, "
        f"{len(answer.unresolved)} this release names no homolog for, "
        f"{_dropped_summary(answer.dropped_partners)}",
        err=True,
    )
    typer.echo(f"  quality  {_quality_summary(answer.null_quality_scores)}", err=True)
    for link in answer.links:
        written = link.as_json()
        typer.echo("\t".join(_cell(written[column]) for column in _HOMOLOG_COLUMNS))
    for stem in answer.unresolved:
        typer.echo(stem + "\t" * (len(_HOMOLOG_COLUMNS) - 1))


def _cell(value: object) -> str:
    """Render one cell of a **Homology link**, in the publisher's word for *nothing here*."""
    return _NULL_CELL if value is None else str(value)


def _dropped_summary(dropped: tuple[str, ...]) -> str:
    """Say which partners a filter removed — counted *and* named, `0` being an answer."""
    if not dropped:
        return "0 dropped partners"
    return f"{len(dropped)} dropped partners: {', '.join(dropped)}"


def _quality_summary(null_columns: tuple[str, ...]) -> str:
    """Say which quality columns this set holds nothing in, before a filter empties itself."""
    if not null_columns:
        return "the publisher scored this pair; both quality columns carry values"
    return (
        f"{' and '.join(null_columns)} null on every link of this set, so a filter on "
        f"{'either' if len(null_columns) > 1 else 'it'} empties rather than narrowing"
    )
