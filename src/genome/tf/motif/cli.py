"""The ``genome motif`` sub-app — a FASTA in, a **Hit table** out, a summary on stdout.

A thin Typer wrapper over :mod:`genome.tf.motif`: it translates arguments, drives the
scan and chooses an output format. It ships from this package so that what the summary
says and what the **Hit table**'s provenance records change in one place. A motif belongs
to no **Assembly**, so this sub-app hangs off the root app rather than under ``tf``.

Examples
--------
>>> from genome.tf.motif.cli import app
>>> [command.name for command in app.registered_commands]
['scan']
"""

from __future__ import annotations

import json as _json
from collections.abc import Iterator as _Iterator
from pathlib import Path as _Path
from typing import Any as _Any

import typer

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
from genome.workers import resolve_workers as _resolve_workers

#: What a failed motif scan raises. The same three, and each already says what to do: a
#: release this package does not prepare, a threshold that is not a p-value and a file that
#: is not FASTA are ``ValueError``s; a FASTA that is not there is an ``OSError``; and a
#: release that cannot be fetched, a directory holding one left unfinished and a process
#: pool that died under the scan are ``RuntimeError``s.
#: Its own name rather than the assembly sub-app's list reused, because a motif belongs to
#: no assembly and the two lists are alike by coincidence rather than by construction.
_MOTIF_SCAN_ERRORS = (ValueError, OSError, RuntimeError)

app = typer.Typer(
    help="Scan a FASTA with a JASPAR release and write the hits to Parquet.",
    no_args_is_help=True,
)


@app.command("scan")
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
    when the worker count is below 1, when the release is not prepared here and cannot be
    fetched — which is what a compute node with no internet looks like — and when a directory
    holds a release left unfinished.
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
