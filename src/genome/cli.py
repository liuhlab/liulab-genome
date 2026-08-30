"""Command-line interface — the root Typer app, and the sub-app per topic hung off it.

Every command that names an **Assembly**, an **Annotation**, a TF list, an **Xref set**, a
**Homology set** or a **Motif set** ships from the package that owns it, beside the result
type it renders. What is left here is the three commands belonging to no topic —
``version``, ``revcomp`` and ``doctor`` — the :meth:`typer.Typer.add_typer` calls that
mount the six sub-apps, and the deprecated-alias table below.

Examples
--------
>>> from genome.cli import app
>>> [group.name for group in app.registered_groups]
['assembly', 'annotation', 'tf', 'xref', 'homology', 'motif']
"""

from __future__ import annotations

import json as _json
from collections.abc import Callable as _Callable

import typer

from genome import __version__ as _package_version
from genome.annotation import cli as _annotation_cli
from genome.assembly import cli as _assembly_cli
from genome.external import ToolNotFoundError
from genome.external import doctor as _doctor
from genome.homology import cli as _homology_cli
from genome.seq import DNA
from genome.tf import cli as _tf_cli
from genome.tf.motif import cli as _motif_cli
from genome.xref import cli as _xref_cli

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


# --- the topic sub-apps ------------------------------------------------------
#
# Each is named for the module it ships from, which is why the Orthology context's is
# `homology`; and `motif` hangs off the root rather than off `tf` because a **Motif set**
# belongs to no assembly and is usable with no genome open.
app.add_typer(_assembly_cli.app, name="assembly")
app.add_typer(_annotation_cli.app, name="annotation")
app.add_typer(_tf_cli.app, name="tf")
app.add_typer(_xref_cli.app, name="xref")
app.add_typer(_homology_cli.app, name="homology")
app.add_typer(_motif_cli.app, name="motif")


#: The flat spelling each command answered to before it moved under a sub-app, mapped to
#: the very function object the sub-app registered — Typer's decorator returns the function
#: unchanged, so nothing here is a second implementation and the two spellings cannot
#: drift. One mapping and one loop so the release after this one deletes the aliases as a
#: unit: cut this table, cut the loop below it, and nothing else moves.
_DEPRECATED_ALIASES: dict[str, _Callable[..., None]] = {
    "register": _assembly_cli.register,
    "verify": _assembly_cli.verify,
    "table-row": _assembly_cli.table_row,
    "register-annotation": _annotation_cli.register_annotation,
    "register-gtf": _annotation_cli.register_gtf,
    "annotations": _annotation_cli.list_annotations,
    "gene-list": _annotation_cli.gene_list,
    "gene-categories": _annotation_cli.gene_categories,
    "tf-gene-list": _tf_cli.tf_gene_list,
    "tf-cofactor-list": _tf_cli.tf_cofactor_list,
    # `xref` is missing from this table and is the one old spelling that cannot be in it:
    # it is a sub-app's name as well as a command's, and a root app holds one of each under
    # one name. `genome xref SPECIES ID...` is answered by `_XrefGroup` in
    # `genome.xref.cli` instead, which is deleted when this table is.
    "match-symbols": _xref_cli.match_symbols,
    "homologs": _homology_cli.homologs,
    "motif-scan": _motif_cli.motif_scan,
}

# `hidden` keeps every one of them out of `genome --help`, so the tree a reader is shown is
# the new one alone; `deprecated` is what prints the notice naming the command — on stderr,
# which is what leaves `--json` on stdout parseable for a script that has not moved yet.
for _old_spelling, _command in _DEPRECATED_ALIASES.items():
    app.command(_old_spelling, hidden=True, deprecated=True)(_command)
