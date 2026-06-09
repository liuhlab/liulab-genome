"""Command-line interface (thin wrapper over the API).

Real commands land in Phase 3. This stub keeps the ``genome`` entry point importable.
"""

import typer

app = typer.Typer(help="Tools for handling genomic files.", no_args_is_help=True)


@app.command()
def version() -> None:
    """Print the installed package version."""
    from genome import __version__

    typer.echo(__version__)
