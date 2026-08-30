"""What every CLI test module needs: one runner, both streams of a result, and help text.

The commands split across `tests/assembly/`, `tests/annotation/`, `tests/tf/`,
`tests/xref/` and `tests/homology/` the way the source does, and these three are the only
things all of those modules share. Fixtures stay with the tests that need them — a fixture
handing a class setup it never had would be a silent change of what is covered.
"""

from __future__ import annotations

import re

from typer.testing import CliRunner

from genome.cli import app

runner = CliRunner()


def output(result: object) -> str:
    """Return a result's stdout and stderr together, wherever the runner put them."""
    return (getattr(result, "stdout", "") or "") + (getattr(result, "stderr", "") or "")


#: Every ANSI escape sequence rich writes — the colour, the bolding and the resets it
#: emits at each line boundary. Stripped before a help string is asserted against, because
#: whether the runner is drawing colour is a property of the terminal that happens to be
#: attached and never of what the help says: CI colours its output and a local run may not,
#: which is a test that passes on one machine and fails on the other.
_ANSI = re.compile(r"\x1b\[[0-9;]*m")


def help_text(*command: str) -> str:
    """One command's ``--help``, with the drawing taken out of the words.

    Rich draws the options in a bordered table, so a sentence too long for one line arrives
    with a ``│`` and a newline through the middle of it, and a colour reset wherever it
    broke. What a test asserts is what the help *says*, so the escapes and the rules go and
    the whitespace collapses.
    """
    rendered = output(runner.invoke(app, [*command, "--help"]))
    return " ".join(_ANSI.sub("", rendered).replace("│", " ").split())
