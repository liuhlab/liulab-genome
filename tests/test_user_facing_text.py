"""The user-facing text guard: a record number is a citation between agents, never to a user.

`AGENTS.md` sanctions citing a record by number — a record number is permanent where a
path is not, and repo-internal prose should cite decisions that way. The defect this
guards is not the citation, it is the **boundary crossing**: the record tree is kept out
of the built site on purpose, so a number that reaches a reader outside this repository
resolves nowhere they can go — not on the site, not in the installed package, not in the
search index. It reads as a reference and functions as noise.

**So the test is audience, not syntax.** Would a reader outside this repository see this
string? Two places where the answer is yes, and this module holds both:

- a string literal that is not a docstring — an exception message, or a module constant
  whose value is printed;
- the docstring of a function registered as a CLI command, which Typer prints verbatim as
  ``--help`` and `mkdocs-typer2` renders onto the site's CLI page. A command's docstring
  **is** runtime output, so "docstrings are exempt" does not hold at the CLI boundary.

Everything else keeps its numbers, deliberately: comments, and the docstrings of ordinary
public objects. Both are agent-facing. Whether `mkdocstrings` should publish record
numbers onto the reference page at all is a separate, still-open question, and a guard
that flagged ordinary docstrings would answer it here by force.

The invariant is asserted, never the list. Ten messages and one command docstring cited
records when this was written; an eleventh added tomorrow fails here, which is the whole
point — a test naming the ten would go green the moment somebody wrote the eleventh.
"""

from __future__ import annotations

import ast
import inspect
import re
from collections.abc import Callable, Iterator
from pathlib import Path

import typer

import genome
from genome.cli import app

#: What a reader outside this repository cannot resolve.
RECORD_NUMBER = re.compile(r"ADR-\d{4}")

#: The package's sources, found from the installed module rather than from this file's
#: own path, so the guard walks what is imported and not a tree that happens to sit here.
SRC = Path(inspect.getfile(genome)).parent


def _docstring_nodes(tree: ast.Module) -> set[int]:
    """Return the ids of every string node that is a docstring, module, class or function.

    Identity, not text: two docstrings with the same content are two nodes, and a runtime
    string that happens to repeat a docstring's text is not one of them.
    """
    holders = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    ids: set[int] = set()
    for node in ast.walk(tree):
        if not isinstance(node, holders):
            continue
        first = node.body[0] if node.body else None
        if (
            isinstance(first, ast.Expr)
            and isinstance(first.value, ast.Constant)
            and isinstance(first.value.value, str)
        ):
            ids.add(id(first.value))
    return ids


def runtime_record_numbers(source: str) -> list[int]:
    """Return the line of every string literal in ``source`` that cites a record and is not a docstring.

    A comment is not a string literal and never reaches this: `ast` drops comments, which
    is exactly the discrimination wanted rather than an accident worth working around.
    """
    tree = ast.parse(source)
    docstrings = _docstring_nodes(tree)
    return [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in docstrings
        and RECORD_NUMBER.search(node.value)
    ]


def _commands(
    group: typer.Typer, prefix: str = "genome"
) -> Iterator[tuple[str, Callable[..., object]]]:
    """Yield ``("genome xref symbols", callback)`` for every command under ``group``, recursively.

    The sub-apps are walked rather than listed, so a topic mounted tomorrow is covered
    without this being edited.
    """
    for command in group.registered_commands:
        if command.callback is None:
            continue
        name = command.name or command.callback.__name__.replace("_", "-")
        yield f"{prefix} {name}", command.callback
    for sub in group.registered_groups:
        if sub.typer_instance is None:
            continue
        name = sub.name or ""
        yield from _commands(sub.typer_instance, f"{prefix} {name}".strip())


def test_no_message_this_package_can_print_cites_a_record_number() -> None:
    """A string literal that is not a docstring is text a user can be handed."""
    offenders = [
        f"{path.relative_to(SRC.parent)}:{line}"
        for path in sorted(SRC.rglob("*.py"))
        for line in runtime_record_numbers(path.read_text())
    ]
    assert not offenders, (
        "these strings can be printed to a user and cite a record they cannot look up; "
        f"name the constraint in words instead: {', '.join(offenders)}"
    )


def test_no_cli_command_help_cites_a_record_number() -> None:
    """A command's docstring is its ``--help`` and its page on the site."""
    offenders = [
        name for name, callback in _commands(app) if RECORD_NUMBER.search(callback.__doc__ or "")
    ]
    assert not offenders, (
        "`--help` for these commands cites a record the reader cannot look up; name the "
        f"constraint in words instead: {', '.join(offenders)}"
    )


def test_the_guard_reads_the_boundary_and_not_the_pattern() -> None:
    """Comments and ordinary docstrings keep their numbers; only the printable string is caught.

    The guard above is only correct if it discriminates, so the discrimination is tested
    directly rather than inferred from a green suite: a module that cites the same record
    three ways is reported once, on the line that a user could read.
    """
    source = '''\
"""A module docstring citing ADR-0001."""


def helper() -> str:
    """A function docstring citing ADR-0001."""
    # A comment citing ADR-0001.
    return "a message citing ADR-0001"
'''
    assert runtime_record_numbers(source) == [7]
