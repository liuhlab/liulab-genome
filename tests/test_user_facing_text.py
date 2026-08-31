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
- the docstring of a function Typer prints as ``--help`` — a registered command's, or a
  group's own callback — which `mkdocs-typer2` also renders onto the site's CLI page. A
  command's docstring **is** runtime output, so "docstrings are exempt" does not hold at
  the CLI boundary.

Everything else keeps its numbers, deliberately: comments, and the docstrings of ordinary
public objects. Both are agent-facing. Whether `mkdocstrings` should publish record
numbers onto the reference page at all is a separate, still-open question, and a guard
that flagged ordinary docstrings would answer it here by force.

The invariant is asserted, never the list. Ten messages and one command docstring cited
records when this was written; an eleventh added tomorrow fails here, which is the whole
point — a test naming the ten would go green the moment somebody wrote the eleventh.

**What the walk cannot see**, so that a green suite is not read as more than it is: it
reads string literals in the package's own ``.py`` sources. A message assembled at runtime
from a shipped table, or interpolated from a value read off disk, is outside it.
"""

from __future__ import annotations

import ast
import re
from collections.abc import Callable, Iterator

import typer

from genome.cli import app

from ._sources import PACKAGE, parse, sources

#: What a reader outside this repository cannot resolve.
RECORD_NUMBER = re.compile(r"ADR-\d{4}")


def _docstring_node_ids(tree: ast.Module) -> set[int]:
    """Return the ``id()`` of every string node holding a module, class or function docstring.

    Identity, not text: two docstrings with the same content are two nodes, and a message
    that happens to repeat a docstring word for word is not one of them.
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


def runtime_record_citation_lines(tree: ast.Module) -> list[int]:
    """Return the line of every string literal in ``tree`` that cites a record and is not a docstring.

    A comment is not a string literal and never reaches this: `ast` drops comments, which
    is exactly the discrimination wanted rather than an accident worth working around.
    """
    docstrings = _docstring_node_ids(tree)
    return [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, str)
        and id(node) not in docstrings
        and RECORD_NUMBER.search(node.value)
    ]


def _help_texts(
    group: typer.Typer, prefix: str = "genome"
) -> Iterator[tuple[str, Callable[..., object]]]:
    """Yield ``("genome xref symbols", callback)`` for every docstring Typer prints under ``group``.

    Three kinds, because Typer prints three. A **command**'s docstring is its own
    ``--help``. A **group callback**'s docstring is the group's ``--help`` — the page a
    reader lands on before any command — and there is none today, which is exactly why it
    is walked: the eleventh site being prevented is the one nobody thought to look at.
    The sub-apps are recursed rather than listed, so a topic mounted tomorrow is covered
    without this being edited.
    """
    if group.registered_callback is not None and group.registered_callback.callback is not None:
        yield prefix, group.registered_callback.callback
    for command in group.registered_commands:
        if command.callback is None:
            continue
        name = command.name or command.callback.__name__.replace("_", "-")
        yield f"{prefix} {name}", command.callback
    for sub in group.registered_groups:
        name = f"{prefix} {sub.name or ''}".strip()
        if sub.callback is not None:
            yield name, sub.callback
        if sub.typer_instance is not None:
            yield from _help_texts(sub.typer_instance, name)


def test_no_message_this_package_can_print_cites_a_record_number() -> None:
    """A string literal that is not a docstring is text a user can be handed."""
    offenders = [
        f"{path.relative_to(PACKAGE.parent)}:{line}"
        for path in sources()
        for line in runtime_record_citation_lines(parse(path))
    ]
    assert not offenders, (
        "these strings can be printed to a user and cite a record they cannot look up; "
        f"name the constraint in words instead: {', '.join(offenders)}"
    )


def test_no_cli_help_text_cites_a_record_number() -> None:
    """A command's or a group's docstring is its ``--help`` and its page on the site."""
    offenders = [
        name for name, callback in _help_texts(app) if RECORD_NUMBER.search(callback.__doc__ or "")
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
    assert runtime_record_citation_lines(ast.parse(source)) == [7]
