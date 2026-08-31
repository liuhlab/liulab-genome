"""Record citations: agent-facing in `src/`, and absent from every page the site publishes.

A record number is a citation between agents. It resolves for anyone reading the source,
because the record tree is in the repository beside it, and for nobody reading the built
site, because `mkdocs.yml` keeps that tree out on purpose. `mkdocstrings` renders `src/`
docstrings onto the reference page and so carries citations across that line; the extension
under test removes them as griffe hands each docstring over, leaving the source untouched.

Two guards, at two altitudes. Here, that the citation's *form* is the one the extension can
remove — which is what makes the removal total rather than a best effort. In
`scripts/check_site_has_no_record_numbers.sh`, run by `pixi run docs-build`, that no number
reached the built artifact by any route at all. This one is in the fast lane on purpose: a
docstring that would leak should fail beside the docstring, not only in the docs job.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from mkdocs_record_citations import CITATION, StripRecordCitations

from ._sources import PACKAGE, parse, sources

#: The record tree, which is in the repository and deliberately not in the built site.
RECORDS = Path(__file__).resolve().parents[1] / "docs" / "adr"

#: A record number wherever it appears, in any shape. `CITATION` matches only the shape the
#: extension removes, so the two together are what says *every* number is in that shape.
ANY_NUMBER = re.compile(r"ADR-\d{4}")

#: The house form: a trailing parenthetical, one record or several. A citation woven into a
#: sentence as its subject cannot be deleted without leaving the prose ungrammatical, which
#: is why the form is a rule and not a habit.
TRAILING_PARENTHETICAL = re.compile(r"\(ADR-\d{4}(?:,\s+ADR-\d{4})*\)")


@pytest.fixture
def strip() -> StripRecordCitations:
    """The extension exactly as `mkdocs.yml` configures it, which is with nothing."""
    return StripRecordCitations()


def test_a_citation_is_removed_and_the_sentence_closes_up(strip: StripRecordCitations) -> None:
    """The space in front of the parenthetical goes with it, or a gap is published instead."""
    assert strip.render("nesting is forbidden by the model (ADR-0008).") == (
        "nesting is forbidden by the model."
    )


def test_several_records_in_one_parenthetical_go_together(strip: StripRecordCitations) -> None:
    """Three docstrings cite a pair at once; half a citation is worse than none."""
    assert strip.render("mixing builds is an error (ADR-0003, ADR-0005), never a warning.") == (
        "mixing builds is an error, never a warning."
    )


def test_a_citation_wrapped_onto_its_own_line_is_removed_with_the_line_break(
    strip: StripRecordCitations,
) -> None:
    """Docstrings wrap, so a citation is often the first thing on a line.

    Taking the newline with it rejoins the sentence exactly as the renderer would have,
    rather than leaving a leading space where the parenthetical stood.
    """
    assert strip.render("a free-form local key\n(ADR-0003), not necessarily a UCSC one.") == (
        "a free-form local key, not necessarily a UCSC one."
    )


def test_prose_that_cites_nothing_is_returned_unchanged(strip: StripRecordCitations) -> None:
    """The extension runs over every docstring in the package, so it must be inert."""
    prose = "Return the reverse complement, which is an involution."

    assert strip.render(prose) == prose


def test_the_rewrite_does_not_reach_the_module_a_caller_imports() -> None:
    """The site's copy loses the citation; ``help()`` and ``__doc__`` keep it.

    This is what keeps the two readers separate rather than trading one for the other. A
    docstring is agent-facing, so the number belongs in it; the page is not, so it does
    not belong there. The extension edits the docstring griffe parsed from the file, and
    the live ``__doc__`` is a different string that nothing here touches.
    """
    from genome import Genome

    docstring = Genome.__doc__ or ""
    assert ANY_NUMBER.search(docstring), "picked an example that no longer cites a record"

    StripRecordCitations().render(docstring)

    assert ANY_NUMBER.search(Genome.__doc__ or "")


def test_every_citation_in_the_package_is_one_the_extension_can_remove() -> None:
    """The form invariant, and the reason the strip is total rather than best-effort.

    Seven citations were once the subject of their own sentence — "the guess ADR-0003
    exists to forbid" — and deleting the number from one leaves a hole in the grammar.
    Naming the idea and citing in a trailing parenthetical costs nothing and makes the
    removal a rule the renderer can keep.
    """
    # `src/**/*.md` is out of scope and tracked as #201: the shipped ATTRIBUTION files reach
    # a user with no record tree, but no render step can strip them and the fix is in the
    # prose. `sources()` walks `.py` only, which is what draws that line here.
    off_form: list[str] = []
    for path in sources():
        # Wrapping is the renderer's business, not the citation's: a parenthetical split
        # across two lines is the same citation, so flatten before asking about form.
        flat = re.sub(r"\s*\n\s*(?:#\s*)?", " ", path.read_text(encoding="utf-8"))
        covered = {i for m in TRAILING_PARENTHETICAL.finditer(flat) for i in range(*m.span())}
        off_form += [
            f"{path.relative_to(PACKAGE.parent)}: …{flat[max(0, m.start() - 55) : m.end() + 20]}…"
            for m in ANY_NUMBER.finditer(flat)
            if m.start() not in covered
        ]

    assert not off_form, (
        "cite a record as a trailing parenthetical — '(ADR-0006)' — and name the idea in "
        f"the sentence, or the site cannot drop the number cleanly: {off_form}"
    )


def test_every_record_the_package_cites_exists() -> None:
    """A renumbered record fails here, in the fast lane, and not only in the docs build.

    The site no longer publishes these numbers, so nothing downstream would notice a
    dangling one. That makes this guard more necessary than it was, not less: the reader
    it now serves entirely is the agent, and a number pointing at no record misleads one.
    """
    numbered = {path.name[:4] for path in RECORDS.glob("[0-9][0-9][0-9][0-9]-*.md")}
    assert numbered, f"no records under {RECORDS}"

    dangling = sorted(
        f"{path.relative_to(PACKAGE.parent)}: {number}"
        for path in sources()
        for number in ANY_NUMBER.findall(path.read_text(encoding="utf-8"))
        if number.removeprefix("ADR-") not in numbered
    )

    assert not dangling, f"these cite a record that does not exist: {', '.join(dangling)}"


def test_no_citation_sits_inside_an_example_the_extension_would_edit() -> None:
    """A citation in a doctest is code, and deleting from code changes what the example says.

    The removal runs over the whole docstring, fences included, so a ``>>>`` line citing a
    record would be published showing a call nobody could have written. Nothing is in that
    position today, and this keeps it that way.
    """
    holders = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    fenced: list[str] = []
    for path in sources():
        for node in ast.walk(parse(path)):
            if not isinstance(node, holders):
                continue
            inside = False
            for line in (ast.get_docstring(node, clean=False) or "").splitlines():
                if line.strip().startswith("```"):
                    inside = not inside
                if ANY_NUMBER.search(line) and (inside or line.lstrip().startswith((">>>", "..."))):
                    fenced.append(f"{path.relative_to(PACKAGE.parent)}: {line.strip()[:60]}")

    assert not fenced, f"a citation inside an example would be edited by the strip: {fenced}"


def test_the_reference_page_publishes_no_source_listing() -> None:
    """The other half of how the page reaches zero, and the easier half to undo by accident.

    A listing is published verbatim, so griffe's rewrite cannot reach the citations inside
    one — thirty-five numbers rode onto the page that way. Turning listings back on returns
    them, and the docs-build guard would catch it, but only in the docs job.
    """
    config = (Path(__file__).resolve().parents[1] / "mkdocs.yml").read_text(encoding="utf-8")

    assert re.search(r"^\s*show_source:\s*false\s*$", config, re.MULTILINE), (
        "show_source must stay false, or source listings republish the citations in them"
    )


def test_the_extension_mkdocs_loads_is_the_one_under_test() -> None:
    """A rename that missed `mkdocs.yml` would leave these tests passing against nothing."""
    config = (Path(__file__).resolve().parents[1] / "mkdocs.yml").read_text(encoding="utf-8")

    assert "scripts/mkdocs_record_citations.py:StripRecordCitations" in config
    assert CITATION.pattern  # the shape the config's extension removes, imported from it
