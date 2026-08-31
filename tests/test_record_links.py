"""The record-link extension: a number the reference page publishes has to resolve.

`mkdocs.yml` keeps the record tree out of the built site, and `mkdocstrings` publishes
`src/` docstrings that cite records by number. Those two facts together put 76 citations on
the reference page that a reader outside the repository could not follow. The extension
under test renders each one as a link to the record on GitHub, so the source reader keeps
the number and the site reader gets somewhere to go.

Tested here rather than left to `mkdocs build --strict`, which only runs in the docs job:
a renumbered record should fail in the fast lane, next to the docstring that cites it.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

from mkdocs_record_links import CITATION, RecordLinks

from ._sources import PACKAGE, parse, sources

BASE_URL = "https://example.invalid/blob/main"

#: The record tree, which is in the repository and deliberately not in the built site.
RECORDS = Path(__file__).resolve().parents[1] / "docs" / "adr"


@pytest.fixture
def links() -> RecordLinks:
    """The extension as `mkdocs.yml` configures it, but pointed at a throwaway host."""
    return RecordLinks(records="docs/adr", base_url=BASE_URL)


def test_a_citation_becomes_a_link_to_that_records_own_file(links: RecordLinks) -> None:
    """The four digits are the record's permanent name; the filename around them is not."""
    rendered = links.render("nesting is forbidden by the model (ADR-0008).")

    assert rendered == (
        "nesting is forbidden by the model ([ADR-0008]"
        f"({BASE_URL}/docs/adr/"
        "0008-a-chimera-is-an-assembly-identified-by-its-component-set.md))."
    )


def test_prose_that_cites_nothing_is_returned_unchanged(links: RecordLinks) -> None:
    """The extension runs over every docstring in the package, so it must be inert."""
    prose = "Return the reverse complement, which is an involution."

    assert links.render(prose) == prose


def test_a_citation_naming_no_record_fails_the_build(links: RecordLinks) -> None:
    """A link into nowhere is worse than the bare number it replaced, so it raises.

    The usual cause is a record that was renumbered while a docstring kept the old
    number — which is exactly the drift the citation-by-number convention exists to
    survive, and the one thing it cannot survive on its own.
    """
    with pytest.raises(KeyError, match="no record is numbered 9999"):
        links.render("as decided in ADR-9999.")


def test_pointing_it_at_a_tree_with_no_records_is_refused(tmp_path: Path) -> None:
    """Otherwise it would link nothing, quietly, and the page would look fine."""
    (tmp_path / "docs").mkdir()

    with pytest.raises(FileNotFoundError, match="no records found"):
        RecordLinks(records=str(tmp_path.name), base_url=BASE_URL)


def test_the_rewrite_does_not_reach_the_module_a_caller_imports() -> None:
    """The site's copy is rewritten; ``help()`` and the source listing beside it are not.

    This is what keeps the two readers separate rather than trading one for the other. The
    extension edits the docstring griffe parsed for rendering, and griffe parsed it from
    the file — so the live ``__doc__`` is a different string that nothing here touches.
    """
    from genome.assembly import chimera_build

    docstring = chimera_build._check_not_nested.__doc__ or ""
    RecordLinks(records="docs/adr", base_url=BASE_URL).render(docstring)

    assert "](" not in (chimera_build._check_not_nested.__doc__ or "")


def test_every_record_the_package_cites_exists() -> None:
    """A renumbered record fails here, in the fast lane, and not only in the docs build.

    The guard is over the whole package rather than over the docstrings mkdocstrings
    happens to publish today: a module promoted onto the reference page tomorrow brings
    its citations with it, and this should already have caught a dangling one.
    """
    numbered = {path.name[:4] for path in RECORDS.glob("[0-9][0-9][0-9][0-9]-*.md")}
    assert numbered, f"no records under {RECORDS}"

    dangling = sorted(
        f"{path.relative_to(PACKAGE.parent)}: ADR-{number}"
        for path in sources()
        for number in CITATION.findall(path.read_text(encoding="utf-8"))
        if number not in numbered
    )

    assert not dangling, f"these cite a record that does not exist: {', '.join(dangling)}"


def test_no_docstring_hides_a_citation_where_the_link_cannot_render() -> None:
    """A citation inside a code example would be published as literal Markdown.

    Found the hard way: a stray paragraph after a ``Returns`` section was parsed as a
    return *type* and rendered as code, so its link came out as `[ADR-0008](…)` on the
    page. Nothing here is inside a fence today, and this keeps it that way.
    """
    holders = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
    fenced: list[str] = []
    for path in sources():
        for node in ast.walk(parse(path)):
            if not isinstance(node, holders):
                continue
            doc = ast.get_docstring(node, clean=False) or ""
            inside = False
            for line in doc.splitlines():
                if line.strip().startswith("```"):
                    inside = not inside
                if CITATION.search(line) and (inside or line.lstrip().startswith((">>>", "..."))):
                    fenced.append(f"{path.relative_to(PACKAGE.parent)}: {line.strip()[:60]}")

    assert not fenced, f"a citation in a code example renders as literal Markdown: {fenced}"


def test_the_configured_base_url_is_the_repository_and_not_this_site() -> None:
    """Linking to a page here would mean publishing the record, which mkdocs.yml declines."""
    config = (Path(__file__).resolve().parents[1] / "mkdocs.yml").read_text(encoding="utf-8")
    configured = re.search(r"base_url:\s*(\S+)", config)

    assert configured is not None
    assert configured.group(1).startswith("https://github.com/liuhlab/liulab-genome/blob/")
