"""Make a record number the reference page publishes resolve, without publishing the record.

`mkdocs.yml` keeps `adr/`, `agents/`, `context/` and `research/` out of the built site on
purpose: the site is the human layer on top of those trees, not a rendering of them.
`mkdocstrings` does not respect that boundary, because it renders `src/` docstrings — and
those cite records by number, which is exactly how repo-internal prose should cite a
decision. Seventy-six citations reached the reference page that way, and a reader outside
the repository could resolve none of them.

Two readers, one docstring. In the source the number resolves, because the tree is right
there. On the site it did not. Stripping it would have fixed the site reader by taking the
pointer away from the source reader, who is the one it was written for; publishing the tree
would have made records a site surface and changed what writing one costs. So the citation
is rendered differently per destination instead: the source keeps the number, and the site
turns it into a link to the record on GitHub. Nothing is stripped and no record becomes a
page.

**This edits griffe's model, never the module.** The rewrite lands on the docstring griffe
parsed for rendering; `genome.__doc__` at runtime and the "Source code in …" block on the
page beside it are both untouched, which is what keeps the source listing a source listing.

A citation naming a record that does not exist fails the build. `mkdocs build --strict`
runs in CI, so a number invented in a docstring — or one left behind by a record that was
renumbered — is caught where it is cheap, rather than shipping as a link into nowhere.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from griffe import Extension, Object

#: A citation as repo-internal prose spells it. The four digits are the record's permanent
#: name; the file around them is not, which is why the number is what a docstring carries.
CITATION = re.compile(r"ADR-(\d{4})")


class RecordLinks(Extension):
    """Rewrite ``ADR-0006`` into a Markdown link to that record, in every rendered docstring.

    Parameters
    ----------
    records : str
        The record tree, relative to the repository root.
    base_url : str
        Where the repository is browsable, up to and including the branch — the records are
        linked there rather than on this site, since publishing them here is the thing
        `mkdocs.yml` declines to do.

    Examples
    --------
    >>> links = RecordLinks(records="docs/adr", base_url="https://example.invalid/blob/main")
    >>> links.render("checksums are taken over unpacked content (ADR-0006).")
    'checksums are taken over unpacked content ([ADR-0006](https://example.invalid/blob/main/docs/adr/0006-checksums-are-taken-over-unpacked-content.md)).'
    """

    def __init__(self, *, records: str = "docs/adr", base_url: str = "") -> None:
        self._base_url = base_url.rstrip("/")
        self._records = records.strip("/")
        root = Path(__file__).resolve().parent.parent / self._records
        self._files = {path.name[:4]: path.name for path in root.glob("[0-9]" * 4 + "-*.md")}
        if not self._files:
            raise FileNotFoundError(
                f"no records found under {root}; RecordLinks would silently link nothing. "
                f"Point `records` at the record tree, relative to the repository root."
            )

    def render(self, text: str) -> str:
        """Return ``text`` with every record citation in it turned into a Markdown link.

        Raises
        ------
        KeyError
            If a citation names a record with no file. The build fails rather than
            publishing a link into nowhere — a renumbered record is the usual cause.
        """

        def link(match: re.Match[str]) -> str:
            number = match.group(1)
            if number not in self._files:
                raise KeyError(
                    f"{match.group(0)} is cited but no record is numbered {number} under "
                    f"{self._records}. Cite a record that exists, or renumber the file."
                )
            return f"[{match.group(0)}]({self._base_url}/{self._records}/{self._files[number]})"

        return CITATION.sub(link, text)

    def on_object(self, *, obj: Object, **kwargs: Any) -> None:
        """Rewrite one object's docstring as griffe hands it over, before anything renders it."""
        if obj.docstring is not None and CITATION.search(obj.docstring.value):
            obj.docstring.value = self.render(obj.docstring.value)
