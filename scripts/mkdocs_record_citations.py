"""Keep record numbers off the built site, without taking them out of the source.

`mkdocs.yml` keeps `adr/`, `agents/`, `context/` and `research/` out of the built site on
purpose: the site is the human layer on top of those trees, not a rendering of them. A
record number is a citation between agents — it belongs in a docstring, a comment and the
repository's own prose, and it resolves for every reader of those because the tree is right
there. It does not belong on a page someone reaches without the repository, where it names
a file they cannot open.

`mkdocstrings` renders `src/` docstrings onto the reference page and so carries those
citations across that line. This removes them as griffe hands each docstring over, which is
why the two readers can be served by one docstring: the source keeps the number, the page
never shows it.

**This edits griffe's model, never the module.** The rewrite lands on the docstring griffe
parsed for rendering; ``genome.__doc__`` and ``help()`` are untouched, so a citation is
still there for whoever reads the source — which is the audience it was written for.

The deletion is total because the citation's shape is fixed: a record is cited as a trailing
parenthetical and never as a noun in a sentence, an invariant `tests/test_record_citations.py`
holds over `src/`. Prose that cited a record as its subject would be left ungrammatical by
this, which is why the form is a rule rather than a habit.
"""

from __future__ import annotations

import re
from typing import Any

from griffe import Extension, Object

#: A citation as repo-internal prose spells it, with the whitespace in front of it so the
#: sentence closes up rather than keeping a gap where the parenthetical stood. The wrapping
#: is what makes ``\s`` rather than a literal space right: a citation is often the first
#: thing on a line, and eating that newline rejoins the sentence exactly as a renderer would.
CITATION = re.compile(r"\s*\((ADR-\d{4}(?:,\s*ADR-\d{4})*)\)")


class StripRecordCitations(Extension):
    """Remove every record citation from the docstrings rendered onto the reference page.

    Examples
    --------
    >>> strip = StripRecordCitations()
    >>> strip.render("checksums are taken over unpacked content (ADR-0006).")
    'checksums are taken over unpacked content.'
    >>> strip.render("mixing builds is an error (ADR-0003, ADR-0005), never a warning.")
    'mixing builds is an error, never a warning.'
    """

    def render(self, text: str) -> str:
        """Return ``text`` with every record citation in it removed."""
        return CITATION.sub("", text)

    def on_object(self, *, obj: Object, **kwargs: Any) -> None:
        """Rewrite one object's docstring as griffe hands it over, before anything renders it."""
        if obj.docstring is not None and "ADR-" in obj.docstring.value:
            obj.docstring.value = self.render(obj.docstring.value)
