"""The glossary guard: a definition that has grown into an explanation, caught by a machine.

The context map and the glossaries take no *file* cap, on the stated ground that a glossary
is measured per entry — how long the file runs reports how many terms the domain has, which
says nothing about the writing. That waiver is declared twice, in `.vale.ini` and under
``[tool.liulab.agent-docs]``, and until this module existed the cap it points at was a claim
with no mechanism: vale checked no entry, and the shared conformance script checked none
either.

**It cannot, and should not be taught to.** Its per-entry rule keys on a `## Glossary`
section this repository does not have and recognises an entry only as a `###` heading, where
this repository writes each term as a bolded label opening a paragraph, grouped under
thematic headings. Teaching a general checker one repository's document idiom would push a
local shape into a shared tool, so this repository carries the cost of its own shape here
(ADR-0024). The cap and the arithmetic are the shared rule's — 200 words, counted as raw
whitespace tokens, markup included — so an entry measured here and an entry measured there
mean the same thing. Vale's own count is smaller and is not what this compares against.

**The invariant is asserted, never the list.** The files are discovered rather than named, so
a tenth glossary added tomorrow is checked without this being edited, and every entry in them
is measured rather than the five that were over the cap when this was written.

**A file that parses to zero entries fails.** That is the half this guard exists for. The
label shape is a convention, not a syntax, and a convention drifts: rename the bold to a
`###` heading, or put the definition on the label's own line, and a naive guard walks the
file, matches nothing, and reports green having checked nothing at all. A guard that cannot
tell *no entry is too long* from *I found no entries* is the same defect it was written to
close, one level up.
"""

from __future__ import annotations

import re
from pathlib import Path

#: The repository root: this file's parent is `tests/`.
REPO = Path(__file__).resolve().parents[1]

#: The shared kernel, at the root, and one glossary per bounded context.
CONTEXT_MAP = REPO / "CONTEXT-MAP.md"
GLOSSARY_DIR = REPO / "docs" / "context"

#: The shared conformance script's number, so the two mean the same thing.
CAP = 200

#: An entry opens with its term in bold, followed by a colon: ``**Gene id stem**:``. Anchored,
#: because the relationship bullets in the context map are ``- **Assembly → Sequence**:`` and
#: are prose about two contexts rather than a definition of one term.
LABEL = re.compile(r"^\*\*(?P<term>[^*]+)\*\*:")

#: Entries are grouped under thematic headings, so a heading ends the entry above it.
HEADING = re.compile(r"^#{1,6} ")


def glossaries() -> list[Path]:
    """Return the context map and every glossary under `docs/context/`, in path order.

    Discovered, not listed, and empty is refused rather than returned: a guard that walked
    no file would pass, and would go on passing.
    """
    found = [CONTEXT_MAP, *sorted(GLOSSARY_DIR.glob("*.md"))]
    missing = [path for path in found if not path.is_file()]
    assert not missing, f"no glossary at {missing}, so nothing there would be checked"
    assert len(found) > 1, f"no glossary files under {GLOSSARY_DIR}, so only the map is checked"
    return found


def entries(text: str) -> list[tuple[str, int]]:
    """Return ``(term, words)`` for each glossary entry in ``text``, in the order written.

    An entry runs from its own label to the next label or the next heading, whichever comes
    first, and its word count includes the label line and the trailing ``Avoid`` line — the
    whole of what a reader reads under that term. Prose before the first label belongs to no
    entry and is not counted anywhere.
    """
    found: list[tuple[str, int]] = []
    term: str | None = None
    words = 0
    for line in text.splitlines():
        label = LABEL.match(line)
        if label or HEADING.match(line):
            if term is not None:
                found.append((term, words))
            term, words = (label["term"].strip(), 0) if label else (None, 0)
        if term is not None:
            words += len(line.split())
    if term is not None:
        found.append((term, words))
    return found


def test_every_glossary_file_yields_entries() -> None:
    """Zero entries is a parser that stopped matching, not a glossary that is clean.

    Every file walked here has terms in it; a file that suddenly has none has changed its
    label convention, and the measurement below has quietly stopped measuring.
    """
    empty = [str(path.relative_to(REPO)) for path in glossaries() if not entries(path.read_text())]

    assert not empty, (
        "these glossary files parse to no entries at all, so the per-entry cap is being "
        "checked against nothing:\n" + "\n".join(f"  {name}" for name in empty) + "\n"
        "An entry is a bolded term label opening a line — `**Term**:`. If that convention "
        "has changed, change LABEL here with it; do not leave the guard reporting green."
    )


def test_no_glossary_entry_runs_past_the_cap() -> None:
    """A glossary says what a word MEANS here. Past the cap it is an explanation.

    The fix is never to raise the number: cut the entry to its definition and move what is
    left to the page that should hold it, or to a record if what grew was a decision.
    """
    over = [
        f"{path.relative_to(REPO)}: {term} is {words} words, over {CAP}"
        for path in glossaries()
        for term, words in entries(path.read_text())
        if words > CAP
    ]

    assert not over, (
        "these glossary entries have grown from a definition into an explanation:\n"
        + "\n".join(f"  {line}" for line in over)
        + "\nCut each to what the word means in this repository and move the explanation to "
        "a document the glossary can point at, or to a record if it is a decision."
    )


def test_an_entry_ends_at_the_next_label_or_the_next_heading() -> None:
    """The split is tested directly, because the count above is only as good as it is.

    Four things at once, and each has been wrong in a hand-rolled splitter: the intro prose
    belongs to no entry; a heading closes the entry above it without opening one; a
    relationship bullet is not a label, though it carries a bolded term and a colon; and the
    last entry in a file is closed by the end of the file.
    """
    glossary = """\
# Context

Intro prose that belongs to no entry at all.

- **Assembly → Sequence**: a relationship bullet, not a definition.

## Language

**Region**:
Three words here.
_Avoid_: locus

**Strand**:
Two words.

## Relationships

**Chimera**:
One word.
"""

    assert entries(glossary) == [("Region", 6), ("Strand", 3), ("Chimera", 3)]
