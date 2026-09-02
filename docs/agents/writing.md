# Writing rules

Three rules. Two are checked by `vale`; the third is on you. Run the checker with
`pixi run vale`, and `pixi run markdownlint` beside it — that one reads structure where vale
reads language.

## Be concise

Shorter beats longer — in documents, issues, commit messages, and replies. If a sentence
survives deletion without loss, delete it. The caps below are ceilings, not targets.

## Agent-facing documents have word caps

| Files | Cap |
| --- | --- |
| `AGENTS.md`, and `CLAUDE.md` which is a symlink to it | 1000 (`Lab.LengthDoc`) |
| `docs/agents/` | 1000 (`Lab.LengthDoc`) |
| `docs/adr/` | 400 (`Lab.LengthAdr`) |
| `docs/research/` | none — a measurement is as long as its method |
| `CONTEXT-MAP.md` and `docs/context/` | none — capped per entry instead, see below |

The same six paths are declared under `[tool.liulab.agent-docs]` in `pyproject.toml`, and
`pixi run conformance` pairs each to the `.vale.ini` section that carries its rule. Add a
document to one of these trees and it is capped; add a tree and it is capped by nothing until
both files name it.

**The context map and the glossaries take no file cap.** A glossary grows one term at a time,
so a file cap on one passes or fails according to how many terms the domain happens to have —
here, a shared kernel plus eight per-context glossaries — and that count is a fact about the
domain, not a writing-quality signal. The cap that measures the writing is per entry, 200 words,
and `tests/test_glossary_entries.py` holds it — including the half that matters, which fails a
glossary file parsing to no entries at all, so a convention change cannot report green by
checking nothing. That is what a cap on a glossary exists to catch, a definition that has grown
into an explanation, and the file total never sees it (ADR-0024).

Vale counts words its own way, skipping table markup, code spans and link syntax, so `wc -w`
overstates — measured at 18% to 40% across this repo's agent docs, and up to 65% elsewhere,
worst on the tables and command references most likely to sit near a cap. Believe
`pixi run vale`, and never size a document, or plan a trim, on `wc`'s say-so.

The caps are dials with tight defaults. Raising one is a one-line diff in `styles/Lab/` — do
that deliberately, and say why in the commit message. Do not raise a cap because a document
ran long; that is the cap working.

## Human-facing prose avoids jargon and stays readable

`README.md`, `CHANGELOG.md` and every `docs/` page the site publishes: no terms from the lab
jargon list (`Lab.Jargon` — architecture-speak and Latinate verbs), and reading grade 11 or
below (`Lab.Readability`). No length cap — a tutorial is as long as the task.

The grade is one number for the whole file, so no single line owns the failure. Shorten the
longest sentences, and never shorten a domain term to buy grades — that is the trade the two
rules exist to prevent. `CHANGELOG.md` sits in this set and is the file the arithmetic bites:
it only ever grows, every entry moves the one number, and a long-sentence entry spends
headroom the entries after it will need.

**This rule also covers what you say, not only what you write.** Nothing checks a chat reply,
so it is on you: when a human asks, answer in plain language and explain the term you would
otherwise reach for. An unexplained term in conversation is the same failure as one in the
README.

## Which is which

Agent-facing means written for a machine that has to act: capped, exempt from the jargon and
readability rules, and kept out of the built site — by `exclude_docs` in `mkdocs.yml`, which
is why no page here carries the search-exclusion front matter the shared conformance rule asks
for. Human-facing means written for a person reading the published site: checked for jargon
and readability, uncapped. A file is one or the other — if you are adding a document and
cannot tell, it is human-facing.
