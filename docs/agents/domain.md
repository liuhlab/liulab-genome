# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the
codebase.

## Before exploring, read these

- **`CONTEXT-MAP.md`** at the repo root — the shared kernel of words every context uses, plus the
  index of the per-context glossaries.
- **`docs/context/<context>.md`** — one glossary per bounded context. Read the ones relevant to the
  topic; they are not meant to be read as a set.
- **`docs/adr/`** — the numbered decision records. Read the ones touching the area you are about to
  work in.
- **`docs/research/`** — dated measurements. A measurement is not a decision: the number lives here
  and whatever it *decided* lives in a record.

## File structure

This repo is **multi-context**, but not in the shape the default template assumes: the glossaries do
not sit beside the code they describe, and there is one ADR directory rather than one per context.

```
/
├── AGENTS.md          ← the router; CLAUDE.md is a symlink to it, one canonical copy
├── CONTEXT-MAP.md     ← shared kernel + index of the eight glossaries
├── docs/
│   ├── adr/           ← every decision, system-wide, numbered
│   ├── agents/        ← this file and its siblings
│   ├── context/       ← one glossary per bounded context
│   │   ├── annotation.md   assembly.md   index.md   motif.md
│   │   └── orthology.md    sequence.md   tf.md      xref.md
│   └── research/      ← dated measurements
└── src/genome/
```

**Do not create `src/<context>/CONTEXT.md` or `src/<context>/docs/adr/`.** That layout was
considered and rejected: ADR-0004 put the glossaries under `docs/context/` because they are read
together, and because TF straddles `tf/gene/`, `tf/cofactor/` and two modules beside them, so
co-locating would put that one somewhere arbitrary.

There is no `CONTEXT.md` anywhere in this repo, and its absence is not a gap to fill.

## Use the glossary's vocabulary

When your output names a domain concept — an issue title, a refactor proposal, a hypothesis, a test
name — use the term as defined, not a synonym an entry lists under *Avoid*. A concept defined
nowhere is a signal either way: usually it is language the project does not use, occasionally a real
gap worth adding.

Two further rules this repo's map sets, which the generic guidance does not:

- **Two vocabularies, and they do not mix.** Domain terms come from the glossaries. Architecture
  terms — module, interface, depth, seam, adapter, leverage, locality — are fixed, and "component",
  "service", "API" and "boundary" are not substitutes for them. One narrowing: "component" is a
  domain term in the Assembly context, so it is banned only as a substitute for *module*.
- **A term marked _(decided, not built — ADR-N)_** names something a record settled and the code
  does not have yet. Use the word; do not call the API it describes.

## The bar for adding a term or a record

Held at review, and nothing mechanises it. Before adding either, the answer must not be readable
from the code. A glossary entry is one or two sentences; a record is one paragraph, twelve lines the
ceiling, and clears all three of hard-to-reverse, surprising-without-context, and a real trade-off.
Prefer editing an existing entry to adding a neighbour.

## Name the idea, never the number

A record number is permanent and citable. The rule numbers in `AGENTS.md` are positional and
re-point silently when a row is inserted, so cite a rule by its idea and never by its number.

## Flag ADR conflicts

If your output contradicts an existing record, surface it rather than silently overriding:

> _Contradicts ADR-0007 (a broken registration raises and names its repair) — but worth reopening
> because…_
