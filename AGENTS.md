# AGENTS.md

The router — `CLAUDE.md` is a symlink to this file: one canonical copy, no fork. Read it in full;
everything else is looked up. Terms are defined once in [`CONTEXT-MAP.md`](CONTEXT-MAP.md) — a
shared kernel of the words every context uses, plus one glossary per bounded context under
[`docs/context/`](docs/context/) — and decisions once in [`docs/adr/`](docs/adr/). This file points
at both and restates neither.

## What this is

`genome` owns **reference-genome files**. Name an assembly and it prepares that assembly on disk —
the FASTA and everything derived from it (`.fai`, `.2bit`, `chrom.sizes`) — then answers sequence
queries by region, registers GTF annotations against the assembly, and builds aligner indexes. The
work itself belongs to native binaries: `samtools`, `bedtools`, STAR, chromap. This package locates
them, drives them and reports what they did. It reimplements none of them.

## Non-negotiable rules: R1–R9

Imperatives only. Each row cites the records that decided it **by number, never by path** — a record
moves between directories and its number does not; `—` means no record cites that rule today.

| # | Rule | Records |
|---|---|---|
| R1 | **pixi only.** Never `pip`, `uv`, `conda`, `venv` or `poetry`. The manifest is `[tool.pixi.*]` in `pyproject.toml`; `pixi.lock` is committed; channels are `conda-forge` then `bioconda` and the order is priority. Never commit secrets or large data files. | ADR-0001 |
| R2 | **Shell out, never reimplement.** `samtools`, `bedtools`, STAR and chromap are **External tool**s: locate on `PATH`, fail naming the install command. | — |
| R3 | **Never read a whole genomic file into memory.** | — |
| R4 | **0-based half-open internally.** Convert to 1-based-inclusive (VCF, GFF/GTF, SAM) only at the I/O boundary, and say so in the docstring. | — |
| R5 | **Never assume assembly, coordinate system, strand or species.** The assembly id travels with the data and is never inferred; mixing builds is an error, not a warning; strand `.` is unknown and never silently `+`; a gene id never crosses species implicitly — homology is answered when a caller asks for it and never applied on their behalf. | ADR-0003, ADR-0019 |
| R6 | **Default to private.** Only `__init__.py` re-exports and the CLI surface are public; promotion is a one-line refactor. | — |
| R7 | **The CLI is a thin client.** Logic lives in the API so `import genome` and the CLI hit one code path; every command emits `--json`; non-zero exit on failure, with errors that name the next action. | — |
| R8 | **A feature without tests is not done.** | — |
| R9 | **Full type annotations on every public function and method.** pyright runs `basic`, so nothing catches a missing one. | — |

## Where to read next

- **A term, or a synonym to avoid** — [`CONTEXT-MAP.md`](CONTEXT-MAP.md) and the eight per-context
  glossaries it lists under [`docs/context/`](docs/context/).
- **A decision — why it is this way, and what lost** — [`docs/adr/`](docs/adr/).
- **A measurement — a number, and the method that produced it** — [`docs/research/`](docs/research/),
  dated. A measurement is not a decision: it goes here, and whatever it *decided* goes to a record.

**A term or a record is the exception.** Before adding either, the answer must not be readable from
the code. A glossary entry is one or two sentences; a record is one paragraph, twelve lines the
ceiling, and clears all three of hard-to-reverse, surprising-without-context, and a real trade-off —
prefer editing an existing entry to adding a neighbour. **This bar is held at review and nothing
mechanises it.**

## Working here

- **Python support points at both ends of the range on purpose.** Tests run on 3.13 only; the 3.12
  floor is held by `ruff target-version` and `pyright pythonVersion`, not by a test lane. Do not
  narrow `requires-python` to match what is tested, and do not raise the language level to match it.
- **Tests.** `pixi run check` is the gate — lint, fmt-check, typecheck, test, run concurrently.
  Write the test from the spec first, then implement until green. pytest, and hypothesis for parsers,
  coordinate conversions and sequence transforms: assert invariants over generated inputs —
  `reverse_complement` is an involution, length is preserved, round-trips are identity. Tests mirror
  `src/`. Fixtures are small real files under `tests/data/`, subsampled; never a large genomic file.
  Coverage is a signal, not a target.
- **The suite is two lanes, and together they are all of it.** `-m aligner` against its negation is
  total, so no test falls between them. A test needing `samtools` or `bedtools` runs inside the pixi
  env — gate it with a skip when the binary is absent. A test needing **STAR or chromap** is
  different: those live in optional features, so gate it with `_needs` in
  `tests/aligner/test_aligner.py`, which applies the marker and the skip under one name. Never the skip alone — that test would sit in
  the unit lane and skip itself green, which is what the split exists to end.
- **Docstrings: NumPy structure is mechanised** — ruff selects `D` with the numpy convention — so the
  bar is what ruff cannot check. At least one runnable example on a public object, and *runnable* is
  mechanised too: the unit lane runs every example, so one that drifts from the code fails the build.
  A line that cannot run offline carries `# doctest: +SKIP` **itself**: a directive does not carry
  to the line below it, and a trailing comment that is prose reads exactly like one and is not.
  Together that is how the only line in the package that would download a genome ended up the only
  line not skipped. Lines that do run stay unskipped beside it — most such blocks are mixed. A subclass docstring describes
  only what differs from what it overrides, never the shared prose again. A short one-liner is enough
  for a small `_helper`.
- **Errors are actionable:** say what was wrong and what the caller should do, with a specific
  exception type, never a bare `Exception`.
- **Side effects live at each context's own edge** — `assembly/`, `annotation/`, `store/`,
  `external.py`, `cli.py` — never in the pure core, which is `region.py` and `seq.py`.
- **Name the idea, never the number.** No rule number and no document section number anywhere outside
  the Records column above. A record number is permanent; an R-number is positional and re-points
  silently when a row is inserted. Held at review; no guard test.
- **Docs are hand-authored Markdown** under [`docs/`](docs/) (MkDocs Material) — the prose page is
  primary, `mkdocstrings` secondary — and stale committed docs are broken code. `docs/adr/`,
  `docs/agents/`, `docs/context/` and `docs/research/` are agent-facing and excluded from the
  built site, and **no page the site publishes carries a record number** — it would name a file
  that reader cannot open. **Keep citing records by number in a docstring**, which is agent-facing:
  cite as a trailing parenthetical, never as a noun in the sentence, and `mkdocstrings` drops the
  citation as it renders. A number naming no record fails the build.
- **Versioning: CalVer `YYYY.MM.MICRO`, from the git tag.** A tag is a release; never hand-edit a
  version. Update the Unreleased section of `CHANGELOG.md`.
- **Git.** Small, single-purpose PRs — one issue, one focused diff. Conventional Commits. Green CI is
  necessary, never sufficient.

## Agent skills

Configuration the engineering skills read before they act. Each file is the contract; this
section points at it and restates none of it.

### Issue tracker

GitHub issues on `liuhlab/liulab-genome`, driven by `gh`. See
[`docs/agents/issue-tracker.md`](docs/agents/issue-tracker.md).

### Triage labels

The five canonical roles, each label string equal to its name. See
[`docs/agents/triage-labels.md`](docs/agents/triage-labels.md).

### Domain docs

Multi-context, and not in the shape the default assumes: `CONTEXT-MAP.md` and eight glossaries
under `docs/context/`, against one `docs/adr/`. There is no `CONTEXT.md` and its absence is not
a gap. See [`docs/agents/domain.md`](docs/agents/domain.md).
