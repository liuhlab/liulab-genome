# CLAUDE.md

> Context and working agreement for any agent (and human) contributing to this repository.
> This file is **load-bearing**: Claude Code reads it locally, and the GitHub Action respects it
> when opening PRs. Keep it concrete and current. When a correction has to be made twice,
> encode the rule here so it never has to be made a third time.

---

## 1. Fill these in before first use

- **Package name:** `liulab-genome` (import name = `genome`)
- **One-line purpose:** handling genomic files — metadata management, file processing, feature extraction.
- **Supported Python:** 3.12, 3.13, 3.14
- **Distribution:** published to PyPI via trusted publishing on a CalVer git tag.
- **Reference genome assumptions:** `<e.g. GRCh38 by default; assembly must be explicit in metadata>`

---

## 2. Toolchain — single source of truth

This project has **both Python and native (non-Python) dependencies** (e.g. `samtools`, `bedtools`),
so the environment manager is **pixi** with **conda-forge + bioconda** channels. Do NOT use `uv`,
`pip`, `python -m venv`, `poetry`, or `conda` directly. The manifest lives in `pyproject.toml`
under `[tool.pixi.*]` tables (single standards-compliant file); `pixi.lock` is committed.

Channels are configured in `[tool.pixi.workspace]` in this order (priority matters):

```toml
channels = ["conda-forge", "bioconda"]
```

All operations go through pixi tasks defined in `[tool.pixi.tasks]`:

| Task | Command | Underlying |
|------|---------|------------|
| Install env (locked) | `pixi install --locked` | resolve from `pixi.lock` |
| Add native/conda dep | `pixi add <pkg>` | conda-forge/bioconda |
| Add PyPI dep | `pixi add --pypi <pkg>` | PyPI |
| Activate shell | `pixi shell` | |
| Lint | `pixi run lint` | `ruff check .` |
| Format | `pixi run fmt` | `ruff format .` |
| Format check | `pixi run fmt-check` | `ruff format --check .` |
| Type check | `pixi run typecheck` | `pyright` |
| Tests | `pixi run test` | `pytest` |
| Tests + coverage | `pixi run test-cov` | `pytest --cov=<import_name> --cov-report=term-missing` |
| All gates | `pixi run check` | lint + fmt-check + typecheck + test |
| Build wheel/sdist | `pixi run build` | `python -m build` (hatchling backend) |
| Docs (live) | `pixi run docs` | `mkdocs serve` |
| Tool doctor | `pixi run -- <import_name> doctor` | verifies samtools/bedtools |

**Before claiming a task is done, `pixi run check` must pass.** Those are the exact gates CI runs.

### Native dependencies (bioconda)
- `samtools` and `bedtools` are runtime conda dependencies and must resolve from bioconda.
- The package shells out to these binaries; never reimplement what they do. Locate them on `PATH`
  (resolved by pixi), and fail with an actionable error if missing — see the `doctor` command.

### Python version matrix
Defined via pixi features + environments, used for CI:

```toml
[tool.pixi.feature.py312.dependencies] python = "3.12.*"
[tool.pixi.feature.py313.dependencies] python = "3.13.*"
[tool.pixi.feature.py314.dependencies] python = "3.14.*"
[tool.pixi.environments]
default   = { features = ["py314", "dev"], solve-group = "default" }
test-py312 = ["py312", "test"]
test-py313 = ["py313", "test"]
test-py314 = ["py314", "test"]
```

---

## 3. Repository layout

```
<repo>/
├── src/<import_name>/        # src layout — always. Code is tested as installed, not from CWD.
│   ├── __init__.py           # exports public API: DNA, RNA, Protein, ...
│   ├── py.typed              # ship type info to consumers
│   ├── seq.py                # DNA / RNA / Protein sequence classes (first feature)
│   ├── core/                 # typed Python API (real logic)
│   ├── io/                   # file readers/writers (I/O boundary)
│   ├── features/             # feature extraction
│   └── cli.py                # thin Typer wrapper over the API — NO logic here
├── tests/                    # mirrors src/ structure; hypothesis-backed where it counts
├── docs/                     # MkDocs Material; hand-authored Markdown pages
├── skills/<import_name>/     # SKILL.md teaching agents to USE this package
├── .github/workflows/        # ci.yml + release.yml + claude.yml
├── pyproject.toml            # project metadata + [tool.pixi.*] manifest + tool config
├── pixi.lock                 # committed lockfile
├── CLAUDE.md                 # this file
└── AGENTS.md                 # tool-neutral pointer to this file
```

---

## 4. Coding standards

- **Full type annotations on every public function and method.** Treat a missing annotation as a bug.
- **Docstring style: NumPy.** Every public function, class, and module gets a docstring with
  Parameters / Returns / Raises and at least one runnable example.
- Prefer pure functions and explicit data flow. Side effects and I/O live at the edges
  (`io/`, `external.py`, `cli.py`), not in `core/` / `features/` / `seq.py`.
- Errors must be actionable: say what was wrong and what the caller should do. Raise specific
  exception types, not bare `Exception`.
- Optimize for the reviewer, not the character count.

---

## 5. Domain invariants — READ BEFORE TOUCHING GENOMIC LOGIC

> These are the landmines. Edit to match real scope, but keep them explicit.

- **Coordinate systems:** internally store all intervals as **0-based, half-open** `[start, end)`.
  Convert to/from 1-based-inclusive (VCF, GFF/GTF, SAM) **only at I/O boundaries**. Any function
  crossing this boundary must say so in its docstring.
- **Assembly/build is never assumed.** Reference build (GRCh37 vs GRCh38) is carried in metadata and
  validated; mixing builds is an error, not a warning.
- **Chromosome naming:** normalize `chr1` vs `1` at ingest; document the canonical internal form.
- **Strand** is `+` / `-` / `.` (unknown). Never default silently to `+`.
- **Large files stream.** Do not read whole genomic files into memory; assume bgzipped + indexed.
- **Metadata is first-class:** sample IDs, assembly, provenance travel *with* the data; dropping
  metadata across a transform is a bug.

---

## 6. Testing — non-negotiable

- **Every new feature ships with tests in the same PR.** A feature without tests is not done.
- Write the test from the spec **first**, confirm intent, then implement until green.
- Use **pytest**. Use **hypothesis** for parsers, coordinate conversions, sequence transforms, and
  feature extractors — assert invariants over generated inputs (e.g. `reverse_complement` is an
  involution; length is preserved; round-trips are identity).
- Keep small real fixtures under `tests/data/`; never commit large genomic files — subsample.
- Tests that require `samtools`/`bedtools` run inside the pixi env (binaries present). If a test is
  written to also run outside pixi, gate it with a skip when the binary is absent.
- Coverage is a signal, not a target; new/changed lines should be covered without hollow tests.

---

## 7. Documentation rules — Markdown-first

- Docs are **hand-authored Markdown** under `docs/` (MkDocs Material), so they are easy to read and
  edit directly. Each feature gets a narrative page with runnable examples.
- `mkdocstrings` may auto-render an API reference from docstrings, but the **primary, editable doc is
  the Markdown prose page** — don't rely solely on autogen.
- A change to any public signature **requires** updating both its docstring and the relevant Markdown
  page in the same change. Stale docs are treated as broken code.
- Keep `skills/<import_name>/SKILL.md` and the usage reference current when the CLI or API changes.

---

## 8. Versioning & releases — CalVer

- Scheme: **`YYYY.MM.MICRO`** (e.g. `2026.6.0`, then `2026.6.1`).
- Version derived from the git tag (via `hatch-vcs`); **a tag is a release.** Never hand-edit a version.
- Tagging `vYYYY.MM.MICRO` triggers `release.yml`: build with hatchling, publish to PyPI via OIDC
  trusted publishing (no API tokens).
- Maintain `CHANGELOG.md` (Keep a Changelog); update the Unreleased section in feature PRs.

---

## 9. CLI & Skill conventions (agent-usability)

- `cli.py` is a **thin** Typer wrapper; logic stays in the API so `import <import_name>` and the CLI
  hit the same code paths.
- Every command supports `--json` for machine-readable output; human-readable is the default.
- Non-zero exit codes on failure, with errors that name the next action.
- A `doctor` command reports whether `samtools`/`bedtools` are discoverable and their versions.
- The Skill (`skills/<import_name>/SKILL.md`) teaches agents *when* and *how* to use the package via
  the CLI, with worked examples and the domain invariants above. Keep its `description` frontmatter
  specific so it triggers on the right phrasing.

---

## 10. Git, PRs, and scope

- **Small, single-purpose PRs.** One issue → one focused PR. If a diff sprawls, the issue was
  underspecified — split it.
- Conventional Commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`).
- Never commit secrets or large data files. Commit `pixi.lock`; respect `.gitignore` (ignore `.pixi/`).
- Mergeable only when CI is green **and** a human has reviewed intent and domain invariants — green
  CI is necessary, never sufficient.

---

## 11. Don't do this

- Don't bypass pixi or add a `pip install` / `conda install` step.
- Don't reimplement samtools/bedtools functionality in Python; shell out to the binaries.
- Don't read whole genomic files into memory.
- Don't add a dependency without justifying it in the PR; prefer stdlib and existing deps.
- Don't perform unrequested refactors or reformat unrelated files.
- Don't weaken or delete a failing test to make CI pass; fix the cause.
- Don't assume coordinate system, assembly, or strand. Make them explicit.
