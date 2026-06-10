# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Calendar Versioning](https://calver.org/) using
`YYYY.MM.MICRO` (e.g. `2026.6.0`).

## [Unreleased]

### Added

- Initial pixi scaffold: `pyproject.toml` with `[tool.pixi.*]` manifest, conda-forge +
  bioconda channels, `samtools`/`bedtools` runtime deps, py312/py313/py314 environments,
  and standard tasks (`lint`, `fmt`, `typecheck`, `test`, `check`, `build`, `docs`).
- Package skeleton: `src/genome/{__init__.py, py.typed, cli.py}` with a stub Typer app.
- MIT license, README, AGENTS.md pointer.
- Quality gates: ruff rule set (E, W, F, I, UP, B, C4, SIM, PT, PTH, N, D, RUF) with
  numpy docstring convention; pyright in basic mode targeting py3.12 minimum; pytest
  with `--strict-config`, `xfail_strict`, and warnings-as-errors; branch coverage config.
- `.pre-commit-config.yaml`: ruff hooks (pinned) + pyright as a local pixi-backed hook.
- Typed biological sequence classes `genome.DNA`, `genome.RNA`, `genome.Protein` —
  `str` subclasses that preserve the value verbatim (case preserved), with typed
  slicing, biological transforms (`complement`, `reverse_complement`, `transcribe`,
  `back_transcribe`), and a `gc_content` property. The alphabet is documentation only,
  not enforced at construction (too costly on large sequences); validate at the I/O
  boundary instead. IUPAC ambiguity codes are out of scope.
- `genome.external` — native tool discovery (`samtools`, `bedtools`) via `shutil.which`
  + `subprocess`, with a `doctor()` function and `ToolNotFoundError` carrying an
  actionable message.
- `genome.io` package (the I/O boundary):
  - `genome.io.download` — `Downloader`, a pooch-backed cache for large files, and
    `UCSCGenomeDownloader`, which fetches and decompresses a reference-genome FASTA
    (`<assembly>.fa.gz`) from the UCSC golden path given an assembly name (e.g. `hg38`).
    `UCSCGenomeDownloader` stores each assembly under `<LIULAB_DATA>/genome/<assembly>/`
    (configurable via the `LIULAB_DATA` env var, default `~/liulab_data`) so all
    reference files for a build are co-located; exposes `liulab_data_dir` and
    `assembly_data_dir` helpers. `UCSCGenomeDownloader.fetch_genome` runs the full
    pipeline end to end — download `<assembly>.fa.gz`, decompress, `samtools faidx`,
    `faToTwoBit`, and `twoBitInfo` — returning a `GenomeFiles` record (the 2bit is
    self-indexed, so its internal index is surfaced as `chrom.sizes`).
    `UCSCGenomeDownloader.validate_assembly` (exposed via `assembly_url`) checks the
    assembly is a real golden-path directory with an HTTP `HEAD` before downloading;
    `fetch_fasta`/`fetch_genome` run it by default (`validate=True`) so a mistyped
    assembly fails fast with an actionable `ValueError` instead of an opaque FASTA 404.
  - `genome.io.fasta` — `prepare_fasta`, which indexes a FASTA (`samtools faidx`),
    converts it to 2bit (`faToTwoBit`), and writes `chrom.sizes` (`twoBitInfo`),
    returning a `GenomeFiles` record; plus the individual `faidx`, `fasta_to_2bit`, and
    `twobit_to_chrom_sizes` steps. Each step is cached at the command-running step: the
    native tool is skipped when its output already exists and is newer than its input
    (a `make`-style freshness check), making re-runs cheap and idempotent. An
    `overwrite=True` flag (on every step, `prepare_fasta`, and `fetch_genome`) forces
    regeneration; the download/decompression cache is handled independently by pooch.
- Runtime dependencies: `pooch`, `requests`, and `tqdm` (PyPI); `ucsc-fatotwobit` and
  `ucsc-twobitinfo` (bioconda) native tools.
- CLI commands `genome revcomp <seq> [--json]` and `genome doctor [--json]`; the existing
  `genome version` stub remains.
- `docs/sequences.md` — hand-authored narrative documentation for the sequence classes.
- Tests: cases covering construction (alphabet not enforced), typed slicing, transforms,
  GC content, plus hypothesis property tests (involution, round-trip, slice closure,
  verbatim construction) and CLI/Typer integration. Tests requiring native binaries skip
  cleanly when they are absent.
- GitHub Actions:
  - `ci.yml` — lint/typecheck in the `default` env (Python-version-independent),
    pytest matrix over `test-py312/313/314`, and a wheel/sdist build job, all
    bootstrapped via `prefix-dev/setup-pixi@v0.8.14` with `locked: true`.
  - `release.yml` — on `v*` tag, build with `pixi run build` and publish to PyPI via
    OIDC trusted publishing (no API tokens; uses `pypa/gh-action-pypi-publish`).
  - `claude.yml` — responds to `@claude` mentions in issues / PR comments / PR reviews,
    with a pixi setup step so the agent can run `pixi run check` in CI.
- MkDocs Material documentation site (`mkdocs.yml`, `docs/{index,usage,reference}.md`)
  with the `mkdocstrings-python` plugin pulling NumPy-style docstrings into the
  secondary API reference page. Live preview via `pixi run docs`.
- `.github/workflows/docs.yml` — builds the MkDocs site and deploys it to
  GitHub Pages on every push to `main` (and on manual dispatch). Uses the
  official `actions/upload-pages-artifact` + `actions/deploy-pages` flow with
  a single `pages` concurrency group. `ci.yml` gains a `docs` job that runs
  `mkdocs build --strict` on PRs so doc breakage is caught before merge.
- `ruff` added to the `docs` feature so mkdocstrings can format rendered signatures.
- `skills/genome/SKILL.md` — agent-usability skill teaching coding assistants
  when to reach for `liulab-genome`, the CLI surface (`revcomp`, `doctor`),
  the typed sequence API, what does *not* stay typed (gotcha), and the
  project's domain invariants.
- `.github/ISSUE_TEMPLATE/feature_request.yml` — structured form prompting for
  the public signature, behavior, edge cases, applicable CLAUDE.md domain
  invariants, worked examples, out-of-scope notes, and a required
  definition-of-done checklist (types + NumPy docstring + tests + Markdown
  doc + skill update + CHANGELOG + green CI). Designed so issues are tight
  enough to hand directly to `@claude`.
