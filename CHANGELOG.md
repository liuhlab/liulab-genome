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
  `str` subclasses validated in `__new__` (case-insensitive, case preserved), with typed
  slicing, biological transforms (`complement`, `reverse_complement`, `transcribe`,
  `back_transcribe`), and a `gc_content` property. IUPAC ambiguity codes are out of scope.
- `genome.external` — native tool discovery (`samtools`, `bedtools`) via `shutil.which`
  + `subprocess`, with a `doctor()` function and `ToolNotFoundError` carrying an
  actionable message.
- CLI commands `genome revcomp <seq> [--json]` and `genome doctor [--json]`; the existing
  `genome version` stub remains.
- `docs/sequences.md` — hand-authored narrative documentation for the sequence classes.
- Tests: 52 cases covering construction, validation, typed slicing, transforms, GC
  content, plus hypothesis property tests (involution, round-trip, slice closure,
  alphabet rejection) and CLI/Typer integration. Tests requiring native binaries skip
  cleanly when they are absent.
