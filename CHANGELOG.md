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
