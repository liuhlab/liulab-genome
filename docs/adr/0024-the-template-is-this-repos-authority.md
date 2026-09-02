# 24. The template is this repository's authority, and a divergence needs a record

`liulab-repo-template` holds the shared conventions and conformance is the default; a difference is
legal only where a record says why. Authority binds the rule, not the file bytes: the conformance
script states rules, not contents, so a legitimate divergence stays green while the conventions stay
checked; where this repository is ahead, the template moves. Standing: mkdocs-material, the
template's own preference; `bioconda` behind `conda-forge` and `osx-64`, for the external tools and
Intel Macs; `docs` inside `default`, plus `aligners` and `test`; a context map and eight glossaries,
no root `CONTEXT.md` (ADR-0004); `conftest.py` at the root, since a `tests/` one cannot reach
`src/genome` doctests; the search-exclusion rule waived, `exclude_docs` keeping those trees off the
site; and the glossary cap per entry not per file, measured by this repository's own suite, since
the shared rule reads an entry as a heading, this one a bolded label. Migrating, adopt pyright
before that script: a rule function missing its trailing return is a typed error under pyright, an
unexplained `AttributeError` without it, and uncovered by a setup skill for new repositories.
