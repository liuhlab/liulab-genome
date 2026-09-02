# 24. The template is this repository's authority, and a divergence needs a record

`liulab-repo-template` holds the shared conventions and conformance is the default; a difference is
legal only where a record says why. Authority binds the rule, not the file bytes — the template's
conformance script is a pull model, stating the rule and not the contents, so a repository that
legitimately diverges stays green while the conventions stay checked; where this repository is
ahead, the fix is to move the template forward, not to regress here. Standing: mkdocs-material, the
mature generator the template's own docs-site record prefers for a production repo; `bioconda`
behind `conda-forge` and `osx-64`, for the external tools and Intel Macs; the `docs` feature inside
`default`, plus the extra `aligners` and `test`, since the site build rides an existing job; a
context map and eight glossaries, no root `CONTEXT.md` (ADR-0004); `conftest.py` at the root, since
`testpaths` reaches `src/genome` for doctests and a `tests/` conftest cannot; the glossary cap per
entry not per file, checked here, not taught to the shared checker; and the search-exclusion rule
waived, `exclude_docs` keeping those trees off the site, so its front matter is inert.
