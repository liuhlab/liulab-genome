#!/usr/bin/env bash
# The outcome guard: no record number reaches a page a reader can open without the
# repository. It reads the built artifact rather than the sources that feed it, so it holds
# whatever the rendering path is — a docstring, a hand-written page, a plugin added later.
# Runs against `site/`, the search index included: a number indexed but not displayed is
# still published.
set -euo pipefail

site="${1:-site}"

if [[ ! -d "$site" ]]; then
  echo "no built site at '$site'. Run \`pixi run docs-build\`, which builds and then checks." >&2
  exit 1
fi

if hits=$(grep -rlE 'ADR-[0-9]{4}' "$site" 2>/dev/null) && [[ -n "$hits" ]]; then
  echo "A record number reached the built site. Records are agent-facing and are excluded" >&2
  echo "from it, so a number here names a file the reader cannot open. Found in:" >&2
  printf '  %s\n' $hits >&2
  echo >&2
  echo "In a docstring, cite the record as a trailing parenthetical — '(ADR-0006)' — which" >&2
  echo "scripts/mkdocs_record_citations.py removes when rendering. In a page under docs/," >&2
  echo "name the idea and drop the number." >&2
  exit 1
fi

echo "no record number in $site"
