#!/usr/bin/env bash
# `pixi run test-aligner` runs this FIRST, and the whole point is that it runs first.
#
# STAR and chromap are the only binaries this package does not ship with itself: they
# live in the optional `star`/`chromap` features, and only the `aligners` environment
# carries them. The tests that need them carry the `aligner` marker AND a skip (see
# `_needs` in tests/test_aligner.py) — so in an environment without the binaries, the
# aligner lane selects three tests and skips all three, exits 0, and reports green
# having built no index at all. A skip is green; that is the whole failure.
#
# So the binaries are proved to ANSWER rather than merely to resolve. `which` says a
# file is on PATH; `--version` says it executes, and on a half-installed environment
# those are different answers. CI runs the same verb, so it is protected by the same
# probe rather than by a second copy of this list.
#
# `set -e` is absent deliberately: a probe that fails must be RECORDED, not fatal, so
# the message below can name every missing binary at once instead of one per run. And
# nothing here may use an empty array — macOS ships bash 3.2, where expanding one under
# `set -u` is an unbound-variable error.
set -uo pipefail

missing=""
for binary in STAR chromap; do
    if ! "$binary" --version >/dev/null 2>&1; then
        missing="$missing $binary"
    fi
done

if [ -n "$missing" ]; then
    cat >&2 <<EOF
aligner lane: these binaries do not answer:$missing

They come from the \`aligners\` environment, which exists to hold them:

    pixi run -e aligners test-aligner

Refusing to select. A missing binary here is not a smaller test run — it is a green one
that built no index, which is the failure this lane exists to end.
EOF
    exit 1
fi
