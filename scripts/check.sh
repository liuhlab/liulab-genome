#!/usr/bin/env bash
# `pixi run check` — the pre-PR gate, with its steps run CONCURRENTLY.
#
# It used to be `{ depends-on = [...] }`, which pixi executes sequentially: lint, then
# fmt-check, then pyright, then pytest, with nothing overlapping. Measured at 9.4 s,
# of which 3.7 s was the static preamble — ruff and pyright share no inputs with each
# other or with pytest, so none of that waiting bought anything.
#
# Each step's output is captured to its own file and printed whole, in a fixed order,
# after every step has finished. Attribution therefore gets better, not worse: a failure
# is one contiguous labelled block instead of four tools interleaved by whoever flushed
# first. A green step is worth three lines; a red one is worth all of them.
#
# Steps are invoked as `pixi run <task>` rather than by spelling their command lines
# again here, so the task table stays the one owner of what each step runs.
#
# Takes its steps as arguments so any gate can share this runner.
#
# Exit status: 1 when a step ran and failed, 2 when the gate itself could not run — a
# usage error, or a capture directory it could not make. `scripts/conformance.py` spells
# its own two the same way, so "could not run" reaches the reader as itself rather than
# as a rule nobody broke.
#
# `-m` puts each step in its own process group, which is what lets the cleanup below
# reach past `pixi` to the pytest workers underneath it. `set -e` is deliberately
# absent: the runner must collect every step's status before it reports.
set -muo pipefail

STEPS=("$@")
[ ${#STEPS[@]} -gt 0 ] || { echo "usage: check.sh <task>..." >&2; exit 2; }

out=$(mktemp -d) || exit 2
# Kill, then delete — in that order, and never one without the other. The steps write
# into $out, so removing it while any of them is alive leaves live processes writing at
# a path that is gone. Each leftover job leads its own process group, so the negated pid
# takes its children with it.
cleanup() {
    local leftover
    leftover=$(jobs -p)
    # shellcheck disable=SC2086  # jobs -p emits one bare pid per line; none may be quoted as one
    for pgid in $leftover; do kill -- "-$pgid" 2>/dev/null; done
    wait 2>/dev/null
    rm -rf "$out"
}
trap cleanup EXIT
# A signal kills the shell without running the EXIT trap, so turn the two that reach an
# abandoned gate into an ordinary exit. Otherwise Ctrl-C leaves the whole gate running.
trap 'exit 130' INT
trap 'exit 143' TERM

start=$SECONDS
pids=()
for step in "${STEPS[@]}"; do
    # shellcheck disable=SC2086  # each step is a single bare task name, never a command
    pixi run --no-progress $step >"$out/$step.log" 2>&1 &
    pids+=("$!")
done

rc=0
# Indexed by position, parallel to STEPS. Not an associative array: macOS ships bash 3.2
# as /bin/bash and 3.2 has none, and a failed `declare -A` under `set -u` would take this
# script down a path where it reports green having verified nothing.
status=()
for i in "${!STEPS[@]}"; do
    if wait "${pids[$i]}"; then
        status[i]="ok"
    else
        status[i]="FAILED"
        rc=1
    fi
done

for i in "${!STEPS[@]}"; do
    printf '\n\033[1m=== %s: %s ===\033[0m\n' "${STEPS[$i]}" "${status[$i]}"
    if [ "${status[$i]}" = "ok" ]; then
        tail -n 3 "$out/${STEPS[$i]}.log"
    else
        cat "$out/${STEPS[$i]}.log"
    fi
done

printf '\n\033[1m=== gate: '
for i in "${!STEPS[@]}"; do printf '%s=%s ' "${STEPS[$i]}" "${status[$i]}"; done
printf 'in %ss ===\033[0m\n' "$((SECONDS - start))"
exit "$rc"
