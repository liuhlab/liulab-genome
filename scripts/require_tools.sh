#!/usr/bin/env bash
# `pixi run check` and `pixi run test` run this FIRST, and the whole point is that it
# runs first.
#
# samtools, faToTwoBit and twoBitInfo are what preparing an assembly costs. Every test
# that writes a real FASTA, `.2bit`, chrom.sizes or annotation database needs all three,
# and every one of those tests carries a `shutil.which` skip of its own — in the chimera
# fixture, and in seven test modules. Each skip is right on its own and stays. Together
# they are the failure: measured, the unit lane with these three off PATH reports 136
# skipped, 694 passed, exit 0 — 16% of the lane, reported as a pass. A skip is green.
# `test-aligner` already had this half of the pattern for its seven tests (see
# require_aligners.sh); the lane that runs the other 830 did not, so "check is green"
# meant 694 of them.
#
# The list is not a second copy of anything: it is `REQUIRED_TOOLS` in
# `src/genome/external.py` — the **External tool**s every run of this package needs,
# which is what `genome doctor` checks and what tests/test_external.py holds equal to
# `PREPARATION_TOOLS`. STAR and chromap are deliberately absent: they are optional
# features and the aligner lane's business. So is bedtools, and for a different reason —
# it is declared in `[tool.pixi.dependencies]` and the External tool table knows what
# installs it, but nothing here shells out to it yet, so no test can fail for its
# absence and requiring it would fail the gate over a binary the suite never runs.
#
# So the binaries are proved to ANSWER rather than merely to resolve; `which` says a file
# is on PATH, running it says it executes, and on a half-installed environment those are
# different answers. CI runs the same verbs, so it is protected by the same probe rather
# than by a second copy of this list.
#
# The probe is not one command, because the tools do not agree on one. samtools answers
# `--version` with exit 0, as STAR and chromap do. The two UCSC binaries reject every
# flag they are offered — `--version`, `-h` and `-help` each print the usage block and
# exit 255 — so no invocation of them exits 0 and "exit 0" cannot be the test. What can:
# run one bare and read who answered. The shell returns 127 for a command it could not
# find and 126 for one it could not execute; any other status came from the program
# itself. Measured on this pair: 255 installed, 127 absent. A bare run is also free of
# consequences — with no arguments both print their usage and exit at once, reading no
# stdin and writing no file.
#
# `set -e` is absent deliberately: a probe that fails must be RECORDED, not fatal, so the
# message below can name every missing binary at once instead of one per run. And nothing
# here may use an empty array — macOS ships bash 3.2, where expanding one under `set -u`
# is an unbound-variable error.
set -uo pipefail

missing=""

if ! samtools --version >/dev/null 2>&1; then
    missing="$missing samtools"
fi

for binary in faToTwoBit twoBitInfo; do
    "$binary" >/dev/null 2>&1
    status=$?
    if [ "$status" -eq 126 ] || [ "$status" -eq 127 ]; then
        missing="$missing $binary"
    fi
done

if [ -n "$missing" ]; then
    cat >&2 <<EOF
unit lane: these binaries do not answer:$missing

They prepare an assembly. Every pixi environment declares them, so one of two things is
true: this environment is out of date with the manifest, or the gate was run beside the
environment rather than inside it.

    pixi install             # re-solve this environment against pyproject.toml
    pixi run check           # run the gate inside it

Refusing to select. A missing binary here is not a smaller test run — it is a green one
that never wrote a FASTA, which is the failure this guard exists to end.
EOF
    exit 1
fi
