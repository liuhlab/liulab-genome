# What CI actually spends, and what the suite costs after being cut to a third

Measured 2026-08-30. Supersedes the worker-count table in
[`test-suite-parallelism-2026-08-14.md`](test-suite-parallelism-2026-08-14.md), whose numbers were
taken when the suite held 517 tests; it held 2425 two weeks later.

**Three findings. (1) CI's largest single item was not a test at all — it was numba recompiling
memelite from source on every run, 27.5 s, because the JIT cache lands inside the pixi env and
setup-pixi saves that env *before* any test has run. Redirecting the cache costs no coverage and
removes 47% of the test lane's wall. (2) The suite was cut 2425 → 933 tests with coverage
unchanged at 98%, which says most of what was removed was re-execution, not checking. (3) A
handful of tests spawn their own process pools, and distributing them by load made the lane's wall
*bimodal* — 13.5 s or 16.1 s on the same commit, decided by scheduling alone. `--dist=loadgroup`
pins each such module to one worker and is faster on both box sizes.**

This is a measurement, not a decision.

## Method

Wall-clock around the whole `pixi run` (so environment activation is inside every number, as it is
for a person typing the command), on the maintainer's box: 14 logical cores, 10 performance.
Repeated runs per configuration; where a configuration is unstable the **maximum** is reported,
because a limit is only met if the worst run meets it. Caches warm unless stated — the first run
of a session is materially slower and is called out where it matters. CI numbers are read from
`gh api` job and step timings over six consecutive runs.

## 1. Where CI's 80 seconds went

Five jobs run concurrently, so the wall is the slowest job plus a 2–4 s queue.

| Job | Wall (median) | setup-pixi | Real work |
|---|---|---|---|
| **test** | **76 s** | 6–9 s | pytest **61 s** |
| lint + typecheck + docs | 35 s | 9–12 s | ruff <1 s, pyright 11 s, mkdocs 9 s |
| test (STAR + chromap) | 28 s | 11–12 s | pytest 9 s |
| build wheel + sdist | 25 s | 11–12 s | build 5 s |

The `test` job was 76 s of an 80 s wall, and 61 s of that was pytest. Inside it, three tests were
half the lane:

| Test | CI |
|---|---|
| `test_compare.py::TestWhatIsAccepted::test_one_motif` | 27.50 s |
| `test_parquet.py::…::test_a_large_result_is_written_without_complaint` | 16.17 s |
| `test_compare.py::…::test_a_de_novo_motif_that_no_release_published` | 7.79 s |

## 2. The 27.5 s was a cold JIT, not a slow test

The same test costs ~0 s locally on a warm box, which is what gave it away. `memelite`'s scan and
compare engines are `@numba.njit(cache=True)`, and numba writes that cache beside the package it
compiled — inside `site-packages`, i.e. inside the pixi env. `setup-pixi` saves its env tarball
directly after `pixi install`, which is before any test has run, so the compiled artifacts were
never in it. Every CI run recompiled from source.

Isolated, with `NUMBA_CACHE_DIR` pointed at an empty directory and then at the same directory once
populated:

| | import | first `tomtom` call | second call |
|---|---|---|---|
| cold cache | 0.64 s | **6.97 s** | 0.00 s |
| warm cache | 0.60 s | **0.02 s** | 0.00 s |

350×, and the cache is 936 KB in 8 files. On CI the compile is larger and the cores slower, which
is why it shows up there as 27.5 s rather than 7 s. Fixed by redirecting `NUMBA_CACHE_DIR` to a
workspace path and caching it.

**The cache key must rotate, and a fixed key is a trap.** numba caches per type signature, so any
one run compiles only the specialisations it happened to reach, and `actions/cache` refuses to
overwrite an existing key — so a first run that saved a partial cache freezes it permanently and
every later run pays the missing compiles forever. This was measured the hard way: two runs of the
same commit, one restoring a 437 KB cache and paying no JIT at all, the other restoring a 474 KB
cache and *still* recompiling four tests at 24–39 s each. The larger cache was the worse one. The
key therefore carries `github.run_id`, with the lock-only prefix as a `restore-key`, so every run
saves what it learned and the cache converges on the full signature set.

A serial warm-up step was tried first and rejected: seeding by one `tomtom` call produced a 936 KB
cache against 1448 KB from a parallel test run, and the suite was *slower* afterwards (8.35 s
against 2.90 s) because the warm-up reached fewer signatures than the tests do. Seeding from the
real run is what fills the cache; the key rotation is what lets it accumulate.

## 3. The suite, cut to a third

| | before | after |
|---|---|---|
| tests collected | 2425 | 933 |
| coverage (statements) | 98% (60 missed of 4149) | 98% (61 missed of 4149) |
| test classes | 169 | 169 |
| distinct `pytest.raises` messages | 154 | 154 |

Coverage moving by one statement and one branch across a 61% cut is the finding: what was removed
was overwhelmingly *re-execution of already-covered lines* — parametrize tables walking one code
path with fifteen rows of data, and classes of six one-assert tests sharing a fixture. Class count
and refusal-message count were checked as whole-tree set differences against `main`, not per file,
so a claim that merely moved between files is not counted as lost.

## 4. Worker count and distribution, after the cut

`--dist=load`, which was right for the old suite, is no longer:

| workers | `--dist=load` min–max | `--dist=loadgroup` min–max |
|---|---|---|
| 4 | 15.91 – 16.84 s | 9.43 – 9.61 s |
| 6 | 14.44 – 15.19 s | 12.33 – 12.78 s |
| 8 | 13.50 – 16.13 s | 8.37 – 12.82 s |
| **10** | 13.35 s (best-of-2) | **8.59 – 8.71 s** |
| 12 | — | 8.92 – 9.19 s |

Two things the table says that a single best-case number would hide:

- **`load` is bimodal at 8 workers: 13.50 s and 16.13 s on the same commit.** Nine tests spawn
  their own `ProcessPoolExecutor`, so when several land on different workers at once the box runs
  ~24 processes on 14 cores. Which mode you get is decided by scheduling, not by the code. That is
  why the maximum is the number reported here: a limit met only on lucky runs is not met.
- **`loadgroup` is faster on a small runner too** — 9.4 s against 12.3 s at four workers, the
  shape a CI runner has. It was adopted for stability and kept for speed. Grouping is by module:
  `test_parallel.py` in one group, `test_motif_mixin.py` and `test_cli.py::TestMotifScan` in
  another. One group for all three was tried and is *worse* (16.3 s) — it serialises 7.4 s of work
  onto one worker while seven idle.

The fixed floor, for reference: spinning up workers and collecting without running anything is
3.77 s at `-n 8` and 1.90 s serial. Bare `pixi run` activation is 0.09 s.

## 5. Where it landed

| | before | after | limit |
|---|---|---|---|
| serial | 56.3 s | **24.0 s** | < 60 s |
| parallel (`pixi run test`) | 17–19 s | **8.4 s** | < 15 s |
| gate (`pixi run check`) | 24 s | **10 s** | — |
| tests | 2425 | **933** | < 1000 |

`--maxprocesses 10` rather than a bare `-n 10`: on a 4-core CI runner `auto` finds four and the cap
does not bind, so the number is a ceiling for a big box, never a floor for a small one.

## What was not adopted

- **A `slow` marker.** A hand-maintained list keyed on a number that changes whenever you optimise.
  The `--durations=10` already printed in every CI log makes the floor visible without a list that
  can go stale.
- **A per-test wall-clock budget asserted in the suite.** The same unchanged test measured 8.4 s
  alone and 12.8 s under contention here, so the assertion would be keyed on hardware and load
  rather than on the test.
