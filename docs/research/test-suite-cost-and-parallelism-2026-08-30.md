# What CI actually spends, and what the suite costs after being cut to a third

Measured 2026-08-30. Supersedes the worker-count table in
[`test-suite-parallelism-2026-08-14.md`](test-suite-parallelism-2026-08-14.md), whose numbers were
taken when the suite held 517 tests; it held 2425 two weeks later.

**Three findings. (1) CI's largest single item was not a test at all — it was numba recompiling
memelite from source on every run, 27.5 s, because the JIT cache lands inside the pixi env and
setup-pixi saves that env *before* any test has run. Redirecting the cache costs no coverage and
removes 47% of the test lane's wall — but redirecting it only makes the artifacts cacheable, and
took two further fixes to make them reliably *usable*: rotating the key so a partial cache cannot
freeze, and pinning the compile target so the key stops varying with the runner's CPU model, the
second measured only after the first had shipped and kept working intermittently.
(2) The suite was cut 2425 → 933 tests with coverage
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

One set of numbers is not from that box and says so where it appears: the codegen comparison under
"What the pin costs" needed x86-64, because that is what CI runs and what the question is about.
Those were taken on a lab compute node — Intel Xeon Platinum 8468, 48 cores — inside a throwaway
pixi env pinning the same memelite 0.4.0 and numba 0.67.0, and they time library calls directly
rather than a pytest lane, so they are not comparable to the wall-clock figures elsewhere here.

## 1. Where CI's 80 seconds went

Five jobs run concurrently, so the wall is the slowest job plus a 2–4 s queue.

| Job | Wall (median) | setup-pixi | Real work |
| --- | --- | --- | --- |
| **test** | **76 s** | 6–9 s | pytest **61 s** |
| lint + typecheck + docs | 35 s | 9–12 s | ruff <1 s, pyright 11 s, mkdocs 9 s |
| test (STAR + chromap) | 28 s | 11–12 s | pytest 9 s |
| build wheel + sdist | 25 s | 11–12 s | build 5 s |

The `test` job was 76 s of an 80 s wall, and 61 s of that was pytest. Inside it, three tests were
half the lane:

| Test | CI |
| --- | --- |
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
| --- | --- | --- | --- |
| cold cache | 0.64 s | **6.97 s** | 0.00 s |
| warm cache | 0.60 s | **0.02 s** | 0.00 s |

350×, and the cache is 936 KB in 8 files. On CI the compile is larger and the cores slower, which
is why it shows up there as 27.5 s rather than 7 s. Redirecting `NUMBA_CACHE_DIR` to a workspace
path and caching it is what makes the artifacts *cacheable* — it is not, on its own, what makes
them *usable*, which took two further findings: the key rotation below, and the target pin under
"The cache was reused only when the runner's CPU model happened to match".

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

## The cache was reused only when the runner's CPU model happened to match

Measured after the redirect above shipped, because the win kept arriving intermittently and the
miss was silent.
A `pull_request` run took the test lane to 32 s with a 2.28 s slowest test; a push to `main`
logged `Cache restored from key: …` and then recompiled anyway — 38.96 s, 26.06 s, 23.63 s and
23.60 s, leaving the lane at 72 s. Deleting every cache and seeding a clean one changed nothing.

`actions/cache` was never the thing failing. numba keys each compiled overload on a tuple whose
middle element is the target description, read here out of a `.nbi` written on an M4:

```text
((float64, float64),
 ('arm64-apple-darwin25.6.0', 'apple-m4', ''),
 ('65ebbeab…', 'e3b0c442…'))
```

**The host CPU model is part of the cache key.** GitHub's runners are heterogeneous and which one
a job lands on is not something the workflow chooses, so whether the restored cache applied was
luck. There is nothing special about the ref; the earlier guess that this was branch scoping was
wrong.

The failure mode is worse than a plain miss because the two caches disagree about what a hit is.
Seeding a cache under one CPU name and re-running under another, counting with numba's own
`Dispatcher.stats`:

| run | numba's counter | index after |
| --- | --- | --- |
| seed as `apple-m4` | miss ×1 | one entry |
| again as `apple-m4` | **hit ×1** | one entry |
| as `apple-m1` | **miss ×1** | *two* entries |

Every file present and readable in the third row, and numba recompiles regardless, adding an entry
beside the old one rather than reporting anything. `actions/cache` reports on the files; numba
decides on the contents; only `--durations` says which happened.

**Fixed by pinning `NUMBA_CPU_NAME=generic` on the `test` job**, which collapses the key's middle
element to `('…', 'generic', '')` — stable across machines of one architecture.

### What the pin costs

Warm runs of `tests/tf/test_compare.py` and `tests/tf/test_scan.py`, 59 tests, three repeats each
after one cold run to fill the cache:

| | cold | warm | cache |
| --- | --- | --- | --- |
| `generic` | 10.82 s | 3.97 / 4.47 / 3.97 s | 1388 KB |
| host (`apple-m4`) | 10.73 s | 4.02 / 4.01 / 4.02 s | 1452 KB |

Indistinguishable, and all 59 pass either way, so the pin changes no result. These numbers replace
an earlier pass over the same two files that read 6.93 s generic against 7.29 s host: same
conclusion, different absolute times because that pass timed the files without the suite's
fixtures warm. The cache is also *not* quite the same size either way, as that pass reported —
1388 KB against 1452 KB, host being the larger, which is what host-specific codegen should do.

**That measurement cannot speak for CI, which is why it is not the one the decision rests on.**
It was taken on arm64, where `generic` still includes NEON — mandatory in ARMv8-A — so
generic-versus-host is nearly free *by construction*. CI is x86-64, where `generic` means the
SSE2 baseline and gives up AVX2 and AVX-512 entirely. That is the architecture the pin actually
runs on, and the gap there could have been real.

Re-measured on an Intel Xeon Platinum 8468 (Sapphire Rapids), which numba resolves to
`sapphirerapids` with the full AVX-512 set on — `avx512f`, `avx512bw`, `avx512dq`, `avx512vnni`,
`avx512fp16` and the rest — so this is the *largest* gap the pin can cost anywhere on x86-64.
memelite's two engines, one numba thread throughout so a varying thread count cannot swamp the
codegen difference, warm calls only, three repeats:

| workload | host (`sapphirerapids`) | `generic` |
| --- | --- | --- |
| `fimo`, 10 motifs × 600 bp (the suite's shape) | 12.6 / 12.2 / 12.3 ms | 12.1 / 11.9 / 15.0 ms |
| `fimo`, 20 motifs × 200 kb (333× the suite) | 30.2 / 30.0 / 29.9 ms | 30.1 / 29.8 / 29.8 ms |
| `tomtom`, 60 × 60 motifs | 69.6 / 70.4 / 70.2 ms | 67.9 / 68.5 / 68.9 ms |
| `tomtom`, 60 × 60 motifs (scan run) | 69.7 / 69.8 / 69.3 ms | 72.5 / 72.0 / 72.2 ms |

**AVX-512 buys these engines nothing measurable.** `fimo` is a dead heat at both scales and does
not diverge as the sequence grows 333×; `tomtom` lands within ±4% and falls on both sides of host
across the two runs, which is noise, not a trend. So the honest framing is no longer *no cost to
CI, unknown cost to a production scan* — it is **no cost to CI, and none detected up to 200 kb on
the widest-vector x86-64 available**, which is the case that was supposed to pay for AVX-512.

What remains unmeasured is a genuinely genome-scale scan — hg38 against a full vertebrate release.
The finding above is evidence that the ratio does not move with scale rather than proof it never
will, and a production question deserves its own dated measurement rather than an extrapolation
from this one. But nothing here suggests host codegen is worth the reintroduced cache bug.

The pin is set on the CI job and deliberately not in the pixi task: a developer's cache is warm
from ordinary use and their CPU does not change, so there is nothing to fix locally.

### The other candidate, and why it is not this one

The index also carries a source stamp, `(st_mtime, st_size)` of the `.py` numba compiled, and the
pixi env is restored from a tarball every run — so mtime drift would invalidate the cache too,
independently of the CPU, with the same silent symptom. It is not ruled out, and under
`NUMBA_CACHE_DIR` it is a sharper risk than it looks: the locator in play is
`UserProvidedCacheLocator`, which stamps at full float precision — the locator that floors the
mtime is a different one, for in-tree caches — so sub-second drift counts.

The two are distinguishable without instrumenting anything, because they lose different amounts.
A stamp mismatch makes `_load_index` return `{}` and discards the **whole** index at once; a CPU
mismatch loses only the overloads compiled elsewhere. The bad run recompiled four compare tests
while the rest of the numba-touching suite stayed fast, which is the CPU's signature and not the
stamp's.

That is an inference from the shape of the loss, not a measurement, and it is worth being clear
that the stamp is **argued down rather than ruled out**. What settles it is two runs, because a
stamp is only unstable relative to another restore of the same env: `scripts/show_numba_cache_key.py`
prints the stamp on every run, so comparing the line across two runs on one lockfile answers it
directly — identical stamps mean the tarball preserves mtimes and this candidate is dead; drifting
stamps mean it is real and the CPU pin alone was never going to be enough. Until a second run
exists to compare against, the honest status is unresolved, and the instrumentation is there so
that resolving it costs nothing but reading.

## 3. The suite, cut to a third

| | before | after |
| --- | --- | --- |
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
| --- | --- | --- |
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
| --- | --- | --- | --- |
| serial | 56.3 s | **24.0 s** | < 60 s |
| parallel (`pixi run test`) | 17–19 s | **8.4 s** | < 15 s |
| gate (`pixi run check`) | 24 s | **10 s** | — |
| tests | 2425 | **933** | < 1000 |

`--maxprocesses 10` rather than a bare `-n 10`: on a 4-core CI runner `auto` finds four and the cap
does not bind, so the number is a ceiling for a big box, never a floor for a small one.

## What the docstring examples cost, and the deadline they surfaced

Measured 2026-08-31, same box and same method, when `--doctest-modules` was turned on over
`src/genome` and 468 doctest items joined the unit lane.

**The lane's wall did not move.** Eight runs each way on one tree, `pixi run test` reporting:

| | items collected | wall |
| --- | --- | --- |
| without the examples | 1090 | 8.76 – 9.33 s |
| with them | 1560 | 8.75 – 9.34 s |

Which is what the bimodal-`load` finding above predicts: the wall of a parallel lane is set by its
longest indivisible item, and 1,283 examples spread over ten workers add none. Serially they are
2.46 s. So the "every example on the critical path of CI" concern prices out at nothing — but
*total work at any instant* went up by 43%, and that is not free.

**What the added contention surfaced is hypothesis's deadline** — the per-test wall-clock budget
this document already declined to assert, arriving by default at 200 ms rather than by choice.
`test_a_lifted_hit_maps_back_to_where_the_scan_found_it` took 300.30 ms on one run and 11.79 ms on
the immediate re-run hypothesis does to check, which it then reports as an *unreliable* test rather
than a slow one:

| tree | runs | occurrences |
| --- | --- | --- |
| examples collected | 16 | 3 |
| examples not collected | 8 | 0 |
| examples collected, deadline off | 12 | 0 |

Same tree, same commit, the collection the only difference — so the items are what did it, and the
test is not what changed. This is the bimodal-`load` finding again in a different costume: a
wall-clock number taken on a ten-worker box where several modules spawn process pools of their own
reports the scheduler, not the code. Turned off for the suite, in one place, as a registered
hypothesis profile — which test is unlucky is decided by scheduling and not by anything about the
test, so per-test suppression would be a list that grows by luck. `--durations=10` in every CI log
still shows what a genuinely slow test costs, which is the same reason a hand-maintained `slow`
marker was declined.

## What was not adopted

- **A `slow` marker.** A hand-maintained list keyed on a number that changes whenever you optimise.
  The `--durations=10` already printed in every CI log makes the floor visible without a list that
  can go stale.
- **A per-test wall-clock budget asserted in the suite.** The same unchanged test measured 8.4 s
  alone and 12.8 s under contention here, so the assertion would be keyed on hardware and load
  rather than on the test.
