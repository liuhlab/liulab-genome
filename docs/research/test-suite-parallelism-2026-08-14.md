# How many pytest workers, and does the gate gain from running concurrently?

> **Superseded for the worker-count and distribution tables, 2026-08-30.** These numbers were
> taken at 517 tests; the suite reached 2425 and was then cut to 933. `--dist=load` and
> `--maxprocesses 8` were both correct here and are both wrong now — see
> [`test-suite-cost-and-parallelism-2026-08-30.md`](test-suite-cost-and-parallelism-2026-08-30.md).
> Section 2 (BLAS thread pinning) still stands; nothing has changed about it.

Measured 2026-08-14 on the maintainer's box (14 logical cores, 10 performance). **Eight workers is
the floor of the curve — 2.8 s against 5.9 s serial. `auto` uncapped is *worse* than eight (3.4 s),
because it takes every core and the workers then contend. Distributing by test (`--dist=load`)
beats distributing by file (`--dist=loadfile`) by ~0.4 s, the suite being lopsided across modules.
Pinning the BLAS thread pools, which the sibling repo needs, makes no measurable difference here
and was not adopted.** Separately, the four-step gate run concurrently is 4 s against 9.4 s
sequential.

This is a measurement, not a decision.

## Method

`pytest -q` under `pixi run -e test`, wall-clock around the whole `pixi run` (so environment
activation is inside every number, as it is for a person typing the command). Two runs per
configuration, the better kept; three for the thread-pinning pair. 517 tests, 4 skipped, at the
time of measurement.

## 1. Worker count

| Configuration | Wall |
| --- | --- |
| serial | 5.92 s |
| `-n 2 --dist=loadfile` | 4.25 s |
| `-n 4 --dist=loadfile` | 3.32 s |
| `-n 6 --dist=loadfile` | 3.13 s |
| `-n 8 --dist=loadfile` | 3.18 s |
| `-n 10 --dist=loadfile` | 3.38 s |
| **`-n 8 --dist=load`** | **2.77 s** |
| `-n 10 --dist=load` | 2.80 s |
| `-n auto --dist=load` | 3.36 s |

Two things the table says that the headline number does not:

- **The curve turns up after eight.** Ten workers is slower than eight and `auto` (fourteen) is
  slower than ten. Each worker is a process that imports pandas, gffutils and the package before it
  runs anything, so past the point where that startup stops being hidden by real work, more workers
  is just more startup.
- **`load` beats `loadfile`.** `loadfile` keeps a module on one worker, and this suite's modules are
  lopsided — `test_gtf`, `test_cli` and `test_download` are most of the work — so a file-grained
  split leaves workers idle waiting for the biggest module. No test here shares mutable state with
  another (every one takes `tmp_path`, and the module-level state that is monkeypatched is
  process-local), so there is nothing `loadfile` was protecting.

`auto --maxprocesses 8` is what the task table spells rather than a bare `-n 8`: it reproduces the
best number here, and on a 2-core CI runner `auto` finds two and the cap does not bind. A bare `8`
would oversubscribe that runner by four.

## 2. BLAS thread pools

| Configuration | Wall |
| --- | --- |
| `-n auto --maxprocesses 8`, no pinning | 2.83 s |
| `-n auto --maxprocesses 8`, `OMP`/`OPENBLAS`/`MKL_NUM_THREADS=1` | 2.82 s |

Nothing. The pathology being guarded against is real — eight worker processes each sizing a thread
pool to the whole box, then fighting for it — but this package's numpy and pandas use is a
`chrom_sizes` Series and nothing else, so no BLAS kernel is ever entered. Declared nowhere, on the
grounds that a knob with no measured effect is a thing to explain later. Revisit if array work
arrives.

## 3. The gate, sequential against concurrent

| Step | Wall, alone |
| --- | --- |
| `lint` | 0.73 s |
| `fmt-check` | 0.07 s |
| `typecheck` | 2.86 s |
| `test` (serial, as it was) | 5.6 s |
| **`check`, sequential `depends-on`** | **9.38 s** |
| **`check`, concurrent via `scripts/check.sh`** | **4 s** |

The static preamble was 3.7 s of the sequential 9.4, spent on tools that share no inputs with each
other or with pytest. Concurrently the gate is bounded by its slowest step rather than by their sum,
and with the suite itself now parallel the two wins compound: 9.4 s to 4 s, most of which is pyright
and pytest overlapping.

The concurrent gate is not a strict superset of CI and is not meant to be: `build` and the strict
docs build have no step here. Each answers a question a code change does not ask, and this gate's
justification is that it is cheap enough to run without thinking about it.
