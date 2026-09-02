# Which motif scan engine, MOODS or memelite's fimo?

Measured 2026-08-28 on the maintainer's box (14 logical cores, 10 performance; macOS arm64),
re-run 2026-08-29 with identical numbers. **The decision was correctness, not speed. At the
default p = 1e-4, `fimo` returns zero hits for 97 of the 879 JASPAR 2024 CORE vertebrate
motifs and says nothing about it, and it never scans the last window of a sequence — a real
site at the final position is missed, and found again the moment one base is appended after
it. Both are properties of the code, not of the input. Speed then agreed: with thresholds
already prepared, MOODS scans 1 Mb against all 879 motifs on both strands in 0.90 s against
`fimo`'s 8.71 s on one core, 9.7×. Single-core MOODS also beats `fimo` given all fourteen
cores (1.51 s).** The two engines agree on only 72.6% of the union of their hit sets, so
they are not interchangeable back-ends.

This is a measurement, not a decision.

**Platform.** The timings are macOS arm64 and are indicative of **ordering only** — a re-run
on a Linux compute node is worth doing before any number here is quoted elsewhere. The
correctness findings are platform-independent: they follow from the arithmetic and the loop
bounds, and section 1 derives one of them without scanning anything at all.

## Method

MOODS 1.9.4.2 (bioconda, `py313` build), memelite 0.4.0 (PyPI), numba 0.67.0, numpy 2.5.0,
Python 3.13.14 — the versions this branch's lockfile pins.

Matrices: the JASPAR 2024 CORE vertebrates non-redundant transfac file, 879 motifs, parsed by
a throwaway parser written for this benchmark (the package's own is not written yet). Both
engines were given the same count matrices, converted to probabilities with a pseudocount of
0.01, a zero-order uniform background, both strands, and a per-position threshold of p = 1e-4.
`fimo` was left to convert those probabilities its own way, with its own `eps` of 1e-4 and its
own bin size of 0.1, because that is what a caller gets.

Targets, both already on this machine and neither downloaded for this: `tests/data/tiny.fa`,
30 kb of real sacCer3 in three records; and 1,000,000 bp of real `ce11` chromosome I
(1,000,001–2,000,000), upper-cased, as the target large enough for a zero hit count to mean
something. Section 3 also uses 10 Mb of the same chromosome as ten records.

Wall clock by `time.perf_counter` inside the process, three repetitions at 1 Mb and five at
30 kb, the best kept. `fimo` is numba-parallel, so it was pinned with `numba.set_num_threads(1)`
for the per-core numbers and its JIT warmed on a throwaway 400 bp call first, so no timing here
contains compilation. Memory is peak RSS from `getrusage`, sampled before and after the scan,
one engine and one threshold per process.

## 1. The 97 motifs that can never produce a hit

`fimo` builds a score-to-p-value table per motif by dynamic programming, then picks the
threshold as the first bin whose p-value falls below the requested one. The table is allocated
longer than the motif can actually score — `_pwm_to_mapping` pads its upper bound by the motif
length — and the padding is never written, so it reads as p = 0. When a motif cannot reach the
requested p-value at all, the search for the first bin below the threshold runs off the end of
the reachable scores and settles in that padding. The threshold comes back above the motif's
best possible score, and `_fast_hits` keeps a window only when `score > thresh`. The motif is
silently unscannable.

This needs no sequence to measure. Replaying `fimo`'s own `_all_pwm_to_mapping` over the 879
matrices and comparing each threshold against that matrix's best attainable score:

| | Motifs |
| --- | --- |
| in the release | 879 |
| **thresholded above their own best possible score** | **97** |
| — of length 5 | 2 |
| — of length 6 | 92 |
| — of length 7 | 3 |
| shorter than 7 bp in the release | 98 |
| shorter than 7 bp yet thresholded as reachable | 4 |

Scanning confirms it exactly. On the 1 Mb target, 97 motifs return zero `fimo` hits, and MOODS
returns hits for every one of those 97 — and for all 879.

| Target | fimo hits | MOODS hits | Zero-hit motifs, fimo | Zero-hit motifs, MOODS | Zero in fimo, non-zero in MOODS |
| --- | --- | --- | --- | --- | --- |
| `tiny.fa`, 30 kb | 6,002 | 7,953 | 148 | 46 | 103 |
| `ce11` chr I, 1 Mb | 214,190 | 284,441 | **97** | **0** | **97** |

Two things the table says that the headline does not:

- **30 kb is too small to see this cleanly.** At 1 Mb the 97 stand alone; at 30 kb they are
  buried among motifs that simply have no site in 30 kb of yeast. The 103 there are the 97
  plus six long motifs that MOODS happened to find once and `fimo` did not. Use a megabase.
- **The cause is mostly, but not only, short motifs.** A 6-mer has 4,096 words, so its best
  possible p-value is 2.44e-4 and it cannot reach 1e-4 by any threshold — 94 of the 97 are
  that. The other three are 7 bp. And the binning cuts the other way too: four motifs under
  7 bp are rated *reachable* by `fimo`'s coarse table when in truth they are not, so the
  table is wrong in both directions near the tail. Neither engine is right here on its own:
  MOODS returned hits for all 97, which is what clamping to the best attainable cutoff looks
  like — calling those motifs at a laxer p than was asked for. That is why the floor of 7 bp
  is a rule of the package rather than a question delegated to either engine.

## 2. The last window of every sequence is never scanned

`_fast_hits` walks `range(end - start - n)` for a motif of length `n` over a record of length
`L`. There are `L - n + 1` windows; it visits `L - n`. The final one is dropped.

Planting the consensus of `MA0139.2` (CTCF, 15 bp) so that it occupies exactly the last window,
after 100 bp of `ACGT` filler:

| Sequence | Length | Site at | fimo starts | MOODS starts |
| --- | --- | --- | --- | --- |
| filler + consensus | 115 | 100 | `[98]` | `[98, 100]` |
| filler + consensus + `A` | 116 | 100 | `[98, 100]` | `[98, 100]` |
| filler + consensus + filler | 215 | 100 | `[98, 100]` | `[98, 100]` |

One appended base is the whole difference. The same shows up at scale without planting
anything: over 879 motifs on the 1 Mb target, the largest end coordinate `fimo` reports is
999,999, never 1,000,000, though the shortest matrix scanned is 5 bp.

For a package that will scan **Region**s and shards of long sequences, this is worse than one
lost site per record: every shard boundary would silently drop a window, and the loss would
depend on how the work was divided.

## 3. Memory, and whether the file has to be resident

`fimo`'s FASTA path concatenates every record into one array before scanning and materialises
every hit as DataFrames before returning any of them; there is no incremental entry point.
MOODS takes one sequence per `scan` call, so a caller can hand it a record — or a shard — and
drain the hits before the next.

Peak RSS attributable to the scan, 1 Mb, four thresholds, marginal cost per hit taken between
adjacent rows:

| p | fimo hits | fimo Δpeak | MOODS hits | MOODS Δpeak |
| --- | --- | --- | --- | --- |
| 1e-5 | 23,449 | 21 MB | 107,467 | 22 MB |
| 1e-4 | 214,190 | 86 MB | 284,441 | 70 MB |
| 5e-4 | 1,114,743 | 412 MB | 1,175,305 | 271 MB |
| 1e-3 | 2,187,306 | 817 MB | 2,295,444 | 541 MB |

Marginally, **`fimo` costs about 380–396 bytes per hit and MOODS about 243–253** when the
caller keeps every match object, which the package will not. The design note's figure of
roughly 270 bytes per hit is close to what MOODS costs, not what `fimo` costs; `fimo` is
around 1.5× worse than that. Consider the 270 superseded by the table above.

The bound that matters is not per-hit but whether it is per-file. Over 10 Mb as ten records:

| Driver | Hits | Δpeak |
| --- | --- | --- |
| `fimo`, one call | 2,025,451 | 736 MB |
| MOODS, record by record, hits drained | 2,830,844 | 111 MB |

MOODS' peak is set by the largest record, `fimo`'s by the whole file, and no argument to
`fimo` changes that. That is the finding against the rule that no whole genomic file may be
read into memory: MOODS can be driven to honour it and `fimo` cannot.

## 4. Speed, per core

879 motifs, both strands, p = 1e-4, one core. Setup is threshold preparation; scan is
everything after it.

| | Setup | Scan, 1 Mb | Scan, 30 kb |
| --- | --- | --- | --- |
| MOODS | 6.60 s | **0.90 s** | **0.03 s** |
| fimo | 0.26 s | 8.71 s | 0.68 s |

**MOODS scans 9.7× faster per core** — 1.11 Mb/s against 0.115 Mb/s. Asking `fimo` for
DataFrames rather than counts adds about half a second, making it 10.2×. The design note's
7.4–8.5× is the right ordering and an understatement of the margin on this box; record 9.7×.

The two setups are the mirror image, and it is the reason the design caches thresholds:

- **MOODS' one slow step is its exact per-matrix threshold**, 6.60 s for 1,758 matrices
  (879 doubled for the reverse complements). `fimo` pays 0.26 s for a binned approximation —
  the same approximation that produces section 1. Paying MOODS' setup cold on every call, the
  crossover is at about 0.8 Mb of sequence; cached, as designed, it never has to be paid twice
  for the same (matrices, background, p).
- **Threads do not rescue `fimo`.** It is numba-parallel across motifs, and all fourteen cores
  bring 1 Mb from 8.97 s to 1.51 s — 5.9× from 14 cores. Single-core MOODS, at 0.90 s, is
  still faster than `fimo` with the whole machine.

## 5. They are not two back-ends for one thing

On the 1 Mb target the two call 288,907 distinct sites between them and share 209,724 — **72.6%
of the union**. MOODS finds 74,717 the other misses (the 97 unscannable motifs are most of it),
`fimo` finds 4,466 MOODS does not. Their score columns are on different scales besides. Offering
both and letting a caller pick would mean two hit tables that cannot be reconciled, which is the
problem the motif subpackage exists to end.
