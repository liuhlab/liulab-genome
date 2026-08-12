# Can a `.2bit` be had without `faToTwoBit`?

Measured 2026-08-11 for [#15](https://github.com/liuhlab/liulab-genome/issues/15). **A prebuilt
`.2bit` can be downloaded from UCSC for the major assemblies, but nothing in Python can *write*
one — every Python 2bit library in existence is read-only, and on bioconda the only FASTA→2bit
converter is `ucsc-fatotwobit` itself. Since this package must accept arbitrary user FASTA, the
real alternative is not a different writer but not needing a `.2bit` at all: `pyfaidx` is already
in the environment for free, agrees with `py2bit` byte-for-byte including soft-masking, and is
*faster* than the repo's current configuration on large queries.**

This is a measurement, not a decision. The decision is
[#16](https://github.com/liuhlab/liulab-genome/issues/16).

## The constraint that frames all of it

Stated by the maintainer while this ticket was being worked:

> Even [if] UCSC ship[s] pre-built 2bit file, we still need the ability from FAST\[A] to 2bit,
> because we have to support other source[s] of FASTA.

`Genome(assembly, path_or_url=...)` seeds from a FASTA UCSC has never seen. No prebuilt artifact
can exist for it. So question 1 below can only ever produce a *fast path* for UCSC assemblies; it
cannot retire `faToTwoBit`. Only question 3 can.

## 1. UCSC does publish a prebuilt `.2bit`

URL shape — same `bigZips` directory the FASTA already comes from:

```
https://hgdownload.soe.ucsc.edu/goldenPath/<db>/bigZips/<db>.2bit
```

| assembly | HTTP | size |
|---|---|---:|
| hg38 | 200 | 835,393,456 B (835 MB) |
| hg19 | 200 | 816,241,703 B (816 MB) |
| mm39 | 200 | 714,181,470 B (714 MB) |
| mm10 | 200 | 714,784,109 B (715 MB) |
| dm6 | 200 | 36,969,050 B (37 MB) |
| sacCer3 | 200 | 3,039,745 B (3.0 MB) |

`md5sum.txt` sits in the same directory and **covers the `.2bit`**, so `pooch`'s `known_hash`
can be populated from it exactly as it could for the FASTA:

```
dcc3ea27079aa6dc3f9deccd7275e0f8  hg38.2bit
1c9dcaddfa41027f17cd8f7a82c7293b  hg38.fa.gz
```

**Caveat — six assemblies were checked, not all ~250.** Every one of the major assemblies has it;
coverage for obscure ones is unverified. A download path would need a 404 fallback regardless.

## 2. No Python library writes 2bit. None.

Every Python 2bit package is a reader:

| package | version | write support |
|---|---|---|
| `py2bit` (lib2bit) | 1.0.1 | no — `py2bit.open` is the entire module surface |
| `twobitreader` | 4.0.0 | no — `TwoBitFile`'s docstring says verbatim *"(note: no writing support)"* |
| `bx-python` `bx.seq.twobit` | 0.14.0 | no — `TwoBitFile` / `TwoBitSequence`, reader methods only |
| `pytwobit` | 0.3.1 | no — *"A fast reader for local or remote UCSC twobit sequence files"* |
| `Bio.SeqIO.TwoBitIO` | current | no — `TwoBitIterator`, a parser; there is no writer class |

Searching bioconda for `twobit` returns the complete set of packaged tools. Exactly one converts
FASTA to 2bit:

```
bioconda/ucsc-fatotwobit    Convert DNA from fasta to 2bit format.     <- the only writer
bioconda/ucsc-twobitinfo    Get information about sequences in a .2bit file.
bioconda/ucsc-twobitmask    Apply masking to a .2bit file, creating a new .2bit file.
bioconda/ucsc-twobittofa    Convert all or part of .2bit file to fasta.
bioconda/ucsc-twobitdup     check to see if a twobit file has any identical sequences in it
bioconda/pytwobit           (reader)
bioconda/twobitreader       (reader)
```

`ucsc-twobitmask` rewrites an *existing* 2bit; it is not a converter. The only other writer
implementation found anywhere is [weng-lab/TwoBit](https://github.com/weng-lab/TwoBit), a C++
reader/writer that is **not** packaged on PyPI or on any conda channel.

So: if the package must emit a `.2bit` from arbitrary FASTA, `faToTwoBit` is the only option that
exists. (Writing one by hand is not absurd — the
[format spec](https://genome.ucsc.edu/goldenpath/help/twoBit.html) is short — but that is new
code to own and test, encoding 3.1 Gbp in Python.)

## 3. The FASTA-only path costs zero new dependencies

### `pyfaidx` is already here, for free

`pyfaidx` 0.9.0.4 is in **every solved environment** in `pixi.lock` — pulled in as a transitive
dependency of `gffutils`, which is a direct dependency:

```
gffutils -> ['pyfaidx>=0.5.5.2']
```

`pysam` is *not* needed for this. It would be a genuinely new dependency (~3.5–4 MB, all four
platforms), and it buys nothing `pyfaidx` does not already provide here.

### `chrom.sizes` needs neither 2bit nor `twoBitInfo`

Columns 1–2 of the `.fai`, read in **0.12 ms**, byte-identical to `py2bit.chroms()`:

```
from .fai: {'chr1': 40000000, 'chr2': 18000000, 'chrM': 16569}
py2bit    : {'chr1': 40000000, 'chr2': 18000000, 'chrM': 16569}
identical : True
```

This was already established while charting; it is confirmed here.

### Sequence: `pyfaidx` and `py2bit` agree byte-for-byte

Over 300 random 500 bp probes, `pyfaidx` and `py2bit(storeMasked=True)` returned **identical**
strings. 160 of the 300 contained soft-masked lowercase, so the mask round-trips through both.

Note the trap this exposed: `py2bit.open(path)` **upper-cases everything** unless `storeMasked`
is passed. `genome.io.twobit.TwoBit` already passes it (`masked: bool = True`), so the repo is
correct — but it means the honest comparison is against the masked configuration, which is the
slower one.

### Throughput

58 Mbp soft-masked FASTA with N runs, M-series Mac, warm cache:

| path | 200 bp × 20,000 | 100 kb × 500 |
|---|---:|---:|
| `py2bit` `storeMasked=True` — **the repo's actual setting** | 0.4 µs/q | 78.7 µs/q |
| `py2bit` `storeMasked=False` | 0.2 µs/q | 21.6 µs/q |
| `pyfaidx`, plain `.fa` + `.fai`, case preserved | 2.1 µs/q | **72.4 µs/q** |
| `pyfaidx`, bgzf `.fa.gz` + `.fai`/`.gzi` | 121.7 µs/q | 395.5 µs/q |
| `samtools faidx` **subprocess per query** | 3,407 µs/q | 3,697 µs/q |

Three things fall out:

1. **On short queries `py2bit` is ~5× faster, and it does not matter.** Both are microseconds.
   At 20,000 lookups the whole difference is 34 ms.
2. **On large queries `pyfaidx` is *faster* than the repo's configuration** (72.4 vs 78.7 µs).
   Applying the soft-mask is what costs: 2bit stores the mask as separate blocks it must re-apply
   per query, while FASTA carries case inline and free. The 2bit speed advantage is real only if
   you give up soft-masking, which this package deliberately does not.
3. **Shelling out per query is ruled out** — ~3.4 ms/query, ~1,700× slower than in-process
   `pyfaidx`. If the 2bit goes, sequence must be served in-process. `samtools faidx` as a
   subprocess is not a candidate for the read path (it remains right for building the index once).

### Disk

Same 58 Mbp genome:

| file | size | % of raw |
|---|---:|---:|
| `.fa` | 58.98 MB | 100% |
| `.fa.gz` (bgzip) | 16.67 MB | 28.3% |
| `.2bit` | 14.51 MB | 24.6% |
| `.fai` / `.gzi` / `.chrom.sizes` | < 10 KB each | — |

**2bit is barely smaller than a bgzipped FASTA.** The compaction argument for 2bit is weak. But
note the bind: bgzf is the compact FASTA and also the *slow* one (121.7 µs/q vs 7.8 µs plain), so
"drop the 2bit, keep only the bgzipped FASTA" trades away query speed. Plain `.fa` + `.fai` is the
fast FASTA configuration and it is the largest of the three.

Today the cache holds `.fa.gz` **and** `.fa` **and** `.fai` **and** `.2bit` **and**
`.chrom.sizes`. Dropping the 2bit saves ~25% of raw genome size per assembly — ~780 MB for hg38.

### Prep time

58 Mbp, and the ~53× extrapolation to hg38's 3.1 Gbp:

| step | measured | hg38 (extrapolated) |
|---|---:|---:|
| `samtools faidx` | 0.56 s | ~30 s |
| `faToTwoBit` | 1.99 s | ~105 s |
| `twoBitInfo` | 0.27 s | ~14 s |
| `bgzip -@4` | 1.80 s | ~95 s |

`faToTwoBit` adds roughly **two minutes of one-time prep per assembly** — against a multi-hundred-MB
download that dominates it. Prep time is not the argument against keeping it.

## 4. `ucsc-fatotwobit` costs very little

`ucsc-fatotwobit` and `ucsc-twobitinfo`, both at version 482, have **identical** platform coverage
and dependency closures:

| platform | package size | builds at v482 |
|---|---:|---:|
| linux-64 | 0.38 MB | 27 versions available |
| linux-aarch64 | 0.40 MB | 9 |
| osx-64 | 2.37 MB | 18 |
| **osx-arm64** | 2.64 MB | **1** |

**osx-arm64 is supported** — the platform question that prompted this. But with a single build at
a single version, it is newly added and has no older version to fall back to if a solve ever needs
one.

The package itself is trivially small; the weight is the closure it drags in:

```
linux-64 / linux-aarch64:  bzip2, libgcc, liblzma, libopenssl-static, libpng,
                           libuuid, libzlib, mysql-connector-c
osx-64:                    same, minus libgcc
osx-arm64:                 same, minus libgcc AND minus mysql-connector-c
```

A UCSC genome-browser utility pulling `mysql-connector-c` and `libpng` to convert a FASTA is
inelegant, but it is small and it resolves on every platform this project targets.

**One thing that constrains the option space:** both binaries carry the *same* closure, so dropping
only `twoBitInfo` — the replaceable one — **saves nothing on install**. The closure disappears only
if both go.

## What this leaves for [#16](https://github.com/liuhlab/liulab-genome/issues/16)

The constraint kills "prebuilt only". Three live options remain:

- **A. Status quo.** Keep `faToTwoBit`. Optionally drop `twoBitInfo` for `.fai` columns 1–2 —
  correct and free, but saves no install weight since the closure stays.
- **B. Keep 2bit, add a UCSC fast path.** Download `<db>.2bit` (md5-verified) for UCSC assemblies,
  fall back to `faToTwoBit` for user FASTA. Saves ~2 min/assembly of prep; keeps the dependency;
  adds a second code path and a 404 fallback to maintain.
- **C. Drop 2bit entirely.** Serve sequence from `.fa` + `.fai` via `pyfaidx` (already present),
  `chrom.sizes` from the `.fai`. Removes `ucsc-fatotwobit`, `ucsc-twobitinfo`, `py2bit`, their
  closure, and ~25% of per-assembly disk. Costs: `genome.io.twobit` and its tests are rewritten,
  short-query lookups get ~5× slower in relative terms (µs either way), and the package can no
  longer hand a `.2bit` to any external tool that wants one.

The last clause in C is the part this ticket cannot answer: whether anything downstream consumes
the `.2bit` as a *file* rather than as this package's internal sequence store.

## Method

Measurements are reproducible from `docs/research/`-adjacent throwaway scripts; the synthetic FASTA
was 3 sequences / 58 Mbp (`chr1` 40 Mbp, `chr2` 18 Mbp, `chrM` 16,569 bp) with telomeric and
assembly-gap N runs and ~50% of blocks soft-masked to lowercase, seeded for reproducibility. Timings
are warm-cache, single run, on an M-series Mac — good to a factor of ~2, which is all any of the
conclusions above rest on. Package facts come from the anaconda.org and PyPI JSON APIs and from
reading library sources directly; UCSC facts from HTTP `HEAD` against `hgdownload.soe.ucsc.edu`.
