---
name: genome
description: |
  Use when handling genomic data — reference genomes and assemblies, fetching the
  sequence of a region, DNA/RNA/protein sequences, reading or writing
  BAM/VCF/BED/GTF, or invoking samtools/bedtools. The `liulab-genome` package
  (import name `genome`) provides the `Genome` class (download an assembly by name
  and query its sequence in 0-based half-open coordinates), a `Region` coordinate
  primitive, typed sequence classes (`DNA`, `RNA`, `Protein`) with biological
  transforms (complement, reverse_complement, transcribe, back_transcribe,
  gc_content), and a `genome` CLI (`revcomp`, `doctor`, `version`).
  TRIGGER when: user wants the sequence of a locus / region (e.g. "chrIV:0-10"),
  loads a reference genome or assembly (hg38, mm39, sacCer3, ...), reads a `.2bit`,
  or needs chromosome sizes; user mentions reverse complement, transcription, GC
  content, or a DNA/RNA/protein string; user imports `genome` or runs
  `genome <subcommand>`; user asks about samtools/bedtools versions or "is the
  tooling installed"; user is working in a project whose pyproject lists
  `liulab-genome`.
  SKIP for: general string manipulation unrelated to biology; IUPAC ambiguity
  codes (intentionally out of scope); raw samtools/bedtools tasks where the user
  has explicitly chosen to shell out themselves.
---

# `liulab-genome` — agent guide

This skill teaches you how to use the `liulab-genome` package (import name
`genome`) correctly. The package wraps mature native tools (`samtools`,
`bedtools`) and provides typed Python primitives for biological sequences.

## When to reach for this package

| Task | Use this package | Don't use this package |
|------|------------------|------------------------|
| Fetch the sequence of a region from a reference | ✅ `Genome("hg38").fetch_sequence("chr1:0-100")` | ❌ hand-rolled FASTA slicing |
| Load a reference assembly by name (download + prepare) | ✅ `Genome("sacCer3")` | ❌ manual `wget` + `samtools faidx` |
| Read sequence from a `.2bit` file | ✅ `TwoBit(path).sequence(...)` | ❌ parsing 2bit bytes by hand |
| Get chromosome sizes / names | ✅ `Genome(...).chrom_sizes` / `.chromosomes` | |
| Reverse-complement a DNA/RNA string | ✅ `DNA(s).reverse_complement()` or `genome revcomp` | ❌ hand-written `str.translate` |
| Transcribe DNA → RNA | ✅ `DNA(s).transcribe()` | |
| Validate a sequence's alphabet | ✅ `DNA(s)` raises `ValueError` listing offending characters | |
| Compute GC fraction | ✅ `DNA(s).gc_content` (defined as 0.0 for empty) | |
| Verify samtools/bedtools are installed | ✅ `genome doctor` or `genome.external.doctor()` | |
| Parse a VCF / GTF / BED file | _(planned — not yet implemented)_ | use samtools/bedtools directly for now |
| Sequences containing `N` or other IUPAC codes | ❌ — out of scope by design | use Biopython or pyfaidx |
| Anything that needs whole-genome streaming | ❌ — these classes hold the full sequence in memory | use `samtools view`/`bedtools intersect` directly |

## The `Genome` class — main entry point

`Genome` is what most users want. Name a UCSC assembly; on construction it
downloads the FASTA (if not cached) and prepares the `.fai`/`.2bit`/`chrom.sizes`,
then serves sequence queries.

```python
from genome import Genome

g = Genome("sacCer3")                 # download + prepare on first use; cached after
g.fetch_sequence("chrIV:0-10")        # DNA('ACACCACACC')
g["chrIV:0-10"]                       # indexing is sugar for fetch_sequence
g.fetch_sequence("chrM")              # bare chromosome name -> whole sequence
g.chromosomes                         # ['chrI', 'chrII', ..., 'chrM']
g.chrom_sizes["chrIV"]                # 1531933  (pandas Series, copy)
```

**Coordinates are 0-based, half-open** (`chrIV:0-10` = the first ten bases). This
matches the package's internal convention — there is **no** 1-based conversion.
Do not "fix" a user's locus by adding 1; if they hand you a 1-based coordinate
(from IGV, VCF, GTF, a paper), subtract 1 from the start yourself before calling.

Key behaviors to rely on (and to tell the user about):

- The result is a [`DNA`](#biological-transforms), so transforms chain directly:
  `g.fetch_sequence("chrIV:0-100").gc_content`.
- **Soft-masking is preserved** (lower-case repeat-masked bases). The reference
  may contain `N` runs — `DNA` stores them verbatim (it does not validate).
- A `Region` with strand `"-"` is returned **reverse-complemented**; a bare
  string is always forward-strand.
- **Out-of-range coordinates raise `ValueError`** — an over-long `end` is never
  silently clamped. `end == chromosome length` is valid; `end >` it is an error.
- `Genome` holds the `.2bit` open; use it as a context manager
  (`with Genome("hg38") as g: ...`) or call `g.close()` to release the handle.

### Coordinates: `Region` and `parse_region`

`genome.region` is the shared coordinate layer (the primitive future features
build on):

```python
from genome import Region
from genome.region import parse_region

Region("chr1", 0, 10)                 # frozen, 0-based half-open, strand '.'
Region("chr1", 0, 10, "-")            # strand is explicit; never default to '+'
len(Region("chr1", 0, 10))            # 10
parse_region("chr1:1,000-2,000")      # ('chr1', 1000, 2000) — separators tolerated
parse_region("chrM")                  # ('chrM', None, None) — bare chromosome
```

### Reading a `.2bit` directly: `TwoBit`

When you already have a `.2bit` file (and don't need the download/prepare flow):

```python
from genome.io.twobit import TwoBit

with TwoBit("sacCer3.2bit") as tb:    # holds the handle open; reuse for many queries
    tb.chroms()                       # {'chrI': 230218, ...}
    tb.sequence("chrIV", 0, 10)       # 'ACACCACACC' (0-based, half-open, bounds-checked)
```

## CLI

The `genome` CLI is a thin Typer wrapper. Logic lives in the Python API, so
the same code paths run whether you call from Python or the shell.

### `genome revcomp <SEQUENCE> [--json]`

Reverse-complement a DNA string. Case is preserved (lowercase is meaningful —
it's the conventional soft-masking marker, not noise).

```bash
$ genome revcomp ATCG
CGAT

$ genome revcomp aTcG
CgAt                                    # case preserved: complement then reverse

$ genome revcomp ATCG --json
{"input": "ATCG", "reverse_complement": "CGAT"}

$ genome revcomp ATCN
error: DNA contains characters outside alphabet {ACGT}: ['N']
$ echo $?
2                                       # exit code 2 on invalid input
```

### `genome doctor [--json]`

Verify that `samtools` and `bedtools` are on `PATH` and report their versions.
Use this when troubleshooting a fresh checkout or before running anything that
shells out.

```bash
$ genome doctor
samtools: samtools 1.21 ...
bedtools: bedtools v2.31.1

$ genome doctor --json
{"samtools": "samtools 1.21 ...", "bedtools": "bedtools v2.31.1"}
```

Exit code `1` if any required tool is missing; the stderr message tells the
user how to fix it (`pixi add` / `pixi shell`).

### `genome version`

Prints the installed package version (derived from the git tag via `hatch-vcs`).

## Python API

Import the typed sequence classes from the top level:

```python
from genome import DNA, RNA, Protein
```

### Construction and validation

Validation runs in `__new__` (because `str` is immutable). It is
**case-insensitive** but **preserves case** in the stored value.

```python
DNA("ATCG")            # DNA('ATCG')
DNA("aTcG")            # DNA('aTcG') — case preserved
DNA("")                # DNA('') — empty is valid
DNA("ATCX")            # raises ValueError naming "X"
RNA("ATCG")            # raises ValueError — RNA alphabet is ACGU, no T
Protein("MKTAYN")      # raises ValueError — N is not a standard amino acid
```

The 20-standard-amino-acid alphabet for `Protein` is `ACDEFGHIKLMNPQRSTVWY`.
**IUPAC ambiguity codes (`N`, `R`, `Y`, `B`, `O`, `U`, `Z`, …) are intentionally
out of scope.** If a user asks for them, do not silently extend the alphabet —
flag it as a deliberate design decision and ask whether to scope new work to
add them.

### Biological transforms

```python
DNA("ATCG").complement()              # DNA('TAGC')   A↔T, C↔G
DNA("ATCG").reverse_complement()      # DNA('CGAT')   complement then reverse
DNA("aTcG").reverse_complement()      # DNA('CgAt')   case preserved
DNA("ATCG").transcribe()              # RNA('AUCG')   T→U, t→u
RNA("AUCG").back_transcribe()         # DNA('ATCG')   U→T, u→t
DNA("GGCC").gc_content                # 1.0 (float in [0.0, 1.0])
DNA("").gc_content                    # 0.0 (not a ZeroDivisionError)
```

### Typed slicing

`__getitem__` is the only inherited `str` method that's overridden — slicing
and indexing return the same subclass:

```python
DNA("ATCG")[1:3]       # DNA('TC')
DNA("ATCG")[0]         # DNA('A')
```

### What does NOT stay typed (gotcha)

All other inherited `str` methods (`upper`, `lower`, `replace`, `strip`,
`split`, `+`, …) return **plain `str`** by design. If a user needs a typed
variant, wrap explicitly — do not silently override these methods in user
code, that would surprise other readers.

```python
type(DNA("aTcG").upper())              # <class 'str'>  (not DNA!)
DNA(DNA("aTcG").upper())               # DNA('ATCG')    — explicit re-wrap
```

### Tool discovery from Python

```python
from genome.external import doctor, tool_version, ToolNotFoundError

doctor()                       # {'samtools': '...', 'bedtools': '...'}
tool_version("samtools")       # first line of `samtools --version`

try:
    tool_version("does-not-exist")
except ToolNotFoundError as err:
    # err message is actionable: names the tool and explains how to install it
    ...
```

## Domain invariants — non-negotiable

When working in this codebase or building features on top, respect these
invariants (they are landmines if violated):

1. **Coordinate systems**: store all intervals internally as **0-based,
   half-open `[start, end)`**. Convert to/from 1-based-inclusive (VCF,
   GFF/GTF, SAM) **only at I/O boundaries**, and document that conversion in
   the function's docstring. Never assume.
2. **Reference assembly is never implicit**. Carry the assembly (GRCh37 vs
   GRCh38, mm10 vs mm39, etc.) in metadata. Mixing assemblies is an error
   that should raise, not a warning.
3. **Chromosome naming**: normalize `chr1` vs `1` at ingest. Document which
   form is canonical internally before relying on equality.
4. **Strand** is `+` / `-` / `.` (unknown). Never default silently to `+`.
5. **Large files stream**. Do not read a whole genomic file into memory.
   Assume bgzipped + indexed input. The sequence classes are for in-memory
   strings (e.g. a single read, an exon, a peak summit ±100 bp), not whole
   chromosomes.
6. **Metadata is first-class**. Sample IDs, assembly, provenance travel
   *with* the data. Dropping metadata across a transform is a bug.

These come from `CLAUDE.md §5`; check the project's `CLAUDE.md` for any
updates before adding new features.

## Adding to the package

If the user asks you to extend `liulab-genome`, follow `CLAUDE.md` strictly:

- src layout; package code lives under `src/genome/`, tests under `tests/`.
- **Full type annotations** on every public function.
- **NumPy-style docstrings** with Parameters / Returns / Raises / Examples on
  every public surface.
- **Tests in the same change** as the feature — pytest + hypothesis for
  parsers, coordinate conversions, and any transform that has an obvious
  invariant (involution, round-trip, length preservation).
- Use `pixi run check` before claiming done — it must be green.
- Never reimplement what `samtools`/`bedtools` already do; shell out via
  `genome.external` and treat that module as the I/O boundary.
- Never `pip install` / `conda install` / `uv add` — only `pixi add` (conda
  dep) or `pixi add --pypi` (PyPI dep). Channels are `conda-forge` then
  `bioconda` in that priority order.
