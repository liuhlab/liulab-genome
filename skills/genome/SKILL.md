---
name: genome
description: |
  Use when handling genomic data — DNA/RNA/protein sequences, reading or writing
  BAM/VCF/BED/GTF, or invoking samtools/bedtools. The `liulab-genome` package
  (import name `genome`) provides typed sequence classes (`DNA`, `RNA`, `Protein`)
  with biological transforms (complement, reverse_complement, transcribe,
  back_transcribe, gc_content) and a `genome` CLI (`revcomp`, `doctor`, `version`).
  TRIGGER when: user mentions reverse complement, transcription, GC content,
  validating a DNA/RNA/protein string, computing alphabet membership; user imports
  `genome` or runs `genome <subcommand>`; user asks about samtools/bedtools
  versions or "is the tooling installed"; user is working in a project whose
  pyproject lists `liulab-genome`.
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
| Reverse-complement a DNA/RNA string | ✅ `DNA(s).reverse_complement()` or `genome revcomp` | ❌ hand-written `str.translate` |
| Transcribe DNA → RNA | ✅ `DNA(s).transcribe()` | |
| Validate a sequence's alphabet | ✅ `DNA(s)` raises `ValueError` listing offending characters | |
| Compute GC fraction | ✅ `DNA(s).gc_content` (defined as 0.0 for empty) | |
| Verify samtools/bedtools are installed | ✅ `genome doctor` or `genome.external.doctor()` | |
| Parse a VCF / GTF / BED file | _(planned — not yet implemented)_ | use samtools/bedtools directly for now |
| Sequences containing `N` or other IUPAC codes | ❌ — out of scope by design | use Biopython or pyfaidx |
| Anything that needs whole-genome streaming | ❌ — these classes hold the full sequence in memory | use `samtools view`/`bedtools intersect` directly |

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
