# The `Genome` class

`Genome` is the package's main entry point. You name an assembly, and it makes
sure every reference file exists locally — then you query it for sequence.
Everything downstream (download, indexing, 2bit conversion, the open file
handle) is handled behind the scenes.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")              # download + prepare on first use (cached after)
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')
```

## Coordinates are 0-based, half-open

**Every coordinate in this package is 0-based and half-open** — `[start, end)`,
the BED convention. `chrIV:0-10` is the first ten bases (positions 0–9);
`chrIV:10-20` is the next ten, with no overlap. This is the same convention used
internally throughout the package (see [`genome.region`](#regions) below); there
is no hidden 1-based conversion.

```python
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')  — bases 0..9
len(sacCer3.fetch_sequence("chrIV:0-10"))  # 10
```

## Constructing a genome

```python
Genome("sacCer3")                        # yeast — small, good for examples
Genome("hg38")                           # human
Genome("mm39", cache_dir="/data/ref")    # override where files are stored
```

On construction `Genome`:

1. validates the assembly name against UCSC (a typo fails fast),
2. downloads `<assembly>.fa.gz` from the UCSC golden path if not already cached,
3. prepares the `.fai` index, `.2bit` encoding, and `chrom.sizes`,
4. opens the `.2bit` for reading.

All of this is cached under `<LIULAB_DATA>/genome/<assembly>/` (configurable via
the `LIULAB_DATA` environment variable; default `~/liulab_data`), so the second
construction is cheap and works offline. The underlying machinery is documented
in [Downloading and preparing genomes](genome-files.md); `Genome` is the
high-level front door to it.

## Fetching sequence

`fetch_sequence` accepts a locus string, a bare chromosome name, or a
[`Region`](#regions), and returns a [`DNA`](sequences.md):

```python
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')
sacCer3["chrIV:0-10"]                     # same thing — indexing is sugar
sacCer3.fetch_sequence("chrM")           # bare name -> the whole chromosome
```

Because the result is a `DNA`, the sequence transforms are right there:

```python
sacCer3.fetch_sequence("chrIV:0-10").reverse_complement()
sacCer3.fetch_sequence("chrIV:0-1000").gc_content
```

### Soft-masking is preserved

Lower-case bases (repeat soft-masking) are kept verbatim — they carry meaning,
so they are not silently upper-cased:

```python
sacCer3.fetch_sequence("chrIV:0-20")     # e.g. DNA('ACACCACACCacacccacac')
```

### Strand

A bare string is always read on the forward strand. Pass a `Region` with strand
`"-"` to get the reverse complement of the interval:

```python
from genome import Region

sacCer3.fetch_sequence(Region("chrIV", 0, 10, "+"))   # forward
sacCer3.fetch_sequence(Region("chrIV", 0, 10, "-"))   # reverse complement
```

### Out-of-range coordinates raise

An `end` past the chromosome length is an error, not a silent truncation:

```python
sacCer3.fetch_sequence("chrIV:0-999999999")
# ValueError: region chrIV:0-999999999: end (999999999) exceeds chrIV length (1531933).
# Coordinates are 0-based half-open, so the maximum valid end is 1531933.
```

`end == length` is valid (it selects through the last base); `end > length`
raises. Unknown chromosomes and malformed loci raise `ValueError` too.

## Inspecting the assembly

```python
sacCer3.assembly                 # 'sacCer3'
sacCer3.chromosomes              # ['chrI', 'chrII', ..., 'chrM'] in reference order
sacCer3.chrom_sizes              # pandas Series: lengths indexed by chromosome name
sacCer3.chrom_sizes["chrIV"]     # 1531933
sacCer3.files                    # GenomeFiles: paths to fasta/.fai/.2bit/chrom.sizes
```

`chrom_sizes` is a pandas `Series` (integer lengths, indexed by chromosome name,
in reference order). It is returned as a copy, so mutating it never corrupts the
genome's own view.

## Gene annotations (GTF)

Beyond sequence, a `Genome` can carry one or more gene annotations. Register a
GTF under a name and it is placed alongside the assembly's files with a gffutils
database built from it:

```python
sacCer3.register_gtf("sacCer3.ensGene.gtf", name="ensembl")
sacCer3.annotations              # ['ensembl'] — registered names
sacCer3.get_gtf_path("ensembl")  # Path to the placed .gtf
sacCer3.default_gtf              # 'ensembl' when it is the only one
```

Annotations are the basis for building aligner indexes (a STAR index is built
against a specific GTF). Registration, the default-annotation rules, and
`Genome.build_star_index(gtf=...)` are covered in
[Annotations & aligner indexes](aligner.md).

## Releasing the file handle

`Genome` holds the `.2bit` file open so repeated queries are fast. Use it as a
context manager (or call `close()`) when you want the handle released
deterministically:

```python
with Genome("sacCer3") as sacCer3:
    seq = sacCer3.fetch_sequence("chrIV:0-100")
# handle closed here
```

## Regions

`genome.region.Region` is the shared coordinate primitive that later features
build on. It is a frozen, validated, 0-based half-open interval with an explicit
strand:

```python
from genome import Region
from genome.region import parse_region

r = Region("chrIV", 0, 10)       # Region(chrom='chrIV', start=0, end=10, strand='.')
len(r)                           # 10
str(r)                           # 'chrIV:0-10'
Region("chrIV", 0, 10, "-")      # strand is explicit; never defaulted to '+'

parse_region("chrIV:0-10")       # ('chrIV', 0, 10)
parse_region("chrIV:1,000-2,000")  # ('chrIV', 1000, 2000) — separators tolerated
parse_region("chrM")             # ('chrM', None, None) — a bare chromosome name
```

Construction enforces the invariants: `start >= 0`, `end >= start`, and `strand`
in `{"+", "-", "."}`.

## Reading sequence directly: `TwoBit`

`Genome` reads sequence through `genome.io.twobit.TwoBit`, a thin wrapper over
`py2bit` that holds an open 2bit handle. You can use it on its own against any
`.2bit` file:

```python
from genome.io.twobit import TwoBit

with TwoBit("sacCer3.2bit") as tb:
    tb.chroms()                  # {'chrI': 230218, ..., 'chrM': 85779}
    tb.sequence("chrIV", 0, 10)  # 'ACACCACACC'  (0-based, half-open)
```

Like `Genome`, `TwoBit.sequence` bounds-checks coordinates — an over-long `end`
raises a `ValueError` instead of being silently clamped (py2bit's default).
`masked=True` (the default) preserves soft-masking; pass `masked=False` to
upper-case everything.
