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

1. looks the assembly up in the curated metadata table,
2. downloads `<assembly>.fa.gz` — from the URL that row pins, or, for an assembly
   the table does not list, from the UCSC golden path after validating the name
   against UCSC (a typo fails fast),
3. checks the unpacked FASTA against the row's checksum, when it pins one,
4. prepares the `.fai` index, `.2bit` encoding, and `chrom.sizes`,
5. writes the registration record that says all of that finished,
6. opens the `.2bit` for reading.

Everything lands under `<LIULAB_DATA>/genome/<assembly>/` (configurable via the
`LIULAB_DATA` environment variable; default `~/liulab_data`). A second construction
reads that assembly's [registration record](genome-files.md#the-registration-record),
confirms every file it claims is present and the right size — no file contents are read
— and opens the `.2bit`, so it is instant and works offline. Nothing is downloaded
twice, and the compressed download is deleted once the record is written. The underlying
machinery is documented in [Downloading and preparing genomes](genome-files.md);
`Genome` is the high-level front door to it.

### A registration that cannot be trusted stops you

If that directory holds files but no record (an interrupted preparation), or a record
that disagrees with what is on disk (a file deleted or truncated afterwards), construction
raises a `RegistrationError` naming the file and the command that fixes it — rather than
rebuilding quietly, or handing back a genome that answers queries from a partial file:

```python
Genome("hg38")
# RegistrationMismatchError: /data/genome/hg38 disagrees with its .completion.json:
# hg38.2bit: recorded 841756144 bytes, found 0. Something changed these files after they
# were registered. Re-register it with `genome register hg38 --force`.
```

An absent or empty directory is not this — that is a fresh registration and proceeds
normally. `genome verify <assembly>` re-reads and re-checksums on demand when you suspect
a problem but nothing has raised. See
[When a registration cannot be trusted](genome-files.md#when-a-registration-cannot-be-trusted)
for the repair and the one trade-off it carries.

### Where the bytes came from

An assembly the lab officially supports carries its source and, once someone has
pinned it, the sha256 of its **unpacked** FASTA. Both are readable off the genome:

```python
sacCer3.source_url   # 'https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz'
sacCer3.sha256       # '6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3'
```

A FASTA that does not match a pinned checksum raises `ChecksumMismatchError` naming
both values, rather than being prepared and quietly used. `None` on either property
means the table pins nothing for this assembly — which is legal, and takes nothing
away: preparation proceeds exactly as it does for an assembly with no row at all.
See [Pinned sources and checksums](genome-files.md#pinned-sources-and-checksums).

### Seeding from your own FASTA (offline, mirrors, custom references)

When UCSC is unreachable (a firewalled compute node, a proxy-only network) or
you have a custom reference that isn't on the golden path, pass `path_or_url=`
to seed the assembly from a FASTA you provide instead of downloading from UCSC:

```python
# a local file — copied into the assembly's directory, then prepared
Genome("ce11", path_or_url="/data/ce11.fa.gz")

# a non-UCSC URL (e.g. a UCSC mirror)
Genome("ce11", path_or_url="https://hgdownload-euro.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz")
```

In this mode **UCSC is never contacted** — there is no assembly-name validation,
so the name only labels the directory and the prepared files, and no pinned
source or checksum from the metadata table is consulted. A gzipped
(`.gz`) source is decompressed automatically. Everything else is identical to the
UCSC path: the `.fai`, `.2bit`, and `chrom.sizes` are prepared under
`<LIULAB_DATA>/genome/<assembly>/` and the same registration record is written — with
the path you gave as its source and the digest of what arrived — so later plain
`Genome("ce11")` calls reuse them. See
[Seeding from a local file or URL](genome-files.md#seeding-from-a-local-file-or-url)
for the underlying `fetch_genome_from`.

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

Beyond sequence, a `Genome` can carry one or more gene annotations. Name one the
annotation table lists for this assembly and it is fetched, checked against the
checksum that table pins, placed alongside the assembly's files and built into a
gffutils database:

```python
sacCer3.register_annotation("ensgene_v101")   # fetch + verify + build + record
sacCer3.annotations                  # ['ensgene_v101'] — registered names
sacCer3.get_gtf_path("ensgene_v101") # Path to the placed .gtf
sacCer3.default_gtf                  # 'ensgene_v101' when it is the only one
```

For a GTF the table does not list, hand over the path instead:

```python
sacCer3.register_gtf("custom.gtf", name="custom")
sacCer3.register_gtf("custom.gtf.gz", name="custom")  # .gz is decompressed for you
```

Either way the registration ends with the same record the assembly writes, and
that record is the only thing that ever says the annotation is finished — never
the database file's existence, which is equally true of a build killed half-way.
Re-registering a name whose record is valid returns it silently: nothing is
fetched and nothing is rebuilt. A directory holding files without a valid record
raises instead, naming `genome register-annotation <assembly> <name> --force`,
which is also what repairs it.

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
