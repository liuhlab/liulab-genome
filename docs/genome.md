# Genome

`Genome` is the entry point. Name an assembly, and every reference file for it is on
disk by the time the constructor returns — then you query sequence, register
annotations, and build aligner indexes off the same object.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')
```

## Preparing an assembly

The first construction fetches the FASTA, checks it against the checksum the metadata
table pins, derives the `.fai`, `.2bit` and `chrom.sizes`, and records that it finished.
Everything lands in one directory per assembly:

```
<LIULAB_DATA>/genome/sacCer3/
├── sacCer3.fa
├── sacCer3.fa.fai
├── sacCer3.2bit
├── sacCer3.chrom.sizes
└── .completion.json
```

Later constructions read that record and open the `.2bit` — instant, and offline.
Nothing is downloaded twice.

`LIULAB_DATA` sets the data root. Unset, the well-known lab paths are tried in order and
`~/liulab_data` is the fallback. Pass `cache_dir=` to place a single assembly elsewhere:

```python
Genome("hg38")
Genome("mm39", cache_dir="/data/ref")
```

Preparing a large assembly takes a while, so it is often worth doing once from a shell
before a pipeline runs: `genome register hg38`. See the [CLI](cli.md).

### Using your own FASTA

When the download source is unreachable, or the reference is one no public site carries,
pass `path_or_url=` — a local file or a URL. A `.gz` source is decompressed for you:

```python
Genome("ce11", path_or_url="/data/ce11.fa.gz")
Genome("ce11", path_or_url="https://hgdownload-euro.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz")
```

Everything else is identical: the same derived files, the same directory, the same
record. Later plain `Genome("ce11")` calls reuse them. Nothing is downloaded, and the
assembly name only labels the directory.

## Fetching sequence

`fetch_sequence` takes a locus string, a bare chromosome name, or a
[`Region`](#regions), and returns a [`DNA`](sequences.md):

```python
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')
sacCer3["chrIV:0-10"]                    # same thing — indexing is sugar
sacCer3.fetch_sequence("chrM")           # a bare name is the whole chromosome
```

Because the result is a `DNA`, the transforms are right there:

```python
sacCer3.fetch_sequence("chrIV:0-10").reverse_complement()
sacCer3.fetch_sequence("chrIV:0-1000").gc_content
```

Lower-case bases (repeat soft-masking) come back verbatim — they carry meaning and are
never upper-cased for you.

A bare string is read on the forward strand. Pass a `Region` with strand `"-"` to get the
reverse complement of the interval:

```python
from genome import Region

sacCer3.fetch_sequence(Region("chrIV", 0, 10, "-"))
```

An `end` past the chromosome length raises rather than truncating silently; `end` equal
to the length is valid. Unknown chromosomes and malformed loci raise `ValueError` too.

```python
sacCer3.fetch_sequence("chrIV:0-999999999")
# ValueError: region chrIV:0-999999999: end (999999999) exceeds chrIV length (1531933).
```

## Inspecting an assembly

```python
sacCer3.assembly                 # 'sacCer3'
sacCer3.chromosomes              # ['chrI', 'chrII', ..., 'chrM'] in reference order
sacCer3.chrom_sizes              # pandas Series of lengths, indexed by chromosome
sacCer3.chrom_sizes["chrIV"]     # 1531933
sacCer3.files                    # GenomeFiles: fasta / fai / twobit / chrom_sizes paths
sacCer3.metadata.species         # 'Saccharomyces cerevisiae'
sacCer3.metadata.ncbi_name       # 'R64-1-1'
sacCer3.metadata.source_url      # where the FASTA came from
sacCer3.metadata.sha256          # digest of the unpacked FASTA (None if unpinned)
```

`chrom_sizes` is returned as a copy, so mutating it never corrupts the genome's own view.
`twobit_path`, `fasta_path` and `chrom_sizes_path` are there to hand to other tools.

`metadata` is always a record, so you read a field off it without checking first. For an
assembly the curated table does not list, every identifier is `None` and `assembly_name`
is the name you opened — the table is a cross-reference, not an allow-list.

## Releasing the file handle

`Genome` keeps the `.2bit` open so repeated queries are fast. Use it as a context manager
(or call `close()`) when you want the handle released at a known point:

```python
with Genome("sacCer3") as sacCer3:
    seq = sacCer3.fetch_sequence("chrIV:0-100")
```

## Gene annotations

A genome can carry several annotations, each under a short name. Name one the annotation
table lists for this assembly and it is fetched, verified, and built into a
[gffutils](https://gffutils.readthedocs.io/) database beside the assembly's files.
Everything about the *collection* of them is `genome.annotations`, one registry object
rather than a handful of methods on `Genome`:

```python
sacCer3.annotations.register("ensgene_v101")   # fetch + verify + build + record
sacCer3.annotations.registered                 # ['ensgene_v101'] — registered here
sacCer3.annotations.offered                    # what the table offers for this assembly
sacCer3.annotations.broken                     # here, and not to be trusted
sacCer3.annotations.path("ensgene_v101")       # Path to the placed .gtf
sacCer3.default_gtf                            # the annotation used when you name none
```

`registered` and `offered` answer different questions: what this machine has, and what
the lab supports. `genome annotations <assembly>` prints one against the other.

The registry is deliberately not a list. It settles a four-way state — registered,
broken, offered, nothing — so there is no `len()`, no `in` and no iterating over it:
name the set you mean. Note that `if sacCer3.annotations:` is therefore always true and
asks nothing; `if sacCer3.annotations.registered:` is the question.

For a GTF the table does not list, hand over the path — a `.gz` is decompressed:

```python
sacCer3.annotations.register_path("custom.gtf", "custom")
```

Each annotation gets its own directory, `<assembly dir>/gtf/<name>/`, holding the GTF,
the `.db`, and the record. Re-registering a name that is already registered does
nothing — no download, no rebuild.

Two things to know when a registration is refused:

- **Chromosome names must match.** Every sequence the GTF names has to be one the
  assembly carries, so an Ensembl-spelled GTF (`1`, `MT`) against a UCSC-spelled assembly
  (`chr1`, `chrM`) raises `ChromosomeMismatchError` instead of building an annotation
  that lines up with nothing. The check is one streaming pass, before the slow database
  build. If you have looked at the mismatch and accept it, pass
  `check_chromosomes=False`. The reverse is fine: an assembly may carry scaffolds the
  annotation never mentions.
- **Gene and transcript inference is off.** GENCODE, Ensembl and RefSeq GTFs declare
  those features already, and inferring them is gffutils' slow path. Turn it on
  (`disable_infer_genes=False`) only for a bare exon-level GTF, which otherwise registers
  as a database of exons and nothing else.

### The default annotation

`default_gtf` is what a caller gets when they name none: an explicit `default_gtf=` at
construction, else the one the table flags for this assembly, else the sole registered
annotation, else `None`.

```python
sacCer3.default_gtf                        # 'ensgene_v101'
Genome("hg38", default_gtf="refseq_2023")  # …unless you say otherwise
```

It names an annotation without locating one — opening a genome never starts a download —
so on a fresh machine the default is typically not registered yet. Asking for
`default_gtf_path` is what surfaces that, naming the command that closes the gap.

## Aligner indexes

Two aligners ship, and they differ in whether an annotation is involved:

| Aligner | Maps | Annotation | Index |
|---------|------|------------|-------|
| [STAR](https://github.com/alexdobin/STAR) | RNA-seq (splice-aware) | **one** — the default unless you name another | `index/star_<gtf>/` — a directory |
| [chromap](https://github.com/haowenz/chromap) | ATAC-seq, ChIP-seq, Hi-C | none | `index/chromap/chromap.index` — one file |

```python
sacCer3.build_star_index(threads=8)                    # against default_gtf
sacCer3.build_star_index("ensgene_v101", threads=8)    # ...or name one
sacCer3.build_chromap_index()
```

The `gtf` key becomes part of the index directory name, so indexes for different
annotations never collide. Omit it and the **default annotation** is used and names the
directory; a genome with no default at all raises rather than guessing, naming both ways
to give it one. A finished index is cached and reused; pass `overwrite=True` to rebuild.

Commonly tuned options are named and everything else is forwarded to the tool as a raw
flag:

```python
sacCer3.build_star_index(
    gtf="ensgene_v101",
    sjdb_overhang=99,        # --sjdbOverhang; ideally read_length - 1
    threads=8,               # --runThreadN
    genomeSAindexNbases=11,  # any other genomeGenerate flag, sans the leading --
)

sacCer3.build_chromap_index(
    kmer=20,                 # -k/--kmer
    window=10,               # -w/--window
    min_frag_length=30,      # any other --build-index flag; underscores become hyphens
)
```

`genomeSAindexNbases` and `genomeChrBinNbits` are sized from the assembly unless you set
them. chromap's index build is single-threaded, so there is no `threads`.

Building an index and using one are separate jobs. The `get_*` methods return the path a
mapping run needs and build nothing:

```python
sacCer3.get_star_index("ensgene_v101")   # -> STAR --genomeDir
sacCer3.get_chromap_index()              # -> chromap -x
```

Neither aligner is in the default environment — install what you need with `pixi add star`
/ `pixi add chromap`, or use the project's `aligners` environment
(`pixi run -e aligners ...`). A missing binary fails fast with the install command.

## Chimera assemblies

A **chimera** is one reference concatenated from assemblies you have already prepared — a
worm and the bacterium it eats — so a library carrying reads from both takes one alignment
pass instead of two.

```python
chimera = Genome.chimera(Genome("ce11"), Genome("ecHT115"))
chimera.assembly              # 'ce11_ecHT115'
```

The name is the component names sorted and joined by `_`, so you never choose it and
either order builds and reopens the same `ce11_ecHT115`. From a shell, naming it is
building it — `genome register ce11_ecHT115`. Nothing is downloaded: every component has
to be registered here already.

Every chromosome carries the component it came from, `<chromosome>__<component>`, and a
**bare name does not resolve**:

```python
chimera["I__ce11:0-10"]       # chromosome I of ce11, first ten bases
chimera["I:0-10"]
# ValueError: unknown chromosome 'I'; ce11_ecHT115 is a chimera, and every chromosome
# name in one carries the component it came from, so a bare name never resolves
# (ADR-0009). It carries 'I' as: I__ce11. Ask for the one you meant.
```

Two accessors read that split back:

```python
chimera.components                    # ['ce11', 'ecHT115'] — None for any other assembly
chimera.chrom_components["I__ce11"]   # 'ce11' — one entry per chromosome
```

The components' own default annotations are merged and registered by the same build, so a
chimera arrives annotated — and aligner indexes are built over it like any other assembly:

```python
chimera.default_gtf                   # 'wormbase_ws298+refseq_rs_2025_06_26'
chimera.build_star_index(threads=8)   # against that merged annotation
```

Its name is what went into it, so it changes when a component's default annotation does.
Rebuilding — `Genome.chimera(..., force=True)`, or `genome register <name> --force` —
registers the new one and **removes the one it replaces**, so a chimera never ends up
carrying two merged annotations with no default between them. An annotation you registered
by hand is never touched.

## Regions

`Region` is the shared coordinate type: frozen, validated, 0-based half-open, with an
explicit strand that is never defaulted to `+`.

```python
from genome import Region
from genome.region import parse_region

r = Region("chrIV", 0, 10)         # Region(chrom='chrIV', start=0, end=10, strand='.')
len(r)                             # 10
str(r)                             # 'chrIV:0-10'
Region("chrIV", 0, 10, "-")

parse_region("chrIV:1,000-2,000")  # ('chrIV', 1000, 2000) — separators tolerated
parse_region("chrM")               # ('chrM', None, None) — a bare chromosome name
```

Construction enforces `start >= 0`, `end >= start`, and `strand` in `{"+", "-", "."}`.

## When a directory cannot be trusted

Registration is finished when its record says so, never because the files look present.
A directory holding files with no record (a run that was interrupted) or a record that
disagrees with what is on disk (a file deleted or truncated afterwards) raises rather
than being rebuilt quietly:

```python
Genome("hg38")
# RegistrationMismatchError: /data/genome/hg38 disagrees with its .completion.json:
# hg38.2bit: recorded 841756144 bytes, found 0. Something changed these files after they
# were registered. Re-register it with `genome register hg38 --force`.
```

Every such message names its own repair, and the repair is always the same shape:
re-register with `--force` (or `overwrite=True` for an index). It keeps whatever is
provably good — an unpacked FASTA whose checksum still matches is reused, and only the
derived files are rebuilt.

An empty or absent directory is not this: that is a fresh registration and proceeds
normally. A broken *annotation* never blocks opening the genome; it is reported instead,
in `annotations.broken`, each entry carrying the command that repairs it.

To re-check integrity when nothing has raised but you suspect a problem, `genome verify
<assembly>` re-reads the FASTA and re-computes its digest.
