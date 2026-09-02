# Assembly

An **assembly** is one published version of one species' reference genome. `sacCer3` is a
build of the yeast genome, `hg38` a build of the human one. Name an assembly and its files
are prepared on this machine. What comes back is an object that reads them.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.chromosomes[:4]    # ['chrI', 'chrII', 'chrIII', 'chrIV']
len(sacCer3.chromosomes)   # 17
```

Reading bases out of a prepared assembly is on [Sequences and regions](sequences.md).

## Which assemblies you can name

The shipped metadata table lists eight assemblies. Each row pins the URL its FASTA comes
from and the sha256 that download is checked against. `genome assembly list` prints them
beside what this machine already holds, and downloads nothing to do it:

```console
$ genome assembly list
assemblies in /Users/hanqing/liulab_data/genome
  hg38          offered, not registered  Homo sapiens GRCh38
  hg19          offered, not registered  Homo sapiens GRCh37
  mm39          offered, not registered  Mus musculus GRCm39
  mm10          offered, not registered  Mus musculus GRCm38
  sacCer3       registered               Saccharomyces cerevisiae R64-1-1
  ce11          registered               Caenorhabditis elegans WBcel235
  ecHT115       offered, not registered  Escherichia coli HT115 ASM435494v1
  ce11_ecHT115  offered, not registered
registered here: 2 — prepare another with `genome assembly register <name>`, or re-check one with `genome assembly verify <name>`
an assembly the table does not list registers too, from the UCSC golden path — with no pinned checksum behind it
```

`registered` means a record on this machine vouches for the four files. A directory in the
tree with no record beside it reads as `here, not registered` — nothing vouches for what
is in it — and one prepared under a name no row lists reads as `registered, not offered`.
Whether a registration's files are still what the record claims is a separate question, and
[`genome assembly verify`](#preparing-an-assembly) is what answers it.

What the listing treats as an assembly is the layout's own rule: one directory per
assembly, named for it, directly under `<LIULAB_DATA>/genome/`. So a directory there counts
unless its name begins with a dot, and a file there counts for nothing. Whether a directory
that counts is registered is then its record's business and never its name's.

The same report is `assembly_status()` from Python, and `assembly_table()` is the table by
itself. `ce11_ecHT115` is a chimera, concatenated from two of the others rather than
downloaded. `lookup_assembly` reads one row by name and returns `None` for a name the table
does not carry:

```python
from genome.assembly import lookup_assembly

lookup_assembly("sacCer3").ncbi_name    # 'R64-1-1'
lookup_assembly("danRer11") is None     # True
```

**A name the table does not list still registers.** Any UCSC assembly name resolves to its
golden-path FASTA, so `genome assembly register danRer11` fetches
`https://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz` and prepares it
the same way. It does not get a pinned checksum, so nothing independent confirms what
arrived; `genome assembly table-row danRer11` computes the row to add if you want one. A
name UCSC does not have fails on the first request, before any bytes are downloaded.

## Preparing an assembly

The first `Genome("sacCer3")` downloads the FASTA, checks it against the checksum the
metadata table pins, derives three companion files from it, and records that all of it
finished. Later constructions read that record and open the `.2bit`, which takes well
under a second and needs no network.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

Four files, in one directory named for the assembly:

| File | What it holds |
| --- | --- |
| `sacCer3.fa` | The reference sequence as text: a header line naming each chromosome, followed by its bases. |
| `sacCer3.fa.fai` | Byte offsets into that text, so a tool can jump straight to a chromosome instead of reading from the top. |
| `sacCer3.2bit` | The same sequence packed two bits to a base. This is what a locus is read out of. |
| `sacCer3.chrom.sizes` | Two columns, chromosome name and length. |

Annotations and aligner indexes land in subdirectories of the same place.

The same work runs from a shell. Do it there before a pipeline starts:

```console
$ genome assembly register sacCer3
registered sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

Running it a second time prints the same thing and downloads nothing. The first run takes
as long as the download does. Yeast is 3.8 MB compressed and finishes in a couple of
seconds; human is 940 MB.

The first use downloads. The lab's CPU compute nodes have no internet, so run `genome assembly register hg38` once on a login node before submitting a job that needs it.

### Using your own FASTA

Pass `path_or_url=` to seed an assembly from a file you already have, or from a URL that
is not the one the table pins. A `.gz` source is decompressed on the way in.

```python
Genome("ce11", path_or_url="/data/ce11.fa.gz")
Genome("ce11", path_or_url="https://hgdownload-euro.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz")
```

The shell spelling is `--source`:

```bash
genome assembly register ce11 --source /data/ce11.fa.gz
```

Everything after that is identical: the same derived files, the same directory, the same
record. A later plain `Genome("ce11")` reuses them. No checksum is compared, because a
seeded FASTA is whatever you handed over, and the assembly name only labels the directory
it lands in.

## Inspecting an assembly

`chromosomes` lists the sequences in the order the FASTA carries them, not sorted.
`chrom_sizes` gives their lengths as a pandas Series indexed by chromosome:

```python
sacCer3.chrom_sizes.head(3)
# chrom
# chrI      230218
# chrII     813184
# chrIII    316620
# Name: length, dtype: int64
```

The Series is a copy, so mutating it never changes the genome's own view of itself.

`files` carries the four paths, under the names `fasta`, `fai`, `twobit` and
`chrom_sizes`, ready to hand to another tool:

```python
sacCer3.files.twobit
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/sacCer3.2bit')
```

### Metadata

`metadata` says which reference this is and how other databases name it. It is always a
record, so a field can be read off it without checking first.

```python
sacCer3.metadata.species      # 'Saccharomyces cerevisiae'
sacCer3.metadata.ncbi_name    # 'R64-1-1'
sacCer3.metadata.ncbi_taxid   # 559292
```

`source_url` says where the FASTA was fetched from and `sha256` is the digest it was
checked against. The rest of the record is in the [API reference](../reference.md).

For an assembly the curated table does not list, every identifier is `None` and
`assembly_name` is the name you opened. The table is a cross-reference, not a list of what
you are allowed to prepare.

## Chimeras

A **chimera** is one reference concatenated from assemblies that are already prepared
here, such as a worm and the bacterium it eats. A library carrying reads from both then
takes one alignment pass instead of two.

```python
chimera = Genome.chimera(Genome("ce11"), Genome("ecHT115"))
chimera.assembly       # 'ce11_ecHT115'
chimera.components     # ['ce11', 'ecHT115']
```

The name is the component names sorted and joined by `_`. You never choose it, and either
order builds and reopens the one `ce11_ecHT115`. From a shell, naming a chimera is
building it:

```bash
genome assembly register ce11_ecHT115
```

Nothing is downloaded either way. **Every component has to be registered on this machine
already**, and one that is not stops the build; the message names the command that
prepares it.

Every chromosome in a chimera carries the component it came from, and a bare name does not
resolve:

```python
chimera.chrom_components["I__ce11"]    # 'ce11'
```

Ask for `I__ce11`, not `I`. The error you get for the bare name lists the spellings this
chimera does carry.

The components' own default annotations are merged and registered by the same build, so a
chimera arrives annotated:

```python
chimera.default_gtf    # 'wormbase_ws298+refseq_rs_2025_06_26'
```

That name is what went into it, so it changes when a component's default annotation
changes. Rebuild with `Genome.chimera(..., force=True)` or `genome assembly register
ce11_ecHT115 --force`, which registers the new merged annotation and removes the one it
replaces. An annotation you registered by hand is never touched.

## When a directory cannot be trusted

A registration is finished when its record says so, not because the files look present.
Opening a genome whose files changed after they were registered raises:

```python
Genome("sacCer3")
# RegistrationMismatchError: /Users/hanqing/liulab_data_scratch/genome/sacCer3 disagrees
# with its .completion.json: sacCer3.2bit: recorded 3039745 bytes, found 0. Something
# changed these files after they were registered. Re-register it with `genome assembly
# register sacCer3 --force`.
```

Run the command the message names. It keeps whatever is provably good: a FASTA whose
checksum still matches is reused, and only the derived files are rebuilt.

A directory holding files but no record is a run that was interrupted, and raises
`UnfinishedRegistrationError`, which names the same repair. An empty or absent directory
is neither of those: that is a fresh registration and proceeds normally.

To branch on it rather than read it, both import from `genome.store`:

```python
from genome.store import RegistrationMismatchError, UnfinishedRegistrationError

try:
    yeast = Genome("sacCer3")
except RegistrationMismatchError as changed:
    print(changed)      # files changed under a good record; re-register with --force
except UnfinishedRegistrationError as partial:
    print(partial)      # a build that never finished; the same repair applies
```

Both messages name the exact command to run. To treat the two the same, catch their
parent instead — `from genome.store import RegistrationError` covers either.

To re-check integrity when nothing has raised but you suspect a problem, `genome assembly
verify sacCer3` re-reads the whole FASTA and recomputes its digest:

```console
$ genome assembly verify sacCer3
/Users/hanqing/liulab_data/genome/sacCer3/sacCer3.fa: sha256 6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3 matches the digest the metadata table pins for it
```

## Releasing the file handle

A `Genome` keeps its `.2bit` open so repeated queries are fast. Use it as a context
manager, or call `close()`, when the handle has to be released at a known point.

```python
with Genome("sacCer3") as sacCer3:
    sizes = sacCer3.chrom_sizes
```
