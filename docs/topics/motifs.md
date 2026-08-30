# Motifs

A **motif** is the short stretch of DNA a transcription factor recognises, written down as a
count of how often each base appeared at each position across the sites the factor was
measured binding. JASPAR publishes them in sets. Name a set and it is prepared on this
machine, and what comes back scans sequence for you.

```python
from genome.tf.motif import JasparDatabase

jaspar = JasparDatabase()
len(jaspar)          # 1019
jaspar.release       # '2026'
```

Which motifs stand for which gene is on
[Transcription factors](transcription-factors.md).

## Getting a release

A **release** is one dated publication of JASPAR's collection. This package prepares `2024`
and `2026`, and uses `2026` when you name none. A **tax group** is the slice of the release
the file covers: `vertebrates`, `plants`, `insects`, `nematodes`, `fungi`, `urochordates`,
`diatoms`, or `all` for the union of the other seven. The default is `vertebrates`.

```python
len(JasparDatabase())                       # 1019
len(JasparDatabase("2024"))                 # 879
len(JasparDatabase("2026", "nematodes"))    # 103
```

The tax group picks which file is downloaded rather than filtering one afterwards, so a worm
scan never pays for a thousand plant matrices.

Constructing a `JasparDatabase` prepares it. The first construction of a given release and
tax group fetches the file, checks that it holds the number of motifs that release is known
to hold, and records that it finished. Every construction after that re-reads what is there
and needs no network.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

The first use downloads. The lab's CPU compute nodes have no internet, so run `genome motif scan` once on a login node before submitting a job that needs it.

Naming a release this package does not prepare raises `ValueError`, and the message lists
the ones it does.

## Looking at a motif

A motif is addressed by its id. Indexing tries the full id, then the id without its version,
then the motif name, and returns exactly one motif:

```python
ctcf = jaspar["MA0139.2"]
ctcf.motif_name       # 'CTCF'
len(ctcf)             # 15
ctcf.consensus        # DNA('GCCACCAGGGGGCGC')
```

The annotation JASPAR publishes beside the matrix comes with it. Four of the fields hold a
tuple, because a dimer has one value per half:

```python
ctcf.tf_class         # ('C2H2 zinc finger factors',)
ctcf.uniprot_ids      # ('P49711',)
ctcf.data_type        # 'ChIP-seq'
ctcf.counts.shape     # (4, 15)
```

`counts` is the matrix itself, four rows in `A`, `C`, `G`, `T` order and one column per
position. It is read-only, so nothing can change a motif's matrix while it keeps the same
id.

**A motif name is a label and several motifs can carry it**, which is why indexing by one
that is not unique refuses rather than picking:

```python
jaspar["CTCF"]
# AmbiguousMotifNameError: the motif name 'CTCF' labels 3 motifs here, so it addresses
# none of them: MA0139.2, MA1929.2, MA1930.2. Index by the motif id of the one you mean,
# or call set.by_name('CTCF') for all of them — a name is a label and only a motif id is
# a key.
```

`by_name` always answers with a tuple, whether the name labels one motif or four:

```python
[motif.motif_id for motif in jaspar.by_name("CTCF")]
# ['MA0139.2', 'MA1929.2', 'MA1930.2']
```

`filter` narrows the set by annotation or by a predicate. The four prose fields match on a
case-insensitive substring, so one spelling finds every spelling:

```python
zinc = jaspar.filter(tf_class="zinc finger")
len(zinc)     # 385
```

What comes back is a plain `MotifSet`, not a database, and it no longer claims to be the
release it was cut from. Hits scanned with it record `None` for release and tax group.

## Scanning regions of a genome

`scan_regions` takes a motif set and one or more
[`Region`](../genome/sequences.md)s, fetches each region's bases, and scans them in one pass.
These two regions gave 227 hits:

```python
from genome import Genome, Region

sacCer3 = Genome("sacCer3")
peaks = [Region("chrIV", 1000, 1500, "+"), Region("chrIV", 9000, 9500, "-")]
hits = sacCer3.scan_regions(jaspar, peaks)
hits.head(3)
#    motif_id motif_name sequence_name  start   end strand      score
# 0  MA0069.1       PAX6         chrIV   1188  1202      +  11.570312
# 1  MA0069.1       PAX6         chrIV   1186  1200      -  11.031250
# 2  MA0069.1       PAX6         chrIV   1447  1461      -  12.843750
```

**The hits come back in the assembly's own coordinates**, not in positions counted from the
start of each region. Both strands are scanned whatever the region's strand is. A `-` strand
region is fetched reverse-complemented, and its hits are flipped back into the forward frame
before you see them, so one interval scanned as `+` and as `-` gives the same rows.

Regions have to be `Region` objects. A locus string is refused, because it carries no strand,
and the strand is what decides whether the flip happens.

Passing `output=` here raises `TypeError`: a scan that streams to a file answers with a path,
and a path holds no coordinates to lift into the assembly's frame.

## The hit table

Every scan in this package answers with the same table, one row per hit.

| Column | What it holds |
|---|---|
| `motif_id`, `motif_name` | Which motif matched. |
| `sequence_name` | The chromosome or FASTA record the hit sits on. |
| `start`, `end` | The bases the matrix scored, 0-based and half-open. |
| `strand` | `+` or `-`, never `.`. A scan knows which of the two it scored. |
| `score` | The log-odds score in bits. Not a p-value. |

What the scan was is on `frame.attrs`, so a table always says where it came from:

```python
hits.attrs["threshold"]              # 0.0001
hits.attrs["background"]             # (0.25, 0.25, 0.25, 0.25)
hits.attrs["release"]                # '2026'
hits.attrs["tax_group"]              # 'vertebrates'
hits.attrs["assembly"]               # 'sacCer3'
len(hits.attrs["motifs_scanned"])    # 904
len(hits.attrs["motifs_skipped"])    # 115
```

`assembly` is there only on a table from `scan_regions`. The other scans answer in the frame
of whatever they were handed and name no assembly.

## Scanning a whole FASTA

`scan_fasta` reads one record at a time, so the file is never held whole. Given `output=` it
streams the hits to Parquet and answers with the path it wrote:

```python
written = jaspar.scan_fasta(
    sacCer3.files.fasta, output="/tmp/sacCer3_hits.parquet", workers=None
)
# PosixPath('/tmp/sacCer3_hits.parquet')
```

A record is named by its header up to the first whitespace, which is what STAR and chromap
write into an alignment made from the same file, so the two join with nobody renaming
anything.

`workers` shards the scan across processes and defaults to one; `None` asks for every core
the allocation granted. More than one produces the identical table, row for row. **A worker
re-imports the script that started it**, so a script scanning with more than one worker needs
its top-level code under `if __name__ == "__main__":`.

The same scan runs from a shell, which defaults to every core instead of one. The skipped
list is trimmed here:

```console
$ genome motif scan /Users/hanqing/liulab_data/genome/sacCer3/sacCer3.fa hits.parquet
scanned 17 sequences with 904 motifs from JASPAR 2026 vertebrates
  background  0.303, 0.194, 0.199, 0.304
  threshold   0.0001
  skipped     115 under 7 positions, so not scanned: MA0004.1, MA0130.1, MA0151.1, ...
  workers     14
  hits        2302168 -> hits.parquet
```

The hits go to the file and the summary to standard output, so `--json` is never mixed with
table data.

## Reading hits back

`read_hits` restores the column dtypes and the provenance both. `pandas.read_parquet` on its
own gives the rows and drops what the scan was.

```python
from genome.tf.motif import hit_count, read_hits

hits = read_hits("/tmp/sacCer3_hits.parquet")
hits.attrs["background"]    # (0.303, 0.194, 0.199, 0.304)
```

`hit_count` counts the rows off the file's footer without reading one of them:

```python
hit_count("/tmp/sacCer3_hits.parquet")    # 2302168
```

`provenance_of` does the same for `attrs` alone. Both cost the same on a 2.3-million-row file
as on an empty one. A path that does not exist raises `FileNotFoundError`.

## Thresholds and backgrounds

`threshold` is one per-position p-value, `1e-4` by default. It is converted per motif into the
score that motif has to clear, so one number means the same stringency for a short matrix and
a long one. A smaller number returns fewer hits:

```python
region = Region("chrIV", 0, 20000, "+")
len(sacCer3.scan_regions(jaspar, region))                    # 4122
len(sacCer3.scan_regions(jaspar, region, threshold=1e-5))    # 853
```

Motifs shorter than 7 positions cannot reach the default threshold at all. They are left out
of the scan and named in `motifs_skipped`, rather than being called at some looser cutoff
without saying so.

`background` is the base composition the scores are taken against, and it moves the answer
more than anything else on this page. Left alone it is chosen from the input: derived from it
when the input holds at least 10,000 unambiguous bases, uniform below that, since a
composition estimated from fewer would distort the cutoffs it sets. **Whichever it turned out
to be is recorded on the result**, so a table always says what its scores were measured
against:

```python
short = Region("chrIV", 1000, 1500, "+")
sacCer3.scan_regions(jaspar, short).attrs["background"]
# (0.25, 0.25, 0.25, 0.25)
sacCer3.scan_regions(jaspar, region).attrs["background"]
# (0.307, 0.203, 0.181, 0.309)
```

`"uniform"` pins it, `"derive"` derives whatever the input holds however little that is, and
four frequencies of your own are accepted in place of a mode. Comparing an unnamed matrix
against a release with `MotifComparison`, reading a transfac file yourself with
`parse_transfac`, and the sharding knobs behind `workers` are in the
[API reference](../reference.md).

## Where the motifs come from

The matrices are JASPAR's. Every database says which release it is and which file it read:

```python
jaspar.source_url
# 'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt'
jaspar.path
# PosixPath('/Users/hanqing/liulab_data/motif/jaspar/2026/vertebrates/JASPAR2026_CORE_vertebrates_non-redundant_pfms_transfac.txt')
```

Cite the JASPAR release you scanned with. The release and tax group travel on every hit table
and into every Parquet file, so a result opened months later still names what to cite.
