# Downloading and preparing genomes

Two modules in `genome.io` cover getting a reference genome onto disk and turning it
into the index/companion files most tools expect:

- [`genome.io.download`](#downloading) — cache large files from the network with
  [pooch](https://www.fatiando.org/pooch/), with a UCSC-aware subclass.
- [`genome.io.fasta`](#preparing-a-fasta) — shell out to `samtools`, `faToTwoBit`, and
  `twoBitInfo` to index a FASTA and derive its `.2bit` and `chrom.sizes`.

Both live at the **I/O boundary**: they touch the network and the filesystem and invoke
native binaries (managed by pixi). They never reimplement what those binaries do.

!!! tip
    Most users don't call these directly — the [`Genome`](genome.md) class wraps the whole
    download-and-prepare flow and then lets you query sequence. Reach for the modules here
    when you need a specific file on disk, a custom cache location, or one preparation step
    on its own.

## Downloading

Every download in the package goes through one fetch step, `fetch_url`, which is the only
place pooch is called. `Downloader` binds that step to a cache directory: it downloads a
URL once and caches it, so asking for the same file again is served from disk.

```python
from genome.io.download import Downloader

dl = Downloader()                      # caches under pooch.os_cache("genome")
path = dl.fetch("https://example.org/annotation.bed.gz")
```

Pass `cache_dir=` to redirect the cache (e.g. shared lab scratch). Supply
`known_hash="md5:…"` to verify the **downloaded bytes**; when omitted, pooch logs the
computed hash so you can pin it next time.

!!! warning "`known_hash` hashes the archive, not the genome"
    `known_hash` is checked by pooch against exactly what it downloaded — for a genome
    that is the `.fa.gz`. Gzip bytes change whenever a file is recompressed, while the
    FASTA inside does not, so an assembly's *pinned* checksum is a different thing: it
    covers the **unpacked** `<assembly>.fa` and is checked after decompression. See
    [Pinned sources and checksums](#pinned-sources-and-checksums).

### Genome FASTA from UCSC

`UCSCGenomeDownloader` fetches an assembly's FASTA: from the URL its metadata row pins,
or — for an assembly no row lists — from a URL built from the golden-path layout. Either
way it downloads `<assembly>.fa.gz` and (by default) decompresses it to `<assembly>.fa`:

```python
from genome.io.download import UCSCGenomeDownloader

dl = UCSCGenomeDownloader("hg38")
dl.fasta_url
# 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

fasta = dl.fetch_fasta()               # -> Path to hg38.fa (cached, multi-GB)
```

Use `decompress=False` to keep the `.fa.gz`, or `progressbar=False` to silence the bar.
When the URL had to be derived, a mistyped assembly fails fast on a `HEAD` request with a
clear error naming the unknown name; when a row pins the URL there is nothing to guess,
so that check is skipped.

#### Pinned sources and checksums

The curated metadata table that ships in the package (`genome/data/assembly_metadata.tsv`)
carries two columns beyond an assembly's identifiers: `source_url`, the URL its FASTA is
fetched from, and `sha256`, the digest of the **unpacked** FASTA that URL yields. Both may
be blank — the table is a cross-reference, not an allow-list, and an assembly with no row
at all keeps working exactly as before.

```python
from genome.metadata import lookup_assembly

lookup_assembly("sacCer3").source_url
# 'https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz'
lookup_assembly("sacCer3").sha256
# '6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3'
```

`fetch_genome` checks the unpacked FASTA against that digest and raises
`ChecksumMismatchError` — naming the file, the expected value and the actual one — when
they disagree. A row with a blank `sha256` has nothing to compare against, so it prepares
normally; `verify_fasta()` returns the computed digest either way.

The digest deliberately covers the decompressed file rather than the download, so a FASTA
recompressed at another level, or copied from a mirror or by hand, still matches the
official row. Filling the column in is a copy-paste:

```bash
$ genome table-row sacCer3          # downloads, unpacks, hashes, prints the row
$ genome table-row sacCer3 --json   # the same row as a JSON object
```

#### One-step pipeline: `fetch_genome`

`fetch_fasta` only gets the FASTA onto disk. `fetch_genome` runs the **whole**
download-and-prepare pipeline in a single call — download, decompress, index, and derive
the companion files — and returns the same `GenomeFiles` record as
[`prepare_fasta`](#preparing-a-fasta):

```python
from genome.io.download import UCSCGenomeDownloader

dl = UCSCGenomeDownloader("hg38")
files = dl.fetch_genome()      # download -> decompress -> faidx -> 2bit -> chrom.sizes

files.fasta        # hg38.fa
files.fai          # hg38.fa.fai
files.twobit       # hg38.2bit
files.chrom_sizes  # hg38.chrom.sizes
```

Everything lands under the assembly's reference directory
(`<LIULAB_DATA>/genome/hg38/`), and the gzipped download is kept alongside the outputs.
The unpacked FASTA is checked against the assembly's pinned checksum before anything is
derived from it (see [Pinned sources and checksums](#pinned-sources-and-checksums)); pass
`known_hash="md5:…"` to additionally have pooch verify the compressed download.

!!! note "Indexing the 2bit"
    The 2bit format is **self-indexed** — it carries an internal per-sequence index,
    which is what `twoBitInfo` reads to write `chrom.sizes`. There is no separate sidecar
    index for it the way `.fai` is for the FASTA, so the pipeline indexes the FASTA
    (`.fai`) and materializes the 2bit's internal index as `chrom.sizes`.

#### Seeding from a local file or URL

When the UCSC golden path is unreachable (firewall/proxy) or you have a custom
reference, `fetch_genome_from` prepares the genome from a FASTA **you** provide
instead of downloading from UCSC. The source is either a local path (copied into
the cache) or an `http(s)`/`ftp`/`sftp` URL (fetched with pooch); a gzipped (`.gz`)
source is decompressed. UCSC is never contacted — there is no assembly-name
validation — and neither the pinned source nor the pinned checksum of any table row
is consulted: what you hand over is what you get.

```python
from genome.io.download import UCSCGenomeDownloader

dl = UCSCGenomeDownloader("ce11")

# from a local file
files = dl.fetch_genome_from("/data/ce11.fa.gz")

# or from a non-UCSC URL, e.g. a mirror
files = dl.fetch_genome_from(
    "https://hgdownload-euro.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.fa.gz"
)
```

The source is normalized to `<assembly>.fa` in the cache and then run through the
same `prepare_fasta` pipeline as `fetch_genome`, so it returns the identical
`GenomeFiles` record and the outputs land in the same per-assembly directory. A
fresh `<assembly>.fa` is reused on later calls unless you pass `overwrite=True`.
The [`Genome`](genome.md#seeding-from-your-own-fasta-offline-mirrors-custom-references)
constructor's `path_or_url=` argument is the high-level front door to this.

!!! note "`sftp://` sources need `paramiko`"
    pooch handles `http(s)` and `ftp` out of the box. An `sftp://` source additionally
    needs `paramiko`, which is not a dependency of this package — install it yourself if
    you need that scheme.

#### Where files are stored

Unlike `Downloader` (which uses pooch's per-user cache), `UCSCGenomeDownloader` stores
each assembly under a **per-assembly reference directory** so the FASTA, its indexes,
and any other reference files for that build live together:

```
<LIULAB_DATA>/genome/<assembly>/
```

The lab data root comes from the `LIULAB_DATA` environment variable and falls back to
`~/liulab_data` when it is unset:

```python
from genome.io.download import liulab_data_dir, assembly_data_dir

liulab_data_dir()            # $LIULAB_DATA, or ~/liulab_data
assembly_data_dir("hg38")    # <LIULAB_DATA>/genome/hg38
```

Pass `cache_dir=` to `UCSCGenomeDownloader` to override this default for a single
downloader.

## Preparing a FASTA

`prepare_fasta` runs the three standard preparation steps and returns a `GenomeFiles`
record of every path it wrote:

```python
from genome.io.fasta import prepare_fasta

files = prepare_fasta("hg38.fa")
files.fai          # hg38.fa.fai      (samtools faidx index)
files.twobit       # hg38.2bit        (faToTwoBit)
files.chrom_sizes  # hg38.chrom.sizes (twoBitInfo: name<TAB>length per sequence)
```

| Step | Tool | Output (default) |
|------|------|------------------|
| Random-access index | `samtools faidx` | `<fasta>.fai` |
| Compact encoding | `faToTwoBit` | `<fasta>.2bit` |
| Chromosome sizes | `twoBitInfo` | `<fasta>.chrom.sizes` |

The individual steps are also exposed — `faidx`, `fasta_to_2bit`, and
`twobit_to_chrom_sizes` — each accepting an explicit destination path if you don't want
the sibling-of-the-input default. They raise `FileNotFoundError` for a missing input and
`RuntimeError` (carrying the tool's stderr) if the binary fails.

### Caching

Every step is cached at the command-running step. Before invoking a binary, the step
checks whether its output already exists and is **newer than its input** (the same
freshness rule `make` uses); if so, the tool is skipped and the existing file reused.
This makes re-running `prepare_fasta` (or `fetch_genome`) cheap and idempotent, while a
regenerated input correctly invalidates downstream outputs.

```python
prepare_fasta("hg38.fa")                  # first call: runs all three tools
prepare_fasta("hg38.fa")                  # again: every step served from cache
prepare_fasta("hg38.fa", overwrite=True)  # force regeneration
```

Pass `overwrite=True` to any of `faidx`, `fasta_to_2bit`, `twobit_to_chrom_sizes`,
`prepare_fasta`, or `UCSCGenomeDownloader.fetch_genome` to rebuild unconditionally. The
download and decompression are cached independently by pooch and are unaffected by
`overwrite`.

!!! note "Binaries come from pixi"
    `samtools`, `faToTwoBit`, and `twoBitInfo` are conda/bioconda runtime dependencies.
    Run inside `pixi shell` (or via `pixi run`) so they are on `PATH`; otherwise these
    functions raise `ToolNotFoundError` with a `pixi add` hint.
