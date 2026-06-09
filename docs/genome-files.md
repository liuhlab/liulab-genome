# Downloading and preparing genomes

Two modules in `genome.io` cover getting a reference genome onto disk and turning it
into the index/companion files most tools expect:

- [`genome.io.download`](#downloading) — cache large files from the network with
  [pooch](https://www.fatiando.org/pooch/), with a UCSC-aware subclass.
- [`genome.io.fasta`](#preparing-a-fasta) — shell out to `samtools`, `faToTwoBit`, and
  `twoBitInfo` to index a FASTA and derive its `.2bit` and `chrom.sizes`.

Both live at the **I/O boundary**: they touch the network and the filesystem and invoke
native binaries (managed by pixi). They never reimplement what those binaries do.

## Downloading

`Downloader` is a thin wrapper over `pooch.retrieve`. It downloads a URL once and caches
it, so asking for the same file again is served from disk.

```python
from genome.io.download import Downloader

dl = Downloader()                      # caches under pooch.os_cache("genome")
path = dl.fetch("https://example.org/annotation.bed.gz")
```

Pass `cache_dir=` to redirect the cache (e.g. shared lab scratch). Supply
`known_hash="md5:…"` to verify the download; when omitted, pooch logs the computed hash
so you can pin it next time.

### Genome FASTA from UCSC

`UCSCGenomeDownloader` knows the golden-path layout. Give it an assembly name and it
builds the URL, downloads `<assembly>.fa.gz`, and (by default) decompresses it to
`<assembly>.fa`:

```python
from genome.io.download import UCSCGenomeDownloader

dl = UCSCGenomeDownloader("hg38")
dl.fasta_url
# 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

fasta = dl.fetch_fasta()               # -> Path to hg38.fa (cached, multi-GB)
```

Use `decompress=False` to keep the `.fa.gz`, or `progressbar=False` to silence the bar.

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
Pass `known_hash="md5:…"` to verify the download.

!!! note "Indexing the 2bit"
    The 2bit format is **self-indexed** — it carries an internal per-sequence index,
    which is what `twoBitInfo` reads to write `chrom.sizes`. There is no separate sidecar
    index for it the way `.fai` is for the FASTA, so the pipeline indexes the FASTA
    (`.fai`) and materializes the 2bit's internal index as `chrom.sizes`.

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
