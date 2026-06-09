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
`<assembly>.fa` in the cache:

```python
from genome.io.download import UCSCGenomeDownloader

dl = UCSCGenomeDownloader("hg38")
dl.fasta_url
# 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'

fasta = dl.fetch_fasta()               # -> Path to hg38.fa (cached, multi-GB)
```

Use `decompress=False` to keep the `.fa.gz`, or `progressbar=False` to silence the bar.

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

!!! note "Binaries come from pixi"
    `samtools`, `faToTwoBit`, and `twoBitInfo` are conda/bioconda runtime dependencies.
    Run inside `pixi shell` (or via `pixi run`) so they are on `PATH`; otherwise these
    functions raise `ToolNotFoundError` with a `pixi add` hint.
