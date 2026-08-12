# Downloading and preparing genomes

Three modules in `genome.io` cover getting a reference genome onto disk and turning it
into the index/companion files most tools expect:

- [`genome.io.download`](#downloading) — fetch large files from the network with
  [pooch](https://www.fatiando.org/pooch/), with a UCSC-aware subclass.
- [`genome.io.fasta`](#preparing-a-fasta) — shell out to `samtools`, `faToTwoBit`, and
  `twoBitInfo` to index a FASTA and derive its `.2bit` and `chrom.sizes`.
- [`genome.io.completion`](#the-registration-record) — the record a finished
  registration writes, which is the only thing that says it finished.

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

fasta = dl.fetch_fasta()               # -> Path to .work/hg38.fa (multi-GB)
```

Both the download and its unpacked form land in the assembly's **working area**, not
beside its prepared files — see [The working area](#the-working-area). Use
`decompress=False` to get the `.fa.gz` back instead, or `progressbar=False` to silence
the bar. When the URL had to be derived, a mistyped assembly fails fast on a `HEAD`
request with a clear error naming the unknown name; when a row pins the URL there is
nothing to guess, so that check is skipped.

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
(`<LIULAB_DATA>/genome/hg38/`). The unpacked FASTA is checked against the assembly's
pinned checksum before it is moved there, so a FASTA that is not the pinned one never
arrives (see [Pinned sources and checksums](#pinned-sources-and-checksums)); pass
`known_hash="md5:…"` to additionally have pooch verify the compressed download.

#### The registration record

`fetch_genome` finishes by writing one **completion record**, `.completion.json`, in the
assembly directory. It holds the URL that was fetched, the sha256 of the unpacked FASTA,
every file it claims with that file's size, the versions of the native tools that
prepared them, the package version, and when it finished:

```python
from genome.io.completion import read_record

record = read_record(assembly_data_dir("hg38"))
record.source_url     # 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
record.sha256         # digest of the unpacked hg38.fa
record.files          # {'hg38.fa': 3099922541, 'hg38.fa.fai': 22262, ...}
record.tool_versions  # {'samtools': 'samtools 1.21'}
record.completed_at   # '2026-08-12T09:14:03+00:00'
```

That record is the **only** thing that answers "is this registered?". It is written
last, after every derived file exists, and atomically, so it is never seen half-written.
Preparing the same assembly again reads it, confirms every file it claims is present and
the right size — presence and size only, no file contents — and returns without
fetching anything. Reopening a prepared human genome therefore costs four `stat` calls
rather than a pass over three gigabytes.

Paths in the record are relative to the assembly directory, so the whole directory can be
moved without invalidating it.

#### When a registration cannot be trusted

A directory that contradicts its record is an error rather than something to rebuild
quietly or trust quietly (ADR-0007). There are four answers, and only one of them raises
in two flavours:

| What is on disk | What happens |
|---|---|
| Nothing, or only `.work/` | A fresh registration — it proceeds normally |
| A record whose every claim holds | Registered — returned without fetching anything |
| Files, but no `.completion.json` | `UnfinishedRegistrationError` — an interrupted or preempted run |
| A record disagreeing with disk | `RegistrationMismatchError`, naming every file that differs and how |

Both errors quote the command that repairs them:

```console
$ genome register hg38
error: /data/genome/hg38 disagrees with its .completion.json: hg38.2bit: recorded
841756144 bytes, found 0. Something changed these files after they were registered.
Re-register it with `genome register hg38 --force`.
```

A forced re-registration keeps whatever is provably good: if the unpacked `<assembly>.fa`
is present and its sha256 is still the one the row pins, it is kept and only the `.fai`,
`.2bit` and `chrom.sizes` are rebuilt, so a preempted cluster job costs seconds rather
than a gigabyte. If that FASTA is missing, if its digest is a different one, or if the row
pins no digest at all — leaving nothing to prove it against — the source is fetched again.

!!! warning "Interrupting a first registration leaves a state that raises"
    While a first registration is running, the directory holds files that no record claims
    yet. A run killed in that window therefore leaves exactly the first error above: the
    next call raises and needs `genome register <assembly> --force`. That is the accepted
    trade-off, chosen over silently resuming — a complete build whose record never landed
    and the wreckage of one killed half-way are indistinguishable on disk, and guessing
    between them risks answering sequence queries from a partial genome. The archive stays
    in `.work/`, so repairing usually downloads nothing.

Deleting `.completion.json` by hand does not reset an assembly, then — it produces the
first of the two errors. `--force` (`overwrite=True` in Python) is the supported way to
register again from scratch.

!!! note "pooch is a downloader; its cache is deliberately not relied on"
    pooch is used here to move bytes — for its retries, progress reporting and
    ftp/sftp support — and for nothing else. Its own cache is not what decides whether
    a file is already here and usable: the completion record owns that judgment, because
    it knows about the unpacked FASTA and the three derived files, and pooch only ever
    knew about the archive it downloaded.

#### By name: `register_assembly` and `verify_assembly`

Two module functions wrap all of the above so that naming an assembly is enough. They are
what `genome register` and `genome verify` call, so a script and the CLI hit one code path:

```python
from genome.io.download import register_assembly, verify_assembly

register_assembly("sacCer3")                            # fetch, verify, prepare, record
register_assembly("sacCer3", force=True)                # repair a directory that raises
register_assembly("ce11", source="/data/ce11.fa.gz")    # seed from your own FASTA

verify_assembly("sacCer3")                              # re-read and re-hash what is registered
verify_assembly("sacCer3", fasta="/tmp/from-a-colleague.fa")
```

`register_assembly` returns the completion record's own fields plus the `directory` they
live in, ready to serialize — that is exactly what `--json` prints.

Verifying is the one operation that reads bytes rather than sizes. Registering and
reopening go by presence and size, which is what makes them instant; this is the
deliberate re-check for when integrity is actually in doubt, and it costs a full pass over
the file. Because the pinned digest covers unpacked content, a FASTA you were handed or
copied from a mirror is checkable against the official row before anything is built on it
— that is what `fasta=` is for, and it needs nothing registered.

#### The working area

Downloads go to `<LIULAB_DATA>/genome/<assembly>/.work/`, a hidden, disposable directory
inside the assembly's own directory:

- it is on the same filesystem as the outputs, so placing the unpacked FASTA is a rename
  rather than a second multi-gigabyte copy;
- nothing in it is ever claimed by the completion record, so it holds working state and
  never a result;
- it survives an interrupted run — the archive stays put, so a preempted cluster job
  repairs from it instead of downloading a genome again;
- the whole thing, archive included, is deleted as soon as the record is written.

Prepared assemblies used to keep their `.fa.gz` beside the FASTA forever. It is plain
gzip rather than bgzip, so no external tool can read it in place; on one developer
machine holding three modest genomes that was 770 MB of pure duplication.

!!! note "Indexing the 2bit"
    The 2bit format is **self-indexed** — it carries an internal per-sequence index,
    which is what `twoBitInfo` reads to write `chrom.sizes`. There is no separate sidecar
    index for it the way `.fai` is for the FASTA, so the pipeline indexes the FASTA
    (`.fai`) and materializes the 2bit's internal index as `chrom.sizes`.

#### Seeding from a local file or URL

When the UCSC golden path is unreachable (firewall/proxy) or you have a custom
reference, `fetch_genome_from` prepares the genome from a FASTA **you** provide
instead of downloading from UCSC. The source is either a local path (copied into
the working area) or an `http(s)`/`ftp`/`sftp` URL (fetched with pooch); a gzipped
(`.gz`) source is decompressed. UCSC is never contacted — there is no assembly-name
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

The source is normalized to `<assembly>.fa` in the assembly directory and then run
through the same `prepare_fasta` pipeline as `fetch_genome`, so it returns the identical
`GenomeFiles` record and the outputs land in the same per-assembly directory. It also
writes the same completion record, with the path or URL you gave as its source and the
digest of what actually arrived — recorded, not compared, since a seeded FASTA is
whatever you handed over. A registered assembly is reused on later calls unless you pass
`overwrite=True`, and a seeded one that needs repairing carries its source into the
message: `genome register ce11 --force --source /data/ce11.fa.gz`, since the plain
command would fetch from a golden path this assembly never came from. The
[`Genome`](genome.md#seeding-from-your-own-fasta-offline-mirrors-custom-references)
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
├── <assembly>.fa               # the unpacked FASTA
├── <assembly>.fa.fai           # samtools faidx index
├── <assembly>.2bit             # faToTwoBit
├── <assembly>.chrom.sizes      # twoBitInfo
├── .completion.json            # the record: written last, read first
└── .work/                      # downloads in progress; gone once the record lands
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
`prepare_fasta`, or `UCSCGenomeDownloader.fetch_genome` to rebuild unconditionally. For
`fetch_genome` that also means ignoring the assembly's
[completion record](#the-registration-record) and registering from scratch; an archive
still sitting in the working area is reused rather than downloaded again.

!!! note "Binaries come from pixi"
    `samtools`, `faToTwoBit`, and `twoBitInfo` are conda/bioconda runtime dependencies.
    Run inside `pixi shell` (or via `pixi run`) so they are on `PATH`; otherwise these
    functions raise `ToolNotFoundError` with a `pixi add` hint.
