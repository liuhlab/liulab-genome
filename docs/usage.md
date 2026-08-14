# Usage

## Installation

The package depends on native tools from bioconda (`samtools`, `bedtools`), so
the recommended path is [pixi](https://pixi.sh):

```bash
git clone https://github.com/lhqing/liulab-genome.git
cd liulab-genome
pixi install --locked        # solve from the committed pixi.lock
pixi shell                   # activate the env in your shell
```

If you only need the Python API and already have `samtools` / `bedtools` on
your `PATH` from somewhere else, a pip install also works:

```bash
pip install liulab-genome
```

You can verify the native toolchain is reachable at any time:

```bash
$ genome doctor
samtools: samtools 1.21 ...
bedtools: bedtools v2.31.1
```

`genome doctor` exits with code `1` and an actionable message if either tool
is missing.

## Command-line interface

The `genome` command is a thin Typer wrapper over the Python API — every
command has a corresponding function in the package, so scripts and notebooks
can call the same code paths without shelling out.

```bash
$ genome --help
$ genome version
$ genome revcomp ATCG
$ genome revcomp aTcG --json
$ genome doctor [--json]
$ genome register hg38 [--source PATH_OR_URL] [--force] [--json]
$ genome register-annotation hg38 gencode_v50 [--force] [--no-check-chromosomes] [--infer-genes] [--infer-transcripts] [--json]
$ genome register-gtf hg38 path/to/annotation.gtf mine [--force] [--no-check-chromosomes] [--infer-genes] [--infer-transcripts] [--json]
$ genome annotations hg38 [--json]
$ genome verify hg38 [--fasta PATH] [--json]
$ genome table-row sacCer3 [--json]
```

`--json` toggles machine-readable output. Errors go to stderr and use non-zero
exit codes (`2` for invalid input, `1` for missing native tools, a failed download,
a checksum that does not match, or a registration that cannot be trusted).

### Registering an assembly

`genome register` prepares an assembly on disk — fetch, verify, index, and write the
record that says it finished — and prints where everything landed:

```console
$ genome register sacCer3
registered sacCer3 in /data/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

Running it again on a registered assembly reads the record and downloads nothing. A
directory that cannot be trusted — files with no record, or a record that disagrees with
what is on disk — exits `1` naming the file and the repair, which is this same command
with `--force`:

```console
$ genome register sacCer3
error: /data/liulab_data/genome/sacCer3 disagrees with its .completion.json:
sacCer3.2bit: recorded 3145728 bytes, found 0. Something changed these files after
they were registered. Re-register it with `genome register sacCer3 --force`.
$ genome register sacCer3 --force
```

`--force` keeps an unpacked FASTA whose checksum still matches and rebuilds only the
derived files; it fetches the source again when it cannot prove that. Interrupting a
**first** registration leaves a directory that raises next time and needs exactly this —
see [When a registration cannot be trusted](genome-files.md#when-a-registration-cannot-be-trusted).

### Registering an annotation

`genome register-annotation` does the same job one level down. The shipped annotation
table lists what each assembly officially supports, so naming one is enough — it is
fetched, the unpacked GTF is checked against the pinned checksum, the gffutils database
is built, and the record that says all of it finished is written last:

```console
$ genome register-annotation sacCer3 ensgene_v101
registered ensgene_v101 for sacCer3 in /data/liulab_data/genome/sacCer3/gtf/ensgene_v101
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
  sha256  d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026
  files   ensgene_v101.db, ensgene_v101.gtf
```

Running it again on a registered annotation reads the record and downloads nothing. It
exits `1` when the table lists no such annotation for that assembly (the message says
what it does list), when the GTF names chromosomes the assembly does not carry, when the
GTF is not the checksum pinned for it, or when the directory cannot be trusted — a
half-built database from an interrupted run, say, which `--force` repairs.

`--no-check-chromosomes` registers one whose chromosome-name mismatch you have looked at
and accept; the record says the names went unchecked, so you can tell months later. See
[The chromosome names have to match](aligner.md#the-chromosome-names-have-to-match).

### Registering a GTF the table does not list

`genome register-gtf` is the escape hatch: you say where the file is, and it is placed,
checked, built and recorded exactly as a listed annotation is. Nothing is downloaded and
no checksum is compared against — an unlisted GTF has none pinned for it — and the name
you give is what addresses it from then on:

```console
$ genome register-gtf sacCer3 ~/annotations/sacCer3.WS298.gtf wormbase_ws298
registered wormbase_ws298 for sacCer3 in /data/liulab_data/genome/sacCer3/gtf/wormbase_ws298
  source  /home/you/annotations/sacCer3.WS298.gtf
  sha256  9e1f0a5c6d2b8e4a1c7f3b0d5a8e2c4f6b9d1e3a7c5f8b2d4e6a0c9f1b3d5e7a
  files   wormbase_ws298.db, wormbase_ws298.gtf
```

A `.gz` source is decompressed on the way in. Naming the assembly is what lets the
chromosome check find its `chrom.sizes` without being told where it is, so an unlisted
GTF is held to the same standard as a listed one — `--no-check-chromosomes` stands the
check down here too. It exits `1` when the file is not there, when its chromosome names
do not line up, or when the directory cannot be trusted, and `--force` repairs the last
of those.

Both registration commands take `--infer-genes` and `--infer-transcripts`, which
reconstruct those features from exon lines. Leave them off for a GENCODE, Ensembl or
RefSeq GTF, which declares both already; turn them on for a **bare exon-level GTF**, one
whose only lines are exons, which otherwise registers as a database of exons and nothing
else. See [Registering an annotation](aligner.md#registering-an-annotation).

### Listing what an assembly offers against what is registered

`genome annotations` answers two questions side by side — which annotations the lab
supports for an assembly, and which are registered on this machine. It downloads and
prepares nothing, so it works for an assembly you have never registered, which is the
case it most needs to serve:

```console
$ genome annotations hg38
annotations for hg38 in /data/liulab_data/genome/hg38
  gencode_v50  offered, not registered  GENCODE v50
  mine         registered, not offered
default: gencode_v50 — not registered here; register it with `genome register-annotation hg38 gencode_v50`
```

The last line names the assembly's [default annotation](aligner.md#the-default-annotation)
— the one the table flags, else the sole registered one — and, when it is not registered
here, the command that closes the gap.

An annotation whose directory is here but cannot be trusted — files left by an
interrupted run, or a record that no longer matches what is on disk — reads as `broken`
rather than as one nobody has fetched, whether or not any table row lists it. The line
under it says what is wrong and names the command that repairs it, which is a command
you can paste:

```console
$ genome annotations hg38
annotations for hg38 in /data/liulab_data/genome/hg38
  gencode_v50  broken                   GENCODE v50
      /data/liulab_data/genome/hg38/gtf/gencode_v50 holds files but no .completion.json, so a previous run was interrupted before it finished: gencode_v50.db. Nothing here can be trusted as complete. Re-register it with `genome register-annotation hg38 gencode_v50 --force`.
default: gencode_v50 — broken here; repair it with `genome register-annotation hg38 gencode_v50 --force`
```

This is where such a thing is discovered, so it is reported and not raised over: the
exit code is still `0`, and one broken annotation never hides the ones beside it. With
`--json`, the row carries `broken`, the `problem` and the `repair`. An unlisted one is
repaired from the GTF it was built from, which its record remembers — and when no record
survives to say, the command is printed with `<path>` where that file goes rather than a
path that would not be there.

### Verifying one

`genome verify` re-reads a FASTA and checks its sha256 against the assembly's pinned
digest, exiting `1` when they differ. Point `--fasta` at any file to check a copy you
were handed against the official row before building anything on it:

```console
$ genome verify sacCer3
/data/liulab_data/genome/sacCer3/sacCer3.fa: sha256 6ff72f07… matches the digest pinned for it

$ genome verify sacCer3 --fasta /tmp/from-a-colleague.fa
error: sha256 mismatch for /tmp/from-a-colleague.fa: expected 6ff72f07…, got 9316629b….
```

### Filling in the metadata table

`genome table-row` downloads an assembly, unpacks it, hashes the unpacked FASTA and
prints the line to paste into the shipped metadata table — how the checksum column
gets filled in without hashing anything by hand:

```bash
$ genome table-row sacCer3
sacCer3	Saccharomyces cerevisiae	sacCer3	R64-1-1	GCF_000146045.2	559292	https://…/sacCer3.fa.gz	6ff72f07…
```

## Python API quickstart

```python
from genome import DNA, RNA, Protein

# Construction preserves the value verbatim (including case). The alphabet is
# NOT enforced — scanning every char is too costly on large sequences.
seq = DNA("aTcGN")  # accepted and stored as-is; no validation
seq = DNA("aTcG")   # OK

# Biological transforms return the right typed result.
seq.complement()                # DNA('tAgC')
seq.reverse_complement()        # DNA('CgAt')
seq.transcribe()                # RNA('aUcG')
seq.gc_content                  # 0.5

# Slicing stays typed.
seq[1:3]                        # DNA('Tc')

# Other str methods (upper, lower, replace, +) return plain str by design;
# wrap explicitly if you want a typed result back.
DNA(seq.upper())                # DNA('ATCG')
```

See [Sequences](sequences.md) for the full contract, including what happens
on invalid input and which `str` methods do or do not preserve the type.

## Running the gates locally

The repository ships with a single `pixi run check` that runs lint, format
check, type check, and the full test suite — the same commands CI runs:

```bash
pixi run check
pixi run test-cov              # tests + coverage report
pixi run docs                  # live MkDocs server at http://127.0.0.1:8000
```

To enable the pre-commit hooks in your local clone (one-off):

```bash
pixi run -- pre-commit install
```
