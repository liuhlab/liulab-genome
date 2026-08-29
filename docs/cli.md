# CLI

The `genome` command is a thin wrapper over the Python API — every command calls the same
function you would call yourself, so a shell and a notebook take one code path.

Every command takes `--json` for machine-readable output. Errors go to stderr with a
non-zero exit code: `2` for invalid input, `1` for everything else (a missing native
tool, a failed download, a checksum mismatch, a registration that cannot be trusted).

```console
$ genome --help
$ genome version
2026.6.0
$ genome version --json
{"version": "2026.6.0"}
```

## `genome doctor`

Report which native tools are on `PATH`, and at what versions. Exits `1` naming the
install command if one is missing. A tool that is installed but rejects `--version` — the
UCSC binaries do — is listed all the same, since presence is what this answers.

```console
$ genome doctor
samtools: samtools 1.21 ...
faToTwoBit: installed; reports no version
twoBitInfo: installed; reports no version
```

STAR and chromap are not checked here: they are optional, and each one checks for itself
when you ask it to build an index.

## `genome register <assembly>`

Prepare an assembly on disk — fetch, verify, index, record — and print where it landed.

```console
$ genome register sacCer3
registered sacCer3 in /data/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

Running it again on a registered assembly downloads nothing.

| Option | Effect |
|---|---|
| `--source PATH_OR_URL` | Seed from your own FASTA instead of the pinned source |
| `--force` | Register again from scratch — the repair for a directory that raises |

A directory that cannot be trusted exits `1` naming the file and the repair:

```console
$ genome register sacCer3
error: /data/liulab_data/genome/sacCer3 disagrees with its .completion.json:
sacCer3.2bit: recorded 3145728 bytes, found 0. Something changed these files after
they were registered. Re-register it with `genome register sacCer3 --force`.
```

`--force` keeps an unpacked FASTA whose checksum still matches and rebuilds only the
derived files; it fetches again when it cannot prove that. Interrupting a *first*
registration leaves a directory that needs exactly this.

### Naming a chimera builds it

An assembly named after two or more registered assemblies, sorted and joined by `_`, is a
[chimera](genome.md#chimera-assemblies): it is concatenated from those components, merged
annotation included, and nothing is downloaded.

```console
$ genome register ce11_ecHT115
registered ce11_ecHT115 in /data/liulab_data/genome/ce11_ecHT115
  components  ce11, ecHT115
  sha256  0f3d6d5e…
  files   ce11_ecHT115.2bit, ce11_ecHT115.chrom.sizes, ce11_ecHT115.fa, ce11_ecHT115.fa.fai
  annotation  wormbase_ws298+refseq_rs_2025_06_26 — the components' own, merged and registered by this build
```

There is no flag for listing the parts — the name carries them. A component this machine
has not prepared exits `1` naming the command that prepares it, and components typed in
the wrong order exit `1` naming the canonical spelling.

## `genome register-annotation <assembly> <name>`

Register one of the annotations the table lists for an assembly: fetch, check the
checksum and the chromosome names, build the gffutils database, record it.

```console
$ genome register-annotation sacCer3 ensgene_v101
registered ensgene_v101 for sacCer3 in /data/liulab_data/genome/sacCer3/gtf/ensgene_v101
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
  sha256  d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026
  files   ensgene_v101.db, ensgene_v101.gtf
  chromosomes checked — every name the GTF uses is one the assembly carries
```

## `genome register-gtf <assembly> <path> <name>`

The same, for a GTF the table does not list. Nothing is downloaded and no checksum is
compared against; `<name>` is what addresses it from then on. A `.gz` is decompressed.

```console
$ genome register-gtf sacCer3 ~/annotations/sacCer3.WS298.gtf wormbase_ws298
```

Both registration commands share these options:

| Option | Effect |
|---|---|
| `--force` | Register again from scratch; also the repair for a broken directory |
| `--no-check-chromosomes` | Register a GTF whose chromosome-name mismatch you accept |
| `--infer-genes` / `--infer-transcripts` | Reconstruct those features from exon lines — for a bare exon-level GTF only |

The last line of the output always says whether the chromosome names were checked, and
if not, why: because you stood the check down, or because the assembly is not registered
yet and there was no `chrom.sizes` to check against. Only the second is worth acting on.

## `genome annotations <assembly>`

What the table offers against what is registered here. Downloads and prepares nothing, so
it works for an assembly you have never registered.

```console
$ genome annotations hg38
annotations for hg38 in /data/liulab_data/genome/hg38
  gencode_v50  offered, not registered  GENCODE v50
  mine         registered, not offered
default: gencode_v50 — not registered here; register it with `genome register-annotation hg38 gencode_v50`
```

An annotation that is here but cannot be trusted reads as `broken`, with the problem and
the repair on the line below. Exit stays `0` — one broken annotation never hides the
others.

## `genome gene-list <assembly> <category>`

Print the gene ids an annotation puts in one [gene category](genome.md#which-genes-are-in-a-category),
one per line. Only the ids go to stdout, so the output pipes; the heading and the
per-source attribution go to stderr.

```console
$ genome gene-list ce11 rRNA > rrna.txt
rRNA for ce11 / wormbase_ws298
  wormbase_ws298  20
```

`--annotation NAME` asks about a registered annotation other than the assembly's default
one. `--json` carries the same answer with the sources kept apart:

```console
$ genome gene-list ce11_ecHT115 rRNA --json
{"assembly": "ce11_ecHT115", "annotation": "wormbase_ws298+refseq_rs_2025_06_26",
 "category": "rRNA", "gene_ids": ["WBGene00004512", …], "sources": [{"component": "ce11", …}]}
```

Exits `1` when the annotation is not registered here, when no curated gene list ships for
it, and when it declares categories but not this one — three different facts, each with
its own message. None of them prints an empty list of genes.

## `genome gene-categories <assembly>`

Which categories that annotation declares, and how many genes are in each — what
`genome gene-list` may be asked for. A merged annotation shows the per-component split.

```console
$ genome gene-categories ce11_ecHT115
categories for ce11_ecHT115 / wormbase_ws298+refseq_rs_2025_06_26
  rRNA  39  (ce11: 23, ecHT115: 16)
```

`--json` emits every category with its gene ids and its sources — the same answer
`genome gene-list` gives for one of them, for all of them at once.

## `genome verify <assembly>`

Re-read a FASTA and check its sha256 against the digest pinned for the assembly. This is
the deliberate full-file re-check; registering and reopening go by size alone. Exits `1`
on a mismatch.

```console
$ genome verify sacCer3
/data/liulab_data/genome/sacCer3/sacCer3.fa: sha256 6ff72f07… matches the digest pinned for it

$ genome verify sacCer3 --fasta /tmp/from-a-colleague.fa
error: sha256 mismatch for /tmp/from-a-colleague.fa: expected 6ff72f07…, got 9316629b….
```

`--fasta` checks any file against the official row — useful for a copy you were handed,
before you build anything on it.

A chimera is also checked against its components: the closing line reads `unchanged` or,
where a digest was missing on one side, `unknown`. A component that is no longer the one
the chimera was built from exits `1`.

## `genome revcomp <sequence>`

Reverse-complement a DNA sequence. Case is preserved; a character outside the `DNA`
alphabet exits `2` naming it. The alphabet is the type's own — `DNA.ALPHABET`, read at the
boundary rather than spelled again here — so the [`outside_alphabet`](sequences.md#construction)
you would call yourself is the check this command applies.

```console
$ genome revcomp ATCG
CGAT
$ genome revcomp aTcG --json
{"input": "aTcG", "reverse_complement": "CgAt"}
```

## `genome table-row <assembly>`

Download an assembly, unpack it, hash the unpacked FASTA, and print the line to paste
into the shipped metadata table — how the checksum column gets filled in. A checksum the
table already pins is reported, never enforced.

```bash
$ genome table-row sacCer3
sacCer3	Saccharomyces cerevisiae	sacCer3	R64-1-1	GCF_000146045.2	559292	https://…/sacCer3.fa.gz	6ff72f07…
```

A chimera exits `1` before anything is downloaded: its row carries the name and nothing
else, so there is no row here to compute. Check one with `genome verify` instead.

## `genome motif-scan <fasta> <output>`

Scan a FASTA with a JASPAR release, write the hits to Parquet, and print a summary of the
run. The batch case — a shell script, a scheduler job — and the only motif command there
is: listing, plotting and comparing motifs are notebook work.

```console
$ genome motif-scan peaks.fa hits.parquet --release 2024 --tax-group all
scanned 41255 sequences with 2338 motifs from JASPAR 2024 all
  background  0.293, 0.207, 0.207, 0.293
  threshold   0.0001
  skipped     8 under 7 positions, so not scanned: MA0261.1, MA2355.1, …
  workers     16
  hits        2841193 -> hits.parquet
```

**The hits go to the named file and the summary to standard output**, so `--json` is never
corrupted by table data and a pipeline can consume the one while a downstream step reads
the other:

```console
$ genome motif-scan peaks.fa hits.parquet --release 2024 --tax-group all --json
{"release": "2024", "tax_group": "all", "motifs_scanned": 2338,
 "motifs_skipped": ["MA0261.1", "MA2355.1", …], "background": [0.293, 0.207, 0.207, 0.293],
 "threshold": 0.0001, "sequences_scanned": 41255, "hits_written": 2841193,
 "workers": 16, "output": "hits.parquet"}
```

A motif under seven positions is named there rather than scanned: a 6-mer's best possible
word has p = 2.44e-4, so it cannot reach the default threshold at all, and an engine asked
for it anyway would fall back to that matrix's best attainable cutoff and over-call in
silence.

Read the hits back with `genome.tf.motif.read_hits(path)`, which restores the compact
dtypes and the provenance — `pandas.read_parquet` gives the rows and drops what the scan
was.

| Option | Effect |
|---|---|
| `--release` | Which JASPAR release: `2024` or `2026` (default) |
| `--tax-group` | Which taxonomic group's file: `vertebrates` (default) through `all` |
| `--threshold` | The per-position p-value each motif's cutoff is converted from; `1e-4` by default |
| `--background` | `auto` (default), `uniform`, or `derive` — see below |
| `--workers` | How many processes; every core the allocation granted by default |

`--background` decides the answer more than any other option: `auto` derives the base
composition from the input when it holds at least 10 000 unambiguous bases and stays
uniform below that, `uniform` pins it, and `derive` derives whatever the input holds.
Whichever it was is in the summary and on the hits. Four frequencies of your own are a
Python call (`scan_fasta(background=[…])`) rather than a flag — a mistyped one on a
command line would change every cutoff and look like a scan that simply found fewer hits.

**It defaults to every core the allocation granted** — `SLURM_CPUS_PER_TASK`, then
`SLURM_CPUS_ON_NODE`, then this process's CPU affinity, then the machine — where the
library defaults to one worker. A console script is a proper entry point, so the
process-pool hazard behind that library default does not apply here. `--workers 1` scans
serially and produces the identical table.

**Prepare the release from a login node.** Naming a release prepares it, and the first
time that is a download. The lab's CPU cluster compute nodes have no internet, so run this
once on a login node before submitting a job that needs it; every run afterwards reads the
cached file out of `<LIULAB_DATA>/motif/jaspar/`, shared by every project on the machine.

Exits `1` when the FASTA is not there or is not FASTA, when the release or tax group is
not one this package prepares, when the threshold is not a p-value in `(0, 1)`, when the
worker count is below 1, and when the release is not cached here and cannot be fetched —
which is what a compute node with no internet looks like. A `--background` that is not one
of the three modes is refused by the parser with `2`, before anything is fetched or read.
