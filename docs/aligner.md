# Annotations & aligner indexes

A reference assembly is sequence; most analyses also need a **gene annotation**
(a GTF) and an **aligner index**. `Genome` ties these together: you register one
or more GTF annotations on the genome, then build an index with whichever aligner
your pipeline maps with.

Two aligners ship today, and they differ in whether an annotation is involved at
all:

| Aligner | Maps | Annotation | Index |
|---------|------|------------|-------|
| [STAR](https://github.com/alexdobin/STAR) | RNA-seq (splice-aware) | **required** — one index per GTF | `index/star_<gtf>/` — a directory |
| [chromap](https://github.com/haowenz/chromap) | chromatin profiles: ATAC-seq/scATAC-seq, ChIP-seq, Hi-C | none | `index/chromap/chromap.index` — a single file |

```python
from genome import Genome

g = Genome("sacCer3")
g.register_gtf("sacCer3.ensGene.gtf", name="ensembl")   # place + index the GTF
g.build_star_index(gtf="ensembl", threads=8)            # STAR index for that annotation
g.build_chromap_index()                                 # chromap needs no annotation
```

## Registering an annotation

A genome can carry several annotations (GENCODE, Ensembl, RefSeq, …), each under
a unique `name`. The shipped annotation table lists which ones the lab supports
for each assembly, so naming one is enough: it is fetched from the URL that table
pins, its unpacked GTF is checked against the pinned sha256 and against the
assembly's chromosome names, and a
[gffutils](https://gffutils.readthedocs.io/) database is built beside it.

```python
g.register_annotation("ensgene_v101")   # fetch + verify + build + record
g.annotations                           # ['ensgene_v101'] — registered names
g.get_gtf_path("ensgene_v101")          # Path to the placed .gtf
```

`register_gtf` is the escape hatch for a GTF no row lists — you say where the
file is, and it is placed, built and recorded the same way:

```python
g.register_gtf("path/to/annotation.gtf", name="custom")
```

The same from a shell, where the assembly is named rather than opened:

```bash
$ genome register-gtf sacCer3 path/to/annotation.gtf custom
```

Each annotation lives in its own directory next to the sequence files:

```
<LIULAB_DATA>/genome/<assembly>/gtf/<name>/
    <name>.gtf          # the annotation, stored uncompressed (a .gz source is decompressed)
    <name>.db           # the gffutils SQLite database built from it
    .completion.json    # the record saying all of that finished
```

A few things to know:

- **The record is what "registered" means.** A database file's existence proves
  nothing — a build killed half-way leaves one that answers queries with most of
  the genes missing. `g.annotations` lists what has a record that agrees with
  disk; anything else is not registered, however many files are lying there.
- **Re-registering a valid annotation is a silent no-op.** Nothing is fetched,
  nothing is rebuilt, no warning is emitted. A directory that cannot be trusted
  raises instead and names `genome register-annotation <assembly> <name> --force`,
  which is also the repair: it keeps a GTF whose checksum still matches and only
  rebuilds the database.
- **Gene/transcript inference is off by default** (`disable_infer_genes` /
  `disable_infer_transcripts` are `True`; `--infer-genes` / `--infer-transcripts`
  from a shell, on either registration command). Standard annotation GTFs already
  declare `gene`/`transcript` features, and inferring them is the classic
  gffutils slow path. Enable it (`disable_infer_genes=False`) only for a bare
  exon-level GTF that lacks those records — registered with inference off, one
  yields a database of exons and nothing else.

### The chromosome names have to match

An annotation and its assembly are one unit or they are nothing, so every
sequence the GTF names must be one the assembly's `chrom.sizes` carries.
Registering an Ensembl-spelled GTF (`1`, `2`, `MT`) against a UCSC-spelled
assembly (`chr1`, `chr2`, `chrM`) would otherwise build an annotation where every
feature sits on a sequence the assembly has never heard of — nothing lines up,
nothing complains, and every query answers nothing while looking healthy.

```python
g.register_annotation("gencode_v44")
# ChromosomeMismatchError: the GTF for 'gencode_v44' names 25 chromosomes the
# assembly does not carry: 1, 10, 11, 12, 13, 14, 15, 16, 17, 18 (and 15 more).
# An annotation and its assembly must spell chromosomes the same way, and the
# usual cause is a UCSC-versus-Ensembl mismatch ('chr1' against '1', 'chrM'
# against 'MtDNA'). The assembly carries: chr1, chr10, ... Register the
# annotation built for this assembly, or pass check_chromosomes=False —
# --no-check-chromosomes from a shell — to register this one anyway.
```

Three things about it are deliberate:

- **It is strict in one direction only.** The GTF's names must be among the
  assembly's; the reverse is not required. An assembly carrying scaffolds,
  patches or alt contigs the annotation never mentions is completely normal and
  is not an error.
- **It runs before the database build**, and on a GTF that has not been placed
  yet — still in the working area when it was fetched, still where you pointed at
  it when it was not. A mismatch therefore costs one streaming pass over the file
  rather than the many minutes the gffutils build takes, and leaves the annotation
  directory exactly as it found it, so the next call reports the same problem
  rather than an interrupted registration. The GTF is streamed, never loaded: a
  GENCODE GTF is well over a gigabyte unpacked.
- **The override is `check_chromosomes=False`**, on `register_annotation`,
  `register_gtf`, their module-level forms and `--no-check-chromosomes` on either
  registration command. It is for the case where you have looked at the mismatch
  and accept it — one unusual contig should not block a legitimate annotation:

```python
g.register_annotation("gencode_v44", check_chromosomes=False)
```

Whether the names were actually checked is recorded, so you can tell months
later:

```python
import json
record = json.loads((g.get_gtf_path("gencode_v44").parent / ".completion.json").read_text())
record["details"]["chromosomes_checked"]      # True when they were checked
```

That flag is `False` for the override, and also when there was nothing to check
against — registering an annotation for an assembly that has not been prepared
yet means no `chrom.sizes` exists to compare with. `register_gtf` at module level
is given a directory and no assembly name, so it cannot find that file on its own:
pass `chrom_sizes=<path>` to have the names checked, which is what
`Genome.register_gtf` and `genome register-gtf` — both of which are told the
assembly — do for you.

### What the table offers, and what is registered here

These are two questions, and a genome answers them separately — "the lab supports
it" is not "this machine has it":

```python
g.offered_annotations   # [AnnotationMetadata(name='ensgene_v101', provider='UCSC', ...)]
g.annotations           # [] — nothing registered on this machine yet
```

`offered_annotations` is this assembly's rows from the shipped table, in table
order, whether or not anyone has registered them; `annotations` is what has a
valid record on this disk. The same pair from a shell, for an assembly you have
not prepared at all:

```console
$ genome annotations hg38
annotations for hg38 in /data/liulab_data/genome/hg38
  gencode_v50  offered, not registered  GENCODE v50
  mine         registered, not offered
default: gencode_v50 — not registered here; register it with `genome register-annotation hg38 gencode_v50`
```

### The default annotation

The annotation everything falls back to when you name none. Four rules, in order:

1. an explicit `default_gtf=` at construction;
2. the annotation the table flags for this assembly — which is how everyone in
   the lab reaches for the same release without discussing it;
3. the sole registered annotation, when exactly one is registered;
4. otherwise none, because a caller who did not choose between several is asked
   rather than guessed at.

```python
g.default_gtf                            # 'gencode_v50' — the table's flag
Genome("hg38", default_gtf="refseq_2023")  # …unless you say otherwise
```

**`default_gtf` names an annotation; it does not locate one.** On a fresh machine
the table's default is exactly the thing nobody has fetched yet, and opening a
genome must never start a gigabyte download and a database build running many
minutes — so construction succeeds, and asking for the *path* is what says the
gap is there:

```python
g.default_gtf_path
# AnnotationNotRegisteredError: no annotation 'gencode_v50' is registered for 'hg38'.
# Registered here: (none). The annotation table offers it for 'hg38', so register it
# with `genome register-annotation hg38 gencode_v50`.
```

An explicit `default_gtf=` behaves the same way, deliberately: one rule for both,
so naming a default up front and registering it on the next line works. `None`
from `default_gtf_path` means no default was decided at all, which is a different
answer from one that is not registered.

`default_gtf` does **not** change what `build_star_index(gtf=...)` requires — that
argument is always explicit (see below).

## Building a STAR index

[STAR](https://github.com/alexdobin/STAR) is a splice-aware RNA-seq aligner. Its
index (the *genomeDir*) is built from the genome FASTA **plus** a gene
annotation, so an index is always tied to one registered GTF:

```python
g.build_star_index(gtf="ensembl")              # required: which annotation
g.build_star_index(gtf="ensembl", threads=8)   # parallelize the build
```

The `gtf` key is **required**: STAR resolves its path via `get_gtf_path` and
passes it as `--sjdbGTFfile` for splice-junction-aware indexing. Because the
index depends on the annotation, each annotation gets its **own** index
directory, so different GTFs never collide:

```
<LIULAB_DATA>/genome/<assembly>/index/star_<gtf>/
```

For example `gtf="ensembl"` builds under `index/star_ensembl/` and
`gtf="refseq"` under `index/star_refseq/` — both reusable, independent indexes.

### Options

`build_star_index` is a thin pass-through to STAR's `index()`. The commonly
tuned options are named; everything else is forwarded as a raw STAR flag:

```python
g.build_star_index(
    gtf="ensembl",
    sjdb_overhang=99,          # --sjdbOverhang; ideally read_length - 1
    threads=8,                 # --runThreadN
    overwrite=True,            # rebuild even if a finished index exists
    genomeSAindexNbases=11,    # any other genomeGenerate flag, sans the leading --
)
```

- **`sjdb_overhang`** (default `100`) — set to `read_length - 1` for best
  splice-junction sensitivity.
- **`threads`** (default `1`) — build threads.
- **`overwrite`** (default `False`) — a finished index is **cached and reused**;
  pass `overwrite=True` to force a rebuild, and to rebuild over a directory that
  cannot be trusted (see below).
- **arbitrary STAR flags** — pass any `genomeGenerate` option by its STAR name
  without the leading `--` (e.g. `genomeSAindexNbases=11`). For small genomes the
  suffix-array index size (`genomeSAindexNbases`) is auto-reduced unless you set
  it yourself.

The call returns the path to the built genome directory:

```python
index_dir = g.build_star_index(gtf="ensembl")
# .../genome/sacCer3/index/star_ensembl
```

### What a finished index looks like

After a successful build the directory holds STAR's binary files plus one piece of
bookkeeping: `.completion.json`, the same completion record every prepared genome
and registered annotation writes. It is written **last**, once STAR has finished,
so its absence means *unfinished* and never *missing*, and it is the only thing
that is ever asked whether the index is usable — never the presence of `SA` or
`Genome`.

It records what the build claimed and how it ran:

- every file the build left in the directory, with its size, so a file deleted or
  truncated afterwards is caught rather than surfacing as a crash inside STAR;
- the exact command, the resolved parameters (including the `gtf` key and the GTF
  path), the STAR version, and the FASTA consumed — so an index is self-describing
  and can be explained months later;
- the package version and the time it finished.

```python
import json
index_dir = g.get_star_index("ensembl")
record = json.loads((index_dir / ".completion.json").read_text())
record["details"]["command"]        # ['STAR', '--runMode', 'genomeGenerate', ...]
record["details"]["parameters"]     # {'threads': 8, 'gtf': 'ensembl', ...}
record["tool_versions"]             # {'STAR': '2.7.11b'}
```

## Building a chromap index

[chromap](https://github.com/haowenz/chromap) is a fast aligner and preprocessor
for chromatin profiles (ATAC-seq/scATAC-seq, ChIP-seq, Hi-C). Its index is built
from the genome FASTA **alone** — no gene annotation — so there is no `gtf`
argument, and one index serves every use of the assembly:

```python
g.build_chromap_index()
```

Unlike STAR's *genomeDir*, the result is a single file:

```
<LIULAB_DATA>/genome/<assembly>/index/chromap/chromap.index
```

That file is what chromap's `-x/--index` expects at mapping time.

### Options

`build_chromap_index` is a thin pass-through to chromap's `index()`. Only the two
minimizer knobs are named; everything else is forwarded as a raw `--build-index`
flag:

```python
g.build_chromap_index(
    kmer=20,               # -k/--kmer
    window=10,             # -w/--window
    overwrite=True,        # rebuild even if a finished index exists
    min_frag_length=30,    # any other --build-index flag; underscores become hyphens
)
```

- **`kmer`** and **`window`** (both default `None`) — left off the command line
  entirely unless you set them, so chromap's own defaults (17 and 7) apply rather
  than a copy of them frozen here that could drift out of date.
- **`overwrite`** (default `False`) — as with STAR, a finished index is **cached
  and reused**; pass `overwrite=True` to force a rebuild, and to rebuild over a
  directory that cannot be trusted.
- **arbitrary chromap flags** — a keyword maps to chromap's long option with
  underscores turned into hyphens (`min_frag_length=30` → `--min-frag-length 30`);
  a list or tuple expands to several arguments after the flag.
- **there is no `threads` option**, unlike STAR. chromap's `-t/--num-threads` is a
  *mapping* parameter; building the index is single-threaded.

The bookkeeping matches STAR's: the same `.completion.json` gates the index and
records the chromap version, the assembly, the resolved parameters, the exact
command and the one index file it claims.

## Retrieving a built index

Building an index and *using* it are separate jobs — a mapping pipeline needs the
path of an index that was built earlier, possibly by someone else. The `get_*`
methods return that path and build nothing:

```python
g.get_star_index("ensembl")   # .../index/star_ensembl           -> STAR --genomeDir
g.get_chromap_index()         # .../index/chromap/chromap.index  -> chromap -x
```

Both are convenience wrappers over a generic form that takes the aligner name
plus whatever selectors identify the index — the annotation for STAR, nothing for
chromap:

```python
g.get_index("star", gtf="ensembl")
g.get_index("chromap")
```

The name is case-insensitive; an unknown aligner raises `ValueError` listing the
ones that are known. Otherwise the index's completion record is read, and each way
it can fail to vouch for the directory gets its own error — every one of them a
`RuntimeError`, and every one naming the call that puts it right:

**Nothing was ever built there** — `IndexNotBuiltError`:

```python
Genome("sacCer3").get_chromap_index()
# IndexNotBuiltError: no chromap index for 'sacCer3' at .../index/chromap:
# nothing has been built there yet. Build it with `Genome.build_chromap_index()`.
```

**A build was interrupted**, leaving index files nothing vouches for —
`UnfinishedRegistrationError`. This is *not* rebuilt silently: the files may be a
complete index whose record was never written, or the wreckage of one killed
half-way, and nothing on disk tells them apart. Rebuild it deliberately with
`Genome.build_chromap_index(overwrite=True)`.

**A file changed after the build** — deleted, truncated, or replaced —
`RegistrationMismatchError`, naming every file that differs and by how much. Same
repair: rebuild with `overwrite=True`.

The last two are the strict-failure trade-off, and it is worth stating plainly:
during a build the directory briefly holds files with no record, so interrupting
one leaves a state that raises next time and needs a forced rebuild. That is
chosen over silently resuming.

One thing to know: looking a path up still constructs the aligner object, and
that constructor checks the binary. Retrieving an index path therefore requires
the aligner to be installed, even though nothing is built.

## The aligners are optional dependencies

Neither STAR nor chromap is in the default environment — an aligner is only
needed when you actually build or look up an index, and each checks for itself at
that moment. Install whichever you need:

```bash
pixi add star            # from bioconda
pixi add chromap         # from bioconda
```

Both are already provisioned together in the project's `aligners` environment
(`pixi run -e aligners ...`). If the binary is missing, the call fails fast with
install instructions and a `genome.external.ToolNotFoundError`, rather than a
cryptic error deep in the run.

## Domain invariant: the annotation is always explicit

Mirroring the assembly rule, the annotation an index is built against is never
implicit. `build_star_index` requires a named `gtf`, the index directory encodes
that name, and the completion record records it — so you can never silently align
against the wrong annotation or overwrite one index with another. A chromap index
has no annotation to get wrong, and its layout says so: one `index/chromap/` per
assembly, with the record still saying exactly how it was built.
