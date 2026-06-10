# Annotations & aligner indexes

A reference assembly is sequence; most analyses also need a **gene annotation**
(a GTF) and an **aligner index** built from both. `Genome` ties these together:
you register one or more GTF annotations on the genome, then build an aligner
index against a chosen annotation.

```python
from genome import Genome

g = Genome("sacCer3")
g.register_gtf("sacCer3.ensGene.gtf", name="ensembl")   # place + index the GTF
g.build_star_index(gtf="ensembl", threads=8)            # STAR index for that annotation
```

## Registering a GTF annotation

A genome can carry several annotations (GENCODE, Ensembl, RefSeq, …), each under
a unique `name`. `register_gtf` copies the GTF into the assembly's data
directory and builds a [gffutils](https://gffutils.readthedocs.io/) database
beside it:

```python
g.register_gtf("path/to/annotation.gtf", name="ensembl")
g.annotations                      # ['ensembl'] — registered names
g.get_gtf_path("ensembl")          # Path to the placed .gtf
```

Each annotation lives in its own directory next to the sequence files:

```
<LIULAB_DATA>/genome/<assembly>/gtf/<name>/
    <name>.gtf      # the (unzipped) annotation, copied in
    <name>.db       # the gffutils SQLite database built from it
```

A few things to know:

- **The GTF must be unzipped.** A `.gz` path is rejected with an actionable
  error — decompress it first (`gunzip`). (Standard annotations are small enough
  uncompressed that this is rarely a burden.)
- **Re-registering the same name raises** `FileExistsError` unless you pass
  `force=True` — annotations are not silently overwritten.
- **Gene/transcript inference is off by default** (`disable_infer_genes` /
  `disable_infer_transcripts` are `True`). Standard annotation GTFs already
  declare `gene`/`transcript` features, and inferring them is the classic
  gffutils slow path. Enable it (`disable_infer_genes=False`) only for a bare
  exon-level GTF that lacks those records.

### The default annotation

When exactly one annotation is registered it becomes the default; with several,
set it explicitly at construction:

```python
g.default_gtf                      # 'ensembl' when it is the only one, else None
g.default_gtf_path                 # its GTF path, or None

Genome("hg38", default_gtf="gencode")   # pick the default up front
```

`default_gtf` is the annotation other features fall back to when you don't name
one; it does **not** change what `build_star_index(gtf=...)` requires — that
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
  pass `overwrite=True` to force a rebuild.
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

After a successful build the directory holds STAR's binary files plus two
bookkeeping artifacts:

- `.success` — a marker written only when the build completes. A half-built
  index (interrupted run) has no marker and is rebuilt rather than trusted.
- `star.index.json` — a sidecar recording the STAR version, the assembly, the
  resolved parameters (including the `gtf` key and the GTF path), and the exact
  command that was run, so an index is self-describing.

## STAR is an optional dependency

STAR is **not** in the default environment — it is only needed when you actually
build an index. Install it into the project's aligner environment:

```bash
pixi add star            # from bioconda
```

If STAR is missing, the build fails fast with install instructions and a
`genome.external.ToolNotFoundError`, rather than a cryptic error deep in the
run.

## Domain invariant: the annotation is always explicit

Mirroring the assembly rule, the annotation an index is built against is never
implicit. `build_star_index` requires a named `gtf`, the index directory encodes
that name, and the metadata sidecar records it — so you can never silently align
against the wrong annotation or overwrite one index with another.
