# Aligner indexes

A read mapper does not search a FASTA. It loads an **index**, a set of binary files derived
from the reference once and reused by every mapping run. Name an aligner and its index is
built for the assembly you have open. The call returns the path it wrote.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.build_chromap_index()
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/index/chromap/chromap.index')
```

Mapping the reads is not part of this package. Registering the annotation a STAR index is
built against is on [Annotations](annotations.md).

## Which aligner

| Aligner | What it maps |
| --- | --- |
| [STAR](https://github.com/alexdobin/STAR) | RNA-seq. It is splice-aware, so a read spanning an exon junction still aligns. |
| [chromap](https://github.com/haowenz/chromap) | ATAC-seq, ChIP-seq and Hi-C, where a read sits on one contiguous stretch. |

STAR learns the junctions from a gene annotation, so a STAR index is built against exactly
one GTF. chromap reads the FASTA and nothing else, so an assembly needs exactly one chromap
index no matter what you map against it.

## Where an index lands

Indexes go under `index/` inside the assembly's own directory. STAR's is
`index/star_<gtf>/`, a directory that is itself what STAR loads. chromap's is the single
file `index/chromap/chromap.index`.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

The annotation name is part of the STAR directory name, so indexes built against different
annotations sit side by side and never collide. Both kinds are large: 12 Mb of yeast gives
a 123 MB STAR directory and a 70 MB chromap index, and a vertebrate assembly runs to tens
of gigabytes.

## Building an index

Two calls. STAR takes the annotation to build against, chromap takes nothing:

```python
sacCer3.build_star_index("ensgene_v101", threads=4)
sacCer3.build_chromap_index()
```

Omit the annotation and the genome's default is used, which is the everyday call:

```python
sacCer3.default_gtf                          # 'ensgene_v101'
sacCer3.build_star_index(threads=4).name     # 'star_ensgene_v101'
```

A genome with no default annotation raises `ValueError` rather than guessing, and the
message names both ways to supply one.

The tool's own output streams to the terminal while it works. Forcing a rebuild of the
index above, trimmed here:

```python
sacCer3.build_chromap_index(overwrite=True)
# Build index for the reference.
# Kmer length: 17, window size: 7
# Built index successfully in 0.50s.
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/index/chromap/chromap.index')
```

Yeast takes a second or two either way. A vertebrate assembly takes an hour or more, which
is why the output is left streaming rather than captured.

There are no CLI commands for indexes. They are built through the API.

## Tuning the build

A handful of options are named arguments. Anything else is forwarded to the tool as a raw
flag, spelled as the tool spells it minus the leading `--`:

```python
sacCer3.build_star_index(
    gtf="ensgene_v101",
    sjdb_overhang=99,        # --sjdbOverhang, ideally read_length - 1
    threads=4,               # --runThreadN
    genomeSAindexNbases=11,  # any other genomeGenerate option
)
```

chromap's two knobs are the minimizer k-mer length and window size, and its build is
single-threaded, so there is no `threads`. Underscores in a forwarded keyword become
hyphens:

```python
sacCer3.build_chromap_index(
    kmer=20,                 # -k/--kmer
    window=10,               # -w/--window
    min_frag_length=30,      # --min-frag-length
)
```

`genomeSAindexNbases` and `genomeChrBinNbits` are computed from the assembly unless you
pass them, so a small genome and a many-sequence reference both skip the hand adjustment
STAR's manual asks for. The sacCer3 build above ran with `--genomeSAindexNbases 10
--genomeChrBinNbits 18`.

## Reusing and rebuilding

A finished index is reused. Calling `build_chromap_index()` a second time returns the same
path without running chromap, so reading an index back does not need the binary installed
at all. Pass `overwrite=True` to rebuild.

What counts as finished is a record written after the tool exited cleanly, not the presence
of index files. A build that was interrupted leaves a directory nothing vouches for, and
opening it raises instead of being reused. **Re-registering the assembly invalidates every
index under it**, because the sequences and chromosome names the index was built from are
no longer the ones a mapping run would use. Each message names the call that rebuilds it.

## Handing an index to a mapper

Building an index and using one are separate jobs. The `get_*` methods return a path and
build nothing, which is what a mapping run wants. STAR's is the directory `--genomeDir`
expects. chromap's is the file `-x/--index` expects, and takes no argument, since one index
covers the assembly:

```python
sacCer3.get_star_index("ensgene_v101")
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/index/star_ensgene_v101')

sacCer3.get_chromap_index()
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/index/chromap/chromap.index')
```

Unlike the build call, `get_star_index` has no default and the annotation has to be named.
Asking for an index nobody built raises `IndexNotBuiltError`:

```python
Genome("ce11").get_star_index("ensgene_v101")
# IndexNotBuiltError: no star index for 'ce11' at
# /Users/hanqing/liulab_data/genome/ce11/index/star_ensgene_v101: nothing has been
# built there yet. Build it with `Genome.build_star_index(gtf='ensgene_v101')`.
```

`get_index("star", gtf="ensgene_v101")` is the general form both wrap, and the `STAR` and
`Chromap` classes underneath them are in the [API reference](../reference.md).

## Installing the aligners

Neither aligner is in the default environment. Add the one you need with `pixi add star` or
`pixi add chromap`, or run against the project's `aligners` environment, which carries both:

```bash
pixi run -e aligners python build_indexes.py
```

A missing binary raises at the point the build would have started rather than when the
genome is opened, and the message names the command that installs it.
