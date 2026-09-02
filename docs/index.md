# Get started

`liulab-genome` keeps reference genomes on disk so a script does not have to. Name an
assembly and the package fetches its FASTA, derives the companion files other tools expect,
and hands back an object that reads sequence out of it. It also registers GTF annotations
against an assembly, builds aligner indexes, and carries published tables of transcription
factors, motifs, gene identifiers and orthologs. The import name is `genome`.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.fetch_sequence("chrIV:0-10")    # DNA('ACACCACACC')
sacCer3.chrom_sizes["chrIV"]            # 1531933
```

Coordinates in the Python API are 0-based and half-open, so `chrIV:0-10` is the first ten
bases. Every class and function is listed in the [API reference](reference.md).

## Install

The work is done by native binaries, so [pixi](https://pixi.sh) is the supported path. One
lock file brings the Python package, `samtools`, `bedtools`, `faToTwoBit`, `twoBitInfo`,
`gffutils` and the `moods` motif scanner together:

```bash
git clone https://github.com/liuhlab/liulab-genome.git
cd liulab-genome
pixi install --locked
pixi shell
```

STAR and chromap are not in that environment. They are large and most work never touches
them, so they live in a second pixi environment named `aligners`. Reach it with
`pixi run -e aligners ...` or `pixi shell -e aligners`.

`pip install liulab-genome` installs the Python package and its Python dependencies alone.
It brings no native binary, so anything that shells out fails until you install `samtools`
and the rest yourself, and it brings no `gffutils`, so annotations do not work either.

`genome doctor` reports what it found on `PATH`, and exits non-zero when a required tool is
missing:

```console
$ genome doctor
samtools: samtools 1.22.1
faToTwoBit: installed; reports no version
twoBitInfo: installed; reports no version
```

## The data directory

Everything the package downloads lands in one directory per machine. Every project reads
the same copy, so nothing is fetched twice and a prepared assembly reopens instantly, offline.

`LIULAB_DATA` names the root. Unset, the well-known lab paths are tried in order, then
`~/liulab_data`. The roots you will meet under it:

```text
<LIULAB_DATA>/
├── genome/<assembly>/     FASTA, .fai, .2bit, chrom.sizes
│   ├── gtf/<name>/        annotation, and the gffutils database built from it
│   └── index/<name>/      STAR / chromap index
├── motif/                 JASPAR releases and score thresholds
├── xref/                  identifier tables
└── homology/              ortholog tables
```

Ask `assembly_data_dir` where an assembly landed. `genome assembly register` prints it too.

```python
from genome.assembly import assembly_data_dir

assembly_data_dir("sacCer3")
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3')
```

Five things download on first use: an assembly, an annotation, a JASPAR release, an xref set
and a homology set. **The lab's CPU compute nodes have no internet**, so the first use of
each has to happen on a login node. Every job after that reads the prepared copy.

`Genome(..., cache_dir=...)` puts one assembly's directory somewhere else, and its
annotations and indexes go with it, for when the shared root is full or slow.

## A walkthrough

Yeast is the quickest thing to try this on. Its FASTA is 3.8 MB compressed, and the whole
walkthrough below finishes in about a minute.

`genome assembly list` is the first thing to run. It names the eight assemblies the
shipped metadata table lists, `sacCer3` among them, and says which are already prepared on
this machine — including any this table does not list. It downloads nothing:

```console
$ genome assembly list
assemblies in /Users/hanqing/liulab_data/genome
  hg38          offered, not registered  Homo sapiens GRCh38
  hg19          offered, not registered  Homo sapiens GRCh37
  mm39          offered, not registered  Mus musculus GRCm39
  mm10          offered, not registered  Mus musculus GRCm38
  sacCer3       registered               Saccharomyces cerevisiae R64-1-1
  ce11          registered               Caenorhabditis elegans WBcel235
  ecHT115       offered, not registered  Escherichia coli HT115 ASM435494v1
  ce11_ecHT115  offered, not registered
registered here: 2 — prepare another with `genome assembly register <name>`, or re-check one with `genome assembly verify <name>`
an assembly the table does not list registers too, from the UCSC golden path — with no pinned checksum behind it
```

A listed assembly is pinned to a source URL and a checksum. Any other UCSC name registers
too, without a pin, which is how you prepare a non-model organism;
[which assemblies you can name](genome/assembly.md#which-assemblies-you-can-name) has the
detail.

Preparing an assembly downloads the FASTA, checks it against a pinned checksum, and derives
the three companion files. Do it from a shell before a pipeline starts:

```console
$ genome assembly register sacCer3
registered sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

Running it a second time prints the same thing and downloads nothing.

`Genome` opens what that produced. Constructing one prepares the assembly too if it is not
prepared already, so the shell step is a convenience rather than a prerequisite:

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.chromosomes[:4]    # ['chrI', 'chrII', 'chrIII', 'chrIV']
```

Ask for bases by region. `fetch_sequence` takes a `chrom:start-end` string and returns a
`DNA`, which carries the sequence transforms:

```python
sacCer3.fetch_sequence("chrIV:1000-1030")
# DNA('TCTATAGTCATACAGACGCTTTTACTTCAC')
sacCer3["chrIV:1000-1030"].reverse_complement()
# DNA('GTGAAGTAAAAGCGTCTGTATGACTATAGA')
```

An annotation is a GTF registered against the assembly under a short name. `genome
annotation list` gives the annotation names for an assembly, says which are built here, and
marks the default:

```console
$ genome annotation list sacCer3
annotations for sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3
  ensgene_v101  registered  UCSC ensGene.v101
default: ensgene_v101
```

On a fresh machine that row reads `offered, not registered` and the last line names the
command to run; [annotations](genome/annotations.md#what-is-registered-offered-or-broken)
covers the other states. `register` takes the assembly and the name, both required, so
registering the default means typing it out. It downloads the GTF and builds a gffutils
database beside the genome files:

```console
$ genome annotation register sacCer3 ensgene_v101
registered ensgene_v101 for sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3/gtf/ensgene_v101
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
  sha256  d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026
  files   ensgene_v101.db, ensgene_v101.gtf
  chromosomes checked — every name the GTF uses is one the assembly carries
```

A genome can hold several annotations. `registered` lists the ones built here, and
`default_gtf` says which one a call uses when you do not name one:

```python
sacCer3.annotations.registered    # ['ensgene_v101']
sacCer3.default_gtf               # 'ensgene_v101'
```

Aligner indexes are built from Python and need the `aligners` environment. chromap uses no
annotation, so one index serves the whole assembly:

```python
sacCer3.build_chromap_index()
# PosixPath('/Users/hanqing/liulab_data/genome/sacCer3/index/chromap/chromap.index')
```

chromap prints its own progress while it works, and yeast takes about a second. An index
that already finished is reused, so calling this again returns the path and runs nothing.
`build_star_index()` is the STAR equivalent. It builds against one annotation, defaults to
`default_gtf`, and takes `threads=`, which starts to matter once the genome is mammalian.

## Where to go next

| I want to… | Page |
| --- | --- |
| Prepare a reference genome, or concatenate two of them | [Assembly](genome/assembly.md) |
| Read bases out of a region, or work with coordinates | [Sequences and regions](genome/sequences.md) |
| Register a GTF and ask what genes it carries | [Annotations](genome/annotations.md) |
| Build a STAR or chromap index and find it again | [Aligner indexes](genome/aligner.md) |
| List a species' transcription factors and cofactors | [Transcription factors](topics/transcription-factors.md) |
| Scan sequence for JASPAR motifs and read the hits | [Motifs](topics/motifs.md) |
| Convert gene ids between databases, or match symbols | [Gene identifiers](topics/gene-identifiers.md) |
| Find the mouse ortholog of a human gene | [Homology](topics/homology.md) |
| Do any of this from a shell script | [CLI overview](cli/index.md) |
| Look up a signature, an argument name or an attribute | [API reference](reference.md) |
