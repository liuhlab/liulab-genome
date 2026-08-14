# liulab-genome

Reference genomes on disk, ready to query. Name an assembly and `genome` fetches it,
prepares the companion files every tool expects (`.fai`, `.2bit`, `chrom.sizes`), and
answers sequence queries by region. It also registers GTF annotations against an
assembly and builds STAR and chromap indexes.

Import name: `genome`.

## Install

The package drives native tools from bioconda (`samtools`, `faToTwoBit`, `twoBitInfo`),
so [pixi](https://pixi.sh) is the supported path:

```bash
git clone https://github.com/liuhlab/liulab-genome.git
cd liulab-genome
pixi install --locked
pixi shell
```

`pip install liulab-genome` installs the Python API alone — you supply the native tools
yourself, and `gffutils` too if you need annotations.

Check the toolchain at any time:

```bash
$ genome doctor
samtools: samtools 1.21 ...
faToTwoBit: installed; reports no version
twoBitInfo: installed; reports no version
```

## Quickstart

Coordinates are **0-based, half-open** everywhere: `chrIV:0-10` is the first ten bases.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")                   # fetch + prepare on first use, cached after
sacCer3.fetch_sequence("chrIV:0-10")          # DNA('ACACCACACC')
sacCer3["chrIV:0-10"].reverse_complement()    # indexing is sugar; the result is a DNA
sacCer3.chrom_sizes["chrIV"]                  # 1531933
```

The same from a shell:

```bash
$ genome register sacCer3
$ genome revcomp ATCG
CGAT
```

## Where to go next

- [**Genome**](genome.md) — the main guide: prepare an assembly, fetch sequence,
  register annotations, build aligner indexes.
- [**Sequences**](sequences.md) — the typed `DNA` / `RNA` / `Protein` classes.
- [**CLI**](cli.md) — every command, with its options and exit codes.
- [**API reference**](reference.md) — generated from docstrings.
