# liulab-genome

A library for handling reference genomes, including genome assembly, annotation, and aligner index and other artifacts.

## Where to go next

- [**Usage**](usage.md) — installation, environment setup, and the CLI.
- [**Genome**](genome.md) — the `Genome` class, the main entry point: name an
  assembly, query its sequence.
- [**Sequences**](sequences.md) — the typed `DNA` / `RNA` / `Protein` classes,
  with worked examples.
- [**Genome files**](genome-files.md) — the download/prepare machinery behind
  `Genome` (`UCSCGenomeDownloader`, `prepare_fasta`).
- [**API reference**](reference.md) — auto-generated from docstrings (secondary
 to the hand-authored pages above).

## At a glance

The `Genome` class is the main entry point — name an assembly and query it in
**0-based, half-open** coordinates:

```python
from genome import Genome

sacCer3 = Genome("sacCer3")              # download + prepare on first use (cached)
sacCer3.fetch_sequence("chrIV:0-10")     # DNA('ACACCACACC')
sacCer3["chrIV:0-10"].reverse_complement()  # indexing is sugar; result is a DNA
```

The typed sequence classes it returns stand alone too:

```python
from genome import DNA

s = DNA("aTcG")
s.reverse_complement()        # DNA('CgAt') — case preserved
s.transcribe().gc_content     # 0.5
s[1:3]                        # DNA('Tc') — slicing stays typed
```

```bash
$ genome revcomp ATCG
CGAT

$ genome doctor
samtools: samtools 1.21 ...
bedtools: bedtools v2.31.1
```

## Project conventions

This project follows a strict set of domain invariants (0-based half-open
intervals, explicit reference assembly, normalized chromosome names, streaming
I/O, metadata-as-first-class). They are documented in
[`CLAUDE.md`](https://github.com/lhqing/liulab-genome/blob/main/CLAUDE.md);
read it before contributing.
