# liulab-genome

A small library for handling genomic files — metadata management, file processing,
and feature extraction. Import name: `genome`.

The package wraps mature native tools (`samtools`, `bedtools`) with a thin, typed
Python API and CLI. Sequence operations get first-class types
(`DNA` / `RNA` / `Protein`) that keep the right type across slicing and
biological transforms.

## Where to go next

- [**Usage**](usage.md) — installation, environment setup, and the CLI.
- [**Sequences**](sequences.md) — the typed `DNA` / `RNA` / `Protein` classes,
  with worked examples.
- [**API reference**](reference.md) — auto-generated from docstrings (secondary
  to the hand-authored pages above).

## At a glance

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
