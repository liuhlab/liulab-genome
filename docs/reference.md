# API reference

Generated from the source. Each entry carries the signature, the argument types and the
attributes exactly as the code declares them.

The written pages come first. They cover the same calls in the order you would make them:
[Assembly](genome/assembly.md), [Sequences and regions](genome/sequences.md),
[Annotations](genome/annotations.md) and [Aligner indexes](genome/aligner.md) for working
with a genome; [Transcription factors](topics/transcription-factors.md),
[Motifs](topics/motifs.md), [Gene identifiers](topics/gene-identifiers.md) and
[Homology](topics/homology.md) for the tables the package ships; and
[CLI overview](cli/index.md) for the shell. Come here once you know the name you want.

Each section below is one package's public surface: what its `__init__.py` re-exports, and
nothing else. A name you can only reach through a submodule path is internal, and it can
move or disappear between releases.

| Package | What it holds |
|---|---|
| `genome` | `Genome`, `Region`, the sequence types, and the results the common calls return |
| `genome.aligner` | The STAR and chromap index builders |
| `genome.annotation` | Registering a GTF and asking an annotation for gene ids |
| `genome.assembly` | The assembly table, the files on disk, and chimera naming |
| `genome.homology` | Ensembl Compara ortholog sets |
| `genome.tf` | The shipped table linking transcription factors to motifs |
| `genome.tf.motif` | JASPAR matrices, scanning, and reading hit files |
| `genome.xref` | Identifier conversion and symbol matching |

::: genome
    options:
      show_root_full_path: true

::: genome.aligner
    options:
      show_root_full_path: true

::: genome.annotation
    options:
      show_root_full_path: true

::: genome.assembly
    options:
      show_root_full_path: true

::: genome.homology
    options:
      show_root_full_path: true

::: genome.tf
    options:
      show_root_full_path: true

::: genome.tf.motif
    options:
      show_root_full_path: true

::: genome.xref
    options:
      show_root_full_path: true
