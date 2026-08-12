# Annotation

What a set of gene models is once an **Assembly** owns it. This context covers `io/gtf.py` and the
registry a **Genome** keeps over it: a GTF is filed under one assembly, given a name, and addressed
by that name forever after.

Words every context shares — **Assembly**, **Genome**, **Chromosome**, **0-based half-open** and the
rest — are defined once in the repo-root `CONTEXT-MAP.md`.

## Language

### Naming

**Annotation**:
A named set of gene models — genes, transcripts, exons — belonging to exactly one **Assembly** and
filed beneath it at `<assembly dir>/gtf/<name>/`. Registering a **GTF** whose **Chromosome** names
disagree with the assembly's is an error and not a warning: the two are one unit or they are nothing.
_Avoid_: gene set, gene model, features, GFF; and never as a synonym for the **GTF**, which is one of
the annotation's two files

**Registered name**:
The short key an annotation is addressed by everywhere — `"gencode_v44"` — unique within its
**Assembly**, and simultaneously the name of its directory and of both files inside. Callers pass
names, never paths; a path in an argument means the annotation was never registered.
_Avoid_: id, key, label, alias; and never "the GTF" as a way of naming one

**Default annotation**:
The annotation used when a caller names none — set explicitly, or adopted on its own when the
assembly has exactly one. Two annotations and no explicit choice leaves no default at all, because a
caller who did not choose should be asked rather than guessed at.
_Avoid_: primary, main, active, current

### Files

**GTF**:
An annotation's source text, one feature per line, kept decompressed inside the annotation's
directory. Its coordinates are 1-based and inclusive, so they convert to **0-based half-open** at the
I/O boundary and are never seen in that form anywhere else.
_Avoid_: GFF (a different format, not a different spelling), annotation file, gene file

**Annotation database**:
The gffutils SQLite file built from the **GTF** and kept beside it as `<name>.db` — what makes an
annotation queryable rather than merely stored. Its presence is what *registered* means: a directory
holding a GTF and no database is a failed registration, not an annotation.
_Avoid_: index, cache, store; and never "the gffutils db" in the API surface

**Feature inference**:
gffutils reconstructing `gene` and `transcript` rows from exon lines in a **GTF** that declares
neither. Off by default and deliberately: GENCODE, Ensembl and RefSeq all declare those features
already, and inferring them is gffutils' slow path — turn it on only for a bare exon-level GTF.
_Avoid_: gene inference, building the hierarchy, auto-detect
