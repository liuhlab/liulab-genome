# Annotation

What a set of gene models is once an **Assembly** owns it. This context covers `annotation/*` and
the registry a **Genome** keeps over it: a GTF is filed under one assembly, given a name, and
addressed by that name forever after.

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

**Annotation metadata**:
The curated TSV row keyed by **Assembly** plus **Registered name**, saying who publishes that
annotation, which release it is, where its **GTF** is fetched from and the sha256 of the *unpacked*
GTF that source yields, plus whether it is the assembly's **Default annotation**. Naming an
annotation is enough to register it because the row knows the rest. A cross-reference and never an
allow-list — no row means one of three things: the annotation was registered by path, a complete
record handed in at the call site replaced the row wholesale, or it was derived here rather than
fetched, which is a **Merged annotation** and leaves a row nothing to say.
_Avoid_: registry, catalog, manifest; and "the GTF table", which names one of the files rather than
the annotation

**Default annotation**:
The annotation used when a caller names none: an explicit choice, else the one the **Annotation
metadata** table flags for that **Assembly** — which is how everyone in the lab reaches for the same
release without discussing it — else the sole registered annotation, and otherwise none, because a
caller who did not choose between several should be asked rather than guessed at. It names an
annotation without locating one, so a default nobody has registered on this machine is the ordinary
state of a fresh install and only asking for its path is an error.
_Avoid_: primary, main, active, current

**Merged annotation**:
The **Annotation** a **Chimera** build derives rather than fetches — a streaming pass over each
component's own **Default annotation**, registered in the same act that writes the chimera's
**FASTA**. Its **Registered name** is the `+`-join of those names in sorted-component order,
`wormbase_ws298+refseq_rs_2025_06_26`, and it carries no **Annotation metadata** row at all, since
nothing was downloaded for one to describe.
_Avoid_: combined GTF, concatenated annotation, chimera GTF; and never "the merge", which names the
pass rather than the annotation it produced

### Files

**GTF**:
An annotation's source text, one feature per line, kept decompressed inside the annotation's
directory. Its coordinates are 1-based and inclusive, so they convert to **0-based half-open** at the
I/O boundary and are never seen in that form anywhere else.
_Avoid_: GFF (a different format, not a different spelling), annotation file, gene file

**Annotation database**:
The gffutils SQLite file built from the **GTF** and kept beside it as `<name>.db` — what makes an
annotation queryable rather than merely stored. Its presence proves nothing: a build killed half-way
leaves one that answers queries with most of the genes missing, so what *registered* means is a
**Completion marker** that agrees with what is on disk, and that record is the only thing ever asked.
_Avoid_: index, cache, store; and never "the gffutils db" in the API surface

**Feature inference**:
gffutils reconstructing `gene` and `transcript` rows from exon lines in a **GTF** that declares
neither. Off by default and deliberately: GENCODE, Ensembl and RefSeq all declare those features
already, and inferring them is gffutils' slow path — turn it on only for a bare exon-level GTF.
_Avoid_: gene inference, building the hierarchy, auto-detect

### Genes

**Gene category**:
One named group of genes inside one **Annotation** — `rRNA` is the only one today — declared by that
annotation's **Curated gene list** and never by this package: which categories exist differs per
annotation and is data. A category is drawn for counting reads, so it is inclusive and its genes may
overlap: `rRNA` holds everything rRNA-derived, pseudogene copies and mitochondrial genes among them.
One that is declared always holds at least one gene, so an annotation that cannot answer for a
category says so rather than answering with none.
_Avoid_: biotype, gene type, `gene_type`/`gene_biotype` — those name the **GTF** attribute a category
is curated *instead of*; class, group, gene set

**Curated gene list**:
The hand-maintained JSON shipped inside the package, one per **Annotation**, saying which of its
genes are in each **Gene category** and, per category, what membership means there and where it came
from. It is the source of truth rather than the **GTF**'s own biotype attribute, which four
publishers spell two ways over three disagreeing taxonomies and one omits altogether (ADR-0011); an
annotation no list ships for is unanswerable, which is not the same fact as one whose category is
empty.
_Avoid_: biotype table, gene type table, allow-list, whitelist; and never the bare "the gene list",
which names no annotation — the ban is on the bare phrase and not on a qualified compound, so the TF
context's **TF gene list**, which names one, is a legal term
