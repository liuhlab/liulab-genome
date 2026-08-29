# TF gene

Which genes a published census judges transcription factors, and which JASPAR motifs answer for
them. This context covers `tf/gene/*` — the shipped censuses, the loader that reads them, and the
**TF gene list** an **Annotation** answers with — and `tf/link.py` beside it, which holds the
**Motif link** join. It is the gene half of the TF context, keyed by gene where the **Motif** half
is keyed by motif; the join between the two is many-to-many, so `tf/link.py` imports both halves and
neither half imports it.

Nothing here decides what a transcription factor is. Every verdict travels with the census that
reached it — Lambert et al. 2018 for human, AnimalTFDB 4.0 for mouse — so two censuses that classify
one factor differently are two answers rather than a contradiction to resolve, and this package's
job is to say which one is speaking.

Nothing in this half is built at the time of writing. The words are settled first, as the motif
half's were, so that every issue, test name, error message and docstring on the way uses one set —
which is why no entry below carries the *(decided, not built)* marker the context map describes:
every one of them would. Use the words; do not call the API they describe until its module lands.

Words every context shares — **Assembly**, **Genome**, **Data dir** and the rest — are defined once
in the repo-root `CONTEXT-MAP.md`. **Annotation**, **Registered name**, **Annotation database** and
**Curated gene list** belong to the [Annotation](./annotation.md) context; **Motif**, **Motif id**,
**Motif name**, **Release**, **Tax group** and **Information content** to the [Motif](./motif.md)
half of this one.

## Language

### What a census says

**TF gene**:
A gene one published census judges a transcription factor, named in one **Annotation**'s own gene
ids. The verdict is the census's and never this package's, and it is graded rather than binary — see
**TF assessment**. Absence from a census and rejection by one are different answers and stay
distinguishable everywhere: Lambert assessed 2,765 genes and judged 1,639 of them TFs, so a human
gene can be a TF, judged not to be, or never looked at, and only the first two are the census
speaking.
_Avoid_: TF and factor on their own — the **Motif** half already bans both as names for a motif, and
here they blur a gene with the protein it encodes; transcription factor as a bare noun in an API
surface, which does not say *whose* judgement; regulator, DNA-binding protein, TF candidate

**TF assessment**:
The census's own graded verdict on one gene, in the publisher's own spelling — Lambert's spans
`Known motif`, `Likely to be sequence specific TF`, `Inferred motif` and `Unlikely to be sequence
specific TF` among others. It is the column a caller tightens or loosens on: a **TF gene list**
carries the assessed-positive genes, and wanting only `Known motif`, or wanting `Inferred motif`
included, is a re-filter on this rather than a second flag invented here. One census's word, never
compared with another's.
_Avoid_: confidence, quality, evidence level, score — a graded verdict is not a number and nothing
ranks two censuses' grades against each other; status, tier, class; and "is a TF" as though the
answer were one bit

**DBD family**:
The DNA-binding-domain family a census classifies a gene under, in the publisher's own vocabulary:
75 values under Lambert's `DBD`, 72 under AnimalTFDB's `Family`. Uniform across **TF gene table**s
in position and deliberately not in content — the two vocabularies are not crosswalked, so
`ARID/BRIGHT` and `ARID` are not asserted equivalent (ADR-0014). Group by it within a species; never
across two.
_Avoid_: family alone — a **Motif** carries a `family` annotation of JASPAR's own, a different
vocabulary about a different object; class, superclass, TF class, domain; DBD by itself, which names
the domain rather than the grouping

### What ships, and what it answers

**TF gene table**:
One census as shipped data: external published knowledge with per-gene attributes, filed in its own
subdirectory of the package's data and found by enumeration and a filename convention, so adding a
species is dropping in a file. Distinct from a **Curated gene list**, which is a membership set
derived from one **Annotation**'s own attribute (ADR-0011) — a census is keyed by **Gene id stem**,
belongs to a species rather than to an annotation, and beyond four uniform columns its columns are
its publisher's own. Its provenance — publisher, version, PubMed id, source URL, and a checksum over
unpacked content — sits in a metadata table beside it, in the shape the assembly and annotation
metadata tables already use.
_Avoid_: TF list, TF database, TF catalogue; **Curated gene list**, which is the other shape of
shipped file and answers a different question; annotation, which names a registered GTF here

**TF gene list**:
One **Assembly**'s **TF gene**s: the census's **Gene id stem**s resolved into the registered
**Annotation**'s own gene ids, carrying the stems that resolved to nothing so that what the census
holds and this annotation does not is visible rather than dropped. Assessed-positive by default. The
species comes from the assembly's own metadata and is never passed in, so asking for human TFs while
holding a mouse assembly is not expressible (ADR-0003); an assembly whose species has no census
raises and names the species that do.
_Avoid_: TF set, TF panel, TF universe; and the bare phrase "the gene list", which the Annotation
context bans for naming no annotation — the qualifier is exactly what makes this compound legal

**Gene id stem**:
A gene id with its version suffix dropped: `ENSG00000123456` for `ENSG00000123456.7`. A stem may
name more than one gene id inside one **Annotation** — nine do in `gencode_v50lift37`, eight of them
pseudoautosomal-Y — so resolution answers with every id a stem names and never picks one. A stem
carrying no version resolves to itself, which is why an annotation whose ids are not Ensembl-shaped
is unaffected. The mechanism belongs to the annotation registry rather than to this context's code,
and reading gene ids is the first thing that opens the **Annotation database** this package has
always built; TF genes are merely its first caller.
_Avoid_: base id — that already names a **Motif id** with its version dropped, `MA0139` for
`MA0139.2`, and one term must not name two things inside one package; unversioned id, stripped id,
gene id (the versioned thing this is a stem *of*), ENSG (one publisher's spelling)

### The join to motifs

**Motif link**:
One row saying a JASPAR **Motif** answers for a **TF gene**: the **Motif id** and **Motif name**,
the **Role** and the partners the profile also names, the profile's tax ids, its **Cross-species
link** flag, the matrix's total **Information content**, and its rank under **Attribution
specificity**. Links are shipped as gzipped TSVs, one per species per **Release** — data and not a
rule (ADR-0015) — so the mapping is readable in R or a shell by collaborators who never import this
package. Only assessed-positive genes receive links, and a profile naming no gene at all, such as an
oncogenic fusion, stays unlinked rather than asserting one.
_Avoid_: motif assignment, TF-motif mapping (that names the whole table, not one row); annotation,
which is cisTarget's word for this idea and a different, largely inferred, thing; binding site,
target — a link is about attribution and claims nothing about where the factor binds

**Role**:
What a **Motif link**'s profile is a motif *of*: `monomer` where the profile names one gene,
`complex` otherwise, with the other named partners recorded on the row. It exists so that a
heterodimer matrix is never read as a monomer's — `FOS::JUN` links to both genes as a complex and to
neither as a monomer — and so that a gene whose only motifs are complexes, AHR, DDIT3, TAL1 and TLX1
among them, is linked rather than reported motif-less.
_Avoid_: type, kind, link type; homodimer, heterodimer, dimer — a complex may name more than two
genes and the row says which; direct/indirect, which grades evidence rather than saying what the
matrix describes

**Attribution specificity**:
The order one gene's **Motif link**s come back in: **Role** `monomer` before `complex`, then
species-matched before **Cross-species link**, then higher total **Information content**, then
**Motif id** — four keys, so the order is total and stable and "the motif for this factor" means the
same thing on two machines and in two **Release**s. It states what a matrix is attributable to and
explicitly not which motif is better: the canonical AP-1 matrix is a complex and describes JUN's
binding better than any JUN monomer does. No quality score is computed or shipped anywhere here.
JASPAR publishes none, and matrix depth is normalised per assay — SMiLE-seq sits near 1,000 sites
throughout a range that runs from 10 to 322,803 — so ranking on depth ranks the assay. A caller who
disagrees re-sorts on the attributes the row already carries.
_Avoid_: quality, confidence, score, weight; best motif, primary motif, canonical motif — each names
a judgement nothing here makes; rank alone, which is the ordering's output rather than its meaning

**Cross-species link**:
A **Motif link** whose profile was measured on a vertebrate other than the gene's own species — the
mouse `Yy1` matrix answering for human `YY1`. Kept and marked, never excluded (ADR-0013): JASPAR's
CORE **Tax group** `vertebrates` files an orthologous pair's matrix under whichever species was
assayed, so a profile's species is an artefact of the experiment rather than a claim about the
factor. Marked on every row because the flag is the only thing a caller who needs species-matched
profiles can filter on — 732 of mouse's 896 links on the 2026 **Release** are cross-species.
_Avoid_: ortholog link, orthologous motif — this package holds no ortholog table and asserts no
orthology, only that JASPAR filed one vertebrate's assay under one vertebrate's name; foreign,
non-native, imputed, borrowed
