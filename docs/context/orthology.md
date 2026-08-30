# Orthology

Which genes in another species a gene is homologous to, on one publisher's own gene trees. This
context covers `homology/*`: the **Homology set** type and the reader behind it. Like an **Xref
set** it is anchored to a species pair and a pinned **Release** rather than to an **Assembly**, is
downloaded once into the **Data dir** as a sibling of the assembly tree, and answers with no
**Genome** open — and, like one, it reads bulk files fetched once and never a remote API.

**Orthology is served and never consumed** (ADR-0019). This context answers a user's cross-species
question and nothing else: no table this package publishes — no **TF gene table**, no **Cofactor
table**, no list of its own — is derived through homology, and no answer is ever silently
species-mapped. A species with no census or no shipped table raises and names the ones that have
them, rather than being answered with a translated guess.

Every term below is built. `genome.homology` answers for all three pairings among human, mouse and
worm off Ensembl Compara's protein gene-tree dumps at release 116 — the per-species dump that holds
a pair is fetched once, verified against the publisher's own md5, sliced to the pair and read back
locally — and `resolve_homologs` puts an answer into one registered **Annotation**'s own gene ids
through the call the Xref half already uses.

Words every context shares — **Assembly**, **Genome**, **Data dir**, **Completion marker** and the
rest — are defined once in the repo-root `CONTEXT-MAP.md`, **Gene id stem** among them. **TF gene
table**, **Cofactor table** and **Cross-species link** are defined in the [TF](./tf.md) glossary;
**Annotation** belongs to [Annotation](./annotation.md), **Release** to [Motif](./motif.md), and
**Xref set** and **Xref source** to [Xref](./xref.md).

## Language

**Homology set**:
One Ensembl Compara **Release** for one species pair, sliced out of the publisher's per-species dump
and read locally. It answers one question — which genes of the other species a **Gene id stem**'s
gene is homologous to — for the three species the lab works on, and raises for any other rather than
inventing one. Compara's per-species files are a de-duplicated partition *at the pair level*, so
either file of a pair may be the one holding it and which is not promised stable across releases: a
fetch finding zero rows raises naming the other file rather than answering empty, and zero is a
trustworthy signal because a pair is never partially present.
_Avoid_: ortholog table, ortholog map, ortholog set, homology database; Compara on its own, which
names the publisher's whole resource rather than this slice of one release of it; and "the
orthologs", which drops the paralogy the same set carries

**Homology link**:
One row of a **Homology set**: two **Gene id stem**s in two species, the **Homology type** the
publisher assigned, and the confidence fields Compara publishes beside it — its high-confidence
flag, `goc_score` and `wga_coverage`. Those are carried through unchanged so that a caller can
filter on them, and they are null for every human↔worm row, which the answer says out loud rather
than leaving a quality filter to empty itself silently.
_Avoid_: ortholog pair, orthology (the relation, not the row this package holds), mapping, match,
hit; and **Cross-species link**, which is the TF context's flag on a JASPAR profile and asserts no
orthology whatever (ADR-0013)

**Homology type**:
The publisher's own label on a **Homology link**, derived from its gene tree rather than counted off
rows: `ortholog_one2one`, `ortholog_one2many` and `ortholog_many2many` for speciation, and Compara's
paralogy labels for duplication. It is never recomputed — not after a filter, not after resolution
into an **Annotation**, not after a caller slices the answer (ADR-0020) — so it stays a claim about
evolution rather than a count of what happens to be in front of you, and a caller who needs a
one-to-one relationship requires this label instead of inferring one from what came back.
_Avoid_: cardinality on its own, which is what counting would give and this is not counted;
relationship, relation, orthology type (the paralogy labels live in the same field); one-to-one
flag; confidence and quality, which are the separate fields a **Homology link** already carries

**Dropped partner**:
A homologous gene an answer no longer names, counted and reported on that answer rather than quietly
discarded — removed by a **Homology type** filter, or lost when the stems were resolved into a
registered **Annotation** that does not spell them. It is what keeps a link that merely *looks*
one-to-one in your view distinguishable from one the publisher called one-to-one, the same
distinction the unresolved bucket already draws for **Gene id stem**s.
_Avoid_: missing, filtered, excluded — each says something went away without saying it is still
counted; unresolved, which names ids that resolved to nothing rather than partners a filter removed;
loss, attrition

**Paralogy link**:
A **Homology link** whose **Homology type** is one of the publisher's duplication labels rather than
a speciation one — the two genes descend from a gene duplication. Kept and marked, never excluded,
so that *not an ortholog* stays distinguishable from *absent* — the stance ADR-0013 already takes
for a **Cross-species link** — and returned only when asked for, so that the common question stays
the easy one.
_Avoid_: paralog on its own, which names a gene where this names the row relating two; homolog and
homology as loose synonyms for ortholog; in-paralog, out-paralog and co-ortholog, a vocabulary this
package neither ships nor computes
