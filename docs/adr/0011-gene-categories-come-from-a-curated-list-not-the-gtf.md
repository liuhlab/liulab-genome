# 11. A gene category comes from a curated list shipped here, not from the GTF's biotype attribute

An annotation can now be asked which of its genes are in a named category, and the answer comes from
a hand-maintained JSON shipped beside the metadata tables, one per annotation. The obvious reading —
read `gene_type` or `gene_biotype` off the registered GTF — lost on three counts measured across the
annotations the lab registers, each fatal alone: the attribute is spelled two ways across four
publishers; the taxonomies behind those spellings disagree, GENCODE splitting `rRNA`,
`rRNA_pseudogene` and `Mt_rRNA` where WormBase and RefSeq do not, with worm's `rrn-3.56` sitting
inside the chrI 45S repeat and labelled a pseudogene — a four-point swing in one plate's measured
rRNA share, which is the caller's judgment and not ours to bake in; and `sacCer3`'s `ensgene_v101`
carries no biotype attribute at all, so a caller deriving categories reports 0% rRNA for yeast and
never finds out. Curating the lists puts that judgment in one place and makes yeast answerable at
all. The cost is unmitigated: a curated list is pinned to one release of one annotation and goes
stale silently when that annotation is re-released, and the suite proves a file well-formed and
internally consistent, never that it still matches the GTF.
