# 20. A Homology type is the publisher's tree-derived label and is never recomputed

Every **Homology link** carries Compara's own `homology_type` — `ortholog_one2one`,
`ortholog_one2many`, `ortholog_many2many` and its paralogy labels — read from the file and never
recalculated: not after a filter, not after resolution into an **Annotation**, not after a caller
slices the answer. Deriving the cardinality instead, by grouping whatever rows are in hand, is the
reading a table of pairs invites and it lost on what the label would then mean. Compara assigns it
from the gene tree, so it survives slicing; a count over rows is an artefact of which per-species
file was downloaded, which quality filter ran and which **Annotation** the stems resolved against,
and it would change when the annotation changed while the evolution did not. The distinction is
where the disagreement lives: across Compara, NCBI and Alliance, 31.3% of asserted human↔mouse pairs
are claimed by exactly one source, yet the strict one-to-one subsets agree at Jaccard 0.83–0.91. The
cost is that an answer can look one-to-one and still be labelled `ortholog_one2many`, so the
**Dropped partner** count is reported beside it rather than the label being quietly corrected.
