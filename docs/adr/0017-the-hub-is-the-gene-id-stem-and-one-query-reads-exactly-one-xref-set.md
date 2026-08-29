# 17. The hub is the Gene id stem, and one query reads exactly one Xref set

Every **Namespace** an **Xref set** carries is a spoke off one hub, the **Gene id stem**, and the
public surface is two verbs: to the hub and from the hub. There is no foreign-to-foreign conversion,
so a caller wanting Entrez → HGNC makes two calls and owns the join — the hop stays visible in their
code instead of invisible in ours — and a query names exactly one **Xref source**, so the source is
a property of the whole answer rather than a column on every row and merging two publishers is not
expressible rather than merely discouraged. The convenient reading, one table over every source with
the shortest path found for you, lost on what it would have to invent: NCBI and Ensembl agree on
only 57.5% of human gene-level (GeneID, ENSG) pairs, so a merged table decides nearly half its rows
by a rule nobody published, and a composed answer inherits the fan-out of both hops with nothing on
it saying so — Ensembl alone reaches 72 stems for one GeneID. The cost is real and is paid by the
caller: two calls where one would read better, no union, vote or "best" mapping anywhere, and a
species with no Ensembl presence — the registered *E. coli* HT115 assembly — met with an error.
