# Xref

Which foreign identifiers name a gene, and which genes a foreign identifier names. This context
covers `xref/*`: the **Xref set** type, the reader behind each **Xref source**, and the curated
metadata table that says which sets exist and where their bytes come from. It answers those two
questions and no third — it holds no coordinates, opens no **Genome**, and belongs to no
**Assembly**, an identifier being a name and not a place, which is why its files sit beside the
motif subtree rather than inside the assembly tree.

Nothing here decides which identifiers name the same gene. Every pair travels with the **Xref
source** that asserted it, so two publishers who disagree are two answers rather than a
contradiction to resolve: NCBI and Ensembl agree on 57.6% of human gene-level (GeneID, ENSG) pairs,
which makes the choice of publisher nearly half the answer.

Every term below is built. `genome.xref` answers for human, mouse and worm off the Alliance of
Genome Resources, and for human and mouse off Ensembl's per-species TSV dumps as well — a second
source, selectable and pinned to its own numbered **Release**, never merged with the first. Two
further sources carry symbols, which those two do not: HGNC's quarterly archive for human, with
approved, previous and alias spellings each typed, and the Alliance's per-species gene submissions
for mouse and worm, with the current approved spelling alone.

Words every context shares — **Assembly**, **Genome**, **Data dir**, **Completion marker** and the
rest — are defined once in the repo-root `CONTEXT-MAP.md`, **Gene id stem**, the hub every entry
below hangs off, among them. **Annotation**, **Registered name** and **Annotation database** belong
to [Annotation](./annotation.md); **Release**, **Motif set** and **Motif name** to
[Motif](./motif.md), and a release here is the same word it is there — one dated, versioned
publication of the source a downloaded set is pinned to.

## Language

### What a set is, and whose it is

**Xref set**:
One species, one **Xref source**, one pinned **Release**, opened the way a **Motif set** is —
constructing it fetches it once into the **Data dir**, re-reads it on every construction after, and
answers with no **Genome** open. It answers two questions and only two: which **Gene id stem**s a
foreign id names, and which foreign ids a stem names. Gene level only, because a gene, a transcript
and a protein have different keys and different sources and are three objects rather than one table
with a level column — the illegal state that is already sitting in the shipped human census. A query
runs against exactly one set (ADR-0017), which is what makes merging two publishers inexpressible
rather than merely discouraged.
_Avoid_: crosswalk, id map, mapping table, id translation table — each names a flat file that has
forgotten who asserted it and when; xref database, xref table; and the bare plural "the xrefs",
which says neither whose nor from when

**Xref source**:
The publisher whose assertions one **Xref set** carries — named on every answer, and a row in a
shipped metadata table beside a reader rather than a branch in code, so adding one is data. A
publisher is eligible only if it keeps dated releases at stable URLs (ADR-0018): the Alliance of
Genome Resources, Ensembl's per-species dumps, HGNC's quarterly archive and WormBase qualify, while
NCBI Gene, UniProt idmapping and MGI are rebuilt daily and overwritten in place, so a checksum
shipped in the wheel would be wrong within a day. Pinnable is not the same as reachable: WormBase's
own download host answers 403 to an automated client, so its assertions are carried through the copy
the Alliance serves of WormBase's own submission — the same bytes, at a URL a script can fetch.
Two sources are not two grades of one thing and
neither is a fallback for the other: they number their releases differently, carry different
**Namespace**s, and some grade each assertion where others grade none — a filter on that grading is
therefore a capability of the source, and one no row of a release satisfies raises rather than
building a set that can only answer nothing.
_Avoid_: provider, vendor, database, resource; authority on its own — a species authority is one
source among several and the word does not say which; and **Source** unqualified, which the
[Assembly](./assembly.md) context already owns for where an assembly's bytes came from

**Default xref source**:
Per species **and per question**, the **Xref source** a caller who names none is answered by, so
that everyone in the lab reaches for the same one without discussing it — exactly the job **Default
annotation** already does. It is a default and not a recommendation: naming a source is how the
scientific choice gets made deliberately, and this exists so that declining to make it is still
reproducible rather than arbitrary. There are two questions and so two flags in the curated table
(ADR-0021), because the publisher carrying a species' identifiers is usually not the one carrying
its **Symbol match**es: human's identifiers default to the Alliance, which publishes no human symbol
at all, and its symbols to HGNC. The question is named where the source is filled in and nowhere
else — a set already built for one source is never answered from another's bytes, and one that
carries no symbol raises naming the one that does.
_Avoid_: preferred source, canonical source, best source, authoritative source — no such judgement
is made here; primary and fallback, which imply a second source is tried when the first misses, and
none is; and a second term for the symbol half — it is one idea with two flags, not two ideas

### What an identifier is

**Namespace**:
The identifier system one foreign id belongs to — Entrez GeneID, HGNC id, UniProt accession, MGI id,
WormBase gene id, symbol — named explicitly on every call, because the string on its own does not
say. One namespace is the hub: Ensembl's, whose ids reduce to **Gene id stem**s and whose shape is
per-species, `ENSG…`, `ENSMUSG…`, and for worm the WormBase gene id, which *is* the Ensembl stable
gene id. Every other namespace is a spoke reached only through that hub and never directly from
another spoke (ADR-0017), and one an **Xref source** does not carry raises and names the ones it
does.
_Avoid_: id type, id system, id space; database, which is the **Xref source**; prefix — a CURIE
prefix spells a namespace and is not one; and never for a Python module path or for a publisher's
namespaced column names, which is what `genome.tf.cofactor` and a **Cofactor table**'s columns mean

**Symbol match**:
One hit of a gene symbol against an **Xref set**, carrying the kind of spelling that matched —
`approved`, `previous` or `alias`. The kind rides on the match rather than being filtered away on
the way out, and a symbol naming several genes answers with all of them, because ambiguity here is
the return type and not an edge case: a table spelling a gene the way it was spelled five years ago
is otherwise dropped without a word, which is what would have happened to 31 of 801 EpiFactors rows.
Matching is exact by default and case-insensitive only when asked for, and the insensitive path
still answers with every gene matched rather than picking one. Which kinds a set can match is the
**Xref source**'s and never the species' — mouse and worm match approved spellings only, their
authorities' typed spellings belonging to publishers that cannot be pinned or cannot be fetched
(ADR-0018) — and why the others are missing rides back on the answer, since *this gene is absent*
and *this source cannot match that spelling* would otherwise both be silence. The reverse direction
is not symmetric: from a stem the answer is the authority's single current approved symbol, which is
labelling.
_Avoid_: gene name; HGNC symbol, which is one authority's spelling of a general idea; synonym — the
authorities' own word is *alias*, and one idea gets one word; fuzzy match, which nothing here does;
and **Motif name**, which labels a matrix and names no gene at all
