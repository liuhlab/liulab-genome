# Attribution — the Xref sources

**No cross-reference table ships here.** This directory holds one curated row per **Xref set**, and
nothing else: which species, which publisher, which **Release**, where its file is fetched from and
what that file's *unpacked* bytes hash to. The tables themselves are downloaded into the **Data
dir** on first construction, because the smallest qualifying file is 25 MB and because a checksum
that ships in a wheel must still match a year later (ADR-0018). What follows is who is to be cited
for each row, and the facts about each publisher's file that the reader beside it encodes.

**Every assertion an `XrefSet` answers with is its publisher's.** This package asserts no pair of
its own, merges no two publishers and composes no two hops (ADR-0017). Cite the publisher whose set
you used; `XrefMetadata.attribution()` renders the line to print.

## Alliance of Genome Resources — `alliance`, release 9.0.0

Alliance of Genome Resources Consortium. **Updates to the Alliance of Genome Resources central
infrastructure.** *Genetics* 227(1), 2024. PMID 38552170. doi:10.1093/genetics/iyae049 — the
`GENECROSSREFERENCE_COMBINED` gene cross-reference file of release **9.0.0**, from
<https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz>.

One file, gene level only, covering ten organisms; three of them are the lab's, and this package
prepares a per-species slice of each. It is the **Default xref source** for all three.

**How the release was pinned.** The Alliance's file-management API answers
`https://fms.alliancegenome.org/api/datafile/by/GENECROSSREFERENCE/COMBINED?latest=true` with the
record's own `s3Url` and `md5Sum`, and `/api/releaseversion/current` names the current release. Two
things that a URL built from memory gets wrong, both checked against the live API: the artifact's
path carries **9.0.0** while the current release is 9.1.0 — 9.1.0 re-serves the 9.0.0 file rather
than rebuilding it, so a URL assembled from the current release 404s — and the file itself states
`# Alliance Database Version: 9.0.0` in its own header, which is what the `version` column records.
Old releases stay live: the same endpoint under `/api/datafile/by/<release>/GENECROSSREFERENCE/COMBINED`
lists 521 builds back to 3.0.0 in 2019, and files that old still download. That is what makes the
Alliance eligible at all (ADR-0018).

**The published md5 is of the uncompressed TSV.** `md5Sum` is
`2575a9b53d3de2346267c7cb34476585`, which is the digest of the 517 MB TSV *inside* the gzip and not
of the 25 MB `.tsv.gz` that is served — the served file's own md5 is the S3 ETag,
`f8bbb8156c75d2d2a5279b24d914de4a`. That is the convention this package already holds itself to
(ADR-0006), so the row pins the publisher's number unchanged and the reader decompresses before it
hashes. Checking the `.gz` against `md5Sum` fails every time, which is worth knowing before anyone
"fixes" it.

### What each species' slice is cut from

The Alliance keys every row by its own gene curie — `HGNC:…` for human, `MGI:…` for mouse, `WB:…`
for worm — and the **Gene id stem** is the `ENSEMBL:` cross-reference on that gene. A gene the
Alliance lists with no `ENSEMBL:` cross-reference has no hub and appears in no slice: **3,904 of
44,569 human genes and 10,668 of 88,150 mouse genes**, and **none at all** of the 46,926 worm genes.

| Species | Taxon | Genes | With a hub | Namespaces carried |
|---|---|---|---|---|
| *Homo sapiens* | `NCBITaxon:9606` | 44,569 | 40,665 | Ensembl, Entrez, UniProt, HGNC |
| *Mus musculus* | `NCBITaxon:10090` | 88,150 | 77,482 | Ensembl, Entrez, UniProt, MGI |
| *Caenorhabditis elegans* | `NCBITaxon:6239` | 46,926 | 46,926 | Ensembl, Entrez, UniProt, WormBase |

**The three hops have three different shapes, and none is inferable from an id string.** For worm it
is the identity: all 46,926 `WB:WBGene…` genes carry `ENSEMBL:WBGene…`, the same string, with zero
differing — so a WormBase gene id *is* a **Gene id stem**. For mouse it is a real join onto `ENSMUSG…`.
For human it reaches Ensembl only through HGNC, and **2,535 human genes carry more than one
`ENSEMBL:` cross-reference**, so 6.2% of HGNC ids name two stems or more and nothing here picks one.
That spread is why this is an object a caller opens rather than a step performed for them.

### What the reader drops, and why

The file carries cross-references to PANTHER, RGD, ExpressionAtlas, RNAcentral and RefSeq besides
the four namespaces above. None is a gene-level identifier system this package answers in, so none
is carried — a level discriminator in the data would make the illegal state representable. Two
details that would otherwise mislead:

- **A cross-reference's prefix says nothing about species.** `RGD:` is the single most frequent
  cross-reference prefix on *human* rows. Species is the `TaxonID` column and only that column.
- **The species authority's id is the `GeneID` column and never the authority-prefixed
  cross-reference row.** For worm those rows carry *symbols* — `WB:WBGene00000001 → WB:aap-1` — so
  reading the authority off the cross-reference column would key 47,156 worm rows by gene symbol
  while calling them WormBase gene ids.

### The duplication, measured

**2,659,704 rows reduce to 1,811,267 distinct `(GeneID, GlobalCrossReferenceID, TaxonID)`** — 31.9%
redundant, because the Alliance emits one row per web page it links a pair from. A whole-row `uniq`
removes **none** of it, the differing columns being the URL and the page name. Deduplication is on
the key and it happens as the file is read, so an answer counts genes rather than links.

### What ships on disk, and where

A prepared set is `$LIULAB_DATA/xref/alliance/9.0.0/<species slug>/<species slug>.xref_table.tsv.gz`
— a **plain gzipped TSV** of `namespace`, `xref_id`, `gene_id_stem`, sorted and unique, so a
collaborator who does not use Python reads it in R or a shell with no library at all. Beside it a
**Completion marker** records the URL, the publisher's own md5 as provenance and the slice's own
sha256 as the integrity check; two checksums because what is stored is a derived slice rather than
the publisher's bytes. The gzip is written with no modification time, so one release sliced on two
machines produces byte-identical files.

## Adding a source

A row here plus a reader in `genome/xref/`, and nothing else. The row must pin a release that stays
retrievable at a stable URL: a publisher that overwrites its file in place, or that keeps no archive
of the release the row names, cannot be pinned and does not belong here, however good its data
(ADR-0018).
