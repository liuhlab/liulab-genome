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

## Ensembl — `ensembl`, release 116

Dyer, S. C., *et al.* **Ensembl 2025.** *Nucleic Acids Research* 53(D1), 2025. PMID 39656687.
doi:10.1093/nar/gkae1071 — the per-species `entrez` TSV dump of release **116**, from
<https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/Homo_sapiens.GRCh38.116.entrez.tsv.gz> and
<https://ftp.ensembl.org/pub/release-116/tsv/mus_musculus/Mus_musculus.GRCm39.116.entrez.tsv.gz>.

One file per species, and a **second** source rather than a better one. Release 116 is what the
lab's registered `gencode_v50` annotation corresponds to, which is why it is the release pinned.

**How the release was pinned.** `Homo_sapiens.GRCh38.<release>.entrez.tsv.gz` answers 200 for every
release probed from **88** through **116**, so old releases stay live and the source is eligible
(ADR-0018). One thing a URL written from memory gets wrong: **`ftp.ensembl.org/pub/current_tsv/` is
a 404**. There is no `current_` shortcut for these dumps and the release number is mandatory — which
is convenient, since it means every URL in the table above already names the release it pins.

### Ensembl is not the equal of the first source

Measured on human release 116 against NCBI's own `gene2ensembl`, both sides reduced to distinct
(Entrez GeneID, ENSG stem) pairs:

| | Pairs |
|---|---|
| NCBI asserts | 38,577 |
| Ensembl asserts | 36,824 |
| both assert | 27,554 |
| **agreement (intersection over union)** | **57.6%** |

**The cause is method, not release skew.** NCBI's mapping is a sequence match at a published overlap
threshold and comes out all but one-to-one — no GeneID names more than two stems. Ensembl's fans out
by two orders of magnitude at the tail:

| | NCBI | Ensembl |
|---|---|---|
| distinct GeneIDs | 38,546 | 28,433 |
| GeneIDs naming exactly one stem | 38,515 (99.9%) | 25,890 (91.1%) |
| **most stems named by one GeneID** | 2 | **72** (GeneID `79166`) |
| **most GeneIDs naming one stem** | 4 | **208** (`ENSG00000278233`) |

So an answer from here is wider than the Alliance's for the same id — `79166` names two stems in
Alliance 9.0.0 and seventy-two here — and the width is Ensembl's assertion rather than a fault.
Nothing narrows it and nothing reconciles the two: two publishers are two answers (ADR-0017).

### Every row is `DEPENDENT` and none is `DIRECT`

Ensembl grades each cross-reference in an `info_type` column, so the obvious quality filter is
`DIRECT`. Counted over the whole release-116 file:

| `db_name` | Rows | `DIRECT` | `DEPENDENT` |
|---|---|---|---|
| `EntrezGene` (human) | **552,633** | **0** | 552,633 |
| `EntrezGene` (mouse) | **358,853** | **0** | 358,853 |
| `EntrezGene_trans_name` (human) | 644,668 | — | `MISC` |

**Filtering to `DIRECT` yields an empty set rather than a smaller one**, which would answer every
query with nothing and look exactly like a gene list that matched none of it. Asking for it raises
and names what the release actually carries, rather than building a set that can only say nothing.

### The published checksum covers the *served* bytes — the opposite of the Alliance's

Each directory ships a `CHECKSUMS` file in BSD `sum` format, computed over the `.tsv.gz` **exactly
as downloaded**: `28782 5952` for release 116's human dump, verified against the file. The Alliance
publishes the other convention, an md5 of the TSV *inside* its gzip. Neither can be assumed of the
other, and a BSD `sum` is a 16-bit checksum — no integrity check at all for a 6 MB file. So the rows
above pin an **md5 of the unpacked bytes** computed here on the pinned release, which is the one
convention every row of this table uses (ADR-0006), and Ensembl's own `sum` values are recorded here
as the publisher's cross-check rather than in the table.

### What the reader keeps, and what it drops

Only the `EntrezGene` rows are cross-references. The dump also carries `EntrezGene_trans_name` rows
whose `xref` column holds a **transcript name** — `KU-MEL-3-201`, not a GeneID — because, as the
file's own `README_entrez.tsv` warns, it "contains all Ensembl external database names which started
with entrez so duplication of hits is possible". Reading them would key the Entrez namespace by
transcript labels. They are a different assertion and are dropped.

The file is transcript- and protein-grained, so 552,633 human rows collapse to 36,824 gene-level
pairs; deduplication is on the pair and happens as the file is read.

| Species | Assembly | Gene-level pairs | Namespaces carried |
|---|---|---|---|
| *Homo sapiens* | GRCh38 | 36,824 | Ensembl, Entrez |
| *Mus musculus* | GRCm39 | 28,681 | Ensembl, Entrez |

**No worm row, and not by oversight.** Ensembl files *C. elegans* under Ensembl Genomes' own
numbering: release-116's worm directory holds `Caenorhabditis_elegans.WBcel235.63.entrez.tsv.gz`, so
"116" would name a file that does not exist. Worm is answered by the Alliance, where the hop is the
identity function, and no worm row is invented here to make the table look symmetrical.

## HGNC — `hgnc`, the quarterly archive file of 2026-07-07

Seal RL, Braschi B, Gray K, McClay J, Tweedie S, Bruford EA. **Genenames.org: the HGNC and PGNC
resources in 2026.** *Nucleic Acids Research* 54(D1):D1098-D1107, 2026. PMID 41287213.
doi:10.1093/nar/gkaf1229 — the quarterly archive file dated **2026-07-07**, from
<https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-07-07.txt>.

Human only, and the **only source here that publishes previous and alias spellings typed**. That is
what it is for: the other three carry ids, and the two that carry symbols at all carry the current
one. HGNC's `prev_symbol` and `alias_symbol` columns say which kind a spelling is, so a **Symbol
match** can say it too rather than this package guessing.

**The failure it prevents, measured on this repo's own shipped data.** Of EpiFactors v2.0's 801
human rows, all 801 carry an HGNC id and all 801 resolve — and **31 still spell the gene by a symbol
HGNC has retired**: `ARNTL` for `BMAL1`, `C11orf30` for `EMSY`, `ACINUS` for `ACIN1`. A symbol join
that knows approved spellings only mis-keys or drops exactly those 31, and says nothing while doing
it. The same measurement is recorded beside the cofactor table it was made on.

**How the release was pinned, and why not from a date.** The archive's file names are dated and the
dates are **irregular** — `2024-07-02`, `2025-01-06`, `2025-10-07`, and both `2026-07-03` and
`2026-07-07` in one quarter — so a URL assembled from *the first of the quarter* is a 404 about half
the time. The pin names a file read out of the bucket's own listing:
`https://storage.googleapis.com/storage/v1/b/public-download-files/o?prefix=hgnc/archive/archive/quarterly/tsv`,
which returned **27 quarterly snapshots** back to `hgnc_complete_set_2020-07-01.txt` when this row
was written. Old files stay live, which is what makes HGNC eligible (ADR-0018).

**The path has a doubled segment, and the remembered one is dead.**
`ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt` is a live 404, and so is the
`…/new/tsv/…` variant. HGNC serves from a Google Cloud Storage bucket whose archive path repeats
itself — `…/hgnc/archive/archive/quarterly/tsv/…`. Both `archive` segments are real.

**The published checksum covers the served bytes, which here are also the unpacked ones.** The file
is plain text and not gzipped, so the bucket listing's own `md5Hash` —
`cd41d33955722de9ac0e14a2557ef5fc` for 16,913,890 bytes, verified against the download — is the
digest of exactly what this table's one convention pins (ADR-0006). The two conventions that sit
side by side in this directory coincide for once; **do not carry that across** to a source whose
file is compressed, where the Alliance's md5 is of the unpacked bytes and Ensembl's `sum` is of the
served ones.

### The schema has drifted, so the reader parses by header name

The header row is **52 columns** wide in `hgnc_complete_set_2020-07-01.txt` and **54** in
`…_2023-01-01.txt` and every file since. A reader indexing by position reads the wrong column on the
older snapshots and reports nothing wrong, so this one finds its seven columns by name — `hgnc_id`,
`symbol`, `prev_symbol`, `alias_symbol`, `ensembl_gene_id`, `entrez_id`, `uniprot_ids` — and raises
naming the missing one if the publisher renames it. Reordering or adding a column cannot change the
answer, and a named test reverses the whole header to prove it.

**A multi-valued cell is quoted and pipe-separated; a single value is bare.** `alias_symbol` reads
`"MOP3|JAP3|PASD3|bHLHe5|ARNTL1"` with the quotes in the file, and `AC3` without them. Splitting
before stripping the quotes keys the namespace by `"MOP3` and `ARNTL1"`, which nobody types.

| Species | Rows | With a hub | Namespaces carried |
|---|---|---|---|
| *Homo sapiens* | 45,019 | 42,337 | Ensembl, Entrez, UniProt, HGNC, symbol (approved, previous, alias) |

**2,682 rows carry no `ensembl_gene_id`** — mostly pseudogenes and non-coding RNAs Ensembl does not
carry — and a row with no hub has nothing to hang a namespace off, so it appears in no slice. That
is the same silence the Alliance reader keeps for a gene listed with no `ENSEMBL:` cross-reference.

## The Alliance's per-species gene submissions — `alliance_bgi`, release 9.0.0

Alliance of Genome Resources Consortium. **Updates to the Alliance of Genome Resources central
infrastructure.** *Genetics* 227(1), 2024. PMID 38552170. doi:10.1093/genetics/iyae049 — the
`BGI` (Basic Gene Information) submissions of release **9.0.0**, from
<https://download.alliancegenome.org/9.0.0/BGI/MGI/1.0.2.5_BGI_MGI_0.json.gz> and
<https://download.alliancegenome.org/8.3.0/BGI/WB/1.0.2.5_BGI_WB_4.json.gz>.

**Each file is a species authority's own submission**, and the row records who curated it as well as
who serves it: MGI's states `"release": "MGI 6.27 2026-04-21"` in its own header and WormBase's
states `"release": "WS298"`. Cite the Alliance, whose file it is; the `publisher` column names both.

**It is here because neither authority can be pinned or fetched directly.** MGI keeps no dated
archive — its `downloads/reports/` has one `archive/` subdirectory and that holds only `iphone/`,
and `MGI_Gene_Model_Coord.rpt` opens with a column header and no build or date stamp — so MGI is not
an eligible **Xref source** (ADR-0018). `downloads.wormbase.org` answers **403 to plain `curl`**,
the release directories and not only the blog, measured from two networks on two days. The Alliance
holds a dated, immutable, checksummed copy of both submissions, and is in any case the ongoing
publisher for *C. elegans* now that WS298 is WormBase's final release.

**The worm file's path names 8.3.0 and its release is 9.0.0**, which is the same wrinkle the
cross-reference file has one directory up: the Alliance re-serves an unchanged submission under a
later release rather than rebuilding it, and its file-management API lists 9.0.0 among that file's
`releaseVersions`. A URL assembled from the release number alone would 404. Read the `s3Url` from
`https://fms.alliancegenome.org/api/datafile/by/<release>/BGI/<MGI|WB>?latest=true`.

**The published md5 is of the unpacked JSON**, as the cross-reference file's is of the unpacked TSV
and unlike Ensembl's. Verified by decompressing both and hashing:

| Submission | Served bytes | md5 of the `.gz` | md5 of the unpacked JSON — **what the row pins** |
|---|---|---|---|
| MGI | 6,754,968 | `39c0fe0a8db6aa3009dd2109f0127fa8` | `a834098c9505ec7fb4a0151480a90734` |
| WB | 5,265,513 | `d7064a579923b3ee0e8cdbf00b1b0a6f` | `4a45ce6beb26dd0dc8c053e5b2e1a835` |

The second column of each row is the file-management API's `md5Sum`, matched only after `gzip -d`.

### Approved spellings only, and that is a decision

Each record carries one `symbol` and one `synonyms` list, and **the list is undifferentiated**.
WormBase's `daf-16` record files `daf-17` — a name the gene genuinely went by — in the same list as
the sequence names `R13H8.1` and `CELE_R13H8.1`, and nothing in the file says which is which.
Reading them would mean labelling each one `previous` or `alias` on this package's own authority,
which is the claim it never makes (ADR-0017). So this source publishes `approved` alone, and every
answer it produces carries the sentence saying so — the explanation is behaviour rather than a
comment, because a spelling MGI retired otherwise comes back unresolved and looks exactly like a
gene that is absent.

| Species | Records | With a hub | With a symbol | Namespaces carried |
|---|---|---|---|---|
| *Mus musculus* | 90,776 | 77,476 | **all 90,776** | Ensembl, Entrez, UniProt, MGI, symbol (approved) |
| *Caenorhabditis elegans* | 48,769 | 46,926 | **all 48,769** | Ensembl, Entrez, UniProt, WormBase, symbol (approved) |

**Every record carries a symbol and not every record carries a hub**, which is the shape that makes
this source worth having: the authority names each of its genes, and the join to Ensembl is where
the losses are. Worm's 46,926 hubs are the same 46,926 the cross-reference file reports, and mouse's
77,476 are within six of its 77,482 — two files, two submissions, one authority, agreeing.

### The file is JSON, and it is read a record at a time

One object with a `metaData` header and a `data` array. The mouse submission unpacks to 72 MB and
the worm one to 74 MB, so neither is decoded whole: the reader finds the array and peels one gene
object off the front of a small rolling buffer. That also makes the two publishers' two different
pretty-printers — two-space for MGI, three-space for WormBase — a non-issue, both being JSON.

### Where a symbol is stored, and under what name

A prepared set is the same plain gzipped TSV every other source writes, and the **kind of spelling**
is a value in its `namespace` column: `symbol`, `previous_symbol`, `alias_symbol`. A fourth column
saying *what sort of row this is* would be the level discriminator this design refuses everywhere
else, and the three spellings read for themselves in a shell:

```
zcat homo_sapiens.xref_table.tsv.gz | awk -F'\t' '$1=="previous_symbol"'
```

Only `symbol` is a **Namespace** a caller may name — a previous and an alias spelling are the same
identifier system with a different standing, not systems of their own.

## Adding a source

A row here plus a reader in `genome/xref/`, and nothing else. The row must pin a release that stays
retrievable at a stable URL: a publisher that overwrites its file in place, or that keeps no archive
of the release the row names, cannot be pinned and does not belong here, however good its data
(ADR-0018).

Two things a second source made explicit that the first alone did not. **A publisher's own checksum
may cover the served bytes or the unpacked ones**, and the two conventions sit side by side in this
directory — check which before pinning, because a row pinned to the wrong scope rejects every
download. And **an evidence filter is a capability of the source rather than of every set**: a file
that grades nothing refuses a filter rather than ignoring it, since a filter silently dropped is a
quality claim the caller believes and nobody made.
