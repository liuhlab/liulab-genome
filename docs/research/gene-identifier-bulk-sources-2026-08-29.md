# Gene identifier bulk sources: what each publisher asserts, and which of them pin

**Date:** 2026-08-29. **Method:** the real published files, fetched and counted. Nothing here is
recalled and nothing is taken from a publisher's prose where a file could be counted instead. Every
URL, release and checksum is in the table below so the counts can be re-derived. Large files were
read on GPU71FM; the rest locally. No package code was involved — this measures the publishers.

This note re-runs a survey whose working files were lost. Where a number here differs from the one
that survey reported, the difference is stated at the point it occurs.

## What was fetched

| File | Retrieved from | Bytes | Checksum |
|---|---|---|---|
| `Homo_sapiens.GRCh38.116.entrez.tsv.gz` | `ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/` | 6,093,982 | md5 `0b2237c4bf392930efa130926ba9ceb1`; `sum` 28782 5952, **matching** the directory's own `CHECKSUMS` |
| `Mus_musculus.GRCm39.116.entrez.tsv.gz` | `ftp.ensembl.org/pub/release-116/tsv/mus_musculus/` | 3,156,296 | md5 `cb143db7fa830bb2b8f67468db1081cf` |
| `gene2ensembl.gz` | `ftp.ncbi.nlm.nih.gov/gene/DATA/` | 313,546,863 | md5 `8bf66a95cf203ca7cdbcd26b3da2bb76`; publisher ships none |
| `GENECROSSREFERENCE_COMBINED_11.tsv.gz` | `download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/` | 25,776,975 | md5 `f8bbb8156c75d2d2a5279b24d914de4a` |
| `hgnc_complete_set_2026-07-07.txt` | `storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/` | 16,913,890 | md5 `cd41d33955722de9ac0e14a2557ef5fc` |
| `gencode.v50.metadata.EntrezGene.gz` | `ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/` | 2,235,757 | md5 `959aa97f7dd5c17aafd4fccc0d718f45` |
| `1.0.2.5_BGI_MGI_0.json.gz` | `download.alliancegenome.org/9.0.0/BGI/MGI/` | 6,754,968 | md5 `39c0fe0a8db6aa3009dd2109f0127fa8`; unpacked `a834098c9505ec7fb4a0151480a90734`, **the publisher's own** |
| `1.0.2.5_BGI_WB_4.json.gz` | `download.alliancegenome.org/8.3.0/BGI/WB/` | 5,265,513 | md5 `d7064a579923b3ee0e8cdbf00b1b0a6f`; unpacked `4a45ce6beb26dd0dc8c053e5b2e1a835`, **the publisher's own** |

`gene2ensembl.gz` carried `Last-Modified: Sat, 29 Aug 2026 08:05:04 GMT` — the morning of this
measurement. That is the caveat on section 1 and the evidence for section 3, and it is why the NCBI
md5 above is a record of what this session read rather than something a later reader can match.

**The Alliance's published md5 is of the uncompressed TSV, not of the `.gz` that is served.** The
file-management API reports `2575a9b53d3de2346267c7cb34476585` for
`GENECROSSREFERENCE_COMBINED_11.tsv.gz`; that is the md5 of the stream after `gzip -d`, and the md5
of the bytes on the wire is the one in the table. The same holds for
`ORTHOLOGY-ALLIANCE_COMBINED_13.tsv.gz` (published `7d091b84f3227ea576826ed1f02642e9`, uncompressed;
`c684c49f4d6a25f2bb8b9f325b383db1` compressed). Both were checked by decompressing to `md5sum`.

## 1. NCBI and Ensembl agree on 57.6% of the human gene-level pairs they assert between them

Both sides reduced to distinct (Entrez GeneID, ENSG stem) pairs: `gene2ensembl` filtered to
`tax_id == 9606` and read at its `GeneID` and `Ensembl_gene_identifier` columns; Ensembl's
`entrez.tsv` filtered to `db_name == EntrezGene` and read at its `xref` and `gene_stable_id`
columns, with the version suffix stripped. Ensembl's file is transcript- and protein-grained —
552,633 rows collapse to 36,824 gene-level pairs.

| | Pairs |
|---|---|
| NCBI asserts | 38,577 |
| Ensembl asserts | 36,824 |
| **both assert** | **27,554** |
| either asserts | 47,847 |
| **agreement (intersection over union)** | **57.6%** |
| of NCBI's own assertions | 71.4% |
| of Ensembl's own assertions | 74.8% |

**The cause is method, not release skew.** The two files are built by different procedures and the
shape of the disagreement follows the procedures:

| | NCBI | Ensembl |
|---|---|---|
| distinct GeneIDs | 38,546 | 28,433 |
| distinct ENSG stems | 38,295 | 31,492 |
| GeneIDs naming exactly one stem | 38,515 (99.9%) | 25,890 (91.1%) |
| **most stems named by one GeneID** | **2** | **72** (GeneID `79166`) |
| **most GeneIDs naming one stem** | **4** | **208** (`ENSG00000278233`) |

NCBI's is a sequence match at a published overlap threshold and comes out all but one-to-one — no
GeneID names more than two stems. Ensembl's fans out by two orders of magnitude at the tail.

**Caveat, and the one number that moved.** A rolling NCBI build was compared against a pinned
Ensembl release, so a little of the gap is skew. The earlier survey reported 38,545 NCBI pairs,
27,526 agreeing and 57.5%; today's build gives 38,577, 27,554 and 57.6%. The Ensembl side reproduced
**exactly** — 36,824 — which is what a pinned release is for, and which also identifies release 116
as the release that survey read. The drift is entirely NCBI's, at 32 pairs over the interval, and it
cannot be resolved by re-fetching, because the build that produced 38,545 was overwritten in place.

## 2. Every human EntrezGene row in Ensembl's TSV is `DEPENDENT`, and none is `DIRECT`

Counted over the whole file, by `info_type`:

| `db_name` | Rows | `DIRECT` | `DEPENDENT` |
|---|---|---|---|
| `EntrezGene` | 552,633 | **0** | 552,633 |
| `EntrezGene_trans_name` | 644,668 | — | — |

Mouse is the same shape: 358,853 `EntrezGene` rows, all `DEPENDENT`, none `DIRECT` (28,681 distinct
gene-level pairs). A filter to `info_type == DIRECT` on either species returns an empty set, not a
smaller one.

The file's own `README_entrez.tsv` notes that it "Dumps contain all Ensembl external database names
which started with entrez so duplication of hits is possible" — the `EntrezGene_trans_name` rows
above are that duplication, and they are a different assertion (a transcript name) from the
`EntrezGene` rows.

## 3. Version suffixes disagree inside a single NCBI row

Over the 96,205 human rows of `gene2ensembl`, counting the values that are not `-`:

| Column | Bare | Versioned |
|---|---|---|
| `Ensembl_gene_identifier` | **96,205** | 0 |
| `Ensembl_rna_identifier` | 0 | **86,904** |
| `Ensembl_protein_identifier` | 0 | **70,278** |

One row carries a bare gene id beside a versioned transcript id and a versioned protein id. Joining
either of the latter to a bare id returns zero rows and raises nothing.

## 4. Which sources keep dated releases, and what was fetched to establish it

Every row below is a live probe made on 2026-08-29, not a reading of documentation.

| Source | Pins? | What was fetched |
|---|---|---|
| **Alliance of Genome Resources** | yes | `fms.alliancegenome.org/api/datafile/by/GENECROSSREFERENCE?latest=true` returns a versioned `s3Path` (`9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz`), an `md5Sum`, an `uploadDate` and the releases it belongs to. The file header repeats the version and the generation date. |
| **Ensembl per-species TSV** | yes | `Homo_sapiens.GRCh38.<r>.entrez.tsv.gz` answers 200 for every release probed from **88** through **116**. Each directory ships `CHECKSUMS` in `sum` format, verified above. |
| **HGNC quarterly archive** | yes | 27 quarterly snapshots listed under `hgnc/archive/archive/quarterly/tsv/`, from `hgnc_complete_set_2020-07-01.txt` to `…_2026-07-07.txt`. |
| **WormBase** | yes, permanently | see section 6 |
| **NCBI Gene** | no | `gene2ensembl.gz` and `gene_info.gz` both carried a `Last-Modified` of 2026-08-29, hours before the fetch. One path, overwritten; no dated directory exists beside it. |
| **UniProt idmapping** | no | `current_release/knowledgebase/idmapping/` is the only path that serves these files. `previous_releases/release-2026_01/knowledgebase/` holds four tarballs and **no `idmapping/` directory**; the same for every dated release directory listed. |
| **MGI** | no | `downloads/reports/` has one `archive/` subdirectory and it contains only `iphone/`. `MGI_Gene_Model_Coord.rpt` (`Last-Modified` 2026-08-24) opens with a column header and no build or date stamp; `MGI_EntrezGene.rpt` has no header line at all. |

Three details a URL written from memory gets wrong, each verified by probe:

- **`ftp.ensembl.org/pub/current_tsv/` is a 404.** There is no `current_` shortcut for the TSV dumps;
  the release number is mandatory.
- **HGNC's EBI FTP path is gone.** `ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt`
  and the `…/new/tsv/…` variant both 404. The live path is a Google Cloud Storage bucket with a
  **doubled path segment** — `…/hgnc/tsv/tsv/…` for current, `…/hgnc/archive/archive/quarterly/tsv/…`
  for the archive.
- **HGNC's schema has drifted.** The header row is 52 columns wide in `hgnc_complete_set_2020-07-01.txt`
  and 54 in `…_2023-01-01.txt`, `…_2026-04-01.txt` and the current file. A reader indexing by position
  reads the wrong column on the older snapshots.

**One correction to the earlier survey.** It grouped UniProt with NCBI and MGI as "rebuilt daily".
UniProt's idmapping files are not daily: `HUMAN_9606_idmapping.dat.gz` carried `Last-Modified`
2026-06-10, roughly a UniProt release cadence. The verdict is unchanged and rests on a different
fact — the files exist at one path only, are replaced there at each release, and no dated release
directory retains a copy.

## 5. The Alliance file: duplicate rows, and the worm identity

2,659,704 data rows across ten taxa.

| | Rows |
|---|---|
| data rows | 2,659,704 |
| distinct whole rows | 2,659,704 |
| **distinct on (`GeneID`, `GlobalCrossReferenceID`, `TaxonID`)** | **1,811,267** |
| keys appearing more than once | 495,212 |

No two rows are identical, so a naive whole-row `uniq` removes nothing. On the three columns that
carry the assertion, **31.9% of rows are redundant**; the repeats differ only in
`CrossReferenceCompleteURL` and `ResourceDescriptorPage`.

Two per-species counts from the same file:

- **Worm: 46,926 genes carry an `ENSEMBL` cross-reference, and in all 46,926 the WormBase gene id and
  the Ensembl id are the same string** — `WB:WBGene00000001 → ENSEMBL:WBGene00000001`, with zero rows
  where they differ.
- **Human: 40,665 genes carry at least one `ENSEMBL` cross-reference, and 2,535 of them carry more
  than one.**

## 6. WormBase WS298 is WormBase's final release

The announcement is on the WormBase blog, which is Cloudflare-blocked to automated clients. The
Alliance community forum carries the same text and serves it to an unauthenticated client:
<https://community.alliancegenome.org/t/announcing-the-final-release-of-wormbase/8461> — HTTP 200,
21,978 bytes, fetched with plain `curl` and no headers. It states that WS298 "will be the final major
WormBase release", schedules a test release for mid-October and the production release for the end
of November, and says the site will be "anchored to the WS298 database release and entered into an
archival maintenance mode where only critical bug fixes are addressed", with *C. elegans* curation
moving to the Alliance.

**`downloads.wormbase.org` is itself 403 to `curl`** — the release directories, not only the blog.
Anything automated that needs WormBase bytes needs a different host or a browser-shaped client.

## 7. GENCODE's metadata files are Ensembl's xref output, measured rather than asserted

GENCODE v50's `_README.TXT` states "This release corresponds to Ensembl version 116". Comparing
`gencode.v50.metadata.EntrezGene.gz` against Ensembl release-116's `entrez.tsv`, both reduced to
distinct (ENST stem, GeneID) pairs:

| | Pairs |
|---|---|
| GENCODE asserts | 532,987 |
| **also in Ensembl's `entrez.tsv`** | **532,859 (99.98%)** |
| in GENCODE only | 128 |
| in Ensembl only | 19,774 |

GENCODE's file is very nearly a subset of Ensembl's. Two files agreeing at 99.98% are one assertion
counted twice, not two publishers corroborating each other.

A second property of the same files, from the README: **every GENCODE metadata file is keyed by
transcript id** — `metadata.EntrezGene` is (transcript id, Entrez Gene id), `metadata.HGNC` is
(transcript id, symbol, HGNC id), `metadata.RefSeq` and `metadata.SwissProt` likewise. There is no
gene-level metadata file; reading one at gene level requires a join the file does not carry.

## 8. Which of the pinnable sources carries a gene symbol, and of what kind

Added while building symbol matching, on the same day and the same files. Every count below is over
the whole published file, not a sample.

| Source | Carries a symbol? | Typed? |
|---|---|---|
| **HGNC quarterly archive** | yes, human only | **yes** — `symbol`, `prev_symbol` and `alias_symbol` are three columns |
| **Alliance `GENECROSSREFERENCE_COMBINED`** | **worm only** | n/a |
| **Alliance `BGI` per-species submissions** | yes, every contributing database | **no** — one `symbol` and one undifferentiated `synonyms` list |
| **Ensembl per-species TSV dumps** | no | — |

**The Alliance cross-reference file carries no human or mouse symbol at all.** Counted over all
2,659,704 rows, grouped by taxon and cross-reference prefix: worm's `WB:`-prefixed rows include
47,156 under the page `gene/spell` whose value is a symbol (`WB:WBGene00000001 → WB:aap-1`), and
every human `HGNC:`-prefixed row (44,569) and every mouse `MGI:`-prefixed row (246,330 across six
pages) carries the gene's own id and never a name. So this file cannot supply mouse symbols, and
the earlier assumption that it could was wrong.

**`GENE_TSV_COMBINED` would have supplied all three species and is not pinnable.**
`9.0.0/downloads/GENE_TSV_COMBINED.tsv.gz` (54 MB) is one row per gene over nine species with
`GeneSymbol`, `GeneSynonyms` and `GeneCrossReferences` columns — and the same path under `8.3.0/`
and `7.4.0/` **404s**, and `fms.alliancegenome.org/api/datafile/by/GENE/COMBINED` errors. It is a
website download rather than a file-management artifact, so it fails the eligibility bar (ADR-0018).

**The `BGI` submissions do pin, and their md5 covers the unpacked bytes.** The file-management API
answers `by/9.0.0/BGI/MGI` and `by/9.0.0/BGI/WB` with a versioned `s3Path`, an `md5Sum` and the
releases each belongs to. Worm's is `8.3.0/BGI/WB/1.0.2.5_BGI_WB_4.json.gz` re-served under 9.0.0
and 9.1.0 — the same re-serving the cross-reference file does — so a URL built from the release
number alone 404s. Verified by decompressing and hashing:

| Submission | Served | Unpacked | md5 of served | md5 of unpacked — **matches the API's `md5Sum`** |
|---|---|---|---|---|
| MGI | 6,754,968 | 71,819,688 | `39c0fe0a8db6aa3009dd2109f0127fa8` | `a834098c9505ec7fb4a0151480a90734` |
| WB | 5,265,513 | 74,123,056 | `d7064a579923b3ee0e8cdbf00b1b0a6f` | `4a45ce6beb26dd0dc8c053e5b2e1a835` |

| Submission | Own `release` | Records | With an `ENSEMBL:` cross-reference | With a symbol |
|---|---|---|---|---|
| MGI | `MGI 6.27 2026-04-21` | 90,776 | 77,476 | **all 90,776** |
| WB | `WS298` | 48,769 | 46,926 | **all 48,769** |

Worm's 46,926 hubs are exactly the 46,926 section 5 counts in the cross-reference file, and mouse's
77,476 are within six of that file's 77,482 — two Alliance products, one authority each, agreeing.

**`downloads.wormbase.org` 403s from a second network too.** Re-probed from GPU71FM as well as
locally, on `releases/WS298/` and on `ftp.wormbase.org/pub/wormbase/releases/WS298/…`, with and
without a browser `User-Agent`: 403 every time. Section 6's finding holds, and it is why worm's
symbols are read from the Alliance's copy of WormBase's own submission rather than from WormBase.

**HGNC's archive listing, re-read.** 27 quarterly `hgnc_complete_set_*` files from `2020-07-01` to
`2026-07-07`, and the dates are irregular enough that a URL built from *the first of the quarter* is
wrong about half the time: `2024-07-02`, `2025-01-06`, `2025-10-07`, `2026-01-06`, and **both**
`2026-04-01` and `2026-04-07`, **both** `2026-07-03` and `2026-07-07`. Each file is also served
under a `.tsv` name with the identical md5. The file is plain text, so the bucket's `md5Hash` is
both the served and the unpacked digest — a coincidence of this source and not a convention.

## 9. What this note could not verify

- **HCOP.** Neither `ftp.ebi.ac.uk/pub/databases/genenames/hcop/…` nor any `hgnc/hcop` prefix in the
  HGNC bucket returns a file; the bucket listing for that prefix is empty. So the bulk files are not
  fetchable from either the old or the new location today, but the earlier survey's specific claim
  that HCOP has been *stale since 2025-05-16* is **carried, not re-measured** — there is no file left
  to date.
- **OrthoList 2 frozen since 2018** is carried from the earlier survey unverified.
- **The earlier survey's NCBI figures** (38,545 / 27,526 / 57.5%) cannot be reproduced at all, in
  either direction, because the build that produced them no longer exists. Section 1 reports today's
  build.
