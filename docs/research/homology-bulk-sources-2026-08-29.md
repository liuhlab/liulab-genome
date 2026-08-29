# Homology bulk sources: Compara's per-species partition, and how far three publishers agree

**Date:** 2026-08-29. **Method:** Ensembl Compara release 116's per-species protein homology dumps
for human, mouse and worm, downloaded to GPU71FM and counted row by row; NCBI's `gene_orthologs` and
the Alliance's `ORTHOLOGY-ALLIANCE` file read the same way. All three Compara files were verified
against the publisher's own `MD5SUM` before anything was counted. No package code was involved.

This note re-runs a survey whose working files were lost. Section 4 is the one place where the
re-measurement lands somewhere different, and it says so there.

## What was fetched

| File | Retrieved from | Bytes | md5 |
|---|---|---|---|
| `Compara.116.protein_default.homologies.tsv.gz` (human) | `ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/` | 109,478,724 | `59857f48bbbdf6812999d58d7a24ccc4` ✓ |
| the same, mouse | `…/homologies/mus_musculus/` | 111,917,367 | `8f9870f0f12ece5032f8e62117de9924` ✓ |
| the same, worm | `…/homologies/caenorhabditis_elegans/` | 85,130,626 | `2d2b4ef12d7cb8acd84bd560f055daa0` ✓ |
| `gene_orthologs.gz` | `ftp.ncbi.nlm.nih.gov/gene/DATA/` | 128,872,562 | `1b4c83b33d1662407eee7e3c114a9d33` |
| `ORTHOLOGY-ALLIANCE_COMBINED_13.tsv.gz` | `download.alliancegenome.org/9.0.0/ORTHOLOGY-ALLIANCE/COMBINED/` | 15,302,140 | `c684c49f4d6a25f2bb8b9f325b383db1` |
| `GENECROSSREFERENCE_COMBINED_11.tsv.gz` | `download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/` | 25,776,975 | `f8bbb8156c75d2d2a5279b24d914de4a` |

✓ marks an md5 checked against the `MD5SUM` file in the same directory. `gene_orthologs.gz` carried
`Last-Modified: Sat, 29 Aug 2026 08:15:49 GMT`, the morning of the fetch, and ships no checksum. The
Alliance md5s published by its file-management API are of the **uncompressed** TSV, not of the `.gz`
— see the companion note on gene identifier sources.

A first attempt at the mouse file was resumed with `curl -C -` after a timeout and produced a
112,195,895-byte file that `gzip -t` reported as "decompression OK, trailing garbage ignored" and
whose md5 did not match. It was discarded and re-fetched. A resumed download of a gzip stream can
decompress and still be wrong.

## 1. The per-species files partition at the pair level, and Compara says so

Cross-species rows among the lab's three species, counted in each file:

| File | Total rows | → human | → mouse | → worm |
|---|---|---|---|---|
| **human** | 3,878,214 | — | **0** | **23,982** |
| **mouse** | 4,522,853 | **23,764** | — | **25,006** |
| **worm** | 4,517,796 | **0** | **0** | — |

The human file holds no human↔mouse row at all. The human and mouse files together carry all three
pairings among the lab's species, and the worm file carries none of them.

This is documented, not accidental. `README.gene_trees.tsv_dumps.txt` in each directory:

> To eliminate redundancy, each genome-specific homology TSV file contains an arbitrary subset of
> orthologies involving the given genome. To access all available orthologies between two genomes
> […] you will need to download the genome-specific files of both genomes.

"Arbitrary" is the publisher's own word: which file a pair lands in is not a promise. What the counts
add is that the split is total at the pair level — a pair is wholly present or wholly absent, never
partial — so a zero count is a reliable signal that the rows are in the other file.

**Every cross-species row is an ortholog row.** In the human file, by type:

| `homology_type` | Rows | Species |
|---|---|---|
| `ortholog_one2one` | 2,421,012 | cross-species |
| `ortholog_one2many` | 885,783 | cross-species |
| `ortholog_many2many` | 430,246 | cross-species |
| `other_paralog` | 128,020 | same-species |
| `within_species_paralog` | 13,144 | same-species |
| `gene_split` | 9 | same-species |

The split is clean: no `ortholog_*` row is same-species, and no paralogy or `gene_split` row is
cross-species. In release 116 there are **zero between-species paralog rows** for human↔mouse or
human↔worm; the paralogs that arrive in these files are paralogs within one species.

## 2. Cardinality collapses between mouse and worm

19,670 distinct human genes appear in the human file — the human genes Compara places in a protein
gene tree, and the denominator used throughout this section.

| | vs mouse | vs worm |
|---|---|---|
| ortholog rows | 23,764 | 23,982 |
| — `ortholog_one2one` | 16,335 | 2,764 |
| — `ortholog_one2many` | 2,181 | 3,073 |
| — **`ortholog_many2many`** | 5,248 (22.1%) | **18,145 (75.7%)** |
| human genes with ≥1 ortholog | 17,826 | 8,056 |
| human genes with ≥1 `ortholog_one2one` | 16,335 | 2,764 |
| **as a share of all 19,670** | **83.0%** | **14.1%** |
| as a share of those with any ortholog | 91.6% | 34.3% |
| `is_high_confidence` = 1 | 17,233 of 23,764 | 2,979 of 23,982 |

The earlier survey reported 83.1% and 14.1%. The exact fractions are 16,335/19,670 = 83.05% and
2,764/19,670 = 14.05%, so the mouse figure is 83.0% at one decimal place, not 83.1%. Nothing else
about it moved.

**Both quality scores are null in every worm row.**

| Pair | Rows | `goc_score` null | `wga_coverage` null |
|---|---|---|---|
| human↔mouse | 23,764 | 0 | 0 |
| human↔worm | 23,982 | **23,982 (100%)** | **23,982 (100%)** |
| mouse↔worm | 25,006 | **25,006 (100%)** | **25,006 (100%)** |

Human↔mouse has both scores on every row; both worm pairings have neither on any row. A filter that
requires either score empties a worm answer completely and reports nothing about why.

## 3. Compara is the only surveyed source that writes cardinality down

The columns the dumps carry, from the same README: `identity`, `homology_type`,
`homology_identity`, `dn`, `ds`, `goc_score`, `wga_coverage`, `is_high_confidence`, `homology_id`.
`homology_type` is described as "homology type and cardinality (e.g. 'ortholog_one2one')" — it is a
property of the gene tree, present on the row as published.

NCBI's `gene_orthologs` carries a `relationship` column whose only value on the pairs read here is
`Ortholog`; the Alliance's file carries `Algorithms`, `AlgorithmsMatch`, `OutOfAlgorithms`,
`IsBestScore` and `IsBestRevScore`. Neither states a cardinality. A consumer of either derives one by
grouping, which makes it a property of the subset in hand.

## 4. Three publishers on human↔mouse: agreement concentrates in the one-to-one core

**Method, in full, because the number is sensitive to it.** Compara's human↔mouse ortholog rows are
already in Ensembl gene ids. NCBI's `gene_orthologs`, filtered to the 9606/10090 pairs, is in Entrez
GeneIDs; the Alliance's orthology file is in MOD ids (`HGNC:…`, `MGI:…`). Both were reduced to
(ENSG, ENSMUSG) through the Alliance's `GENECROSSREFERENCE_COMBINED` file — `NCBI_Gene:…` → MOD id →
`ENSEMBL:…` for NCBI, MOD id → `ENSEMBL:…` for the Alliance — keeping every published Ensembl
cross-reference and dropping only rows where one side has none. **Strict one-to-one** means Compara's
own `ortholog_one2one` label for Compara, and, for the other two, the subset of that source's pairs
in which the human gene appears in exactly one pair and the mouse gene in exactly one pair, taken in
the source's native id space before mapping.

| | Native pairs | Mapped pairs | Unmappable |
|---|---|---|---|
| Compara | 23,764 | 23,764 | — (already Ensembl ids) |
| NCBI | 17,096 | 18,746 | 191 rows (**1.1%**) |
| Alliance | 24,592 | 29,177 | 259 rows (**1.1%**) |

| | Pairs | Share |
|---|---|---|
| union of all asserted pairs | 32,772 | |
| **claimed by exactly one source** | **10,249** | **31.3%** |
| claimed by exactly two | 6,131 | 18.7% |
| claimed by all three | 16,392 | 50.0% |
| — unique to Compara | 3,474 of its 23,764 | |
| — unique to NCBI | 99 of its 18,746 | |
| — unique to the Alliance | 6,676 of its 29,177 | |

| Jaccard | whole sets | strict one-to-one | one-to-one, on genes all three cover |
|---|---|---|---|
| Compara vs NCBI | 0.629 | 0.842 | 0.974 |
| Compara vs Alliance | 0.620 | 0.833 | 0.982 |
| NCBI vs Alliance | 0.636 | 0.907 | 0.968 |

The last column restricts both sides to the 16,456 human and 16,435 mouse genes that all three
sources place in some pair. On that shared universe the union is 16,496 pairs and **31 of them —
0.2% — are claimed by one source only**. The disagreement is not spread through the data: it is
almost entirely in genes one source covers and another does not, and in the non-one-to-one tail.

**This is where the re-measurement differs from the survey it replaces.** That survey reported 24.6%
of pairs claimed by exactly one source and a strict one-to-one Jaccard of 0.91–0.93. Re-derived here,
the figures are **31.3%** and **0.83–0.91**. The recipe that produced 24.6% is gone and could not be
recovered, so the two cannot be reconciled; what is above is the recipe stated in full and the
numbers it gives. The gap is method, and the method is sensitive — four defensible variants:

| Variant | Exactly one source | One-to-one Jaccard |
|---|---|---|
| A: every Ensembl cross-reference kept (**above**) | 31.3% | 0.83–0.91 |
| B: only genes with exactly one Ensembl cross-reference | 31.8% | 0.87–0.92 |
| C: Alliance gated on `IsBestScore` and `IsBestRevScore` | 30.1% | 0.82–0.92 |
| D: B and C together | 34.0% | 0.86–0.94 |

None reaches 24.6%; the one-to-one Jaccard brackets 0.91–0.93 from below in three of the four. What
survives unchanged is the direction and its size: whole-set agreement is around 0.62–0.64, one-to-one
agreement is 0.83–0.94, and on a shared gene universe one-to-one agreement is 0.97–0.98.

**All of these are lower bounds on agreement.** They lean on a third-party reconciliation of ids —
the Alliance cross-reference file — and 1.1% of NCBI's rows and 1.1% of the Alliance's are dropped
because one side carries no Ensembl cross-reference at all. Under variant B those losses rise to
7.3% and 15.8%. A dropped row cannot agree with anything.

## 5. What pins, and what has stopped moving

**Compara's dumps go back to release-90 and the checksum beside them does not.** Probed per release
for `homo_sapiens`:

| Release | homology file | `MD5SUM` | `CHECKSUMS` |
|---|---|---|---|
| 88, 89 | absent | — | — |
| 90 | present | present | — |
| 95, 99, 100, 105, 110, 111, 112 | present | — | — |
| 113 | present, **uncompressed `.tsv`** | — | present |
| 114, 115 | present | — | present |
| 116 | present | present | present |

So the files pin back to release-90, but "an MD5 beside every file" holds only at release-90 and
release-116; releases 91 through 112 publish neither checksum file, and release-113 ships the dumps
uncompressed, under a different filename from every release around it.

Sources that have stopped moving, each verified from the publisher's own bytes:

- **NCBI HomoloGene is retired.** Its FTP `README`, last updated 2024-02-01, states "The HomoloGene
  website has been retired, replaced by the NCBI Orthologs site" and cites an announcement dated
  2024-01-30. The last build, `build68/`, is dated 2014-04-15 in the directory listing, and
  `build68/homologene.data` carries `Last-Modified: Tue, 06 May 2014 15:12:54 GMT`. The current
  symlink path is a 404.
- **WormBase WS298 is WormBase's final release**, announced on the Alliance community forum at
  <https://community.alliancegenome.org/t/announcing-the-final-release-of-wormbase/8461> — HTTP 200
  to a plain unauthenticated `curl`, 21,978 bytes. The WormBase blog copy is Cloudflare-blocked, and
  so is `downloads.wormbase.org` itself, which answers 403.
- **HCOP's bulk files are not fetchable from either of its locations.** The EBI path 404s and the
  `hgnc/hcop` prefix of HGNC's bucket lists nothing.

## 6. What this note could not verify

- **The 24.6% figure and the 0.91–0.93 Jaccard range** of the earlier survey. Section 4 reports what
  today's re-derivation gives and the four variants around it; the original recipe is lost.
- **HCOP "stale since 2025-05-16"** is carried from the earlier survey. There is no longer a file to
  read a date from.
- **OrthoList 2 frozen since 2018** is carried from the earlier survey, unverified here.
- **Whether the human/mouse split of section 1 is stable across releases** was not tested — only
  release 116 was counted at the pair level. The README's word "arbitrary" is the publisher's own
  and is the reason it was not assumed.
