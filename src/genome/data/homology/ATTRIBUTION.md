# Attribution — the homology provenance table

**No homology data ships here.** This directory holds one small table, `homology_metadata.tsv`, that
says where a species pair's homologies are published and what those bytes hash to — so naming a pair
is enough to prepare it, and a **Homology set** is downloaded rather than redistributed. Everything
the package answers with is the publisher's assertion, fetched from the publisher and read locally.

## Ensembl Compara — the only publisher here

Herrero J, Muffato M, Beal K, Fitzgerald S, Gordon L, Pignatelli M, Vilella AJ, Searle SM, Amode R,
Brent S, Spooner W, Kulesha E, Yates A, Flicek P. **Ensembl comparative genomics resources.**
*Database (Oxford)* 2016:bav096. PMID 26896847. doi:10.1093/database/bav096.

The files read are the per-species gene-tree homology TSV dumps of the `protein_default` collection,
release **116**, from <https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/>. The
ncRNA dumps and the whole-clade collection dumps beside them are deliberately not read: they are
different objects with different membership, not more rows of this one.

| pair | file that holds it | rows | md5 |
|---|---|---|---|
| *C. elegans* ↔ *H. sapiens* | `homo_sapiens` | 23,982 | `59857f48…` |
| *C. elegans* ↔ *M. musculus* | `mus_musculus` | 25,006 | `8f9870f0…` |
| *H. sapiens* ↔ *M. musculus* | `mus_musculus` | 23,764 | `8f9870f0…` |

**Cite Ensembl Compara for every answer.**
`HomologyMetadata.attribution()` renders the line to print.

## Why release 116, and why the URL is written down rather than built

Compara's dumps pin cleanly back to release 90, and they are not uniformly pinnable. Measured
against the published listings: an `MD5SUM` file exists only at releases **90 and 116**; a
`CHECKSUMS` file exists at **113–116**; releases 91–112 publish neither. Release **113 ships these
dumps uncompressed** — `.tsv`, not `.tsv.gz` — so a URL built from a template is a live 404 there.
Release 116 is the newest, has a published checksum and has the compressed naming, which is why it
is the release pinned. Each row carries the URL read off that release's own listing, and a row with
no md5 is refused as it is read: that refusal is what stops a release with no checksum from being
pinned by accident.

The checksum is checked against the bytes **as they are fetched**, and the derived slice's own
sha256 goes into the set's **Completion marker** and is re-checked on every read. Both are
load-bearing: a resumed download of one of these gzips has been seen to pass `gzip -t` —
"decompression OK, trailing garbage ignored" — while its md5 was wrong. Opening cleanly is not
evidence.

## The partition, counted

Compara's per-species files are a de-duplicated partition **at the pair level**, which its own
README states: each file holds "an arbitrary subset of orthologies involving the given genome", and
to get everything between two genomes you must take both files. Counted on release 116:

| file | rows against human | against mouse | against worm |
|---|---|---|---|
| `homo_sapiens` | — | **0** | 23,982 |
| `mus_musculus` | 23,764 | — | 25,006 |
| `caenorhabditis_elegans` | **0** | **0** | — |

So the human and mouse files together cover all three pairings among the lab's species and the worm
file is needed for none of them. Which file holds a pair is arbitrary and **is not promised stable
across releases** — release 110's human file held 16,242 mouse rows where 116's holds none — so the
`holding_species` column is a measurement of one release and is verified again every time a set is
prepared. A pair whose slice comes back empty raises and names the other file; a pair is never
*partially* present, which is what makes zero a trustworthy signal.

## What the release actually contains, where it differs from what one might assume

- **No cross-species paralogy.** Counted over the whole human dump — 4.0 M rows across roughly 200
  partner species — there are **zero** `between_species_paralog` rows, and every `other_paralog`
  (128,020), `within_species_paralog` (13,144) and `gene_split` (9) row relates two genes of *one*
  species. A **Homology link** relates two species, so those rows are not a pair's and are not in
  its set. `paralogs=True` is kept as the place a cross-species duplication label would land, and on
  release 116 it changes nothing for human, mouse and worm.
- **Both quality scores are null for either worm pairing**, not only for human↔worm: `goc_score` and
  `wga_coverage` are `NULL` on 100% of the 23,982 human↔worm rows and 100% of the 25,006 mouse↔worm
  rows. An answer says which score columns its set holds nothing in, measured over the slice, so a
  filter written on one is told rather than left to empty in silence.
- **Cardinality is far from one-to-one for worm.** Of the human↔worm rows, 75.66% are
  `ortholog_many2many` and one-to-one reaches 14.1% of human genes, against 83.0% for mouse. That
  spread is the argument for **Homology type** being a typed fact on every link rather than a
  footnote, and it is why nothing this package publishes is derived through homology (ADR-0019).

## What this package adds, and does not

Nothing. No quality score, no ranking and no "best ortholog" of this package's own exists; the
**Homology type** is Compara's tree-derived label and is never recomputed, not after a filter and
not after resolution into an **Annotation** (ADR-0020). No **TF gene table**, no **Cofactor table**
and no list this package ships is derived through homology, and no answer is ever silently
species-mapped (ADR-0019).
