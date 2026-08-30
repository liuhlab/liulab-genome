# Attribution — the shipped cofactor tables

External published knowledge, redistributed here so that asking whether a gene is a **Transcription
cofactor** needs no network. **Every classification is its publisher's — cite the publisher whose
table you used.** Membership is the publisher's too for mouse and worm; **human membership is this
package's own**, the union of two lists neither of them publishes (ADR-0016), and the section below
says what that costs. The two provenance tables beside the data repeat these facts
machine-readably: `cofactor_metadata.tsv` is keyed by species and carries the taxid, the file and the
sha256 of the **unpacked** TSV inside each gzip (ADR-0006); `cofactor_source_metadata.tsv` is keyed
by species *and* source and carries the publisher, version, PubMed id and URL. Two tables and not
one, because one row cannot describe a table built from several publishers, and joining them
positionally inside a cell is the shape that breaks quietly.

## Homo sapiens — `homo_sapiens.cofactor_table.tsv.gz`

**A union this package publishes, and the only place in `genome.tf` where the verdict is ours.**
Three sources, each cited in its own row of `cofactor_source_metadata.tsv`:

| source | contributes | rows |
|---|---|---|
| **AnimalTFDB 4.0** | membership, family, category | 1,024 genes |
| **EpiFactors v2.0** | membership, function, target, modification, complex | 796 genes, 801 rows |
| **HGNC**, one pinned dated monthly archive | the Ensembl gene id of every EpiFactors row | no membership |

Shen WK, Chen SY, Gan ZQ, Zhang YZ, Yue T, Chen MM, Xue Y, Hu H, Guo AY. **AnimalTFDB 4.0.**
*Nucleic Acids Research* 51(D1):D39-D45, 2023. PMID 36268869 — the human cofactor list, from
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Homo_sapiens_Cof>.

Marakulina D, Vorontsov IE, Kulakovskiy IV, Lennartsson A, Drabløs F, Medvedeva YA. **EpiFactors
2022: expansion and enhancement of a curated database of human epigenetic factors and complexes.**
*Nucleic Acids Research* 51(D1):D564-D570, 2023. PMID 36350659. doi:10.1093/nar/gkac989 — the
**EpiFactors v2.0** main gene table, from
<https://epifactors.autosome.org/public_data/v2.0/EpiGenes_main.csv>.

Seal RL, Braschi B, Gray K, McClay J, Tweedie S, Bruford EA. **Genenames.org: the HGNC and PGNC
resources in 2026.** *Nucleic Acids Research* 54(D1):D1098-D1107, 2026. PMID 41287213.
doi:10.1093/nar/gkaf1229 — the monthly archive dated **2026-08-07**, from
<https://storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/hgnc_complete_set_2026-08-07.txt>.

**1,466 genes: 354 both publishers list, 670 AnimalTFDB alone, 442 EpiFactors alone.** Across the
whole table, 85 AnimalTFDB families in 6 categories, and EpiFactors' three published vocabularies
at 19 functions, 7 targets and 24 modifications once a multi-valued cell is split.

**Membership is unioned; classification is not.** A gene either publisher lists is a row. Its
AnimalTFDB columns are filled only if AnimalTFDB listed it and its EpiFactors columns only if
EpiFactors did, so a blank group is a publisher that never named the gene. Nothing is inferred
across the two in either direction for any pair of values (ADR-0014), and a `source` of `both` says
the two agree that the gene belongs and says nothing about how either of them classified it.

**EpiFactors joins to Ensembl through its HGNC id, never through its symbol.** All 801 of its rows
carry one and all 801 resolve; 31 of them still name the gene by a symbol HGNC has retired —
`ACINU` for `ACIN1`, `ARNTL` for `BMAL1`, `C11orf30` for `EMSY` — so a symbol match would key those
genes wrongly or drop them. The archive is a **pinned dated file** and never the rolling current
one, so that the 442 stems only HGNC can supply are reproducible; the archive's dates are irregular,
so the pin names a file read from its listing rather than one constructed from a date. The `symbol`
column is HGNC's currently approved spelling on every human row, reached through the id.

**Five genes have two EpiFactors rows each and ship as one.** `ALKBH1`, `HSPA1A`, `HSPA1B`, `NAT10`
and `PTBP1`. They are not duplicate rows — one carries a histone-modification annotation and the
other an RNA-modification one — and a table is one row per stem, so their cells are unioned and
deduplicated within a cell as well as across the two rows. **The cost is that for those five the
pairing between a function and its own modification is lost**: `ALKBH1` ships
`Histone modification;RNA modification;DNA modification` beside
`DNA demethylation;RNA demethylation`, and nothing in the row says which demethylation belongs to
which function. Go back to EpiFactors' own two rows if that pairing is what you need.

**The TF list and the cofactor list overlap: 151 human genes are both.** A Lambert-positive **TF
gene** and a **Transcription cofactor** — 57 of them from the AnimalTFDB side, 122 from the
EpiFactors side, 28 from both — so **a caller who unions `tf_gene_list` with `tf_cofactor_list`
double-counts those 151 genes.** `TBP`, `KMT2A` and `DNMT1` are among them. Being a cofactor never
suppresses the motifs a census already reached; the two are independent verdicts about different
questions, and `tests/tf/test_tf_cofactor.py` pins the number so it cannot drift unnoticed.

## Mus musculus — `mus_musculus.cofactor_table.tsv.gz`

Shen WK, Chen SY, Gan ZQ, Zhang YZ, Yue T, Chen MM, Xue Y, Hu H, Guo AY. **AnimalTFDB 4.0: a
comprehensive animal transcription factor database updated with variation and expression
annotations.** *Nucleic Acids Research* 51(D1):D39-D45, 2023. PMID 36268869. doi:10.1093/nar/gkac907
— the **AnimalTFDB 4.0** mouse cofactor list, from
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Mus_musculus_Cof>, keyed by
Ensembl mouse gene ids on GRCm39 (Ensembl 105), unversioned as published.

**970 genes, 84 families, 6 categories.** Every gene listed is one AnimalTFDB accepts: it publishes
no rejected set and says nothing whatever about the genes it left out, so `is_cofactor` reads `yes`
on every row.

## Caenorhabditis elegans — `caenorhabditis_elegans.cofactor_table.tsv.gz`

The same publisher, release and paper — the **AnimalTFDB 4.0** worm cofactor list, from
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Caenorhabditis_elegans_Cof>,
keyed by WormBase `WBGene…` gene ids.

**317 genes, 57 families, 6 categories.**

## Worm has cofactors and no TF census, and that is not a bug

`tf_cofactor_list` answers for a worm assembly while `tf_gene_list` **raises** for the same
assembly. AnimalTFDB assessed *C. elegans* cofactors; no publisher has released a *C. elegans*
transcription-factor census, and this package publishes no census of its own. The asymmetry is the
publishers' shape rather than a gap in the reader — absence is theirs and not ours — and it is
recorded here so that nobody files it as a defect. `tests/tf/test_tf_cofactor.py` pins it.

## What the columns are, and whose vocabulary each one is

Four uniform columns lead every table — `gene_id_stem`, `symbol`, `is_cofactor`, `source` — and
everything after them is one publisher's own column under a namespaced snake_case name, so a table
built from two publishers is more columns and one more provenance row rather than a change of
format. Human carries `animaltfdb_family`, `animaltfdb_category`, `epifactors_function`,
`epifactors_target`, `epifactors_modification` and `epifactors_complex_name`; mouse and worm carry
the two AnimalTFDB columns and nothing else. `source` is a closed vocabulary (`animaltfdb`,
`epifactors`, `both`) validated as the file is read; it asserts agreement on **membership only** and
never on classification, and two publishers' vocabularies are never crosswalked (ADR-0014). A blank
cell is that publisher recording nothing and reads back as `None`. A cell holding more than one
value spells them apart with `;`, the separator every multi-valued cell in this package uses —
`interpro_ids` in a **TF gene table**, a **Motif link**'s partners — and the build refuses a
publisher's value that already contains one rather than writing a cell that splits into values
nobody published.

The provenance table's own `source` column is a *wider* vocabulary than the shipped tables':
`animaltfdb`, `epifactors`, `hgnc`. `both` is a fact about a row of a table and describes no
publisher, so it never appears there; `hgnc` lists no gene and so never appears in a table, but the
stems of 442 human genes exist only because it said so, which is a contribution to cite rather than
an implementation detail to bury.

`is_cofactor` reads `yes` on every row that ships today and is kept anyway. Dropping it would make
presence in the file the verdict, at which point a future source could not record a rejection
without a format change.

## The family and category come from two of the publisher's own files

AnimalTFDB publishes membership and `Family` in the per-species cofactor list and the category each
family sits in in a separate summary,
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/cof_info_summary.tsv> — and its own two files
spell five families differently. The build joins them through five hand-written rules:

| gene-list family | summary family |
|---|---|
| `Lysine methyltransferase` | `Lysine methyltransferase family` |
| `Histone lysine methyltransferase` | `Lysine methyltransferase family` |
| `Other_Co-activator/repressors` | `Other_Co-activator_repressors` |
| `MYB` | `Others` |
| `MYC` | `Others` |

Each rule is proved by the publisher's own arithmetic: mapping every family the gene list uses and
summing the summary's `family_count` over the *distinct* families that result gives exactly the
number of genes the list carries — 1,024 across 82 summary keys for human, 970 across 81 for mouse,
317 across 56 for worm. Human and mouse need all five rules; worm needs only
`Other_Co-activator/repressors` and `MYC`. The build
re-runs that arithmetic and **fails loudly** both when a family survives the map with no category
and when the two numbers stop reconciling, so a release that renames or re-counts a family breaks
the build rather than quietly blanking a column. The shipped `animaltfdb_family` is the gene list's
own spelling, never the summary's: the map exists to find each gene's category and not to re-spell
what the publisher classified it under.

## Rebuilding

`python scripts/build_tf_cofactor.py animaltfdb4_mouse|animaltfdb4_worm <the species cofactor list>
<cof_info_summary.tsv>`, and for human `python scripts/build_tf_cofactor.py union_human
Homo_sapiens_Cof cof_info_summary.tsv --epifactors EpiGenes_main.csv --hgnc
hgnc_complete_set_2026-08-07.txt`. The generator lives outside the wheel — **every curation rule is
in it and none of it is in the package** — takes file paths and downloads nothing, fails loudly when
a publisher re-spells a column, when a family survives the spelling map with no category, and when
an identifier does not resolve, and writes byte-stable output, so a re-release is a re-run and a
reviewable diff. **No test proves a shipped table still matches its publishers, and none can**:
regenerating one needs downloads and CI has no network, the limitation ADR-0011 already accepts. The
counts above are pinned in `tests/tf/test_tf_cofactor.py`, turning silent drift into a loud failure —
the most that is available.
