# Attribution — the shipped cofactor tables

External published knowledge, redistributed here so that asking whether a gene is a **Transcription
cofactor** needs no network. **Membership and classification are the publishers' and none of them
are this package's — cite the publisher whose table you used.** The two provenance tables beside the
data repeat these facts machine-readably: `cofactor_metadata.tsv` is keyed by species and carries the
taxid, the file and the sha256 of the **unpacked** TSV inside each gzip (ADR-0006);
`cofactor_source_metadata.tsv` is keyed by species *and* source and carries the publisher, version,
PubMed id and URL. Two tables and not one, because one row cannot describe a table built from several
publishers, and joining them positionally inside a cell is the shape that breaks quietly.

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
recorded here so that nobody files it as a defect. `tests/test_tf_cofactor.py` pins it.

## What the columns are, and whose vocabulary each one is

Four uniform columns lead every table — `gene_id_stem`, `symbol`, `is_cofactor`, `source` — and
everything after them is one publisher's own column under a namespaced snake_case name, so a table
built from two publishers is more columns and one more provenance row rather than a change of
format. `source` is a closed vocabulary (`animaltfdb`, `epifactors`, `both`) validated as the file
is read; it asserts agreement on **membership only** and never on classification, and two
publishers' vocabularies are never crosswalked (ADR-0014). A blank cell is that publisher recording
nothing and reads back as `None`.

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
number of genes the list carries — 970 across 81 summary keys for mouse, 317 across 56 for worm.
Mouse needs all five rules; worm needs only `Other_Co-activator/repressors` and `MYC`. The build
re-runs that arithmetic and **fails loudly** both when a family survives the map with no category
and when the two numbers stop reconciling, so a release that renames or re-counts a family breaks
the build rather than quietly blanking a column. The shipped `animaltfdb_family` is the gene list's
own spelling, never the summary's: the map exists to find each gene's category and not to re-spell
what the publisher classified it under.

## Rebuilding

`python scripts/build_tf_cofactor.py animaltfdb4_mouse|animaltfdb4_worm <the species cofactor list>
<cof_info_summary.tsv>` — the generator lives outside the wheel, takes file paths and downloads
nothing, fails loudly when a publisher re-spells a column, and writes byte-stable output, so a
re-release is a re-run and a reviewable diff. **No test proves a shipped table still matches
AnimalTFDB, and none can**: regenerating one needs a download and CI has no network, the limitation
ADR-0011 already accepts. The counts above are pinned in `tests/test_tf_cofactor.py`, turning silent
drift into a loud failure — the most that is available.
