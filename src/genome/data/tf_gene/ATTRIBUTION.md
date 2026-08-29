# Attribution — the shipped TF gene tables

External published knowledge, redistributed here so that asking which genes a census judges
transcription factors needs no network. **The verdicts are the publishers' and none of them are this
package's — cite the publisher whose census you used.** `census_metadata.tsv` repeats these facts
machine-readably, with the publisher's own **DBD family** column and the sha256 of the **unpacked**
TSV inside each gzip (ADR-0006).

## Homo sapiens — `homo_sapiens.tf_gene_table.tsv.gz`

Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT.
**The Human Transcription Factors.** *Cell* 172(4):650-665, 2018. PMID 29425488.
doi:10.1016/j.cell.2018.01.029 — database extract **v_1.01**, from
<https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv>.

2,765 genes assessed, 1,639 judged transcription factors. Thirteen of the release's twenty-nine
columns ship; the blobs keyed to 2018 Ensembl protein and transcript ids, the free-text curation
notes, the curator names and the cross-references unusable without them are dropped. Three published
identifiers are UniProt entry names rather than Ensembl gene ids — `ZNF73_HUMAN`, `DUX1_HUMAN`,
`DUX3_HUMAN` — and ship as written, so they resolve against no annotation and are reported.

## Mus musculus — `mus_musculus.tf_gene_table.tsv.gz`

Shen WK, Chen SY, Gan ZQ, Zhang YZ, Yue T, Chen MM, Xue Y, Hu H, Guo AY. **AnimalTFDB 4.0: a
comprehensive animal transcription factor database updated with variation and expression
annotations.** *Nucleic Acids Research* 51(D1):D39-D45, 2023. PMID 36268869. doi:10.1093/nar/gkac907
— the **AnimalTFDB 4.0** mouse TF list, from
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF>, keyed by Ensembl
mouse gene ids on GRCm39 (Ensembl 105), unversioned as published.

1,611 genes, all judged transcription factors. **AnimalTFDB publishes no rejected set**: it lists the
genes it accepts and says nothing whatever about the rest, so mouse has no assessed-negative rows
rather than fabricated ones, and a mouse gene is either in this file or unassessed. Four of its six
columns ship, which is every judgement it makes about a gene; dropped are `Species` (one word on all
1,611 rows) and the cross-references `Protein` and `Entrez_ID`, the latter already `NA` for 109.

## The two family vocabularies are not crosswalked

Each census classifies its genes under its own publisher's **DBD family** vocabulary — 75 values
under Lambert's `DBD`, 72 under AnimalTFDB's `Family` — deliberately left un-harmonised (ADR-0014).
`ARID/BRIGHT` and `ARID` are **not** asserted equivalent, nor is any other near-match: an equivalence
nobody has checked is worse than two vocabularies that say who spelled them. Group within a species.

## Rebuilding, and licensing

`python scripts/build_tf_census.py lambert2018|animaltfdb4_mouse <the publisher's file>` — the
generator lives outside the wheel, fails loudly when a publisher re-spells a column, and writes
byte-stable output, so a re-release is a re-run and a reviewable diff. Neither census states a licence;
both are public downloads published as community resources whose papers ask to be cited.
