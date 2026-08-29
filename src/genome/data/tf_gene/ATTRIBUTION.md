# Attribution — the shipped TF gene tables

Every file in this directory is external published knowledge, redistributed here so
that asking which genes a census judges transcription factors needs no network and no
download step. **The verdicts are the publishers' and none of them are this package's.**
Cite the publisher whose census you used.

`census_metadata.tsv` carries the same facts per census, machine-readable: species,
NCBI taxid, file, publisher, version, PubMed id, which of the publisher's columns the
**DBD family** was taken from, the source URL, and the sha256 of the **unpacked** TSV
inside the gzip (ADR-0006).

## Homo sapiens — `homo_sapiens.tf_gene_table.tsv.gz`

Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR,
Weirauch MT. **The Human Transcription Factors.** *Cell* 172(4):650-665, 2018.
PMID 29425488. doi:10.1016/j.cell.2018.01.029

Database extract v_1.01, downloaded from
<https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv>.

2,765 genes assessed, 1,639 judged transcription factors. Thirteen of the release's
twenty-nine columns ship; the identifier blobs keyed to 2018 Ensembl protein and
transcript ids, the free-text curation notes, the curator names and the 2018
cross-references unusable without them are dropped. Three of the release's identifiers
are UniProt entry names rather than Ensembl gene ids — `ZNF73_HUMAN`, `DUX1_HUMAN`,
`DUX3_HUMAN` — and ship as the publisher wrote them, so they resolve against no
annotation and are reported rather than dropped.

## Mus musculus — `mus_musculus.tf_gene_table.tsv.gz`

Shen WK, Chen SY, Gan ZQ, Zhang YZ, Yue T, Chen MM, Xue Y, Hu H, Guo AY.
**AnimalTFDB 4.0: a comprehensive animal transcription factor database updated with
variation and expression annotations.** *Nucleic Acids Research* 51(D1):D39-D45, 2023.
PMID 36268869. doi:10.1093/nar/gkac907

AnimalTFDB 4.0 mouse TF list, downloaded from
<https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF>.
Keyed by Ensembl mouse gene ids on GRCm39 (Ensembl 105), unversioned as published.

1,611 genes, all of them judged transcription factors. **AnimalTFDB publishes no
rejected set**: it lists the genes it accepts and says nothing whatever about the rest,
so mouse has no assessed-negative rows rather than fabricated ones, and a mouse gene is
either in this file or unassessed. Four of the release's six columns ship, which is
every judgement AnimalTFDB makes about a gene. Dropped: `Species`, the same word on all
1,611 rows and already recorded in `census_metadata.tsv`; and the two identifier
cross-references, `Protein` (semicolon-joined Ensembl 105 protein ids) and `Entrez_ID`
(NCBI's own id, already `NA` for 109 of the genes).

## The two family vocabularies are not crosswalked

Each census classifies its genes under its own publisher's **DBD family** vocabulary —
75 values under Lambert's `DBD`, 72 under AnimalTFDB's `Family` — and the two are
deliberately left un-harmonised (ADR-0014). Lambert's `ARID/BRIGHT` and AnimalTFDB's
`ARID` are **not** asserted equivalent here, and neither is any other near-match: an
equivalence nobody has checked is worse than two vocabularies that each say who spelled
them. Group by family within a species; never across two.

## Rebuilding

`python scripts/build_tf_census.py lambert2018 <the publisher's file>`, and
`python scripts/build_tf_census.py animaltfdb4_mouse <the publisher's file>` — the
generator lives outside the wheel, fails loudly when a publisher re-spells a column, and
writes byte-stable output, so a re-release is a re-run and a reviewable diff.

## Licensing

Neither census states a licence. Both are public downloads published as community
resources whose papers ask to be cited, and both are redistributed here with the
attribution above.
