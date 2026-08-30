# CLI: lookup commands

Six commands that answer from a published table rather than from one genome's files. `tf`
names an assembly, `xref` and `homology` name a species, and `motif` names a FASTA. Every
answer says which publisher and which release produced it.

```console
$ genome xref ids "Homo sapiens" --from-stems symbol ENSG00000141510
gene id stems -> symbol ids for Homo sapiens (hgnc 2026-07-07)
  source  https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-07-07.txt
  1 resolved, 1 symbol ids, 0 this release names none for
ENSG00000141510	TP53
```

The tables land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

Each synopsis below is its group's own `--help`, so it leaves the leading `genome` off:
`xref ids` is run as `genome xref ids`. `homology` and `motif` carry one command each, which
their synopses fold into the group name; they are run as `genome homology links` and `genome
motif scan`. `--json`, exit codes and the split between stdout and stderr are on the
[CLI overview](index.md).

## genome tf

`gene-list` prints the gene ids a published census judges transcription factors, and
`cofactor-list` the ids a publisher lists as cofactors. Which census or table answers is
decided by the species the assembly's own metadata row names, never by anything you pass.
Both resolve the publisher's gene id stems into the ids the annotation actually spells, so
the output joins straight to a counts matrix. Who published the verdict, and what it does
and does not cover, is on [Transcription factors](../topics/transcription-factors.md).

```console
$ genome tf cofactor-list ce11 > cofactors.txt
TF cofactors for ce11 / wormbase_ws298 (Caenorhabditis elegans)
  AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Caenorhabditis_elegans_Cof
  317 cofactors, 317 gene ids, 0 stems this annotation carries no gene for
```

**A worm assembly is answered by `cofactor-list` and refused by `gene-list`**, because
AnimalTFDB assessed worm cofactors and no publisher has released a worm TF census. The
message names the species that do have one. An assembly that is not registered here, and one
whose species nothing names, exit `1` too, each with its own message.

The two commands, with every argument and option:

::: mkdocs-typer2
    :module: genome.tf.cli
    :name: tf
    :engine: native

## genome xref

`ids` converts identifiers to and from gene id stems, which is how a column of Entrez
GeneIDs from a GEO series or UniProt accessions from a mass-spec run reaches the rest of
this package. The direction is named and never inferred: `--to-stems NAMESPACE` reads the
ids as that namespace, `--from-stems NAMESPACE` answers in it, and naming neither or both
exits `2`. What the namespaces are and which source answers by default is on
[Gene identifiers](../topics/gene-identifiers.md).

```console
$ genome xref ids "Homo sapiens" --to-stems hgnc HGNC:11998 HGNC:13666 HGNC:10041 > stems.tsv
hgnc ids -> gene id stems for Homo sapiens (alliance 9.0.0)
  source  https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz
  2 resolved, 3 gene id stems, 1 this release names none for

$ cat stems.tsv
HGNC:11998	ENSG00000141510
HGNC:13666	ENSG00000094914
HGNC:13666	ENSG00000291836
HGNC:10041
```

**Every id you passed gets at least one row.** One naming two genes prints both rather than
whichever came first, and one this release names nothing for gets a row with an empty second
column, so nothing you asked about goes missing from the output.

`symbols` is the separate command for going the other way from a gene symbol. It matches
approved, previous and alias spellings, and each row says which kind matched, which is why
it is not a third direction of `ids`:

```console
$ genome xref symbols "Homo sapiens" ARNTL ADCY3 Brca1 > genes.tsv
gene symbols -> gene id stems for Homo sapiens (hgnc 2026-07-07)
  source   https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-07-07.txt
  columns  asked, symbol, gene_id_stem, kind
  matching exact, on approved, previous, alias spellings
  2 resolved, 3 matches, 1 this release matched nothing for

$ cat genes.tsv
ARNTL	ARNTL	ENSG00000133794	previous
ADCY3	ADCY3	ENSG00000138031	approved
ADCY3	ADCY3	ENSG00000155897	previous
Brca1
```

`ARNTL` is a spelling HGNC retired and it still reaches its gene. Matching is exact by
default, because the species already fixes the authority; `--case-insensitive` folds both
sides and still answers with every gene it matched.

The two commands, with every argument and option:

::: mkdocs-typer2
    :module: genome.xref.cli
    :name: xref
    :engine: native

## genome homology

`links` prints the genes of another species that a gene id stem's gene is homologous to, on
Ensembl Compara's own gene trees. Any pairing among human, mouse and worm answers, either
way round, off one file fetched once and read locally. What the labels mean and how to put
an answer back into an annotation's own gene ids is on [Homology](../topics/homology.md).

```console
$ genome homology links "Caenorhabditis elegans" "Homo sapiens" \
      WBGene00020462 WBGene00008317 WBGene00008352 > homologs.tsv
Caenorhabditis elegans -> Homo sapiens orthologs (Ensembl Compara 116)
  source   Ensembl Compara release 116 (PMID 26896847) — https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz
  columns  gene_id_stem, homolog_gene_id_stem, homology_type, is_ortholog, is_high_confidence, goc_score, wga_coverage
  2 resolved, 2 links, 1 this release names no homolog for, 0 dropped partners
  quality  goc_score and wga_coverage null on every link of this set, so a filter on either empties rather than narrowing

$ cat homologs.tsv
WBGene00020462	ENSG00000177479	ortholog_one2one	True	True	NULL	NULL
WBGene00008317	ENSG00000164074	ortholog_one2one	True	True	NULL	NULL
WBGene00008352
```

Every cell is the publisher's, printed verbatim and never recomputed. **An empty row is not
`NULL`**: an empty row is a gene this release names no homolog for, and `NULL` is Compara's
own word for a cell it recorded nothing in on a link that does exist. Orthologs are the
answer by default and `--paralogs` returns every link the publisher wrote, marked by its own
`homology_type` rather than filtered out.

The command, with every argument and option:

::: mkdocs-typer2
    :module: genome.homology.cli
    :name: homology
    :engine: native

## genome motif

`scan` reads a FASTA, scores every sequence against a JASPAR release, writes the hits to
Parquet and prints a summary of the run. It is the batch case and the only motif command.
Listing, plotting and comparing motifs are notebook work, on
[Motifs](../topics/motifs.md).

```console
$ genome motif scan peaks.fa hits.parquet --release 2024
scanned 500 sequences with 781 motifs from JASPAR 2024 vertebrates
  background  0.298, 0.198, 0.198, 0.306
  threshold   0.0001
  skipped     98 under 7 positions, so not scanned: MA0004.1, MA0130.1, MA0151.1, …
  workers     14
  hits        22251 -> hits.parquet
```

The skipped list is trimmed here; the command names all 98. A motif under seven positions
cannot reach the default threshold at all, so it is named rather than scanned at some looser
cutoff you did not ask for.

The hits go to the named file and the summary to stdout, so `--json` is never corrupted by
table data. Read the file back with `genome.tf.motif.read_hits`, which restores the compact
dtypes and the provenance that says what the scan was; `pandas.read_parquet` gives the rows
and drops both. **`--background` decides the answer more than any other option**, and
whichever mode was used is recorded in the summary and on the hits.

The command, with every argument and option:

::: mkdocs-typer2
    :module: genome.tf.motif.cli
    :name: motif
    :engine: native
