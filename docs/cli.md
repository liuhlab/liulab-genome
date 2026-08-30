# CLI

The `genome` command is a thin wrapper over the Python API — every command calls the same
function you would call yourself, so a shell and a notebook take one code path.

Every command takes `--json` for machine-readable output. Errors go to stderr with a
non-zero exit code: `2` for invalid input, `1` for everything else (a missing native
tool, a failed download, a checksum mismatch, a registration that cannot be trusted).

```console
$ genome --help
$ genome version
2026.6.0
$ genome version --json
{"version": "2026.6.0"}
```

## `genome doctor`

Report which native tools are on `PATH`, and at what versions. Exits `1` naming the
install command if one is missing. A tool that is installed but rejects `--version` — the
UCSC binaries do — is listed all the same, since presence is what this answers.

```console
$ genome doctor
samtools: samtools 1.21 ...
faToTwoBit: installed; reports no version
twoBitInfo: installed; reports no version
```

STAR and chromap are not checked here: they are optional, and each one checks for itself
when you ask it to build an index.

## `genome register <assembly>`

Prepare an assembly on disk — fetch, verify, index, record — and print where it landed.

```console
$ genome register sacCer3
registered sacCer3 in /data/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

Running it again on a registered assembly downloads nothing.

| Option | Effect |
|---|---|
| `--source PATH_OR_URL` | Seed from your own FASTA instead of the pinned source |
| `--force` | Register again from scratch — the repair for a directory that raises |

A directory that cannot be trusted exits `1` naming the file and the repair:

```console
$ genome register sacCer3
error: /data/liulab_data/genome/sacCer3 disagrees with its .completion.json:
sacCer3.2bit: recorded 3145728 bytes, found 0. Something changed these files after
they were registered. Re-register it with `genome register sacCer3 --force`.
```

`--force` keeps an unpacked FASTA whose checksum still matches and rebuilds only the
derived files; it fetches again when it cannot prove that. Interrupting a *first*
registration leaves a directory that needs exactly this.

### Naming a chimera builds it

An assembly named after two or more registered assemblies, sorted and joined by `_`, is a
[chimera](genome.md#chimera-assemblies): it is concatenated from those components, merged
annotation included, and nothing is downloaded.

```console
$ genome register ce11_ecHT115
registered ce11_ecHT115 in /data/liulab_data/genome/ce11_ecHT115
  components  ce11, ecHT115
  sha256  0f3d6d5e…
  files   ce11_ecHT115.2bit, ce11_ecHT115.chrom.sizes, ce11_ecHT115.fa, ce11_ecHT115.fa.fai
  annotation  wormbase_ws298+refseq_rs_2025_06_26 — the components' own, merged and registered by this build
```

There is no flag for listing the parts — the name carries them. A component this machine
has not prepared exits `1` naming the command that prepares it, and components typed in
the wrong order exit `1` naming the canonical spelling.

## `genome register-annotation <assembly> <name>`

Register one of the annotations the table lists for an assembly: fetch, check the
checksum and the chromosome names, build the gffutils database, record it.

```console
$ genome register-annotation sacCer3 ensgene_v101
registered ensgene_v101 for sacCer3 in /data/liulab_data/genome/sacCer3/gtf/ensgene_v101
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz
  sha256  d3f33fbf97deef26e2495f709f1c5bb2e2e1bf1ce71fb80758c2c9de42ad7026
  files   ensgene_v101.db, ensgene_v101.gtf
  chromosomes checked — every name the GTF uses is one the assembly carries
```

## `genome register-gtf <assembly> <path> <name>`

The same, for a GTF the table does not list. Nothing is downloaded and no checksum is
compared against; `<name>` is what addresses it from then on. A `.gz` is decompressed.

```console
$ genome register-gtf sacCer3 ~/annotations/sacCer3.WS298.gtf wormbase_ws298
```

Both registration commands share these options:

| Option | Effect |
|---|---|
| `--force` | Register again from scratch; also the repair for a broken directory |
| `--no-check-chromosomes` | Register a GTF whose chromosome-name mismatch you accept |
| `--infer-genes` / `--infer-transcripts` | Reconstruct those features from exon lines — for a bare exon-level GTF only |

The last line of the output always says whether the chromosome names were checked, and
if not, why: because you stood the check down, or because the assembly is not registered
yet and there was no `chrom.sizes` to check against. Only the second is worth acting on.

## `genome annotations <assembly>`

What the table offers against what is registered here. Downloads and prepares nothing, so
it works for an assembly you have never registered.

```console
$ genome annotations hg38
annotations for hg38 in /data/liulab_data/genome/hg38
  gencode_v50  offered, not registered  GENCODE v50
  mine         registered, not offered
default: gencode_v50 — not registered here; register it with `genome register-annotation hg38 gencode_v50`
```

An annotation that is here but cannot be trusted reads as `broken`, with the problem and
the repair on the line below. Exit stays `0` — one broken annotation never hides the
others.

## `genome gene-list <assembly> <category>`

Print the gene ids an annotation puts in one [gene category](genome.md#which-genes-are-in-a-category),
one per line. Only the ids go to stdout, so the output pipes; the heading and the
per-source attribution go to stderr.

```console
$ genome gene-list ce11 rRNA > rrna.txt
rRNA for ce11 / wormbase_ws298
  wormbase_ws298  20
```

`--annotation NAME` asks about a registered annotation other than the assembly's default
one. `--json` carries the same answer with the sources kept apart:

```console
$ genome gene-list ce11_ecHT115 rRNA --json
{"assembly": "ce11_ecHT115", "annotation": "wormbase_ws298+refseq_rs_2025_06_26",
 "category": "rRNA", "gene_ids": ["WBGene00004512", …], "sources": [{"component": "ce11", …}]}
```

Exits `1` when the annotation is not registered here, when no curated gene list ships for
it, and when it declares categories but not this one — three different facts, each with
its own message. None of them prints an empty list of genes.

## `genome gene-categories <assembly>`

Which categories that annotation declares, and how many genes are in each — what
`genome gene-list` may be asked for. A merged annotation shows the per-component split.

```console
$ genome gene-categories ce11_ecHT115
categories for ce11_ecHT115 / wormbase_ws298+refseq_rs_2025_06_26
  rRNA  39  (ce11: 23, ecHT115: 16)
```

`--json` emits every category with its gene ids and its sources — the same answer
`genome gene-list` gives for one of them, for all of them at once.

## `genome tf-gene-list <assembly>`

Print the gene ids a published census judges transcription factors, one per line — Lambert
et al. 2018 for human, AnimalTFDB 4.0 for mouse, chosen by the species the assembly's own
metadata row names and never by anything you pass. Only the ids go to stdout, so the
output pipes; the heading and the census's attribution go to stderr.

```console
$ genome tf-gene-list hg38 > tf_genes.txt
TF genes for hg38 / gencode_v50 (Homo sapiens)
  Lambert et al. 2018 v_1.01 (PMID 29425488) — https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
```

A third stderr line closes the account: how many genes and gene ids those came to, and how
many gene id stems this annotation carries no gene for. A census is keyed by *stems* —
gene ids with the version suffix dropped — and a registered annotation is not, so each
stem is resolved into the ids that annotation actually spells (one naming two genes prints
both, rather than whichever came first), and every one of Lambert's 1,639 assessed-positive
stems is either a gene on stdout or on that count, never quietly dropped.

Nothing here decides what a transcription factor is. The verdict is the census's, which is
why whose it is prints beside the answer; and the two censuses classify domain families
under their own publishers' vocabularies, which are not crosswalked — group by
`dbd_family` within a species, never across two.

`--annotation NAME` asks a registered annotation other than the assembly's default one.
`--json` carries the whole record: every gene with the census's own assessment and DBD
family, the provenance to cite, and the stems that resolved to nothing.

```console
$ genome tf-gene-list hg38 --json
{"assembly": "hg38", "annotation": "gencode_v50", "species": "Homo sapiens",
 "provenance": {"publisher": "Lambert et al. 2018", "version": "v_1.01", "pubmed_id": 29425488, …},
 "genes": [{"gene_id_stem": "ENSG00000137203", "gene_ids": ["ENSG00000137203.12"],
            "symbol": "TFAP2A", "is_tf": true, "dbd_family": "AP-2",
            "judgements": {"tf_assessment": "Known motif", …}}, …],
 "gene_ids": ["ENSG00000137203.12", …], "unresolved": […]}
```

The mouse census records nothing beyond the four columns every census shares, so a mouse
gene's `judgements` is empty rather than spelled differently — read it with `.get`.

Assessed-positive genes only, and there is no flag to widen it: a bare id list has nowhere
to say which of the two an id is, and a pipeline would read the rejected ones as
transcription factors. `Genome("hg38").tf_gene_list(include_rejected=True)` is where the
genes a census assessed and turned down are expressible, because there each id travels
with the verdict reached on it.

Exits `1` when the annotation is not registered here, when no census ships for the
assembly's species, and when nothing says what species the assembly is — three different
facts, each with its own message. None of them prints an empty list of genes.

## `genome tf-cofactor-list <assembly>`

Print the gene ids a publisher lists as transcription cofactors, one per line — the other
half of the machinery, shaped exactly like `genome tf-gene-list`. A cofactor is a chromatin
remodeller, a histone-modifying enzyme, a Mediator subunit: it recognises no sequence of its
own, so no scan will ever find it, but which genes are cofactors is published. The table is
chosen by the species the assembly's own metadata row names, never by anything you pass.

```console
$ genome tf-cofactor-list mm39 > cofactors.txt
TF cofactors for mm39 / gencode_vM39 (Mus musculus)
  AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Mus_musculus_Cof
```

A third stderr line closes the account as `genome tf-gene-list`'s does: how many cofactors
and gene ids those came to, and how many gene id stems this annotation carries no gene for.
Each of the mouse table's 970 stems is either a gene on stdout or on that count, never
quietly dropped.

`--annotation NAME` asks a registered annotation other than the assembly's default one.
`--json` carries the whole record: every gene with the publisher that listed it and that
publisher's own classification, one provenance entry per publisher to cite, and the stems
that resolved to nothing.

```console
$ genome tf-cofactor-list mm39 --json
{"assembly": "mm39", "annotation": "gencode_vM39", "species": "Mus musculus",
 "provenance": {"species": "Mus musculus", "ncbi_taxid": 10090, …,
                "sources": [{"publisher": "AnimalTFDB", "version": "4.0", "pubmed_id": 36268869, …}]},
 "cofactors": [{"gene_id_stem": "ENSMUSG00000000085", "gene_ids": ["ENSMUSG00000000085.16"],
                "symbol": "Scmh1", "is_cofactor": true, "source": "animaltfdb",
                "classifications": {"animaltfdb_family": "Others",
                                    "animaltfdb_category": "Other Cofactors"}}, …],
 "gene_ids": ["ENSMUSG00000000085.16", …], "unresolved": […]}
```

Nothing here decides what a cofactor is. Membership and classification travel with the
publisher who reached them, which is why whose list it is prints beside the answer, and why
each entry's `source` says who listed the gene — `both` there is agreement that the gene is
a cofactor and nothing about how either publisher classified it. The vocabularies are not
crosswalked: group by `animaltfdb_category` within one publisher, never across two.

`genome tf-cofactor-list ce11` answers while `genome tf-gene-list ce11` exits `1`. AnimalTFDB
assessed worm cofactors and no publisher has released a worm TF census — the publishers'
shape, not a defect here.

Exits `1` when the annotation is not registered here, when no cofactor table ships for the
assembly's species — the message names the species that have one — and when nothing says
what species the assembly is. Three different facts, each with its own message. None of them
prints an empty list of genes.

## `genome xref <species> <ids>...`

Convert identifiers to and from gene id stems against one published xref set — the way a
column of Entrez GeneIDs from a GEO series, UniProt accessions from a mass-spec run or HGNC
ids from a curated resource reaches this package's answers without writing Python, and
without the hand-built join everyone in the lab writes slightly differently. No assembly is
named and no genome is opened: an identifier is a name and not a place.

**The direction is named, never inferred.** `--to-stems NAMESPACE` reads the ids as that
namespace and answers in gene id stems; `--from-stems NAMESPACE` reads them as stems and
answers in that namespace. Exactly one is given — they carry the namespace, so neither can
be half-specified — and naming neither or both exits `2`. A string does not say which system
it belongs to, so `HGNC:11998` asked the wrong way answers *nothing found* rather than
quietly turning around.

```console
$ genome xref "Homo sapiens" --to-stems hgnc HGNC:11998 HGNC:13666 HGNC:10041 > stems.tsv
hgnc ids -> gene id stems for Homo sapiens (alliance 9.0.0)
  source  https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz
  2 resolved, 3 gene id stems, 1 this release names none for

$ cat stems.tsv
HGNC:11998	ENSG00000141510
HGNC:13666	ENSG00000094914
HGNC:13666	ENSG00000291836
HGNC:10041
```

The pairs go to stdout, tab-separated, so the output pipes — `cut -f2` is the answer and
`cut -f1` says what asked for it — and the heading, the publisher's URL and the counts go to
stderr. **Every id you passed gets at least one row.** One naming two genes prints both
rather than whichever came first (6.2% of human HGNC ids name two stems in this release),
and one that resolved to nothing gets a row with an empty second column: what your list
holds and this release does not is the one thing a hand-rolled join drops silently.

`--source NAME` picks the xref source; omitting it answers from the species' default one, so
everyone in the lab reaches for the same one without discussing it. It is a default and not
a recommendation — naming a source is how the scientific choice gets made deliberately, and
NCBI and Ensembl agree on only 57.6% of human gene-level (GeneID, ENSG) pairs. Every answer
names the source and release that produced it either way, so a result is reproducible a year
later. **The default is per species and per question** (ADR-0021): `--from-stems symbol` is
the one question here that is about symbols, and it is answered by the source that carries
them rather than by the identifier default, which for human carries none.

`--json` carries the same answer the API renders, keyed by what was asked about, with the
ids that named nothing under `unresolved`.

```console
$ genome xref "Homo sapiens" --from-stems uniprot ENSG00000141510.18 ENSG00000288541 --json
{"species": "Homo sapiens", "source": "alliance", "release": "9.0.0", "namespace": "uniprot",
 "resolved": {"ENSG00000141510.18": ["P04637"]}, "unresolved": ["ENSG00000288541"],
 "xref_ids": ["P04637"]}
```

There is no third direction: Entrez to HGNC is two calls and the join is yours, which keeps
the hop visible in your pipeline rather than invisible in ours. `genome.xref.XrefSet` in the
[API reference](reference.md) is the same two verbs from Python, on one code path with this.

**`--from-stems symbol` labels; the other way round is a command of its own.**
`--from-stems symbol` gives the authority's single current approved spelling for each stem,
which is what a figure axis wants, and reaches the source that carries symbols — `hgnc` for
human, `alliance_bgi` for mouse and worm — without your naming it. Toward stems is not its
mirror: a symbol also matches
spellings the authority has retired, answers with every gene any of them names, and carries
which kind matched — so it is [`genome match-symbols`](#genome-match-symbols-species-symbols)
that answers it, and `--to-stems symbol` exits `2` naming that command. The `--to-stems` help
lists the namespaces it does convert and no longer offers the one it refuses.

Naming a species prepares its set, which the first time is a download — so run it once on a
login node before submitting a job that needs it, as the lab's compute nodes have no
internet. Exits `1` when no set exists for the species (the message names the ones that do),
when the source is not one this package prepares, when the set is not here and cannot be
fetched (the message names the call to make on a login node), when the namespace is not one
the set carries (the message names the ones it does), and when a directory holds a set an
interrupted download left unfinished.

## `genome match-symbols <species> <symbols>...`

Print the genes each gene symbol names, and which kind of spelling matched — the way a gene
list copied out of a paper becomes usable without first finding its ids, and without the join
that silently drops every row spelling its gene the way the authority used to. No assembly is
named and no genome is opened: a symbol is a name and not a place.

```console
$ genome match-symbols "Homo sapiens" ARNTL ADCY3 Brca1 > genes.tsv
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

**A symbol is matched, never converted**, which is why this is a command rather than a third
direction of `genome xref`. Approved, previous *and* alias spellings are matched, every gene
id stem any of them names comes back, and each match says which kind of spelling it was.
`ARNTL` is a spelling HGNC retired and reaches `BMAL1` anyway: of EpiFactors v2.0's 801 human
rows, **31 spell their gene that way**, and an approved-only match drops exactly those.
`ADCY3` is HGNC's approved symbol for one gene and a symbol it retired from another, so it
answers with both — ambiguity is what you are handed rather than something resolved on your
behalf.

The matches go to stdout, tab-separated, so the output pipes — `cut -f3` is the answer,
`cut -f1` says what asked for it and `cut -f4` says which kind of spelling matched — and the
heading, the URL, the counts and what this source could not have matched go to stderr. The
last three columns are the keys `SymbolMatch.as_json()` writes, in its order, so the text and
`--json` cannot drift apart. **Every symbol you passed gets at least one row**, the ones that
matched nothing getting one with the other columns empty. Column 2 is the authority's own
spelling and column 1 is yours, which differ whenever the answer is worth reading twice.

**Matching is exact by default**, because the species is fixed by the set: `Brca1` asked of a
human set is a mouse spelling asked of the wrong authority, and matching nothing beats half
working. `--case-insensitive` folds both sides and still answers with **every** gene matched
rather than picking one.

```console
$ genome match-symbols "Homo sapiens" brca1 --case-insensitive --json
{"species": "Homo sapiens", "source": "hgnc", "release": "2026-07-07",
 "case_insensitive": true, "kinds": ["approved", "previous", "alias"], "limits": null,
 "resolved": {"brca1": [{"symbol": "BRCA1", "gene_id_stem": "ENSG00000012048",
                         "kind": "approved"}]},
 "unresolved": [], "gene_id_stems": ["ENSG00000012048"]}
```

**What this source could not have matched is printed too.** Only `hgnc` publishes previous
and alias spellings typed; `alliance_bgi` matches current approved symbols alone, because
each of its records files retired names and sequence names in one undifferentiated synonyms
list and typing them would put a kind on a claim no publisher made (ADR-0018). So a mouse
answer says which kinds it could match and, on a `limits` line, why the rest are missing —
without which *this gene is not in the release* and *this source does not publish the
spelling you used* would both be silence.

```console
$ genome match-symbols "Mus musculus" Arntl
gene symbols -> gene id stems for Mus musculus (alliance_bgi 9.0.0)
  source   https://download.alliancegenome.org/9.0.0/BGI/MGI/1.0.2.5_BGI_MGI_0.json.gz
  columns  asked, symbol, gene_id_stem, kind
  matching exact, on approved spellings
  0 resolved, 0 matches, 1 this release matched nothing for
  limits   this source publishes one current approved symbol per gene beside an undifferentiated synonyms list, …
Arntl
```

**Omitting `--source` answers from the species' default source for symbols**, which is not
the row its identifiers default to: human's ids come from `alliance`, whose cross-reference
file carries no human symbol at all, and its symbols from `hgnc`; mouse's and worm's from
`alliance_bgi`. A default is per species and per question for that reason (ADR-0021), and
every answer names the source and release that produced it either way. Naming a source is
still how the scientific choice gets made deliberately and a named one is never swapped, so
`--source alliance` exits `1` saying that set carries no symbol and naming the one that does,
rather than quietly answering out of somebody else's file. `XrefSet.for_symbols(species)` in
the [API reference](reference.md) is the same fill-in, on one code path with this, and
`genome.xref.XrefSet.match_symbols` the same one verb.

Naming a species prepares its set, which the first time is a download — so run it once on a
login node before submitting a job that needs it, as the lab's compute nodes have no
internet. Exits `1` when no set exists for the species (the message names the ones that do),
when the source is not one this package prepares, when a named source carries no symbols (the
message names the one that does), when the set is not here and cannot be fetched (the message
names the call to make on a login node), and when a directory holds a set an interrupted
download left unfinished.

## `genome homologs <species> <other-species> <stems>...`

Print the genes of another species a gene id stem's gene is homologous to, on Ensembl
Compara's own gene trees — the way a hit carries across species without leaving the package
and without the BioMart web API, whose intermittent failures make a pipeline built on it
fail irreproducibly. Everything here is a bulk file fetched once and read locally. No
assembly is named and no genome is opened: a homology set is anchored to a species pair and
a release, not to a build. Any pairing among human, mouse and worm answers, either way
round, off one prepared file per pair.

```console
$ genome homologs "Homo sapiens" "Caenorhabditis elegans" \
      ENSG00000152670 ENSG00000177479 ENSG00000000000 > homologs.tsv
Homo sapiens -> Caenorhabditis elegans orthologs (Ensembl Compara 116)
  source   Ensembl Compara release 116 (PMID 26896847) — https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz
  columns  gene_id_stem, homolog_gene_id_stem, homology_type, is_ortholog, is_high_confidence, goc_score, wga_coverage
  2 resolved, 4 links, 1 this release names no homolog for, 0 dropped partners
  quality  goc_score and wga_coverage null on every link of this set, so a filter on either empties rather than narrowing

$ cat homologs.tsv
ENSG00000152670	WBGene00001598	ortholog_one2many	True	False	NULL	NULL
ENSG00000152670	WBGene00001599	ortholog_one2many	True	False	NULL	NULL
ENSG00000152670	WBGene00001600	ortholog_one2many	True	False	NULL	NULL
ENSG00000177479	WBGene00020462	ortholog_one2one	True	True	NULL	NULL
ENSG00000000000
```

The links go to stdout, tab-separated, so the output pipes — `cut -f2` is the answer and
`cut -f1` says what asked for it — and the heading, the attribution, the counts and the two
qualifications go to stderr. **Every stem you passed gets at least one row.** One with three
homologs prints three rather than whichever came first, and one this release names no
homolog for gets a row with every other column empty: what your list holds and this release
does not is visible rather than dropped. An empty cell there is not `NULL`, which is
Compara's own word for a cell it recorded nothing in on a link that does exist.

**Every cell is the publisher's.** `homology_type` is Compara's tree-derived label printed
verbatim, and it is never recomputed from what came back — an answer can show one partner
and still read `ortholog_one2many`, which is why the label is carried instead of counted
(ADR-0020). The high-confidence flag and both quality scores come through the same way. This
package publishes no score, ranking or "best ortholog" of its own.

Two qualifications ride on every answer. **Dropped partners** — homologous genes a filter
removed — are counted *and* named, so a link that merely looks one-to-one in your view stays
distinguishable from one the publisher called one-to-one. And whichever quality columns the
set holds no value in *anywhere* are named up front: Compara records neither `goc_score` nor
`wga_coverage` on any link of **either** worm pairing, so `awk -F'\t' '$6 > 50'` would empty
itself rather than narrow, and you are told before you write it.

`--paralogs` returns every link the publisher wrote rather than only the ones its label calls
a speciation event. A paralogy link is **marked** by that label in the `homology_type` column
and never excluded, so *not an ortholog* stays distinguishable from *absent* — the stance
ADR-0013 already takes. Measured: release 116 publishes no cross-species paralogy for these
three species at all, so on it the flag changes nothing and the heading says which question
was asked either way.

`--release` names the Compara release, and `--json` carries the same answer the API renders,
keyed by the stems asked about, with the ones that named nothing under `unresolved`.

```console
$ genome homologs "Mus musculus" "Homo sapiens" ENSMUSG00000074698 --json
{"species": "Mus musculus", "other_species": "Homo sapiens", "release": "116",
 "resolved": {"ENSMUSG00000074698": [
   {"gene_id_stem": "ENSMUSG00000074698", "homolog_gene_id_stem": "ENSG00000101266",
    "homology_type": "ortholog_one2many", "is_ortholog": true, "is_high_confidence": true,
    "goc_score": 100, "wga_coverage": 100.0}, …]},
 "unresolved": [], "dropped_partners": [], "null_quality_scores": [],
 "homolog_gene_id_stems": ["ENSG00000101266", "ENSG00000254598"]}
```

Putting an answer into a registered annotation's own gene ids is
`genome.homology.resolve_homologs` from Python, and it is deliberately not a flag here: it
names an assembly, and this command names none. `genome.homology.HomologySet` in the
[API reference](reference.md) is the same one verb, on one code path with this.

Naming a pair prepares its set, which the first time is a download — so run it once on a
login node before submitting a job that needs it, as the lab's compute nodes have no
internet. Exits `1` when no set is pinned for a species (the message names the ones that
are), when the release is not pinned, when both species are the same one, when a stem carries
a version (the message names the stem to pass), when the set is not here and cannot be
fetched (the message names the call to make on a login node), when a directory holds a set an
interrupted download left unfinished, and when the file recorded as holding this pair holds
none of its rows — Compara's per-species dumps are a de-duplicated partition at the pair
level and the assignment moves between releases, so that is an error naming the other file
rather than an empty answer that would read as *these species share no homologs*.

## `genome verify <assembly>`

Re-read a FASTA and check its sha256 against the digest pinned for the assembly. This is
the deliberate full-file re-check; registering and reopening go by size alone. Exits `1`
on a mismatch.

```console
$ genome verify sacCer3
/data/liulab_data/genome/sacCer3/sacCer3.fa: sha256 6ff72f07… matches the digest pinned for it

$ genome verify sacCer3 --fasta /tmp/from-a-colleague.fa
error: sha256 mismatch for /tmp/from-a-colleague.fa: expected 6ff72f07…, got 9316629b….
```

`--fasta` checks any file against the official row — useful for a copy you were handed,
before you build anything on it.

A chimera is also checked against its components: the closing line reads `unchanged` or,
where a digest was missing on one side, `unknown`. A component that is no longer the one
the chimera was built from exits `1`.

## `genome revcomp <sequence>`

Reverse-complement a DNA sequence. Case is preserved; a character outside the `DNA`
alphabet exits `2` naming it. The alphabet is the type's own — `DNA.ALPHABET`, read at the
boundary rather than spelled again here — so the [`outside_alphabet`](sequences.md#construction)
you would call yourself is the check this command applies.

```console
$ genome revcomp ATCG
CGAT
$ genome revcomp aTcG --json
{"input": "aTcG", "reverse_complement": "CgAt"}
```

## `genome table-row <assembly>`

Download an assembly, unpack it, hash the unpacked FASTA, and print the line to paste
into the shipped metadata table — how the checksum column gets filled in. A checksum the
table already pins is reported, never enforced.

```bash
$ genome table-row sacCer3
sacCer3	Saccharomyces cerevisiae	sacCer3	R64-1-1	GCF_000146045.2	559292	https://…/sacCer3.fa.gz	6ff72f07…
```

A chimera exits `1` before anything is downloaded: its row carries the name and nothing
else, so there is no row here to compute. Check one with `genome verify` instead.

## `genome motif-scan <fasta> <output>`

Scan a FASTA with a JASPAR release, write the hits to Parquet, and print a summary of the
run. The batch case — a shell script, a scheduler job — and the only motif command there
is: listing, plotting and comparing motifs are notebook work.

```console
$ genome motif-scan peaks.fa hits.parquet --release 2024 --tax-group all
scanned 41255 sequences with 2338 motifs from JASPAR 2024 all
  background  0.293, 0.207, 0.207, 0.293
  threshold   0.0001
  skipped     8 under 7 positions, so not scanned: MA0261.1, MA2355.1, …
  workers     16
  hits        2841193 -> hits.parquet
```

**The hits go to the named file and the summary to standard output**, so `--json` is never
corrupted by table data and a pipeline can consume the one while a downstream step reads
the other:

```console
$ genome motif-scan peaks.fa hits.parquet --release 2024 --tax-group all --json
{"release": "2024", "tax_group": "all", "motifs_scanned": 2338,
 "motifs_skipped": ["MA0261.1", "MA2355.1", …], "background": [0.293, 0.207, 0.207, 0.293],
 "threshold": 0.0001, "sequences_scanned": 41255, "hits_written": 2841193,
 "workers": 16, "output": "hits.parquet"}
```

A motif under seven positions is named there rather than scanned: a 6-mer's best possible
word has p = 2.44e-4, so it cannot reach the default threshold at all, and an engine asked
for it anyway would fall back to that matrix's best attainable cutoff and over-call in
silence.

Read the hits back with `genome.tf.motif.read_hits(path)`, which restores the compact
dtypes and the provenance — `pandas.read_parquet` gives the rows and drops what the scan
was.

| Option | Effect |
|---|---|
| `--release` | Which JASPAR release: `2024` or `2026` (default) |
| `--tax-group` | Which taxonomic group's file: `vertebrates` (default) through `all` |
| `--threshold` | The per-position p-value each motif's cutoff is converted from; `1e-4` by default |
| `--background` | `auto` (default), `uniform`, or `derive` — see below |
| `--workers` | How many processes; every core the allocation granted by default |

`--background` decides the answer more than any other option: `auto` derives the base
composition from the input when it holds at least 10 000 unambiguous bases and stays
uniform below that, `uniform` pins it, and `derive` derives whatever the input holds.
Whichever it was is in the summary and on the hits. Four frequencies of your own are a
Python call (`scan_fasta(background=[…])`) rather than a flag — a mistyped one on a
command line would change every cutoff and look like a scan that simply found fewer hits.

**It defaults to every core the allocation granted** — `SLURM_CPUS_PER_TASK`, then
`SLURM_CPUS_ON_NODE`, then this process's CPU affinity, then the machine — where the
library defaults to one worker. A console script is a proper entry point, so the
process-pool hazard behind that library default does not apply here. `--workers 1` scans
serially and produces the identical table.

**Prepare the release from a login node.** Naming a release prepares it, and the first
time that is a download. The lab's CPU cluster compute nodes have no internet, so run this
once on a login node before submitting a job that needs it; every run afterwards reads the
cached file out of `<LIULAB_DATA>/motif/jaspar/`, shared by every project on the machine.

Exits `1` when the FASTA is not there or is not FASTA, when the release or tax group is
not one this package prepares, when the threshold is not a p-value in `(0, 1)`, when the
worker count is below 1, and when the release is not cached here and cannot be fetched —
which is what a compute node with no internet looks like. A `--background` that is not one
of the three modes is refused by the parser with `2`, before anything is fetched or read.
