# Test fixtures

Small, subsampled **real** files — never a large genomic file. Every sequence byte here came from
UCSC's `sacCer3` golden path, every motif byte from a published JASPAR release and every
cross-reference row from a published Alliance of Genome Resources release; all three were cut
down on the lab cluster. Three fixtures depart from the source bytes and every departure is named
below: names are this repo's where a fixture needs a spelling `sacCer3` does not have; two stretches
are lower-cased; and `planted_motifs.fa` has **39 bases overwritten** with three JASPAR consensus
words, which is the one place a base here was chosen rather than found. Nothing else is synthesised.

| File | What it is |
| --- | --- |
| `tiny.fa` | The first 10 000 bases of `chrI`, `chrII` and `chrIII` of `sacCer3`, cut with `samtools faidx` and re-headed to the bare chromosome names |
| `tiny.fa.gz` | `tiny.fa`, gzipped — the compressed-source path |
| `tiny.gtf` | The 85 `sacCer3.ensGene` features that fall wholly inside those three windows |
| `tiny.gtf.gz` | `tiny.gtf`, gzipped |
| `ensembl_style.gtf` | `tiny.gtf` with the `chr` prefix stripped (`I`, `II`, `III`) — an Ensembl-spelled annotation for the chromosome-name mismatch case |
| `tiny_jaspar_transfac.txt` | Ten whole records of the JASPAR 2024 `all` union file, in its own order — see below |
| `planted_motifs.fa` | Two 600-base windows of `tiny.fa` with three motif consensus words written into them at known positions and strands — see below |
| `xref/alliance_genecrossreference_tiny.tsv.gz` | Fourteen whole genes of the Alliance of Genome Resources release 9.0.0 cross-reference file, copied verbatim — see below |
| `xref/hgnc_complete_set_tiny.txt` | Ten whole rows of HGNC's quarterly archive file of 2026-07-07, under its own 54-column header, copied verbatim — see below |
| `xref/bgi_mgi_tiny.json.gz` | Seven whole gene records of MGI's Alliance gene submission, each object's own bytes — see below |
| `xref/bgi_wb_tiny.json.gz` | Five whole gene records of WormBase's WS298 Alliance gene submission, likewise — see below |
| `homology/compara116.*.homologies.tsv.gz` | Whole rows of three Ensembl Compara 116 per-species homology dumps, in each file's own order — see below |

Sources:

- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz`
- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz`
- `https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_transfac.txt`
- `https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz`
- `https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/<species>/Compara.116.protein_default.homologies.tsv.gz`

`tiny.gtf` coordinates are 1-based inclusive, as every GTF is; they convert at the I/O boundary and
are never seen in that form inside the package.

## `tiny_jaspar_transfac.txt` — the motif records

Ten records of the JASPAR 2024 `all` union file, copied whole and kept in that file's own order.
Nothing is edited: every count, every annotation value and every blank annotation is what JASPAR
published. Each one is here for a rule it breaks, and the union file is the source because only it
holds every **Tax group** — five are represented, `diatoms` included, which is a group of exactly
one motif in both releases.

| Record | Name | Positions | Tax group | The rule it exists to break |
| --- | --- | --- | --- | --- |
| `MA0119.1` | `NFIC::TLX1` | 14 | vertebrates | A **dimeric name**, and the record behind it carries **two** classes, **two** families and **two** UniProt accessions — semicolon separated, which is why those four annotations are tuples and not strings |
| `MA0789.1` | `POU3F4` | 9 | vertebrates | **Two PubMed ids** on one matrix, so the plural id lists are exercised without a dimer |
| `MA0079.5` | `SP1` | 9 | vertebrates | **Fractional counts** — `1.05485`, which the `.jaspar` serialization rounds away and this one keeps. The only record here whose counts are not whole |
| `MA0139.2` | `CTCF` | 15 | vertebrates | The everyday case, and the **first of three motifs sharing one name** |
| `MA1929.2` | `CTCF` | 31 | vertebrates | The second. Two would make a name ambiguous; three make the error list more than a pair |
| `MA1930.2` | `CTCF` | 33 | vertebrates | The third, the **longest matrix in the release**, and the one with the **least informative flanks** — 0.36 and 0.31 bits, which trim at 0.4. It also carries twelve interior positions under 0.25 bits, so it is real data for the rule that trimming acts **only on the ends** |
| `MA2355.1` | `PK06791.1` | 6 | plants | **Below the minimum scannable length**; its class is `C3H(C),C2HC zinc-fingers like factors`, one value **containing a comma**; and its UniProt list is **empty** |
| `MA0261.1` | `lin-14` | 6 | nematodes | Below the minimum length again, and **both** its class and its family are blank — the source stated nothing, which is common and is not an error |
| `MA0283.1` | `CHA4` | 8 | fungi | Its data type is `PBM, CSA and/or DIP-chip`: **one value containing a comma**, in the other annotation where that happens |
| `MA1407.2` | `bZIP14` | 8 | diatoms | The release's **only diatom motif**, so the degenerate tax group is a case the fixture actually holds |

**The separator is a semicolon and never a comma.** Four records above turn on that: two carry
several values in one annotation, separated by `;`, and two carry a comma *inside* a single value.
Splitting on the comma corrupts about fifty records per release and fails nothing while doing it,
which is what these four are here to catch.

None of this is prose to be trusted either: `tests/tf/test_jaspar.py` asserts every row of the table —
each id, name, length and tax group, each multi-valued and each blank annotation, the fractional
counts, the two below-minimum lengths and the long record's flank bits — against the committed
bytes.

## `planted_motifs.fa` — the sequences a scan is checked against

**This is the one fixture whose bases are not all found bases.** Everything else here is `sacCer3`
as UCSC published it; this file is two 600-base windows of `tiny.fa` — `chrI:1-600` as `plantedI`
and `chrII:1-600` as `plantedII`, 1-based inclusive — with **three consensus words written over 39
of those 1 200 bases**. A scan needs a site it knows the answer for, and yeast does not oblige on
demand. Every base outside the three intervals below is `sacCer3`'s, and `tests/tf/test_scan.py`
asserts exactly that: it puts the source bases back over the three intervals and gets `tiny.fa`'s
own windows out again.

What was planted, where, and on which strand — coordinates **0-based half-open**, as everything
inside this package is:

| Record | Interval | Strand | Motif | Bases written | The rule it exists to break |
| --- | --- | --- | --- | --- | --- |
| `plantedI` | `[100, 115)` | `+` | `MA0139.2` CTCF | `GCCACCAGGGGGCGC` | The everyday case: a site on the forward strand, found at the position it was planted at |
| `plantedI` | `[300, 315)` | `-` | `MA0139.2` CTCF | `GCGCCCCCTGGTGGC` | The **same word reverse-complemented**. It must be reported on `-` over *these* bases — a forward-frame `[300, 315)` and not a position counted from the other end, which is the off-by-one the whole adapter exists to centralise |
| `plantedII` | `[200, 209)` | `+` | `MA0789.1` POU3F4 | `tatgcaaat` | **Lower case**: a soft-masked site, which must be found exactly as its upper-case equivalent is, because a scan discards **Soft-masking** without being asked (ADR-0012) |

Two more things the bytes carry, and neither is decoration:

- **`plantedII`'s header is `>plantedII  sacCer3 chrII:1-600, bases 180-240 soft-masked`.** The
  record name is `plantedII` and everything after the first whitespace is description — which is
  what STAR and chromap write into an alignment made from the same file, so a **Hit table** joins
  against that alignment with nobody renaming anything. `plantedI`'s header carries no description,
  so both shapes are present.
- **`plantedII[180:240)` is lower-cased** — exactly one wrapped line, so it is visible in the file,
  and it holds the planted POU3F4 site with flanks either side of it. The bases are `sacCer3`'s;
  only their case is this repo's.

The motifs are the ones `tiny_jaspar_transfac.txt` already holds, so this file adds no matrices.
That fixture's two below-minimum records, `MA2355.1` and `MA0261.1`, are why nothing needed to be
planted to test the length floor: they can never be scanned, and every scan names them.

None of this is prose to be trusted: `tests/tf/test_scan.py` asserts every row of the table above
against the committed bytes — the two record names and their headers, the 600-base lengths, each
planted word at its offset, that the reverse word really is the forward one flipped, the bounds of
the masked window, the 60-base wrap, and the `sacCer3` backbone underneath all of it.

## `chimera/` — the tiny component assemblies

Four component assemblies, cut from the files above rather than downloaded again, so every base and
every annotation line is still real `sacCer3`. Each carries something the naming contract has an
opinion about and no shipped assembly can demonstrate:

| Component | Chromosomes | Wrap | What it is here for |
| --- | --- | --- | --- |
| `tinyCe` | `I`, `II`, `MtDNA` | 60 | `ce11`'s shape, and the only **soft-masked** component — so a chimera of these is heterogeneously masked. **Collides** with `tinySc` on `I` and `II`, which the shipped pair never does |
| `tinyEc` | `NZ_TINY01000001.1`, `NZ_TINY01000002.1`, `chr1_KI270706v1_random` | 80 | Names that already **hold an underscore**, in both real shapes: `ecHT115`'s accession, an underscore and a dot in one name, so no split may be a first-occurrence one; and the `hg38` name that a *single*-underscore separator could not tell from a suffixed one. Wrapped at a second width, so no build may rewrap |
| `tinySc` | `I`, `II`, `III` | 60 | `sacCer3` spelled the Ensembl way, where each name is a **strict prefix** of the next. Being the third makes N > 2 an ordinary case |
| `tinyEcDub` | `NZ_TINY02000001.1`, `NZ_TINY02__000002.1` | 60 | A name already carrying a **doubled underscore** — the only thing that pushes the separator past `__`, and something no real assembly does. Its own name is also a strict prefix trap, `tinyEc` inside `tinyEcDub`. Ships no GTF: it exercises a naming rule |

`tinyCe`, `tinyEc` and `tinySc` are the everyday set and ask for the `__` separator, as all seven
shipped assemblies do; adding `tinyEcDub` is what makes a chimera derive a longer one.

Each chromosome is one **disjoint** slice of `tiny.fa`, renamed — disjoint so a test can always tell
which component a chromosome came from, and so no gene id is carried by two components, which would
hand the annotation merge an id collision the shipped pair does not have. The GTF beside each FASTA
is every `tiny.gtf` transcript falling wholly inside that slice, with the chromosome name rewritten
and the coordinates shifted to the slice — so each GTF's chromosome names are **set-equal** to its
component's, as both shipped GTFs are against their assemblies:

| Component | Chromosome | Slice of `tiny.fa` (1-based inclusive) | Length |
| --- | --- | --- | --- |
| `tinyCe` | `I` | `chrII:2701-5200` | 2500 |
| `tinyCe` | `II` | `chrIII:6301-8650` | 2350 |
| `tinyCe` | `MtDNA` | `chrI:6901-9200` | 2300 |
| `tinyEc` | `NZ_TINY01000001.1` | `chrIII:2601-6300` | 3700 |
| `tinyEc` | `NZ_TINY01000002.1` | `chrI:2401-4600` | 2200 |
| `tinyEc` | `chr1_KI270706v1_random` | `chrII:5701-7700` | 2000 |
| `tinyEcDub` | `NZ_TINY02000001.1` | `chrII:7701-8300` | 600 |
| `tinyEcDub` | `NZ_TINY02__000002.1` | `chrII:8301-8800` | 500 |
| `tinySc` | `I` | `chrI:1-2400` | 2400 |
| `tinySc` | `II` | `chrII:1-2700` | 2700 |
| `tinySc` | `III` | `chrIII:1-2600` | 2600 |

Every length is distinct, so a chimera's chromosome can be named by its length alone. The one
departure from the source bytes is the first **200 bases of `tinyCe`'s `I`**, lower-cased to stand
in for the repeat masking a real assembly carries — the bases are `sacCer3`'s, only their case is
this repo's.

None of this is prose to be trusted: `tests/assembly/test_chimera_fixtures.py` asserts every slice, length,
wrap width and masked stretch against the committed bytes, and `CHIMERA_COMPONENTS` in
`tests/conftest.py` is the same table as data, with a `chimera_component` fixture that registers one
as an assembly.

## `xref/alliance_genecrossreference_tiny.tsv.gz` — the cross-reference rows

Fourteen whole genes of the Alliance of Genome Resources `GENECROSSREFERENCE_COMBINED` file of
release 9.0.0 — every row those genes have, copied verbatim, under the file's own fourteen comment
lines and its own header. Nothing is edited and nothing is synthesised: every duplicate row, every
cross-reference to a database this package does not carry, and the file's own column order are what
the Alliance published. 2.7 KB gzipped, against the publisher's 25 MB.

Whole genes rather than sampled rows, because a gene's rows are what a slice is built from: the
**Gene id stem** is the `ENSEMBL:` cross-reference on the gene, so half a gene's rows would be a
gene this package could not resolve, which is a case the file genuinely has and one the fixture
should not manufacture.

| Gene | Species | The case it is here for |
| --- | --- | --- |
| `HGNC:11998` (`TP53`) | human | The plain case: one Ensembl id, one Entrez id, one UniProt accession |
| `HGNC:1100` (`BRCA1`) | human | Its `NCBI_Gene:672` pair is listed **twice**, under two different pages — the duplication, on real bytes |
| `HGNC:13666` | human | **Two `ENSEMBL:` cross-references**, so one HGNC id and one Entrez id each name two stems and nothing may pick one |
| `HGNC:7622` | human | The same again, reached from the other side: `NCBI_Gene:4661` names two stems |
| `HGNC:10041` | human | A real human gene the Alliance lists with **no `ENSEMBL:` cross-reference at all**, so it has no hub and appears in no slice |
| `MGI:98834` (`Trp53`) | mouse | Fourteen UniProt accessions on one stem — the reverse verb's many-valued answer |
| `MGI:105105`, `MGI:1921534` | mouse | Two `ENSEMBL:` cross-references each, so the collision is not a human-only artefact |
| `WB:WBGene00000001` (`aap-1`) | worm | The identity hop — its `ENSEMBL:` id **is** its WBGene id — and a `WB:aap-1` **symbol** row under the authority's own prefix, which is the trap the reader must not read as a gene id |
| `WB:WBGene00000912` (`daf-16`) | worm | A second worm gene, so worm is not one row |
| `WB:WBGene00003425/32/49/63` | worm | Four genes sharing `UniProtKB:P05634`, so one foreign id names **four** stems |

Two things the fixture deliberately does **not** show, because the publisher's file does not:

- **No versioned identifier.** Not one of the 43,867 human ids in release 9.0.0 carries a version
  suffix. The ingest-side version stripping is covered instead by rows the test writes itself, in
  `TestVersionsOnTheSourceSide`, and it is kept because other publishers do spell versions — NCBI's
  `gene2ensembl` writes the gene id bare and the transcript id versioned in the same row.
- **No duplicate whole rows.** A whole-row `uniq` removes nothing here, exactly as it removes
  nothing from the full file: the duplication is on `(GeneID, GlobalCrossReferenceID, TaxonID)`,
  where 2,659,704 rows reduce to 1,811,267 distinct.

None of this is prose to be trusted: `TestFixtureBytes` in `tests/xref/test_xref.py` asserts the header,
the three taxa, the duplication's shape, the hub-less gene and the worm symbol row against the
committed bytes, and `FIXTURE_SLICES` pins what each species' slice comes to.

## `homology/` — the Compara homology dumps

Three subsamples of Ensembl Compara release 116's per-species `protein_default` homology dumps, one
per species the package prepares a **Homology set** for. Every line is a whole published row copied
verbatim, in the source file's own order, under the publisher's own header — nothing is edited, and
no row is synthesised. Each source file is 85–110 MB gzipped and holds millions of rows; these hold
tens.

| File | Rows | What it is here for |
| --- | --- | --- |
| `compara116.homo_sapiens.homologies.tsv.gz` | 17 | The human dump: 11 human↔worm rows, 4 human↔human paralogy rows and 2 human↔zebrafish rows — and, as the real file has, **zero human↔mouse rows** |
| `compara116.mus_musculus.homologies.tsv.gz` | 27 | The mouse dump: 10 mouse↔human rows, 11 mouse↔worm rows, 4 mouse↔mouse paralogy rows and 2 mouse↔zebrafish rows — one file holding **two** of the three pairs |
| `compara116.caenorhabditis_elegans.homologies.tsv.gz` | 4 | The worm dump: four worm↔worm paralogy rows and **no row of either pair**, which is the real file's shape too |

What each is in the set *for*:

- **The partition, for real.** The human file holds no human↔mouse row, so preparing that pair from
  it is the published partition trap rather than a staged one, and the worm file holds no pair at
  all. Which file holds which pair is what `homology_metadata.tsv` records, measured.
- **All three speciation Homology types.** Both cross-species pairs carry `ortholog_one2one`,
  `ortholog_one2many` and `ortholog_many2many` rows, and every gene taken was taken with **all** of
  its rows against that species — so a **Homology type** here is the publisher's own and not an
  artefact of the sampling. Two human genes, `ENSG00000163331` and `ENSG00000112977`, are
  `many2many` to the *same* two worm genes, which is what a partner named by two asked stems
  looks like.
- **Both null quality scores.** Every human↔worm and mouse↔worm row carries `NULL` in `goc_score`
  and `wga_coverage`, as 100% of the real rows of both pairings do; the mouse↔human rows carry real
  values in both, and the zebrafish rows carry them in the human file too — so a set that holds
  nothing in a column is a fact about the pair and not about the file it was cut from.
- **Paralogy rows that are not a pair's.** The paralogy rows are `other_paralog` (human, mouse) and
  `within_species_paralog` (worm), and every one of them relates two genes of *one* species, as
  every duplication row in release 116 does — there is no `between_species_paralog` row anywhere in
  the published human dump. They are here to hold the boundary: a **Homology link** relates two
  species, so these are in the file and never in a pair's set.
- **Compara's own null spelling**, `NULL`, and its flag spellings `1`, `0` and `NULL` — the
  high-confidence flag is set on the ortholog rows and null on every paralogy row.

None of this is prose to be trusted: `tests/homology/test_homology.py` asserts every row of the table above
against the committed bytes.

## `xref/ensembl_entrez_*_tiny.tsv.gz` — the second source's rows

Cut from Ensembl release 116's per-species dumps,
`Homo_sapiens.GRCh38.116.entrez.tsv.gz` and `Mus_musculus.GRCm39.116.entrez.tsv.gz`, under each
file's own header. Verbatim: nothing edited, nothing synthesised, every column kept.

**Whole closures rather than sampled rows.** Each fixture is every row the release has that touches
a chosen set of GeneIDs, closed under the (GeneID ↔ stem) relation — so if a GeneID is here, every
stem it names is here, and if a stem is here, every GeneID naming it is here. A fan-out counted on
the fixture is therefore the fan-out the release asserts, not an artefact of where the sample
stopped. 151 human rows against the publisher's 1.2 million; 75 mouse rows.

| Fixture | Rows | What it is here for |
| --- | --- | --- |
| human | 151 | 149 `EntrezGene` rows, **every one `DEPENDENT` and not one `DIRECT`** — the empty-filter trap on real bytes, which is release 116's whole shape in miniature (552,633 rows, zero direct) |
| human | | GeneID **`79166` naming 72 stems** — the fan-out, and the id the two sources disagree about |
| human | | Stem **`ENSG00000173213` named by four GeneIDs** — the fan-out running the other way, which is what makes `from_stems` many-valued |
| human | | GeneID `4661` naming two stems, and `7157` (`TP53`) naming one — the ordinary cases beside the extreme ones |
| human | | Two `EntrezGene_trans_name` rows carrying `KU-MEL-3-201` and `BMS1P4-202` in the `xref` column — **transcript names, not GeneIDs**, and `MISC` rather than `DEPENDENT`, so a reader that did not split on `db_name` would key the Entrez namespace by transcript labels |
| mouse | 75 | The same shape at a second species: 71 rows all `DEPENDENT`, GeneID `100040298` naming three stems, stem `ENSMUSG00000047675` named by two, four `MISC` transcript-name rows, and `22059` (`Trp53`) naming `ENSMUSG00000059552` — which the Alliance fixture also carries |

## `xref/alliance_genecrossreference_disagreeing.tsv.gz` — one gene, two publishers

All eight rows Alliance release 9.0.0 has for `HGNC:15497`, under the file's own comment block and
header, copied verbatim. It exists so that one test can construct **two `XrefSet`s over one species
and ask both about one id**: Entrez GeneID `79166` names **two** stems here — `ENSG00000170858` and
`ENSG00000293273` — and **seventy-two** in Ensembl release 116. Neither answer is merged into the
other and each names its own source (ADR-0017).

The existing 14-gene Alliance fixture cannot show this. Every human id it carries — `672`, `4661`,
`7157`, `8086` — Ensembl 116 agrees with exactly, so a disagreement had to be cut from a gene it
does not hold rather than manufactured in one it does.

The same eight rows also carry, incidentally and for real, two traps the reader already handles: the
`NCBI_Gene:79166` pair listed **twice** under two different pages, and an `RGD:` cross-reference on a
*human* row.

None of this is prose to be trusted: `TestFixtureBytes` in `tests/xref/test_xref_ensembl.py` asserts the
header, the evidence types, the transcript-name rows and both fan-outs against the committed bytes,
and the module's constants pin what each slice comes to.

## `xref/hgnc_complete_set_tiny.txt` — the typed spellings

Ten whole rows of HGNC's quarterly archive file dated **2026-07-07**, under the file's own header,
copied verbatim and sorted by `hgnc_id`. Plain text and not gzipped, because the published file is:
the URL a curated row pins ends `.txt`, and the fake fetch serves the fixture under that name.

The header is **54 columns**, which is the current schema; the same file was **52** in 2020, which
is why the reader finds its seven columns by name. Every row is real, and each is there for a
reason:

| Row | What it is in the fixture for |
| --- | --- |
| `BMAL1` | `prev_symbol` **`ARNTL`** — one of the 31 EpiFactors rows the cofactor work measured as spelling its gene the way HGNC used to. Also five aliases in one quoted, pipe-separated cell |
| `EMSY` | `prev_symbol` **`C11orf30`**, the second named one of those 31, and an **empty** alias cell beside it |
| `ACIN1` | `prev_symbol` `ACINUS`, the third — and a two-value alias cell |
| `ADCY3` and `ADCY8` | **One spelling, two genes, two kinds**: `ADCY3` is HGNC's approved symbol for the first and a symbol it retired from the second, so a match on it answers with both and says which was which |
| `TP53` | The alias **`p53`**, written in lower case by HGNC itself, so that exact and case-insensitive matching differ on a real spelling rather than an invented one |
| `BRCA1` | Whose mouse spelling is `Brca1` — the exact-by-default case, asked of a human set |
| `KDM1A` | Two previous and three alias spellings, the widest row here |
| `MECP2` | Three previous spellings and no alias, the mirror of `EMSY` |
| `AAVS1` | **An empty `ensembl_gene_id`** — no hub, so the row contributes nothing at all, as 2,682 of the archive's 45,019 rows do |

## `xref/bgi_mgi_tiny.json.gz` and `xref/bgi_wb_tiny.json.gz` — current symbols

Seven mouse and five worm gene records from the Alliance's per-species gene submissions —
`9.0.0/BGI/MGI/1.0.2.5_BGI_MGI_0.json.gz` and `8.3.0/BGI/WB/1.0.2.5_BGI_WB_4.json.gz`, the latter
re-served under release 9.0.0. Each record is the publisher's **own bytes**, indentation included,
lifted whole out of the `data` array and re-wrapped with that file's own `metaData` header — so the
two publishers' two different pretty-printers, two-space for MGI and three-space for WormBase, are
both in the suite.

| File | Genes | What it is in the fixture for |
| --- | --- | --- |
| MGI | `Trp53`, `Brca1`, `Bmal1`, `Kdm1a`, `Mecp2`, `Scgb1a1`, `Adcy3` | Mouse spellings of genes the human fixture carries, so that a mouse-cased symbol asked of a human set is a real mouse symbol rather than a typo |
| WB | `aap-1`, `daf-16`, `msp-10`, `lin-12`, `unc-22` | The worm hop, where a WormBase gene id **is** the **Gene id stem**. `daf-16`'s record carries `daf-17`, `R13H8.1` and `CELE_R13H8.1` in **one undifferentiated `synonyms` list** — a genuine former name beside two sequence names, with nothing saying which is which, which is why this source publishes approved spellings only |

The worm file's own header states `"release" : "WS298"` — WormBase's final release, and the one the
lab has registered. `downloads.wormbase.org` answers **403** to an automated client, so this copy,
served by the Alliance at a pinned path with a published md5, is how those bytes are reachable.

None of this is prose to be trusted: `TestFixtureBytes` in `tests/xref/test_xref_symbols.py` asserts the
column count, the quoted multi-valued cell, the hub-less row and the WS298 header against the
committed bytes, and the module's constants pin what each slice comes to.
