# What do `ecHT115` and `ce11` actually contain?

Measured 2026-08-12 for [#39](https://github.com/liuhlab/liulab-genome/issues/39), on **GPU71FM**,
from the two assemblies this package had already prepared there. **`ecHT115` is a contig-level
SPAdes draft — 87 unplaced scaffolds, not a chromosome plus plasmids — so no design may assume a
chimera has few large sequences; `ce11` is 7. Both GTFs' seqnames are set-equal to their assembly's,
so a merged annotation must preserve every spelling or `ChromosomeMismatchError` fires on all of
them. A last-occurrence split on the separator recovers the component for every candidate character;
a first-occurrence split does not. The two differ in FASTA shape — `ecHT115` 80-column, upper-case,
unmasked; `ce11` 60-column and 21.95% soft-masked — so a chimera is heterogeneously masked.**

This is a measurement, not a decision. It corrected a fact the chimera map
([#38](https://github.com/liuhlab/liulab-genome/issues/38)) had been charted on, and fed the
naming, build and merge decisions from
[#41](https://github.com/liuhlab/liulab-genome/issues/41) onward.

## Findings

Read off **GPU71FM**, from the two assemblies this package already prepared there, on
2026-08-12. Data root is `/home/hanliu/liulab_data` (`$LIULAB_DATA` unset in a
non-interactive shell; the completion records show the same tree via its
`/large_storage/zhoulab/hanliu/...` real path). Both were prepared by package version
`0.0.1.dev51+dev.g0f644e123` with `samtools 1.23.1`, and both carry finished
`.completion.json` records whose `sha256` matches the pinned digest in
`src/genome/data/assembly_metadata.tsv`. Cross-checked against the NCBI assembly report
for `GCF_004354945.1`.

**The one thing that is not what the map assumed:** `ecHT115` is **not** a chromosome plus
plasmids and carries **no** `NZ_CP…` accession. `GCF_004354945.1` / `ASM435494v1` is a
**contig-level SPAdes draft**: 87 unplaced WGS contigs, all named `NZ_SMTD010000NN.1`.

---

## 1. `ecHT115` — GCF_004354945.1 / ASM435494v1

**87 sequences, 4,602,240 bp total.**

### What the assembly is

From `GCF_004354945.1_ASM435494v1_assembly_report.txt`:

```text
# Assembly name:  ASM435494v1
# Organism name:  Escherichia coli HT115 (E. coli)
# Infraspecific name:  strain=HT115
# Taxid:          634469
# BioSample:      SAMN11094052
# BioProject:     PRJNA526261
# Submitter:      Universidad Mayor
# Date:           2019-03-18
# Assembly level: Contig
# Genome representation: full
# WGS project:    SMTD01
# Assembly method: SPAdes v. 3.1.0
# Genome coverage: 200.0x
# Sequencing technology: Illumina MiSeq
# RefSeq assembly and GenBank assemblies identical: yes
```

Across all 87 rows of that report:

| column | value, all 87 rows |
| --- | --- |
| `Sequence-Role` | `unplaced-scaffold` |
| `Assigned-Molecule` | `na` |
| `Assigned-Molecule-Location/Type` | `na` |
| `Assembly-Unit` | `Primary Assembly` |
| `UCSC-style-name` | `na` |

**No plasmids are identified, and neither is a chromosome.** Nothing in the assembly is
molecule-assigned, so there is no plasmid naming question to answer: a plasmid, if one is
in there, is an unlabelled contig indistinguishable from the rest. Sequence lengths in the
report sum to 4,602,240 — the same total as the prepared `chrom.sizes`.

The FASTA deflines carry the original SPAdes node names, e.g.

```text
>NZ_SMTD01000010.1 Escherichia coli HT115 NODE_10_length_132761_cov_25.566785, whole genome shotgun sequence
```

Accessions are assigned in size-descending order of the node. Note the accessions are
**not** in one-to-one numeric correspondence with the node numbers: `NODE_61` was dropped,
so `NZ_SMTD01000061.1` is `NODE_62` and the offset persists to the end
(`NZ_SMTD01000087.1` = `NODE_88`). Only the accession reaches `chrom.sizes`; the node name
lives only in the defline comment.

### `ecHT115.chrom.sizes` — verbatim, in file order

This is byte-identical to columns 1–2 of `ecHT115.fa.fai` and to STAR's
`chrNameLength.txt`, i.e. **`chrom.sizes` order is FASTA order** (verified by `diff`).

```text
NZ_SMTD01000010.1 132761
NZ_SMTD01000011.1 121051
NZ_SMTD01000012.1 117739
NZ_SMTD01000013.1 112487
NZ_SMTD01000014.1 107747
NZ_SMTD01000015.1 103890
NZ_SMTD01000016.1 92981
NZ_SMTD01000017.1 88592
NZ_SMTD01000018.1 87159
NZ_SMTD01000019.1 83572
NZ_SMTD01000001.1 327216
NZ_SMTD01000020.1 78749
NZ_SMTD01000021.1 72745
NZ_SMTD01000022.1 70023
NZ_SMTD01000023.1 67479
NZ_SMTD01000024.1 61547
NZ_SMTD01000025.1 58548
NZ_SMTD01000026.1 58531
NZ_SMTD01000027.1 57161
NZ_SMTD01000028.1 55003
NZ_SMTD01000029.1 49636
NZ_SMTD01000002.1 286405
NZ_SMTD01000030.1 43952
NZ_SMTD01000031.1 42626
NZ_SMTD01000032.1 41423
NZ_SMTD01000033.1 41013
NZ_SMTD01000034.1 41006
NZ_SMTD01000035.1 40308
NZ_SMTD01000036.1 35465
NZ_SMTD01000037.1 35284
NZ_SMTD01000038.1 33817
NZ_SMTD01000039.1 31727
NZ_SMTD01000003.1 265311
NZ_SMTD01000040.1 31564
NZ_SMTD01000041.1 29560
NZ_SMTD01000042.1 29519
NZ_SMTD01000043.1 27807
NZ_SMTD01000044.1 26598
NZ_SMTD01000045.1 25380
NZ_SMTD01000046.1 24131
NZ_SMTD01000047.1 19551
NZ_SMTD01000048.1 19279
NZ_SMTD01000049.1 18754
NZ_SMTD01000004.1 264659
NZ_SMTD01000050.1 16037
NZ_SMTD01000051.1 14204
NZ_SMTD01000052.1 13923
NZ_SMTD01000053.1 13222
NZ_SMTD01000054.1 10583
NZ_SMTD01000055.1 9367
NZ_SMTD01000056.1 8230
NZ_SMTD01000057.1 8199
NZ_SMTD01000058.1 6804
NZ_SMTD01000059.1 4367
NZ_SMTD01000005.1 262082
NZ_SMTD01000060.1 4267
NZ_SMTD01000061.1 3317
NZ_SMTD01000062.1 2970
NZ_SMTD01000063.1 2729
NZ_SMTD01000064.1 2447
NZ_SMTD01000065.1 2277
NZ_SMTD01000066.1 2258
NZ_SMTD01000067.1 1894
NZ_SMTD01000068.1 1742
NZ_SMTD01000006.1 213735
NZ_SMTD01000069.1 1735
NZ_SMTD01000070.1 1711
NZ_SMTD01000071.1 1708
NZ_SMTD01000072.1 1708
NZ_SMTD01000073.1 1344
NZ_SMTD01000074.1 1283
NZ_SMTD01000075.1 1255
NZ_SMTD01000076.1 1179
NZ_SMTD01000077.1 1124
NZ_SMTD01000078.1 1061
NZ_SMTD01000007.1 181795
NZ_SMTD01000079.1 1058
NZ_SMTD01000080.1 848
NZ_SMTD01000081.1 754
NZ_SMTD01000082.1 734
NZ_SMTD01000083.1 651
NZ_SMTD01000084.1 574
NZ_SMTD01000085.1 569
NZ_SMTD01000086.1 541
NZ_SMTD01000087.1 533
NZ_SMTD01000008.1 176290
NZ_SMTD01000009.1 159375
```

**That order is neither lexicographic nor size-sorted** — it is the order NCBI ships the
`_genomic.fna.gz` in, and the package preserves it end to end. Two sequences share a
length (`NZ_SMTD01000071.1` and `NZ_SMTD01000072.1`, both 1708), so length is not a key.

### `ecHT115.fa` properties

| property | value |
| --- | --- |
| line width | 80 bases |
| case | **all upper-case** — zero soft-masking |
| alphabet | strictly `ACGT`; **no `N`, no IUPAC ambiguity codes at all** |
| `.fa` size | 4,669,182 B |

### GTF `refseq_rs_2025_06_26`

Source: `.../GCF_004354945.1_ASM435494v1_genomic.gtf.gz`, provider RefSeq, version
`RS_2025_06_26`, registered with `chromosomes_checked: true`.

**87 distinct seqnames — the set is exactly equal to `chrom.sizes` (verified by `diff` of
the sorted sets), spelled identically, accession-with-version and all.** The GTF also emits
its records in the same order as the FASTA. Every contig, including the 533 bp smallest,
carries at least one feature. Unpacked GTF is 7,659,939 B.

So for `ecHT115` the annotation is a **total** match, not a subset: any chimera merge must
keep all 87 spellings intact or `ChromosomeMismatchError` fires for all 87.

---

## 2. `ce11` — WormBase WS298 (PRJNA13758)

**7 sequences, 100,286,401 bp total.** Confirms what the map already had.

### `ce11.chrom.sizes` — verbatim, in file order

```text
V 20924180
X 17718942
IV 17493829
II 15279421
I 15072434
III 13783801
MtDNA 13794
```

Byte-identical to columns 1–2 of `ce11.fa.fai` and to STAR's `chrNameLength.txt`. This
order is **size-descending**, which is how WormBase ships the FASTA. Lexicographic order
would be `I II III IV MtDNA V X` — different.

### `ce11.fa` properties

| property | value |
| --- | --- |
| line width | 60 bases |
| case | **soft-masked** — 22,011,127 lower-case bases, 21.95% of the assembly |
| alphabet | strictly `ACGTacgt`; **no `N`, no IUPAC ambiguity codes** |
| `.fa` size | 101,957,983 B |

### GTF `wormbase_ws298`

Source: `.../c_elegans.PRJNA13758.WS298.canonical_geneset.gtf.gz`, provider WormBase,
version `WS298`, registered with `chromosomes_checked: true`. Unpacked GTF is
182,928,275 B.

**7 distinct seqnames, set-equal to `chrom.sizes`, spelled identically** (bare roman
numerals and `MtDNA`, no `chr` prefix). Feature-line counts:

| seqname | GTF lines |
| --- | ---: |
| `I` | 98,364 |
| `II` | 105,695 |
| `III` | 88,602 |
| `IV` | 159,254 |
| `MtDNA` | 136 |
| `V` | 139,258 |
| `X` | 116,122 |

---

## 3. For the naming contract (#42)

Characters, lengths and orders a suffix separator has to survive.

**Character inventory of the names as they exist today:**

| assembly | distinct characters across all names | name lengths |
| --- | --- | --- |
| `ecHT115` | `. 0 1 2 3 4 5 6 7 8 9 D M N S T Z _` | all exactly 17 |
| `ce11` | `A D I M N V X t` | 1, 2, 3, 5 (`I`, `IV`, `III`, `MtDNA`) |

So the existing corpus already contains `_`, `.`, digits, upper-case and one lower-case
letter (`t` in `MtDNA`). It contains **no** `-`, `|`, `~`, `:`, `@`, `#`, `+`, `;`, `/`,
whitespace, or any other punctuation.

**Component assembly names are clean.** Every `assembly_name` in
`src/genome/data/assembly_metadata.tsv` — `hg38 hg19 mm39 mm10 sacCer3 ce11 ecHT115` — is
strictly alphanumeric. No shipped assembly name contains a separator candidate.

**What that means for parsing.** Suffixing every name with `<sep><assembly>`:

| candidate `<sep>` | already in a chromosome name? | recover assembly by `rsplit(sep, 1)` | recover by `split(sep, 1)` |
| --- | --- | --- | --- |
| `_` | **yes** (`NZ_SMTD…`) | works | **fails** |
| `.` | **yes** (`…000010.1`) | works | **fails** |
| `-` | no | works | works |
| `\|` | no | works | works |
| `~` | no | works | works |
| `@` | no | works | works |
| `:` | no | works | works — but see below |

Because component assembly names are alphanumeric, a *last-occurrence* split recovers the
assembly for every candidate including `_` and `.`. A *first-occurrence* split is only
safe for a separator absent from the corpus. The regex a consumer writes must therefore be
right-anchored (`(?:^|.*)<sep>(ecHT115|ce11)$`), not left-anchored.

**Anchoring matters for `ce11` independently of the separator.** `I` is a prefix of `II`
and `III`; `I`, `V` and `X` are all substrings of other names. Any attribution or lookup
regex over chimera names must anchor both ends — an unanchored `I_ce11` matches inside
`III_ce11`.

**Longest name after suffixing** the current corpus: 25 characters
(`NZ_SMTD01000010.1_ecHT115`). All 94 suffixed names are unique under every candidate
separator.

**SAM legality** — the SAM spec (hts-specs `SAMv1.tex`, "Character set limitations")
allows in a reference name:

```text
[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*
```

i.e. any printable ASCII except backslash, comma, quotation marks, backtick, apostrophe,
`()`, `[]`, `{}`, `<>`, and it may not start with `*` or `=`. Every candidate above is
legal. Two footnoted cautions from that same section: `|` "historically appeared in
reference names derived from NCBI FASTA files", and `:` is legal but collides with
`name:begin-end` region notation — `samtools faidx`, `bedtools` and the `Genome.sequence`
region string would all have to parse it out of ambiguity. `.` carries the analogous
hazard for anything that treats a trailing dot-suffix as a version or file extension.

**Sort order.** Neither assembly's `chrom.sizes` is name-sorted today, and the two use
different orders (`ecHT115` = NCBI source order, which is arbitrary; `ce11` =
size-descending). A suffix changes lexicographic order for `ce11` only in that all seven
gain a common tail — relative order among them is unchanged for any single separator. But
if a chimera's `chrom.sizes` were ever *sorted* by name rather than concatenated in
component order, the two components would **interleave** for a separator that sorts before
alphanumerics (e.g. `-`, ASCII 0x2D) versus staying grouped for one that sorts after
(`_` 0x5F, `~` 0x7E). Concatenating in component order sidesteps that entirely and is what
both source FASTAs already do internally.

---

## 4. For the index ticket (#48)

**Totals, and what the package's own formula yields.** `star.py` computes
`genomeSAindexNbases = min(14, max(2, int(log2(chrom_sizes.sum()) / 2 - 1)))`, matching
STAR's documented `min(14, log2(GenomeLength)/2 - 1)`.

| | sequences | total bp | `genomeSAindexNbases` |
| --- | ---: | ---: | ---: |
| `ecHT115` | 87 | 4,602,240 | **10** (as built) |
| `ce11` | 7 | 100,286,401 | **12** (as built) |
| chimera `ecHT115` + `ce11` | **94** | **104,888,641** | **12** (unchanged from `ce11` alone) |

Both existing indexes were built with STAR 2.7.11b and chromap 0.3.2-r518, `--sjdbOverhang
100`, `--runThreadN 48`; the values above are read back from the committed
`genomeParameters.txt` and `.completion.json`, so the formula is confirmed against what
actually ran.

**`genomeChrBinNbits` is the parameter this map actually disturbs, and the package does
not set it.** STAR's `--help` (2.7.11b, verbatim):

> `genomeChrBinNbits 18` — int: =log2(chrBin), where chrBin is the size of the bins for
> genome storage: each chromosome will occupy an integer number of bins. For a genome with
> large number of contigs, it is recommended to scale this parameter as
> `min(18, log2[max(GenomeLength/NumberOfReferences, ReadLength)])`.

| | `GenomeLength/NumberOfReferences` | recommended `genomeChrBinNbits` | used |
| --- | ---: | ---: | ---: |
| `ecHT115` | 52,899 | **15** | 18 |
| `ce11` | 14,326,629 | 18 | 18 |
| chimera | 1,115,837 | **18** | — |

Two consequences:

1. **The chimera does not need it.** At 94 references and 104.9 Mbp the recommendation is
   still 18, which is STAR's default. Concatenating `ecHT115` into `ce11` does not push
   `genomeChrBinNbits` off its default.
2. **`ecHT115` standalone is already over-binned today, and it is measurable.** At
   `genomeChrBinNbits=18` the 87 contigs occupy 91 bins of 262,144 B — and the built
   `Genome` file is exactly 91 × 262,144 = **23,855,104 B**, i.e. 5.2× the 4.6 Mbp it
   stores, all of it padding. At the recommended 15 it would be 197 bins × 32,768 =
   6,455,296 B. `chrStart.txt` confirms the last bin boundary at 23,855,104. STAR emitted
   no warning about this.

Index sizes as built, for scale:

| | STAR `Genome` | STAR `SA` | STAR `SAindex` | chromap index |
| --- | ---: | ---: | ---: | ---: |
| `ecHT115` | 23,855,104 | 37,968,483 | 6,116,787 | 34,220,408 |
| `ce11` | 124,598,858 | 1,019,544,912 | 97,867,203 | 591,765,632 |

`ecHT115`'s GTF yields **zero splice junctions** (`sjdbList.fromGTF.out.tab` and
`sjdbList.out.tab` are both 0 bytes, `sjdbInfo.txt` is 6 bytes) — expected for a bacterium,
but worth knowing: a merged chimera annotation contributes junctions from the `ce11` half
only.

---

## Method and caveats

- Everything above is read from files already on GPU71FM under
  `/home/hanliu/liulab_data/genome/{ecHT115,ce11}/`; no compute was run on the login node
  beyond `grep`/`cut`/`diff`/`awk` over these files and `STAR --help`.
- The NCBI assembly report was fetched fresh from
  `https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/354/945/GCF_004354945.1_ASM435494v1/GCF_004354945.1_ASM435494v1_assembly_report.txt`.
- The SAM character set is quoted from `samtools/hts-specs` `SAMv1.tex` at `master`.
- Chimera figures (94 sequences, 104,888,641 bp, `genomeSAindexNbases` 12,
  `genomeChrBinNbits` 18) are arithmetic over the two measured assemblies, not something
  built and observed. No chimera was constructed.
- These facts are pinned by checksum: both assemblies and both annotations verified
  against the `sha256` in the shipped tables, so they will not drift unless a table row
  changes.

Written up as `docs/research/echt115-ce11-contents.md` on branch
`research/echt115-ce11-contents` (not opened as a PR).
