# What do STAR and chromap do with genome size and contig count?

Measured 2026-08-13 for [#40](https://github.com/liuhlab/liulab-genome/issues/40). **STAR's
`genomeSAindexNbases` is computed from the assembly and reproduces STAR's own source exactly,
including its truncation; `genomeChrBinNbits` is never passed, so STAR's default of 18 always holds
— which pads the 87-scaffold `ecHT115` draft from 4.6 Mb to 23.9 Mb, and STAR warns on none of it.
chromap has no size- or count-dependent index parameter at all. `_` is legal anywhere in a SAM
reference name and there is no length limit, but both tools truncate a FASTA header at the first
whitespace, so a suffix must ride on the first token.**

This is a measurement, not a decision. The decisions it fed are
[#42](https://github.com/liuhlab/liulab-genome/issues/42), which owns the naming contract, and
[#48](https://github.com/liuhlab/liulab-genome/issues/48), which decided what a chimera changes
about an index build. Line numbers are against `main` at 0f644e1.

Two names below have since moved, and the measurement is left as it was taken rather than rewritten:
`_kwargs_to_flags` is now the one renderer `Aligner._flags`, and `Genome.register_gtf` is now
`genome.annotations.register_path`. What was measured about STAR and chromap is unaffected — the
numbers came from the tools, not from this package.

## TL;DR

- **STAR: `genomeSAindexNbases` is computed, `genomeChrBinNbits` is not passed at all.** The repo's
  formula is `src/genome/aligner/star.py:143-145` and it matches STAR's *own source* exactly,
  including the truncation — but is one lower than the STAR *manual's* worked example, because the
  manual rounds and the code truncates. Both are primary; they disagree with each other.
- **chromap has no genome-size or contig-count parameter at all.** `-k 17 / -w 7` are constants;
  nothing in `chromap -i` scales by size or reference count. The repo passes neither unless asked.
- **For `ce11` + `ecHT115` specifically, both parameters land on the value they already have.**
  The chimera is 94 sequences / 104,888,641 bp → `genomeSAindexNbases` 12 (same as `ce11` alone),
  `genomeChrBinNbits` guidance 18 (= the default). Nothing about *this* chimera moves either knob.
- **But `ecHT115` alone is the case that is wrong today**, and it is wrong already, chimera or not:
  the pinned RefSeq assembly is **87 unplaced WGS scaffolds**, not one chromosome, and at the
  default `genomeChrBinNbits=18` STAR pads it from 4.6 Mb to **23.9 Mb — a 5.18× inflation**.
- **Correction to a charted fact.** The map (#38) records that `ecHT115` sequences are named
  `NZ_CP…`. They are not. `GCF_004354945.1` is a scaffold-level draft; its names are
  `NZ_SMTD01000001.1` … `NZ_SMTD01000087.1`.
- **For #42: `_` is legal anywhere in a SAM reference name, first character included, and so are
  `.` `-` `|` `:` `~` `+` `#` `@` `/`.** Illegal everywhere: `\` `,` `"` `'` `` ` `` `()` `[]` `{}`
  `<>` and all whitespace; `*` and `=` are legal but not as the first character. There is no length
  limit. Both `I_ce11` and `NZ_SMTD01000001.1_ecHT115` are valid and reach `@SQ SN:` byte-for-byte
  through both aligners. **But both aligners truncate the FASTA header at the first whitespace**, so
  the suffix must go on the first token, not after the description.
- **The layout and the completion record need nothing.** Neither is contig-count- or
  name-dependent.

---

## 1. What `src/genome/aligner/star.py` does

**`genomeSAindexNbases` is computed from the assembly**, not hard-coded and not omitted —
`src/genome/aligner/star.py:141-145`:

```python
# STAR requires a reduced suffix-array index size for small genomes:
#   min(14, log2(genomeLength) / 2 - 1). Default it unless overridden.
if "genomeSAindexNbases" not in kwargs:
    genome_length = int(self._genome.chrom_sizes.sum())
    kwargs["genomeSAindexNbases"] = min(14, max(2, int(math.log2(genome_length) / 2 - 1)))
```

Three things about this:

- The input is `Genome.chrom_sizes.sum()` — the sum of true sequence lengths, read from
  `chrom.sizes` (`src/genome/genome.py:169`, `src/genome/io/fasta.py:202-238`). That is **exactly
  the quantity STAR itself uses** for its own recommendation (`nGenomeTrue`, excluding bin padding).
  The basis is right.
- It is a **default, not a floor**: an explicit `genomeSAindexNbases=` kwarg wins, which is how the
  docs describe it (`docs/aligner.md:236-238`).
- It is passed through the generic kwargs path, so it reaches STAR as `--genomeSAindexNbases <n>`
  via `_kwargs_to_flags` (`src/genome/aligner/star.py:176-188`).

**`genomeChrBinNbits` is not mentioned anywhere in the repo.** `grep -rn "ChrBinNbits"` over the
whole tree returns nothing — not in `src/`, `tests/` or `docs/`. STAR's default of `18` is therefore
always in force unless a caller passes it by hand as a kwarg.

**Nothing else about the assembly is read.** The build command is assembled in
`STAR.index` (`src/genome/aligner/star.py:132-173`) from: the FASTA path, the resolved GTF path,
`--sjdbOverhang`, `--runThreadN`, `--genomeDir`, plus kwargs. `src/genome/aligner/aligner.py` and
`src/genome/aligner/mixin.py` add no aligner flags at all — `mixin.py:22-45` is a pass-through to
`STAR(...).index(**kwargs)`, and `aligner.py` contributes only layout, installation checking,
`_run` (`aligner.py:212-236`) and the completion record. **Contig count is never read by anything.**

## 2. What `src/genome/aligner/chromap.py` does

**Nothing is computed and nothing is hard-coded — both knobs are omitted unless the caller asks.**
`src/genome/aligner/chromap.py:124-138`: `kmer` and `window` default to `None` and are only added to
the command when explicitly passed. The command is `--build-index --ref <fasta> --output
<index_dir>/chromap.index` plus whatever the caller supplied.

This is correct, because **chromap has no genome-size or contig-count parameter to compute.** From
`chromap -h`, chromap 0.3.2-r518, the repo's own `aligners` pixi environment — the complete set of
indexing options:

```
 Indexing options:
  -i, --build-index          Build index
      --min-frag-length INT  Min fragment length for choosing k and w automatically [30]
  -k, --kmer INT             Kmer length [17]
  -w, --window INT           Window size [7]
```

Three options, and the only automatic scaling chromap does (`--min-frag-length`) is driven by
expected **fragment length**, not by genome size or reference count. There is no chromap analogue of
`genomeSAindexNbases` or `genomeChrBinNbits`.

## 3. What STAR's documentation actually says

Quoted verbatim. Two primary sources that **do not agree with each other**, which matters for #48.

### `genomeSAindexNbases` (default `14`)

`source/parametersDefault` lines 71-72
(<https://raw.githubusercontent.com/alexdobin/STAR/master/source/parametersDefault>), identical to
the text `STAR --help` prints in our own environment:

> `genomeSAindexNbases         14`
> `    int: length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer
> strings will use much more memory, but allow faster searches. For small genomes, the parameter
> --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).`

The manual, section "Very small genome.",
`extras/doc-latex/STARmanual.tex` lines 157-158
(<https://github.com/alexdobin/STAR/blob/master/extras/doc-latex/STARmanual.tex>) — note
`doc/STARmanual.tex` does not exist; `doc/` holds only the compiled PDF:

> For small genomes, the parameter `--genomeSAindexNbases` **must** to be scaled down, with a
> typical value of `min(14, log2(GenomeLength)/2 - 1)`. For example, for 1 megaBase genome, this is
> equal to 9, for 100 kiloBase genome, this is equal to 7.

("must to be" is the manual's own grammatical slip; "chrosomes" below is likewise its own typo.)

**STAR does not compute this for you — it only warns.** `Genome_genomeGenerate.cpp:167-172`:

```cpp
if (pGe.gSAindexNbases > log2(nGenomeTrue)/2-1) {
    ostringstream warnOut;
    warnOut << "--genomeSAindexNbases " << pGe.gSAindexNbases << " is too large for the genome size=" << nGenomeTrue;
    warnOut << ", which may cause seg-fault at the mapping step. Re-run genome generation with recommended --genomeSAindexNbases " << int(log2(nGenomeTrue)/2-1);
    warningMessage(warnOut.str(),P.inOut->logMain,std::cerr,P);
};
```

That string is present in the binary we ship (`strings $(command -v STAR)` in the `aligners` env
returns `", which may cause seg-fault at the mapping step. Re-run genome generation with recommended
--genomeSAindexNbases "`), and `warningMessage` prefixes it with `!!!!! WARNING: `
(`ErrorWarning.cpp:25-38`). So a too-large value is a warning to stderr and `Log.out`, never an
error — the index builds and may seg-fault at mapping time.

**The manual and the code disagree by one.** `log2(1e6)/2 - 1 = 8.966`. The manual's worked example
says **9** (it rounds); STAR's own recommendation is `int(...)` = **8** (it truncates, and has no
`min(...,14)` cap). **The repo's `int(math.log2(L)/2 - 1)` reproduces STAR's code exactly**, and adds
the `min(14, ...)` cap the code lacks. So the repo is not wrong so much as it has picked one of two
disagreeing primary sources. Divergence between the two only appears for a handful of sizes —
`sacCer3` at 12.16 Mb is one (exact 10.768: manual→11, repo→10).

Edge case worth naming: the repo's `max(2, ...)` floor can exceed STAR's warning threshold, which
fires whenever the value is greater than the exact real-valued `log2(L)/2 - 1`. That threshold drops
below 2 only for L < 64 bp, so the floor is inert for anything real.

### `genomeChrBinNbits` (default `18`)

`source/parametersDefault` lines 68-69, again identical to `STAR --help` here:

> `genomeChrBinNbits           18`
> `    int: =log2(chrBin), where chrBin is the size of the bins for genome storage: each chromosome
> will occupy an integer number of bins. For a genome with large number of contigs, it is
> recommended to scale this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).`

The manual, section "Genome with a large number of references.", `STARmanual.tex` lines 160-161:

> If you are using a genome with a large (>5,000) number of references (chrosomes/scaffolds), you
> may need to reduce the `--genomeChrBinNbits` to reduce RAM consumption. The following scaling is
> recommended: `--genomeChrBinNbits = min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])`.
> For example, for 3 gigaBase genome with 100,000 chromosomes/scaffolds, this is equal to 15.

Note the difference in force: `genomeSAindexNbases` is "**must**", `genomeChrBinNbits` is "you may
need to", and its stated trigger is **>5,000 references**.

**STAR emits no warning at all for `genomeChrBinNbits`** — no warning, no error, no clamping. Unlike
the SA parameter, getting this one wrong is silent.

**The cost is padding, and the manual never quantifies it.** `genomeScanFastaFiles.cpp:50-52` pads
each chromosome up to a whole multiple of `2^genomeChrBinNbits`:

```cpp
if (N>0) {//pad the chromosomes to bins boudnaries
    N = ( (N+1)/mapGen.genomeChrBinNbases+1 )*mapGen.genomeChrBinNbases;
};
```

so the worst case is `N_references × 2^genomeChrBinNbits` wasted bases — 262,144 bases per reference
at the default. (5,000 × 256 kB ≈ 1.3 GB, which is where the manual's ">5,000" threshold comes
from.) STAR logs both figures, so the overhead is directly observable in `Log.out`
(`Genome_genomeGenerate.cpp:150-151`): `Genome sequence total length = ` and
`Genome size with padding = `.

**No memory formula for either parameter exists in any STAR document.** The manual says only "much
more memory" and "to reduce RAM consumption". Anything quoting a formula is secondary.

## 4. The numbers, for the assemblies this map actually names

`ecHT115` is pinned to RefSeq `GCF_004354945.1` / ASM435494v1 (`src/genome/data/assembly_metadata.tsv:8`).
From its NCBI assembly report
(<https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/354/945/GCF_004354945.1_ASM435494v1/GCF_004354945.1_ASM435494v1_assembly_report.txt>):

- **87 sequences, all role `unplaced-scaffold`**, total **4,602,240 bp**; longest 327,216, shortest 533.
- Names are `NZ_SMTD01000001.1` … `NZ_SMTD01000087.1` — **not `NZ_CP…`**. Longest name 17 chars.

`ce11` sizes are the ones already charted on #38 (7 sequences, 100,286,401 bp).

| | sequences | length | `genomeSAindexNbases` exact | repo passes | `genomeChrBinNbits` guidance |
|---|---:|---:|---:|---:|---:|
| `ce11` | 7 | 100,286,401 | 12.290 | **12** | 18.00 → 18 (= default) |
| `ecHT115` | 87 | 4,602,240 | 10.067 | **10** | **15.69 → 15** |
| chimera `ce11`+`ecHT115` | 94 | 104,888,641 | 12.322 | **12** | 18.00 → 18 (= default) |

(`genomeChrBinNbits` guidance computed as `min(18, log2(max(L/N, 100)))`, ReadLength = 100.)

Padding at various `genomeChrBinNbits`, computed as `Σ ceil(len/2^b) · 2^b` over the real sequence
lengths:

| | raw | b=18 (default) | b=15 | inflation at default |
|---|---:|---:|---:|---:|
| `ce11` | 100,286,401 | 101,187,584 | 100,368,384 | **×1.01** |
| `ecHT115` | 4,602,240 | 23,855,104 | 6,455,296 | **×5.18** |
| chimera | 104,888,641 | 125,042,688 | 106,823,680 | **×1.19** |

**Read that table carefully — the chimera is not the problem.** Adding 87 small scaffolds to a
100 Mb nematode costs ~20 MB of padding, ×1.19, and leaves the guidance formula sitting exactly on
the default because the nematode dominates `GenomeLength/NumberOfReferences`. The formula, applied
to the chimera, tells you to change nothing.

The case that is already mis-parameterised today is **`ecHT115` on its own**: 87 scaffolds at the
default bin size inflate a 4.6 Mb genome to 23.9 Mb. That is a bug the package has now, for any
fragmented small assembly, independent of this map. Whether it is worth fixing is #48's call.

Two consequences that follow, for #48 to weigh:

- Applying `min(18, log2(max(L/N, ReadLength)))` naively would need a ReadLength the index build
  does not have. `sjdb_overhang` (default 100, `star.py:87`) is `read_length - 1` by construction
  and is the only read-length-shaped number in scope.
- A chimera *dilutes* the signal, and here is the number: for this chimera's 104,888,641 bp, the
  merged-assembly guidance stays pinned at 18 until **N > 400** sequences. At N = 507 it is 17.66,
  at N = 1,007 it is 16.67, at N = 5,007 it is 14.35. A formula over the merged assembly cannot see
  that one component is 87 scaffolds; a per-component view would. Whether that matters is a
  decision, not a fact.

## 5. Reference-sequence names — written for #42

### The SAM specification's grammar — the answer #42 needs

**`_` is legal anywhere in a reference name, including as the first character. So is `.`, and so is
`-`, `|`, `:`, `~`, `+`, `#`, `@` and `/`.** Both candidate spellings —
`I_ce11` and `NZ_SMTD01000001.1_ecHT115` — are valid `@SQ SN:` values with no caveat.

SAM v1, §1.2.1 "Character set restrictions"
(<https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex>, lines 206-221):

> Reference sequence names may contain any printable ASCII characters in the range `[!-~]` apart
> from backslashes, commas, quotation marks, and brackets — i.e., apart from
> `` \ , " ` ' () [] {} <> `` — and may not start with `*` or `=`.
>
> Thus they match the following regular expression:
>
> `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*`
>
> For clarity, elsewhere in this specification we write this set of allowed characters as a
> character class `[:rname:]` and extend the POSIX regular expression notation to use `^*=` to
> indicate the omission of `*` and `=` from the character class. Thus this regular expression can
> be written more clearly as `[:rname:^*=][:rname:]*`.

A footnote on the same page names two of the separators explicitly as deliberately permitted:

> Characters that are *not* disallowed include `|`, which historically appeared in reference names
> derived from NCBI FASTA files, and `:`, which appears in HLA allele names.

And the `@SQ SN` row itself (lines 286-290):

> Reference sequence name. The `SN` tags and all individual `AN` names in all `@SQ` lines must be
> distinct. […] Regular expression: `[:rname:^*=][:rname:]*`

| character | verdict |
|---|---|
| `_` `.` `-` `\|` `:` `~` `+` `#` `@` `/` `!` `$` `%` `&` `;` `?` `^` and alphanumerics | **legal anywhere, including first character** |
| `=` `*` | legal, but **not as the first character** |
| `\` `,` `"` `'` `` ` `` `(` `)` `[` `]` `{` `}` `<` `>` | **illegal everywhere** (named in the prose exclusion list) |
| space and all whitespace | **illegal everywhere** (outside the printable `[!-~]` range) |

**No length limit.** The regex is unbounded, deliberately: the spec bounds QNAME explicitly
(`[!-?A-~]{1,254}`, line 427) and says of query names "(They are also limited in length.)" — no
such sentence or bound exists for RNAME. In BAM, `l_name` is a `uint32_t` marked *limited*, meaning
"limited by available memory […] well before they are limited by […] Java's signed 32-bit integer
maximum array size" (lines 1060-1061, 1096-1098). At 25 characters this is not a live constraint.

**Number of references: `n_ref` < 2^31** in BAM (line 1057-1058), unbounded in text SAM.

**`@SQ` order is load-bearing** — worth knowing if a chimera's components are ever reordered.
Line 285: "The order of `@SQ` lines defines the alignment sorting order." And the `@HD SO`
definition (lines 263-266): "For coordinate sort, the major sort key is the RNAME field, **with
order defined by the order of `@SQ` lines in the header**." Also, RNAME "(if not `*`) must be
present in one of the `SQ-SN` tag" (line 517-519) — so a name that gets mangled between FASTA and
header is not a cosmetic problem.

One implication for the naming contract that follows straight from the grammar: **whitespace is
illegal, and both STAR and chromap truncate a FASTA header at the first whitespace anyway** (see
below). So the suffix must be attached to the *first token* of the header, not appended after the
description.

### STAR

- **Truncates the FASTA header at the first whitespace.** `genomeScanFastaFiles.cpp:38-45` reads the
  header with `lineInStream.ignore(1,' ')` then `>>` into a `std::string`, so `>chr1 some description`
  becomes `chr1`. **The manual never states this**; it is only readable from the code.
- **The only documented character restriction** is `STARmanual.tex:122`:
  > The tabs are not allowed in chromosomes' names, and spaces are not recommended.

  Nothing anywhere in the manual, `parametersDefault` or the source restricts `_`, `.`, `-` or any
  other character. On those, STAR is **silent** — no validation, no rejection, no escaping,
  no transformation.
- **No length cap at build time** (names are `std::string`), but an **undocumented 999-character cap
  at load time**: `Genome.cpp:150-158` reads `chrName.txt` through a `char chrInChar[1000]` buffer,
  so a name of 1000+ characters silently ends the read loop. Irrelevant at 25 characters, but it is
  the only cap that exists.
- **No reordering, no sorting.** Names are `push_back`-ed in FASTA order, files in command-line
  order. Nothing in the codebase sorts `chrName`. The manual treats the order as load-bearing
  (`STARmanual.tex:131`): you may rename entries in `chrName.txt` "while **keeping the order** of the
  chromosomes in this file".
- **The name reaches `@SQ SN:` intact.** `samHeaders.cpp:28-31`:
  ```cpp
  samHeaderStream << "@SQ\tSN:"<< genomeOut.chrName.at(ii) <<"\tLN:"<<genomeOut.chrLength[ii]<<"\n";
  ```
  Verbatim, in genome order, no transformation. So whatever survives the whitespace truncation is
  exactly what downstream attribution reads back out.
- **The GTF must agree with the FASTA on names.** STAR's own error text says so:
  > Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome
  > naming between GTF and FASTA file.

  For a chimera this means the merged GTF's `seqname` column must carry the *same* suffixes as the
  FASTA. STAR offers `--sjdbGTFchrPrefix` for a uniform prefix, which is no help for a per-component
  suffix.

### chromap

- **Zero validation.** No character restriction, no length cap, no rejection, no deduplication —
  grepping `sequence_batch.cc` / `.h` and `index.cc` for validation turns up only file-open errors.
  Names are never inspected.
- **Truncates at the first whitespace**, same as STAR. `sequence_batch.cc:36-37` swaps in kseq's
  `name` field, and the bundled `kseq.h` splits the header on `KS_SEP_SPACE`; the remainder is kept
  as `comment` and never written anywhere.
- **Names are not stored in the index at all.** `Index::Save()` (`index.cc:91-130`) writes exactly
  four things: `kmer_size_`, `window_size_`, the lookup table and the occurrence table. No names, no
  lengths, no contig count. That is why `-r ref.fa` is required at *mapping* time as well as index
  time — `chromap.h:640-644` re-reads the FASTA to recover the names. **Consequence for a chimera:
  the naming contract is entirely a property of the FASTA, and re-suffixing sequence names would not
  invalidate a chromap index — but it would silently change what the SAM header says while the index
  file stays byte-identical.**
- **`@SQ SN:` gets the name byte-for-byte.** The only `@SQ` in the whole source tree,
  `mapping_writer.cc:311-321`:
  ```cpp
  this->AppendMappingOutput(
      "@SQ\tSN:" + std::string(reference_sequence_name) +
      "\tLN:" + std::to_string(reference_sequence_length) + "\n");
  ```
  (chromap emits no `@HD` and no `@PG` — the header is `@SQ` lines only.)
- **No reordering by default.** rids are assigned in FASTA file order; reordering happens only when
  `--chr-order FILE` is passed, and that file is matched by exact whole-line string equality against
  the truncated name.
- **Hard limit of 2^31 − 1 contigs**, from minimizer hit bit-packing (`minimizer_generator.cc:59`,
  `hit_utils.h`), and any single contig must be < 2^32 bases. Both irrelevant here.
- **The contig-count sensitivity chromap does have is in the mapping path, not the index.**
  `chromap.h:336-357` allocates `2 × num_threads × num_contigs` `std::vector`s before mapping a
  single read. At 94 contigs this is nothing; on a heavily fragmented reference it is not. Out of
  scope for this package, which never maps — but relevant to whoever runs the alignment.

### The repo's own name check, which fires before either aligner does

`Genome.register_gtf(..., check_chromosomes=True)` already refuses a GTF naming any sequence the
assembly does not have — `_reject_unknown_chromosomes` at `src/genome/io/gtf.py:850-864`, raising
`ChromosomeMismatchError` (`gtf.py:88`). It is one-directional: extra assembly sequences the GTF
never mentions are fine, unknown GTF names are not. **So a chimera whose merged GTF is not suffixed
in lockstep with the FASTA fails at annotation registration, before STAR is ever invoked** — which is
the right place for it to fail, and means the index build inherits a guarantee rather than needing
its own check.

## 6. Layout and the completion record

**Both work untouched. Nothing here is contig-count- or name-aware.**

- `index_dir` is `assembly_data_dir(<assembly>) / "index" / <name>` (`aligner.py:102-108`,
  `download.py:85`). A chimera is an assembly with its own name, so it gets its own data dir and its
  own `index/` for free.
- `index/star_<gtf>/` (`star.py:50-54`) keys on the **annotation name**, not on anything about the
  sequences. A chimera carrying one merged annotation registered under one key gets one directory,
  exactly like any other assembly. If the merge instead registered N per-component GTFs, this layout
  would build N indexes — but that is the annotation-merge ticket's problem, not the index's.
- `index/chromap/chromap.index` (`chromap.py:59-62`, `_build_arguments` returns `""` at
  `chromap.py:65-67`) — one per assembly, no selector, unchanged.
- The completion record's `details` are `aligner`, `binary`, `assembly`, `fasta`, `command`,
  `parameters` (`aligner.py:268-285`). **No chromosome names and no contig count are stored.** The
  one thing that would change for a chimera is that `parameters` would record whatever
  `genomeChrBinNbits` value a fix for #48 introduces — which is exactly what the record is for, and
  needs no schema change.
- `_claimed_files` (`aligner.py:238-257`) `rglob`s the whole index dir. STAR's *genomeDir* has a
  fixed file list regardless of contig count — `strings` on our STAR binary gives the complete set
  (`Genome`, `SA`, `SAindex`, `chrName.txt`, `chrNameLength.txt`, `chrLength.txt`, `chrStart.txt`,
  `exonGeTrInfo.tab`, `exonInfo.tab`, `geneInfo.tab`, `sjdbInfo.txt`, `sjdbList.fromGTF.out.tab`,
  `sjdbList.out.tab`, `transcriptInfo.tab`), one file each and never one per contig — so the record
  does not grow with the chimera.

## Method and sources

- Repo facts read directly from the tree at `main` 0f644e1.
- STAR facts: `source/parametersDefault`, `extras/doc-latex/STARmanual.tex`, and the C++ sources
  named inline, all at `alexdobin/STAR` master commit `b1edc1208d91a53bf40ebae8669f71d50b994851`
  (2024-01-25), which is identical to release **2.7.11b**. Cross-checked against `STAR --help` and
  `strings` on the binary in this repo's own `aligners` pixi environment (STAR **2.7.11b**).
- chromap facts: `chromap -h` from the same environment (chromap **0.3.2-r518**), plus the README,
  the `chromap.1` manpage and the `src/` sources named inline at `haowenz/chromap` master (the
  manpage self-identifies as `chromap-0.2.6 (r490)`, 2024-01-25). The three indexing options are
  identical between the manpage and the 0.3.2-r518 `--help` we ship, so the source reading applies.
- SAM grammar from the spec's own LaTeX source,
  <https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex> (the master for
  <https://samtools.github.io/hts-specs/SAMv1.pdf>). The restriction dates from Jan 2019 per the
  spec's changelog.
- `ecHT115` sequence count, names and lengths from the NCBI assembly report, linked above.
- Arithmetic is a short script over those published lengths; no index was built.

