# Test fixtures

Small, subsampled **real** files — never a large genomic file. Every sequence byte here came from
UCSC's `sacCer3` golden path and every motif byte from a published JASPAR release; both were cut
down on the lab cluster, and nothing is synthesised. Names are this repo's where a fixture needs a
spelling `sacCer3` does not have, and one 200-base stretch is lower-cased — both are called out
below, and nothing else departs from the source bytes.

| File | What it is |
|---|---|
| `tiny.fa` | The first 10 000 bases of `chrI`, `chrII` and `chrIII` of `sacCer3`, cut with `samtools faidx` and re-headed to the bare chromosome names |
| `tiny.fa.gz` | `tiny.fa`, gzipped — the compressed-source path |
| `tiny.gtf` | The 85 `sacCer3.ensGene` features that fall wholly inside those three windows |
| `tiny.gtf.gz` | `tiny.gtf`, gzipped |
| `ensembl_style.gtf` | `tiny.gtf` with the `chr` prefix stripped (`I`, `II`, `III`) — an Ensembl-spelled annotation for the chromosome-name mismatch case |
| `tiny_jaspar_transfac.txt` | Ten whole records of the JASPAR 2024 `all` union file, in its own order — see below |

Sources:

- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz`
- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz`
- `https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_transfac.txt`

`tiny.gtf` coordinates are 1-based inclusive, as every GTF is; they convert at the I/O boundary and
are never seen in that form inside the package.

## `tiny_jaspar_transfac.txt` — the motif records

Ten records of the JASPAR 2024 `all` union file, copied whole and kept in that file's own order.
Nothing is edited: every count, every annotation value and every blank annotation is what JASPAR
published. Each one is here for a rule it breaks, and the union file is the source because only it
holds every **Tax group** — five are represented, `diatoms` included, which is a group of exactly
one motif in both releases.

| Record | Name | Positions | Tax group | The rule it exists to break |
|---|---|---|---|---|
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
several values in one annotation, separated by `; `, and two carry a comma *inside* a single value.
Splitting on the comma corrupts about fifty records per release and fails nothing while doing it,
which is what these four are here to catch.

None of this is prose to be trusted either: `tests/test_jaspar.py` asserts every row of the table —
each id, name, length and tax group, each multi-valued and each blank annotation, the fractional
counts, the two below-minimum lengths and the long record's flank bits — against the committed
bytes.

## `chimera/` — the tiny component assemblies

Four component assemblies, cut from the files above rather than downloaded again, so every base and
every annotation line is still real `sacCer3`. Each carries something the naming contract has an
opinion about and no shipped assembly can demonstrate:

| Component | Chromosomes | Wrap | What it is here for |
|---|---|---|---|
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
|---|---|---|---|
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

None of this is prose to be trusted: `tests/test_chimera_fixtures.py` asserts every slice, length,
wrap width and masked stretch against the committed bytes, and `CHIMERA_COMPONENTS` in
`tests/conftest.py` is the same table as data, with a `chimera_component` fixture that registers one
as an assembly.
