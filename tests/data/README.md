# Test fixtures

Small, subsampled **real** files — never a large genomic file. Every byte here came from UCSC's
`sacCer3` golden path and was cut down on the lab cluster; no bases are synthesised. Names are
this repo's where a fixture needs a spelling `sacCer3` does not have, and one 200-base stretch is
lower-cased — both are called out below, and nothing else departs from the source bytes.

| File | What it is |
|---|---|
| `tiny.fa` | The first 10 000 bases of `chrI`, `chrII` and `chrIII` of `sacCer3`, cut with `samtools faidx` and re-headed to the bare chromosome names |
| `tiny.fa.gz` | `tiny.fa`, gzipped — the compressed-source path |
| `tiny.gtf` | The 85 `sacCer3.ensGene` features that fall wholly inside those three windows |
| `tiny.gtf.gz` | `tiny.gtf`, gzipped |
| `ensembl_style.gtf` | `tiny.gtf` with the `chr` prefix stripped (`I`, `II`, `III`) — an Ensembl-spelled annotation for the chromosome-name mismatch case |

Sources:

- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz`
- `https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/genes/sacCer3.ensGene.gtf.gz`

`tiny.gtf` coordinates are 1-based inclusive, as every GTF is; they convert at the I/O boundary and
are never seen in that form inside the package.

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
