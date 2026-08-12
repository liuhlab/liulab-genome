# Test fixtures

Small, subsampled **real** files — never a large genomic file. Every byte here came from UCSC's
`sacCer3` golden path and was cut down on the lab cluster; nothing is synthesised.

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
