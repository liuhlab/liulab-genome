# Research

Dated measurements: a number, and the method that produced it. A measurement is not a
decision — what one *decided* lives in `docs/adr/`. Read this index, then open only the
file you need.

A finding is true as of its date and is never edited to match later code. Where a
measurement has since been acted on, this index says so; the file itself does not.

| Measured | File | Finding |
|---|---|---|
| 2026-08-11 | [2bit-without-fatotwobit](2bit-without-fatotwobit.md) | UCSC publishes prebuilt `.2bit`, but nothing in Python can *write* one, so `faToTwoBit` cannot be retired while arbitrary user FASTA is supported. `pyfaidx` agrees with `py2bit` byte-for-byte including soft-masking, and is faster on large queries. |
| 2026-08-12 | [echt115-ce11-contents](echt115-ce11-contents.md) | `ecHT115` is an 87-scaffold contig-level SPAdes draft, not a chromosome plus plasmids — no design may assume a chimera has few large sequences. Both GTFs' seqnames are set-equal to their assembly's, and the two FASTAs differ in wrapping and masking, so a chimera is heterogeneously masked. |
| 2026-08-13 | [aligner-index-params-and-reference-names](aligner-index-params-and-reference-names.md) | STAR's `genomeSAindexNbases` reproduces STAR's own source exactly; `genomeChrBinNbits` was not passed, and the default of 18 padded `ecHT115` 5.18×. chromap has no size- or count-dependent index parameter. Both tools truncate a FASTA header at the first whitespace, so a suffix must ride on the first token. **Since acted on:** both parameters are now sized from the reference. |
| 2026-08-14 | [test-suite-parallelism](test-suite-parallelism-2026-08-14.md) | Eight pytest workers is the floor of the curve (2.8 s against 5.9 s serial); `auto` uncapped is *worse*, and `--dist=load` beats `--dist=loadfile`. Pinning the BLAS thread pools does nothing here. The four-step gate run concurrently is 4 s against 9.4 s sequential. **Since acted on:** both are now how the tasks are spelled. |
| 2026-08-28 | [motif-scan-engine](motif-scan-engine-2026-08-28.md) | `memelite`'s `fimo` silently returns zero hits for 97 of the 879 JASPAR 2024 vertebrate motifs, and never scans the last window of a sequence. MOODS is also 9.7× faster per core with thresholds prepared, and is the only one of the two that can be driven without holding the whole file. **Since acted on:** MOODS is the engine, in the core dependency table. |
