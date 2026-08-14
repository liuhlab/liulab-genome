# 9. A chimera's chromosome names are suffixed unconditionally, by a separator that announces itself

A chimera's chromosome is spelled `<chromosome>__<component>`, collision or not, so attribution is a
regex over the name and nothing is stored — and that spelling is baked into the FASTA, `.fai`,
`.2bit`, `chrom.sizes`, every merged GTF and both aligner indexes, so changing it invalidates a
built reference rather than a line of code. A single `_` cannot carry it: the real hg38 name
`chr1_KI270706v1_random` is indistinguishable from a suffixed one, so a name off a BAM header cannot
be classified without knowing its reference. The separator is therefore derived per chimera and
recorded — the shortest underscore run, minimum two, longer than any run a component's chromosome
names contain. Resolving bare names too, the obvious kindness, loses because it restores the
ambiguity the suffix abolishes and makes resolution depend on which components are present:
`genome["I:0-100"]` raises and names `I__ce11`. Code written against one component does not drop in.
