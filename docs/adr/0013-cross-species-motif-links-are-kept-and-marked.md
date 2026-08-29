# 13. A motif measured on another vertebrate still links, and every such link is marked

A **Motif link** joins a JASPAR profile to a **TF gene** of one species, and a profile whose matrix
was measured on a different vertebrate is kept and flagged rather than dropped. JASPAR's CORE
`vertebrates` collection files an orthologous pair's matrix under whichever species was assayed, so
a profile's species field is an artefact of the experiment and not a claim about the factor: the
mouse `Yy1` matrix is what JASPAR has to say about YY1 in any vertebrate. Restricting links to
species-matched profiles — the tidy reading, and the one a species column invites — lost on what it
returns: it costs human 108 genes and mouse 552 of its 689, leaving mouse at roughly a fifth of the
coverage JASPAR can support, and it discards real measurements on the strength of a field that does
not mean what it looks like. Every such link therefore carries a **Cross-species link** flag, which
is what a caller whose question demands a species-matched profile filters on — 732 of mouse's 896
links on the 2026 release are cross-species. The cost: nothing here asserts orthology, so such a
link says only that JASPAR filed one vertebrate's assay under one vertebrate's name.
