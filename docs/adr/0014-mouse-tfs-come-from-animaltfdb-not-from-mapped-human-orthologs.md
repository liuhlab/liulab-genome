# 14. Mouse TFs come from AnimalTFDB's own census, not from human TFs mapped through orthologs

Each species is answered by a census published for that species — Lambert et al. 2018 for human,
AnimalTFDB 4.0 for mouse — in that publisher's own vocabulary. Serving mouse by mapping Lambert's
genes through an ortholog table, whose appeal was one harmonised **DBD family** vocabulary and one
file to maintain, lost on coverage and on which genes it drops: it reaches at best 88.6% and loses
236 C2H2 zinc fingers, the family that has diverged most between the two species, so what goes
missing is exactly what a mouse-specific question is about. The price paid instead is two
un-harmonised family vocabularies, 75 values under Lambert's `DBD` and 72 under AnimalTFDB's
`Family`, deliberately not crosswalked: `ARID/BRIGHT` and `ARID` are not asserted equivalent,
because inventing an equivalence nobody has checked is worse than shipping two vocabularies that
each say who spelled them. Two further costs are accepted: no census here is derived through
orthologs, whatever homology the package now serves, and mouse has no assessed-negative genes, since
AnimalTFDB lists only the genes it judges TFs — visible asymmetry rather than a fabricated one.

**Status.** Amended 2026-08-29: the cost line read "this package builds no ortholog or homology
support at all", which ADR-0019 ended by adding a **Homology set** a caller can query. Orthology is
served and still never consumed, so the decision above is untouched — mouse TFs come from
AnimalTFDB's own census, and no table this package publishes is derived across species.
