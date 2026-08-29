# 19. Orthology is served to callers and never consumed by what this package publishes

A **Homology set** answers a user's cross-species question and is used for nothing else here. No
**TF gene table**, no **Cofactor table**, no list of this package's own is derived through homology,
and no answer is silently species-mapped: a species with no census or no shipped table raises and
names the ones that have them. Filling those gaps by translation — a worm TF census from human, a
mouse assessed-negative set from Lambert, a cofactor list carried across species — is the obvious
use of a homology table once one is on disk, and it lost on what the mapping actually reaches.
ADR-0014 measured the easy end and refused it there: human to mouse tops out at 88.6% and loses 236
C2H2 zinc fingers. Worm is not the easy end — strict one-to-one falls from 83.0% of human genes for
mouse to 14.1% for worm, and 76% of human↔worm links are `ortholog_many2many`. The cost is that the
gaps stay open with the table sitting right there: worm still has no TF census, mouse still has no
assessed-negative genes, and all that changes is that a user can now cross the species line
themselves, deliberately, with the publisher's **Homology type** in hand.
