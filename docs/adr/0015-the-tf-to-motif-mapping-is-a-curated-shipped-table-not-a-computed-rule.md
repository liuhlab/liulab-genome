# 15. The TF-to-motif mapping ships as a curated table, not as a rule computed at runtime

The **Motif link**s ship as plain TSVs, one per species per **Release**, CORE vertebrates only,
about 250 KB in all, generated offline by a committed script. Computing the mapping when it is asked
for — upper-case a **Motif name**, split it, match it to a census symbol — lost because it is
roughly right and quietly wrong exactly where it matters: it drops the three human genes JASPAR
renamed after the census was published, it merges `FOS::JUN` into a JUN motif and so reads a
heterodimer matrix as a monomer's, and it either discards every profile measured on another
vertebrate or keeps them with nothing on the row saying so. The tables are plain rather than gzipped
because a curated artifact's value is its reviewable diff, and readable without importing this
package because the mapping is wanted outside it as often as inside. The cost is unmitigated, and is
the one ADR-0011 already accepts for the curated gene lists: no test proves a shipped table still
matches JASPAR and none can, since regenerating one needs a download and CI has no network. Pinned
counts turn silent drift into a loud failure, which is the most that is available.
