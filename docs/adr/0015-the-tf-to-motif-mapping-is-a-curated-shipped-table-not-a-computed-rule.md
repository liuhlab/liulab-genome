# 15. The TF-to-motif mapping ships as a curated table, not as a rule computed at runtime

The **Motif link**s ship as data — one table per species per **Release**, CORE vertebrates only,
generated offline by a committed script and readable without importing this package, since the
mapping is wanted outside here as often as inside. Computing it when it is asked for — upper-case a
**Motif name**, split it, match it to a census symbol — lost because it is roughly right and quietly
wrong exactly where it matters: it drops the three human genes JASPAR renamed after the census was
published, merges `FOS::JUN` into a JUN motif and so reads a heterodimer matrix as a monomer's, and
either discards every profile measured on another vertebrate or keeps it with nothing on the row
saying so. The cost is unmitigated, and is the one ADR-0011 already accepts for the curated gene
lists: no test proves a shipped table still matches JASPAR and none can, since regenerating one
needs a download and CI has no network. Pinned counts turn silent drift into a loud failure, which
is the most that is available.

**Status.** Amended 2026-08-29: a clause calling the tables plain rather than gzipped was struck,
having recorded a storage choice nobody decided. Shipped data follows the repository's convention —
bulk tables gzipped, small metadata tables plain — and no record decides that.
