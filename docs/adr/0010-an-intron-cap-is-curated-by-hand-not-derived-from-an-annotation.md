# 10. The assembly table carries a hand-set intron cap, and nothing here derives one

The table cross-referenced an assembly across naming authorities and pinned where its bytes come
from; it now also carries the longest gap a spliced aligner should take for an intron on that
assembly, and why that number. Neither column is a cross-reference and nothing in this package reads
either — aligner parameters belong to whoever aligns — but the bound is a property of the reference,
so the alternative is each consumer re-typing a per-organism constant into a config of its own.
Deriving it from the registered annotation lost, and not for difficulty: an annotation catalogues
the transcripts someone observed, so its longest intron is a floor on what the organism does rather
than a ceiling on it, and a derivation that gets it wrong fails silently in the *tight* direction —
reads from a real intron stop aligning and nothing says a length rule did it. The values are
therefore deliberately loose round numbers, each with its rationale beside it to tell a convention
from a measurement, and an assembly nothing backs a number for keeps a blank cell, which reads as
*no bound has been chosen* and changes no alignment. The cost is that a wrong cap is unauditable
from here: the suite can prove the cell round-trips, never that the number is right.
