# 7. A broken registration raises and names its repair, rather than being resumed

Nothing on disk told a finished preparation from one killed half-way: the old marker said only
that some run had reached the end, and its absence meant "fetch the whole thing again". The
completion record now answers *is this finished?* on its own, and the two ways a directory can
contradict it are errors — files present with no record, and a record that disagrees with what is
on disk — each naming the offending file and quoting `genome assembly register <assembly>
--force`, which is also the repair: it keeps an unpacked FASTA whose sha256 is still the pinned one, and fetches again
whenever it cannot prove that. The obvious reading, resume silently and rebuild quietly, loses
because the two states are indistinguishable, so guessing risks a genome that is silently wrong
rather than loudly absent. An absent or empty directory is not a broken state; the cost accepted is
that interrupting a *first* registration leaves one that raises next time.
