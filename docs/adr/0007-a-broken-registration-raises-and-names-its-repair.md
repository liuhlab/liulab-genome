# 7. A broken registration raises and names its repair, rather than being resumed

Nothing on disk told a finished preparation from one killed half-way: the old marker said only
that some run had reached the end, and its absence meant "fetch the whole thing again". The
completion record now answers *is this finished?* on its own, and the two ways a directory can
contradict it are errors — files present with no record, and a record that disagrees with what
is on disk — each naming the offending file and quoting `genome register <assembly> --force`.
That command is also the repair: it keeps an unpacked FASTA whose sha256 is still the pinned
one and rebuilds only the derived files, and fetches the source again whenever it cannot prove
that, including when the row pins no digest to prove it against. An absent or empty directory
is not a broken state. The obvious reading — resume silently, rebuild quietly — loses because
the two states are indistinguishable, so guessing risks a genome that is silently wrong rather
than loudly absent. The cost accepted: a first registration briefly holds files with no record,
so interrupting one leaves a state that raises next time and needs a forced re-registration.
