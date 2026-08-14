# 8. A chimera is an Assembly, identified by its component set and not by the order it was given in

A reference concatenated from canonical assemblies is not a new kind of thing: chimera-ness is a
third kind of source — *these components* — so ADR-0003's free-form local key and the existing
aligner paths carry it. Identity is the component *set*: sorted names join with `_`, a name derived
and never given, and the FASTA is its components' bytes verbatim, one token per header extended, so
one set means one directory and one index. That name resolves to its components offline, and four
ordered checks tell a chimeric name from an ordinary local key — the third a split on `_` into two
or more parts, each prepared here or listed in the table — but only a prepared set *builds* one, so
a cold machine raises, naming the missing component. Nothing was fetched, so nothing is pinned: its
row ships identifier columns and `sha256` blank — the table pins what was downloaded, a test what
was derived — and its merged annotation gets no row at all. A parallel `Chimera` type, the obvious
reading, loses: its own route to both aligners works against one alignment pass instead of N, the
first thing a chimera is for; argument order is discarded, and no lab nickname is possible.
