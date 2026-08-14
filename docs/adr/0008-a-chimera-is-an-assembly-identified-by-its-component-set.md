# 8. A chimera is an Assembly, identified by its component set and not by the order it was given in

A reference concatenated from canonical assemblies is not a new kind of thing: chimera-ness is a
third kind of source — *these components*, alongside a pinned URL and a hand-supplied path — so the
existing aligner paths are reached with nothing new. ADR-0003 did the hard part: an assembly id is
already a free-form local key, and validation belongs to the source, not the identifier. Identity is
the component *set* — names sort lexicographically and join with `_`, and the name derives rather
than being given, so one set means exactly one directory and one index. A parallel `Chimera` type,
the obvious reading, loses because it would need its own route to both aligners — working against
one alignment pass instead of N, the first thing a chimera is for. The costs: argument order is
discarded, no lab nickname is possible, and `AssemblyMetadata` inherits a class of row whose
identifier columns are singular — what such a row holds, or whether one is needed, stays open.
