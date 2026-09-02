# 26. A chimera's sequence layout is a published contract, not an artefact of concatenation

A chimera's FASTA holds one contiguous block per **Component**, the blocks in the sorted order its
name spells, each component's own declared order kept inside its block. That arrangement falls out
of the obvious implementation, which is exactly why it needs recording: a consumer filtering one
component's sequences back out of an alignment header recovers a single-assembly header only while
it holds, so the layout is depended on off-repo and any change to it breaks callers silently rather
than loudly. Publishing it costs the freedom to reorder — no interleaving, no streaming components
as they finish, no parallel concatenation that lets blocks race — and that price is the point:
those are precisely the optimisations that would quietly invalidate a consumer's filter. ADR-0008
settles identity and the derived name and ADR-0009 the unconditional suffix; neither says a word
about order, and an unrecorded contract is indistinguishable from a guarantee nobody meant to make.
The alternative — treating layout as an implementation detail and telling consumers to parse the
suffix instead — was rejected because the suffix already answers *which component*, while only the
layout answers *contiguously, in a known order*, which is what makes the filter a slice rather than
a scan.
