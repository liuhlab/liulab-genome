# 23. A page the site publishes carries no record number

A record number is a citation between agents: it names a permanent thing and resolves for anyone
holding the repository, which is why docstrings and this repository's own prose cite that way. The
built site excludes `docs/adr/` on purpose, so the same number on a published page names a file
that reader cannot open — and one hundred and eleven were on the API reference page. The previous
answer linked each to GitHub so neither reader lost anything; it was wrong about the audience, not
the mechanism, because a record is written for whoever changes this package and its number is a
door the reader who came for `Genome.sequence` was never offered. So `mkdocstrings` drops the
citation as it renders and `src/` keeps it, `help()` included. Source listings went with it:
griffe cannot reach a citation inside a verbatim listing, and stripping one would publish source
disagreeing with the file it names. What it costs: that reader has no route to the reasoning
behind a behaviour, and no source beside the signature. What it buys: a record stays cheap.

**Status.** Amends the decision in `34e41ee`, which linked the numbers rather than removing them,
and which shipped in no release.
