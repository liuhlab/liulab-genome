# 4. A context's glossary is `docs/context/<name>.md`

Only one of the four bounded contexts maps to a directory (`aligner/`) — Sequence is one root file,
Annotation is one file inside `io/`, and Assembly straddles the root and `io/` at once — so
co-location, the obvious move, would buy one context and cost an unconditional rule. The glossaries
sit together under `docs/context/` instead, a tree excluded from the built site like `adr/` and
`research/`. What it costs: no glossary sits beside the code it describes.
