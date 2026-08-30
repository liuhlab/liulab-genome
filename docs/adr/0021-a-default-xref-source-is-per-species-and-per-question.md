# 21. A Default xref source is per species and per question

A **Default xref source** is flagged twice in the curated table — once for identifiers, once for
symbols — and the question is named where an unnamed source is filled in: `for_symbols` on
`lookup_xref` and on the `XrefSet` constructor, and so on `genome xref symbols` and on `genome xref
ids --from-stems symbol`. One flag per species made the epic's most common question fail on the
first try, since human's default is the Alliance, which publishes no human symbol at all: a gene
list copied out of a paper needed `--source hgnc` before it worked once, and mouse and worm reach a
third source again. Resolving the source at the *question* instead — `match_symbols` picking a
symbol-carrying set for itself — lost on what it breaks: an **Xref set** holds one publisher's
bytes, and answering from another's is the merge that one query reading exactly one set forbids
(ADR-0017). The cost is a question-named constructor beside the plain one, `XrefSet("Homo sapiens")`
still matching no symbol and raising to name `for_symbols` — the surprise this record exists to
explain.
