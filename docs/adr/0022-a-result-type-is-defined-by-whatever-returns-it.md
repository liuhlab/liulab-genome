# 22. A result type is defined by whatever returns it, and `as_json` respells nothing

Every value the API answers with is a frozen dataclass living in the module that returns it, and
its `as_json` writes the record's own keys in the record's own order, so no surface renames what a
record already calls something. Nineteen of them lived together in `io/results.py` on the grounds
that they change for one reason — what an answer must be able to say — but the only thing they
shared was being returned, and a third of the file belonged to contexts no module under `io/`
imports; it had been opened as a seam because both halves of `io` needed two of them and neither
could import the other, and that reason expired once each type sat beside its producer. The one
type two contexts genuinely share, `ResolvedGeneIds`, is imported from beside `resolve_gene_ids`,
because the rule is where a type is returned and not where it is read. What it costs: nothing
mechanises the convention now that no module gathers it — pyright checks the annotation, not the
key order — so a new result earns its own key-order test, and the rest is held at review.
