# 5. The sequence alphabet is documentation; validation happens at the I/O boundary

`DNA`, `RNA` and `Protein` declare an alphabet and enforce nothing at construction: scanning every
character costs too much on a whole chromosome, and the references this package serves are full of
`N` runs a strict constructor would reject. The type states what its contents mean rather than
proving it, and the check moves to the edge — `genome revcomp` rejects a non-alphabet string and
exits 2, because a CLI argument is short and typed by a human. IUPAC ambiguity codes are excluded
from every declared alphabet on the same reasoning: naming them would imply a promise nothing keeps.
What it cost is that an in-memory `DNA` can hold anything, so a caller who needs the guarantee must
ask for it where the data enters.
