"""Transcription factors, cofactors, and the motifs they bind.

Two halves, kept apart because they are keyed differently. :mod:`genome.tf.gene`
is about the *protein* — which genes are TFs or cofactors, what family a factor
belongs to — and so it is keyed by gene. :mod:`genome.tf.motif` is about the
*sequence a factor recognises* — the count matrix and where it occurs
in an assembly — and so it is keyed by motif. The mapping between the two is
many-to-many, which is why neither half owns it.
"""
