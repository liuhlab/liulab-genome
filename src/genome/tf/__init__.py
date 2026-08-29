"""Transcription factors and the motifs they bind.

Two halves, kept apart because they are keyed differently. :mod:`genome.tf.gene`
is about the *gene* — which genes a published census judges transcription
factors, and the DNA-binding-domain family it classifies each one under — and so
it is keyed by gene. :mod:`genome.tf.motif` is about the *sequence a factor
recognises* — the count matrix and where it occurs in an assembly — and so it is
keyed by motif. The mapping between the two is many-to-many, which is why neither
half owns it. Cofactors are out of scope and nothing here answers for them.
"""
