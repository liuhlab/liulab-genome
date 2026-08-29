# 16. The human cofactor list is this package's own union of two publishers

Human **Transcription cofactor** membership is decided here — the union of AnimalTFDB 4.0 and
EpiFactors v2.0, 1,466 genes: 354 both publishers list, 670 AnimalTFDB alone, 442 EpiFactors alone.
It is the first place this package decides anything rather than relaying a publisher's verdict, and
either source on its own was the alternative it beat: AnimalTFDB does not list the 442 and publishes
none of the function, target or complex vocabulary, while EpiFactors is human-only, so it answers
for neither mouse nor worm, and keys on HGNC ids rather than Ensembl, so it cannot be read against a
**Gene id stem** at all. That crosswalk is HGNC, pinned to one dated monthly archive and never the
rolling current file, so the 442 stems are reproducible. Two costs are accepted rather than smoothed
over: 151 human genes are both a **TF gene** and a cofactor — 57 from the AnimalTFDB side, 122 from
the EpiFactors side — so the two lists overlap and a caller who unions them double-counts; and the
two classification vocabularies are deliberately not crosswalked, for the reason ADR-0014 already
gives.
