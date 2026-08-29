# 18. An Xref source must publish dated releases at stable URLs, and its bytes are downloaded

An **Xref source** is eligible only if its old releases stay retrievable at stable URLs, and an
**Xref set** is fetched into the **Data dir** on first construction rather than shipped in the
wheel. The Alliance of Genome Resources, Ensembl's per-species dumps, HGNC's quarterly archive and
WormBase qualify; NCBI Gene, UniProt idmapping and MGI are rebuilt daily and overwritten in place —
MGI's files carry no date stamp at all — so they are excluded while they stay rolling. Shipping the
tables instead, as every other data file here is shipped, lost on the contract a checksum makes: a
metadata row that travels in the wheel must still match a year later, which a rolling source breaks
within a day, and the smallest qualifying file is 26 MB besides. Two costs are accepted rather than
worked around. Unlike every TF answer this package gives today, neither set answers offline on a
fresh install, so the lab's internet-less compute nodes require the first construction on a login
node. And mouse gets current symbols only: previous and alias spellings are MGI's, MGI publishes no
dated archive, and no eligible publisher has them.
