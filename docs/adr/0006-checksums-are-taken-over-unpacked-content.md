# 6. A checksum is taken over unpacked content, never over the archive it arrived in

Nothing downloaded had ever been verified, and the obvious fix was free: pooch already takes a
`known_hash` and checks the bytes it just wrote. Those bytes are the `.fa.gz`, and gzip bytes are
not the genome — sacCer3's archive digests to `e3a70396…` while the FASTA it unpacks to digests to
`6ff72f07…`, so recompressing at another level, or mirroring through a host that does, breaks a
match that ought to hold. The metadata table therefore pins the sha256 of the **unpacked** file,
computed by streaming it after decompression; pooch's `known_hash` stays available for the archive
and is simply not what the table records. Two consequences are accepted in exchange: verification
happens after unpacking rather than before, so a corrupt download is caught only once it has been
written out in full; and the digest covers line wrapping and headers as well as bases, so it asserts
that two files are identical, not that two references carry the same sequence.
