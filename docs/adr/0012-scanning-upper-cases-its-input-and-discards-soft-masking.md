# 12. Scanning upper-cases its input and discards soft-masking, with no option to honour it

**Soft-masking** is data everywhere else here — it survives fetching, slicing and reverse
complement, and the kernel says it is discarded only by asking — and a scan discards it without
being asked: every sequence is upper-cased on its way to the engine, in one place, and no argument
makes a lower-case base score differently or be skipped. Honouring the mask lost on what it would
return: about half a mammalian genome is repeat-masked, so a mask-aware scan of a fetched **Region**
answers with nothing across most of the genome, silently and for a reason no reader of the **Hit
table** would guess. Offering it as a keyword lost separately — a repeat mask is a property of one
**Assembly**'s annotation pipeline, while a **Motif set** belongs to no assembly and is handed a
`DNA`, a mapping or a FASTA that may carry no mask at all, so one keyword would mean three things.
The upper-casing is the package's own and not the engine's, which happens to fold case as well: the
promise must hold whatever the engine is asked in. What it costs is unmitigated — repeat-driven hits
are called like any other, and a caller who wants them out subtracts a repeat annotation this
package does not hold.
