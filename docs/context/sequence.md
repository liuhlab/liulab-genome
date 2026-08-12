# Sequence

Bases as values. This context covers `seq.py`: `DNA`, `RNA` and `Protein` are `str` subclasses that
carry the biological transforms without ever losing their type or their case. A sequence is a
**leaf** — a `DNA` carries no **Assembly**, no **Region** and no **Strand**, and nothing here can
tell you where its bases came from. That is deliberate: coordinates and provenance belong to the
**Genome** that produced the sequence, and a second copy riding along on the letters would only give
the same question two answers.

Words every context shares — **Assembly**, **Genome**, **Region**, **Strand**, **Soft-masking** and
the rest — are defined once in the repo-root `CONTEXT-MAP.md`.

## Language

### The types

**Sequence**:
Bases or residues and nothing else: a `str` subclass whose type names what its characters mean. Case
is never normalised — lowercase is **Soft-masking**, which is information, so slicing, indexing and
every transform hand it back intact.
_Avoid_: string, seq; and never "read" — a read is a sequencing product, and this package has none

**DNA**:
A **Sequence** over `A`, `C`, `G`, `T` — the class every base-returning call in the package produces.
It is the only type that both complements and transcribes; ambiguity codes (`N` and friends) are
outside its **Alphabet** and so outside the package.
_Avoid_: nucleotide sequence, bases, `str` — a fetched sequence is a `DNA`, never a bare string

**RNA**:
A **Sequence** over `A`, `C`, `G`, `U`, reached by transcribing a **DNA** and returnable to one. The
`U`-for-`T` swap is the whole of what the type encodes.
_Avoid_: transcript (a transcript is an annotation feature with coordinates; an `RNA` has none),
mRNA, cDNA

**Protein**:
A **Sequence** over the twenty standard amino acids. A terminus by design: no complement, no
transcription, and no translation into it — the codon table is not this package's business.
_Avoid_: peptide, amino acid sequence, ORF, residues

**Alphabet**:
The character set a sequence class names — **documentation, not a runtime check**. Construction
validates nothing, deliberately: scanning every character costs too much on a whole chromosome, so
the type states what its contents mean rather than proving it. Reject non-alphabet input at the I/O
boundary, where the scan is paid for once.
_Avoid_: charset, valid characters, validation; IUPAC — the ambiguity codes are excluded on purpose,
not by oversight

### Transforms and measures

**Complement / Reverse complement**:
The base-for-base pairing of a **DNA** or **RNA**, and that pairing read back to front. Reverse
complement is an involution — apply it twice and you have the original, case included — and it is
the only sanctioned way to express a sequence as it reads on the other **Strand**.
_Avoid_: revcomp, rc, flip; and never plain "reverse", which is a string operation with no
biological meaning

**Transcribe / Back-transcribe**:
The `T`-to-`U` rewrite that turns a **DNA** into an **RNA**, and its exact inverse. A type change and
a letter swap: no splicing, no promoter, no biology beyond the characters.
_Avoid_: reverse-transcribe (that names a wet-lab reaction yielding cDNA), translate (translation
makes a **Protein** and is out of scope), convert

**GC content**:
The fraction of a sequence's characters that are `G` or `C`, case-insensitive, `0.0` for the empty
sequence. A measure of the bases in hand — a **Region** or a **Chromosome** has no GC content until
someone fetches it.
_Avoid_: GC%, GC ratio, GC count — it is a fraction in `[0, 1]`, not a percentage and not a tally
