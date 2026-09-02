# Motif

What a transcription factor recognises, and where a scan says it occurs. This context covers
`tf/motif/*`: the **Motif** and **Motif set** types, the JASPAR loader, the scan-engine adapter, and
the mixin that gives a **Genome** one scan method. It is the motif half of the TF context — keyed by
motif, where the [gene-keyed halves](./tf.md), TF genes and cofactors, are keyed by gene. A motif
belongs to no **Assembly**: it is a pattern and not a place, which is why its files are the first
thing under the **Data dir** filed outside the assembly tree.

The package is a scaffold today and its modules land one at a time. The words are settled first, so
that every issue, test name, error message and docstring on the way uses one set.

Words every context shares — **Genome**, **Region**, **Strand**, **0-based half-open**,
**Soft-masking**, **Data dir** and the rest — are defined once in the repo-root `CONTEXT-MAP.md`,
**Prepared set** among them: what a downloaded **Release** on disk *is*, and how one is prepared,
is that entry's and is not restated here.

## Language

### What a motif is

**Motif**:
The pattern one transcription factor recognises, as this package holds it: a **Count matrix**, the
**Motif id** and **Motif name** it is addressed by, and the six annotations JASPAR publishes about
it — class, family, **Tax group**, UniProt accessions, PubMed ids and data type. Frozen, and a leaf
in the sense a **Sequence** is: it names no **Assembly**, no **Region** and no **Strand**, and it
holds no **Background**, so scanning one motif against two backgrounds does not make it two motifs.
It carries a trim offset as well, because trimming uninformative flanks returns a motif with the
same id whose column zero is column `offset` of the untrimmed one — which is what lets a **Motif
hit** found with the trimmed matrix be read back in the full motif's frame.
*Avoid*: consensus — a motif reports one and the method keeps the word, but the vocabulary names no
concept with it: the string is a one-letter-per-position rendering of the **Count matrix**, never
the motif itself and never what a **Motif hit** matched; pattern, TF, factor, "the matrix" alone

**Motif id**:
The key one **Motif** is addressed by — JASPAR's versioned matrix accession, `MA0139.2`. A **Motif
set** resolves the bare base id, `MA0139`, to the same motif, unambiguously, because a non-redundant
**Release** ships exactly one version of each. It is the identifier a lookup can never be ambiguous
about, which is why it is what an error names when a **Motif name** is. Unambiguous is not the same
as stable, though: a base id addresses whichever version the **Release** ships, and that version is
not guaranteed to describe the same factor — or even the same organism — as another version of the
same base id cited elsewhere. Both releases ship `MA0095.4`, mouse `Yy1`, so `db["MA0095"]` answers
with the mouse factor to someone who came looking for the human `MA0095.1`.
*Avoid*: accession, matrix name, key; and **Motif name** as a synonym — one addresses a motif, the
other labels it

**Motif name**:
The factor name JASPAR publishes for a **Motif** — `CTCF`, or `Ptf1a::Rbpj` for a dimer. A label and
not a key: 66 names collide in the 2024 **Release** and 71 in 2026, so a lookup by name always hands
back a tuple, of one where the name is unique, and indexing a **Motif set** by an ambiguous name
raises an error naming every matching **Motif id** rather than silently returning one of four CTCFs.
*Avoid*: gene name, gene symbol, symbol — a dimer's name names no gene and joins to no
**Annotation** on its own; id, accession, key

**Count matrix**:
The 4 × L table of observed base counts behind a **Motif** — one row per base, one column per
position — and the single source of truth in it. Its values are floats and not integers, because
JASPAR's own records carry fractional counts. Probabilities come from normalising a column; log-odds
come from a **Background** and a pseudocount, which are arguments and never fields, so no derived
form is ever stored.
*Avoid*: PWM, PSSM, position weight matrix, scoring matrix, log-odds matrix — this package stores
counts and derives every weighted form on demand, so a name promising weights describes something
that is not on the object; and bare "the matrix" in an API surface, which does not say which of the
three forms you are holding

**Information content**:
How much one column of a **Count matrix** says, in bits: 0 for a column that says nothing, 2 for one
fixed base. It is the y-axis of the sequence logo and the quantity trimming thresholds on, so the
height you see is the number you set.
*Avoid*: entropy — the inverse of this, and one idea gets one word; conservation — the alignment
field's word for agreement between species, a different measurement on different data; signal,
information, bits

### Where motifs come from

**Motif set**:
The container every scan, filter and comparison is a method on: **Motif**s held as one addressable
group, indexed by **Motif id**, by bare base id, or by a **Motif name** that is unique. A
**Release** read from disk is a motif set that also knows its release, its **Tax group** and where
its bytes came from; so is a set of de novo matrices out of a model, which is how everything a
release can do applies to matrices JASPAR never published. Filtering — by annotation or by an
arbitrary predicate — hands back a plain motif set rather than a release, because a filtered release
is no longer that release.
*Avoid*: database, except of a whole **Release** on disk, which is one; collection, which is
JASPAR's word for CORE against the redundant sets this package does not read; library, catalog,
motif list

**Release**:
One dated, versioned publication of the source a downloaded set is pinned to — here JASPAR's
non-redundant CORE data, 2024 and 2026 being the two this package prepares, and the same word for
an **Xref source**'s and for Ensembl Compara's. A release and a **Tax group** together name exactly
one JASPAR file, and one such file is one **Prepared set**, filed under the `motif/` subtree of the
**Data dir** and shared by every project on the machine; both are recorded on every **Hit table**,
so a table saved months ago still says what it was scanned with.
*Avoid*: version — the package has one and a **Motif id** has one, and neither is this; JASPAR
version, database version, edition, snapshot

**Tax group**:
JASPAR's taxonomic split of a **Release**, and the other half of what names a file: `vertebrates` —
the default, the lab's common case being human and mouse — plus `plants`, `insects`, `nematodes`,
`fungi`, `urochordates`, `diatoms`, and `all`, which is exactly the union of the other seven. It
chooses which file is downloaded rather than filtering one after it is loaded, so a worm scan never
pays for a thousand plant and fungal matrices. Each **Motif** carries its own as one of its six
annotations.
*Avoid*: species — a group holds many, and `all` holds every one JASPAR curates; taxon, clade,
lineage, taxonomy; and the bare word "group", which does not say on which axis

### Scanning

**Motif hit**:
One place a **Motif** cleared the score its **Threshold** converts to: a sequence name, a
`[start, end)` interval, a **Strand**, and the log-odds score in bits. Coordinates are **0-based
half-open** and always in the forward frame, whichever strand matched, and the strand is always `+`
or `-` — never `.`, because a scan knows which of the two it scored.
*Avoid*: TFBS, binding site, site — a hit is a matrix clearing a cutoff on some bases, and whether
the factor binds there is a question this package does not ask; match, call, peak

**Hit table**:
The one answer every scan returns, whatever it was handed — a `DNA`, named sequences, a FASTA, or
**Region**s resolved against a **Genome**: one row per **Motif hit**, with a fixed column set and
fixed compact dtypes that are part of the contract rather than an optimisation, so nothing
downstream branches on how the scan was invoked. Its provenance travels on the frame — the
**Background** used, the **Threshold**, the **Release**, the **Tax group**, which motifs were
scanned and which were skipped, and the **Assembly** when the scan was of **Region**s — because two
tables missing any of those cannot be reconciled. A
scan too large to hold streams to Parquet and hands back the path instead.
*Avoid*: results, output, scan result, hits file, BED (the coordinates outlive the format); and
never a bare frame passed on with its metadata dropped, which is the same table with its meaning
removed

**Background**:
The base composition log-odds are scored against — four frequencies, zero-order — and an argument to
a scan rather than anything a **Motif** owns. Automatic by default: derived from the input where it
holds at least 10,000 unambiguous bases, uniform below that, since a composition estimated from
fewer would distort the very cutoffs it sets. It decides the answer more than any other parameter —
moving from uniform to one real chromosome's composition changed the hit count by 2.5% and turned
over 26% of the hits — which is why it is recorded on the **Hit table** rather than left to be
assumed.
*Avoid*: null model, prior, GC correction; base composition, which names the measurement rather than
the role it plays here; and background set or background regions, the enrichment sense of the word —
enrichment is out of scope

**Threshold**:
The stringency a scan is called at, given as one per-position p-value — 1e-4 by default — and
converted per **Motif**, against the **Background**, into the score that motif must clear. A p-value
means the same stringency for a short matrix and a long one, which a score does not. A motif that
cannot reach it is left out and named among the skipped rather than called at something looser, so
every row of a **Hit table** was called at what was asked for; the same arithmetic is why the
shortest scannable motif is 7 positions, the best possible word of a 6-mer being p = 2.44e-4.
*Avoid*: cutoff — that is the per-motif score this converts into, not this; score threshold,
significance, FDR, q-value, since nothing here corrects for multiple testing; and "p-value" for the
**Hit table**'s score column, which is log-odds in bits
