# Context map

`liulab-genome` turns `(an assembly name) + (a lab data directory)` into reference files you can
query — a locus into bases, a GTF into a registered annotation, a FASTA into an aligner index. Its
vocabulary splits by bounded context: each file below is the glossary for one part of the source
tree, and the shared kernel at the bottom holds the words every context uses. The glossaries live
under `docs/context/`, not beside the code — only four of the eight map cleanly onto a directory,
and TF straddles `tf/gene/`, `tf/cofactor/` and a module beside them, so co-locating them would put
four in arbitrary places. Rules live in `CLAUDE.md`; this map and the eight files it lists are a
glossary and nothing else.

**Use these words.** When your output names a domain concept — an issue title, a refactor proposal,
a hypothesis, a test name — use the term as defined, not a synonym an entry lists under *Avoid*. A
concept defined nowhere is a signal either way: usually it is language the project does not use,
occasionally a real gap worth adding.

**A term marked _(decided, not built — ADR-N)_ names something the record settled and the
code does not have yet.** Use the word; do not call the API it describes.

**Two vocabularies, and they do not mix.** Domain terms come from these files. Architecture terms —
module, interface, depth, seam, adapter, leverage, locality — are fixed, and "component", "service",
"API" and "boundary" are not substitutes for them. One narrowing, not a loophole: "component" is a
domain term in the [Assembly](./docs/context/assembly.md) context, so it is banned only as a
substitute for *module*.

## Contexts

- [Sequence](./docs/context/sequence.md) — covers `seq.py`: bases as a typed string, and the
  transforms that keep the type
- [Assembly](./docs/context/assembly.md) — covers `genome.py`, `metadata.py`, `external.py`,
  `chimera.py`, `io/{source,components,fetch,chimera,download,registration,fasta,twobit,utils}.py`:
  which reference this is, where its files live, and how a locus becomes bases
- [Annotation](./docs/context/annotation.md) — covers `io/gtf.py`, `gene_list.py` and the GTF
  registry on `Genome`: what a GTF declares over one assembly, the name it is addressed by, and
  which of its genes are in a category
- [Index](./docs/context/index.md) — covers `aligner/*`: what one external mapper needs built before
  it can map, and how a finished build is told from an abandoned one
- [Motif](./docs/context/motif.md) — covers `tf/motif/*`: what a transcription factor recognises,
  stored as counts and belonging to no assembly, and where a scan says it occurs
- [TF](./docs/context/tf.md) — covers `tf/gene/*`, `tf/cofactor/*` and `tf/link.py`: which genes a
  published census judges transcription factors and which a publisher lists as cofactors of
  transcription, in a registered annotation's own gene ids, and which motifs answer for them
- [Xref](./docs/context/xref.md) — covers `xref/*`: which foreign identifiers name a gene and which
  genes a foreign identifier names, on one named publisher's assertions at one pinned release, with
  no genome open
- [Orthology](./docs/context/orthology.md) — covers `homology/*`: which genes in another species a
  gene is homologous to, and how many-to-many the publisher's own gene tree says that is

`io/results.py` sits in Assembly and Annotation both: what a registration answers with, for either.
It is the return types the API hands back and the CLI renders, so it carries the vocabulary of
whichever context asked.

`cli.py` is covered by no context — the map covers what carries domain vocabulary, not the whole tree.

## Relationships

- **Assembly → Sequence**: `Genome.fetch_sequence(region)` is the only place a **Region** becomes
  bases. It returns a `DNA`, never a `str`, and preserves **soft-masking**. The edge runs one way
  only — Sequence is a leaf, and a `DNA` carries no assembly, no region and no strand.
- **Assembly → Annotation**: an annotation belongs to exactly one assembly, filed under the
  assembly's **Data dir** at `<assembly dir>/gtf/<name>/` and addressed everywhere by its registered
  name. Registering a GTF whose chromosome names disagree with the assembly's is an error, not a
  warning.
- **Assembly + Annotation → Index**: STAR indexes an assembly's FASTA *against* one registered
  annotation, one index per annotation (`index/star_<gtf>/`); chromap indexes the FASTA alone
  (`index/chromap/`). Whether an annotation is an input is a property of the aligner, not of the index.
- **Index → Assembly**: an index lives *inside* the assembly it indexes, at
  `<assembly dir>/index/<name>/`, so the assembly's layout owns where indexes go. Both contexts
  assert a finished build with a **Completion marker**.
- **Genome → Motif**: scanning **Region**s is the one place a motif's region-local positions become
  **Chromosome** coordinates — the interval and the hit **Strand** both flipped for a `-` strand
  region, so nobody does that arithmetic twice. The edge runs one way only: the scan method is mixed
  into `Genome`, and the motif modules import `Genome` under type checking alone, since a motif
  belongs to no **Assembly** and a **Motif set** is usable with no genome open.
- **Motif → Sequence**: a scan consumes bases and nothing else — a `DNA`, a mapping of named
  sequences, or a FASTA on disk — and Sequence stays the leaf it is, so nothing crosses back.
  Scanning upper-cases what it is handed: **Soft-masking** does not survive the crossing, and there
  is no option to honour it.
- **TF gene → Annotation**: a census is keyed by **Gene id stem** and answers in a registered
  **Annotation**'s own gene ids — the registry resolves the stems, hands back every gene id a stem
  names rather than picking one, and carries the stems that resolved to nothing on the answer.
  Resolution is a general mechanism on the registry and the first thing to read the **Annotation
  database** this package has always built; TF genes are its first caller. The species is read from
  the **Assembly**'s own metadata and never passed in, and an assembly whose species has no census
  raises and names the ones that do.
- **TF gene ↔ Motif**: the join is many-to-many — one gene has several matrices, one dimer matrix
  answers for several genes — so neither half owns it. It lives in `tf/link.py`, which imports both
  halves and is imported by neither, and it ships as a plain table per species per **Release** rather
  than being computed from a **Motif name** (ADR-0015). Nothing crosses into the motif half: no field
  on **Motif**, no change to the JASPAR loader or the scan path.
- **Xref → Annotation**: an **Xref set** answers in **Gene id stem**s and stops there. Putting those
  into a registered **Annotation**'s own gene ids is `AnnotationRegistry.resolve_gene_ids`, used
  unchanged, so the hop is one existing call the caller makes and no xref module imports the
  registry. The edge runs one way: a set is anchored to a species and a **Release**, never to an
  **Assembly** or an annotation, so an id conversion never needs a **Genome** open.
- **Xref → TF gene**: a census and a **Cofactor table** are keyed by **Gene id stem** too, so an
  **Xref set** is how an Entrez, HGNC or UniProt column reaches them — the gap that let two UniProt
  entry names ship in the human census unnoticed. One way only: nothing shipped here is keyed by a
  foreign **Namespace**, and the xref half reads no census.
- **Xref ↔ Motif**: shape shared, and nothing else. Both a **Motif set** and an **Xref set** belong
  to no **Assembly**, live outside the assembly tree under the **Data dir**, are pinned to a
  **Release** and are prepared by construction — the same object, learned once. Neither imports the
  other, and a **Motif name** is never fed to a **Symbol match**: it labels a matrix and names no
  gene.
- **Orthology → Annotation**: the same one existing call the Xref edge uses, `resolve_gene_ids`,
  against the same registry, unchanged. A **Homology type** crosses it untouched and a **Dropped
  partner** count says what the crossing removed, so a link that only looks one-to-one in your
  annotation is distinguishable from one the publisher called one-to-one (ADR-0020).
- **Orthology → TF gene**: a prohibition rather than a call. No **TF gene table**, **Cofactor
  table** or list of this package's own is ever derived through homology and no answer is silently
  species-mapped (ADR-0019), so neither half imports the other and a species with no census raises
  and names the ones that have one.
- **Orthology ↔ Motif**: a **Cross-species link** is not a **Homology link**. JASPAR files an
  orthologous pair's matrix under whichever species was assayed and asserts no orthology by doing so
  (ADR-0013), and nothing joins the flag to a **Homology set** — a marked motif link stays a fact
  about where an assay was run.
- **Xref ↔ Orthology**: two sets that share a key and read each other never. Both are keyed by
  **Gene id stem** and both are one publisher at one **Release**, but a homology answer is never
  translated through an **Xref set** and an identifier answer is never carried across species — the
  same refusal to compose two hops into a third that ADR-0017 states within the xref half.

## Shared kernel

Eleven words every context uses. Anything defined here is not redefined in a context file.

**Assembly**:
The identity of one reference — a free-form **local** key that names its directory under the **Data
dir** at `genome/<assembly>/` and is the lookup key into the curated metadata table. UCSC is the
default *source*, not the namespace: validation is a property of the source, so `validate_assembly`
runs only on the UCSC fetch path and never on an assembly seeded from your own FASTA.
_Avoid_: build, genome build, reference, species; calling it a UCSC id — the name is local even when
the bytes came from UCSC; and **Genome**, which is this identifier already opened

**Genome**:
One **Assembly**'s materialized files and open handles, bundled behind the object you query. An
assembly is an identifier; a genome is that assembly opened — which is why a genome closes and an
assembly does not.
_Avoid_: reference, genome object, organism; **Assembly** as a synonym

**Chromosome**:
One named sequence of an **Assembly**, sized and ordered by that assembly's `chrom.sizes`. The
assembly's own spelling is authoritative — a name arriving from anywhere else (`chr1` against `1`) is
reconciled against it at ingest, never assumed to match.
_Avoid_: contig, scaffold, seqname, sequence (that is a Sequence-context word); spelled `chrom` in
code, but write the whole word in prose

**Region**:
A `[start, end)` interval on one **Chromosome** carrying an explicit **Strand** — the package's one
interval type, frozen and self-validating (`src/genome/region.py`). It names no assembly, so it means
nothing until resolved against a **Genome**.
_Avoid_: interval, locus, range, window, BED record; "position" for anything wider than one base

**Strand**:
`+`, `-`, or `.` — forward, reverse, unknown. `.` is a real answer meaning nobody knows, and is never
silently promoted to `+`.
_Avoid_: orientation, direction, sense/antisense; `None`, `""`, or an implicit `+` standing in for
unknown

**0-based half-open**:
The only coordinate convention that exists inside this package: `[start, end)`, so `end - start` is
the length and `chr1:0-10` is the first ten bases. 1-based-inclusive coordinates (VCF, GFF/GTF, SAM)
are converted at the I/O boundary and never travel inward.
_Avoid_: BED coordinates (the convention outlives the format), zero-indexed, exclusive end — name
both halves of the convention or neither

**Soft-masking**:
Lower-case bases marking repeat-masked regions. It is data, not formatting: it survives fetching,
slicing and reverse-complement, and is discarded only by asking, since `TwoBit(masked=True)` is the
default. Scanning is the one exception and discards it unasked — ADR-0012.
_Avoid_: case, lowercase, formatting; bare "masking" — hard-masking writes `N` and is a different
thing

**Gene id stem**:
A gene id with its version suffix dropped: `ENSG00000123456` for `ENSG00000123456.7`. A stem may
name more than one gene id inside one **Annotation** — nine do in `gencode_v50lift37`, eight of them
pseudoautosomal-Y — so resolution answers with every id a stem names and never picks one. A stem
carrying no version resolves to itself, which is why an annotation whose ids are not Ensembl-shaped
is unaffected. The mechanism belongs to the annotation registry rather than to any one context's
code, and reading gene ids is the first thing that opens the **Annotation database** this package
has always built; TF genes are merely its first caller.
_Avoid_: base id — that already names a **Motif id** with its version dropped, `MA0139` for
`MA0139.2`, and one term must not name two things inside one package; unversioned id, stripped id,
gene id (the versioned thing this is a stem *of*), ENSG (one publisher's spelling)

**Data dir**:
The root of all lab reference data, read from `$LIULAB_DATA` (`src/genome/io/registration.py`). Most
of it is the assembly tree at `genome/`, under which each **Assembly** owns exactly one directory,
and that per-assembly directory is the layout most other contexts file into — annotations at
`gtf/<name>/`, indexes at `index/<name>/`. Not all of it: data belonging to no assembly is a sibling
of that tree rather than a tenant of it, and the **Motif** files under `motif/` are the first.
_Avoid_: cache, cache dir (a cache may be evicted; this may not — though it is spelled `cache_dir` in
code), data root, download dir, workdir

**Completion marker**:
The record written only after a multi-step build finished, so its absence means *unfinished*, never
*missing*. One spelling — `.completion.json`, in the directory the build filled — carrying the
provenance as well as the verdict: where the bytes came from, their checksum, every file claimed with
its size, the **External tool** versions, the package version, when it finished, and whatever else its
own kind must be able to explain — for an **Index**, the exact command it ran, the parameters and the
FASTA consumed. Confirming one compares presence and size and reads no contents, so it is the cheap
answer to *is this finished* and the only answer to *how was this made*.
_Avoid_: flag, success flag, sentinel, stamp, lock file; and never an output file's mere existence,
which is what this word exists to distrust

**External tool**:
A binary the package shells out to instead of reimplementing — resolved on `PATH`, version-detected
before use, and failing with the exact command that installs it. One module, `external.py`, serves
every one of them: samtools and the two 2bit tools an assembly is prepared with, and the STAR and
chromap an **Index** is built with. Whether a tool's output is captured or streamed is an argument,
not a second implementation.
_Avoid_: dependency, native dependency (that is the packaging view), subprocess, backend, wrapper,
bare "binary"
