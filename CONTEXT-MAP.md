# Context map

`liulab-genome` turns `(an assembly name) + (a lab data directory)` into reference files you can
query — a locus into bases, a GTF into a registered annotation, a FASTA into an aligner index. Its
vocabulary splits by bounded context: each file below is the glossary for one part of the source
tree, and the shared kernel at the bottom holds the words every context uses. The glossaries live
under `docs/context/`, not beside the code — only one context is a directory, so co-locating them
would put three of the four in arbitrary places. Rules live in `CLAUDE.md`; this map and the four
files it lists are a glossary and nothing else.

**Use these words.** When your output names a domain concept — an issue title, a refactor proposal,
a hypothesis, a test name — use the term as defined, not a synonym an entry lists under *Avoid*. A
concept defined nowhere is a signal either way: usually it is language the project does not use,
occasionally a real gap worth adding.

**Two vocabularies, and they do not mix.** Domain terms come from these files. Architecture terms —
module, interface, depth, seam, adapter, leverage, locality — are fixed, and "component", "service",
"API" and "boundary" are not substitutes for them. One narrowing, not a loophole: "component" is a
domain term in the [Assembly](./docs/context/assembly.md) context, so it is banned only as a
substitute for *module*.

## Contexts

- [Sequence](./docs/context/sequence.md) — covers `seq.py`: bases as a typed string, and the
  transforms that keep the type
- [Assembly](./docs/context/assembly.md) — covers `genome.py`, `metadata.py`, `external.py`,
  `io/{download,registration,fasta,twobit,utils}.py`: which reference this is, where its files live,
  and how a locus becomes bases
- [Annotation](./docs/context/annotation.md) — covers `io/gtf.py` and the GTF registry on `Genome`:
  what a GTF declares over one assembly, and the name it is addressed by
- [Index](./docs/context/index.md) — covers `aligner/*`: what one external mapper needs built before
  it can map, and how a finished build is told from an abandoned one

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

## Shared kernel

Ten words every context uses. Anything defined here is not redefined in a context file.

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
default.
_Avoid_: case, lowercase, formatting; bare "masking" — hard-masking writes `N` and is a different
thing

**Data dir**:
The root of all lab reference data, read from `$LIULAB_DATA` (`src/genome/io/registration.py`), under
which each **Assembly** owns exactly one directory. That per-assembly directory is the layout every
other context files into — annotations at `gtf/<name>/`, indexes at `index/<name>/`.
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
before use, and failing with the exact command that installs it. samtools/bedtools in `external.py`
and STAR/chromap in `aligner/` are one concept implemented twice.
_Avoid_: dependency, native dependency (that is the packaging view), subprocess, backend, wrapper,
bare "binary"
