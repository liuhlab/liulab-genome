# Index

What an aligner must be given before it can map anything, built once and left on disk. This context
covers `aligner/`: STAR and chromap are driven far enough to produce an index and no further — the
package builds indexes and never aligns a single read. One word does double duty in this field, and
this context claims it outright: *index* here is a built aligner index, never a `.fai` or a `.2bit`.

Words every context shares — **Assembly**, **Data dir**, **Completion marker**, **External tool** and
the rest — are defined once in the repo-root `CONTEXT-MAP.md`.

## Language

**Index**:
An **Aligner**'s precomputed search structure over one **Assembly**, built from that assembly's FASTA
and read-only ever after. It exists only once its **Completion marker** is present; a directory of
index files without one is an interrupted build, so rebuilding it is a deliberate act — `overwrite`
— and never something that happens quietly.
_Avoid_: never for `.fai` or `.2bit` — those are the Assembly context's own derived artifacts and are
named by extension, never "index" and never "FASTA index". Also: genomeDir, reference, database

**Aligner**:
The **External tool** that builds an **Index** — today STAR and chromap — together with the class
that knows its binary, its version, its flags and its layout. Naming the tool is not a promise to run
it for mapping: this package builds the index and hands over the **Artifact**, and the alignment
command is the caller's to issue.
_Avoid_: mapper, read mapper; and never let "aligner" imply this package aligns — mapping is out of
scope by decision, not merely unimplemented

**Index dir**:
The one place an index may live: `<assembly dir>/index/<name>/` inside the assembly's **Data dir**,
so the assembly owns the layout and an index can never drift from what it indexes. `<name>` identifies
the index rather than just the **Aligner** — STAR builds one per annotation (`index/star_<gtf>/`),
chromap one for the whole assembly (`index/chromap/`).
_Avoid_: output dir, genomeDir, cache; and never a caller-supplied path — the location is derived from
the assembly, not configured

**Artifact**:
The single path the **Aligner** itself consumes at mapping time — the directory for STAR, one file for
chromap — as distinct from everything else in the **Index dir**. Asking for it asserts the build
finished: with no **Completion marker**, or one the **Index dir** no longer bears out, it raises here
and names the call that builds the index, rather than yielding a path that fails deep inside someone
else's command line.
_Avoid_: output, result, index files (most of them are not it), genomeDir

**Install instructions**:
The text an **Aligner** prints when its binary is absent from `PATH`, naming the bioconda package to
add. Every aligner owes one, and it must name the next command — "STAR not found" on its own is a bug,
not an error message.
_Avoid_: error message, help text, usage; and never an auto-install or a silent fallback — the package
tells you what to run, it does not run it
