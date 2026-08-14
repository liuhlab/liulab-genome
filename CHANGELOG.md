# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Calendar Versioning](https://calver.org/) using
`YYYY.MM.MICRO` (e.g. `2026.6.0`).

## [Unreleased]

Preparing a reference assembly is now reproducible: naming one is enough to fetch it from a pinned
source, verify it against a pinned checksum, and record exactly what was done. A half-finished
preparation is no longer indistinguishable from a finished one.

### Added

- **A completion record for every finished build.** One `.completion.json`, written last and
  atomically, in an assembly's directory, in each annotation's, and in each index's. It records the
  source URL, the checksum, every file claimed with its size, the external-tool versions, the
  package version, and when it finished — so a directory can be explained months later. Reading it
  is the only way anything asks whether a build finished.
- **Two curated tables, shipped inside the package.** The assembly table gains a source URL and a
  sha256; a new annotation table lists what each assembly officially supports, with provider,
  version, URL, checksum and a default flag. Seven assemblies are pinned — `hg38`, `hg19`, `mm39`,
  `mm10`, `sacCer3`, `ce11` and `ecHT115` — each with a checksum computed from real bytes.
- **Registering an annotation by name.** `Genome.register_annotation("gencode_v50")` fetches,
  verifies, checks chromosome names, builds the database and records it. The path-based
  `register_gtf` remains for an unlisted GTF.
- **A chromosome-name check that runs before the database build**, so a GTF spelled for the wrong
  assembly fails in seconds rather than after many minutes. Every sequence name in the GTF must
  exist in the assembly's `chrom.sizes`; the reverse is not required, since an assembly may carry
  scaffolds an annotation never mentions. `check_chromosomes=False`, or
  `--no-check-chromosomes`, overrides it, and the record says which annotations were checked.
- **`Genome.offered_annotations`**, the table's rows for this assembly, answering a different
  question from `Genome.annotations`, which stays "registered on this machine".
- **`Genome.broken_annotations`**, the complement of `Genome.annotations`: between the two, every
  directory under `gtf/` is accounted for as registered, broken, or not begun. A half-built
  annotation is now reported where a reader would look for it — `genome annotations` marks it
  broken, says what is wrong with it, and prints the command that repairs it — rather than only
  when re-registering it. It still never stops a genome opening or hides the annotations beside it.
  `registered` and `broken` are never both true, so `registered` keeps its meaning; the `--json`
  rows gain `broken`, `problem` and `repair`.
- **CLI commands**, each emitting `--json` and exiting non-zero on failure: `genome register`,
  `genome register-annotation`, `genome register-gtf`, `genome verify`, `genome table-row` and
  `genome annotations` — the last listing what the lab offers against what is registered locally,
  without preparing anything.
- **`genome register-gtf <assembly> <gtf> <name>`**, the shell route for a GTF no table row lists —
  a collaborator's, a preprint's, one built in-house. It names the assembly rather than a directory,
  which is what lets it find that assembly's `chrom.sizes` and hold an unlisted GTF to the same
  chromosome-name check a listed one gets. Nothing is downloaded and no checksum is compared
  against, since an unlisted GTF has none pinned for it; the record carries the path it came from
  and the digest of what was placed.
- **`--infer-genes` and `--infer-transcripts`**, on both registration commands. Off by default,
  because GENCODE, Ensembl and RefSeq GTFs declare those features already and inferring them is the
  slow path. They are for a bare exon-level GTF, which otherwise registers as a database of exons
  and nothing else without saying so.
- **Chimera assemblies.** `Genome.chimera(worm, bacterium)` concatenates two or more assemblies
  already prepared on this machine into one reference and hands it back open, so a library carrying
  reads from more than one organism takes one alignment pass instead of N. Its name is derived from
  the component names, sorted and joined by `_`, and is never given — `ce11_ecHT115`, whichever
  order they were listed in — so `genome register ce11_ecHT115` builds the same thing from a shell.
  Nothing is downloaded: a component this machine has not prepared is named, with the command that
  prepares it, rather than fetched on the strength of a typed name. Every chromosome
  is suffixed `<chromosome>__<component>` unconditionally, so a bare name no longer resolves against
  a chimera and the refusal names the spellings that do. Aligner indexes are built over one exactly
  as over any other assembly.
- **`Genome.components` and `Genome.chrom_components`.** The first is the component assembly names,
  and `None` for an assembly that is not a chimera — the single test of which it is. The second says
  which component each chromosome came from, as a Series mirroring `chrom_sizes`, and is total: a
  non-chimera maps every chromosome to its own assembly.
- **A merged annotation, registered by the chimera's own build.** Each component contributes its
  default annotation, and the result is addressed by the `+`-join of their names in
  sorted-component order — `wormbase_ws298+refseq_rs_2025_06_26`. A chimera therefore arrives
  annotated, and `genome register <chimera> --force` repairs the annotation and the FASTA together.
  Components that between them contribute nothing leave the chimera with no annotation rather than
  an empty one.
- **`ce11_ecHT115` ships as a curated row** — the name, and every other column blank, including the
  sha256, since a chimera's bytes are derived here from components that are themselves pinned. It
  carries no annotation row at all. The tables stay a cross-reference and never an allow-list, so a
  chimera no row lists is still legal.
- **Test fixtures under `tests/data/`** — real subsampled `sacCer3` bytes, replacing inline fixtures
  for the work that needs real FASTA and GTF content. `tests/data/chimera/` adds four tiny component
  assemblies cut from those same bytes, between them carrying what no shipped assembly can
  demonstrate: a chromosome-name collision, a name that is a strict prefix of another, names already
  holding an underscore, and one holding a doubled underscore.
- **`DNA.outside_alphabet(text)` and the public `DNA.ALPHABET`** — the alphabet check `genome
  revcomp` applies, now askable from `import genome` and inherited correctly by `RNA` and
  `Protein`. It reports the offending characters and refuses nothing; construction still validates
  nothing.
- **`genome version --json`**, the one command that lacked the flag: `{"version": "..."}`.
- **A metadata lookup takes the table to read.** The four lookups take `table=`, defaulting to the
  shipped rows (`assembly_table()`, `annotation_table()`), and `AssemblyMetadata.from_row` /
  `AnnotationMetadata.from_row` build a record from one row — raising `MetadataRowError`, naming the
  column, rather than half-building one. Curating a row no longer means reaching past the API.

### Changed

- **`Genome.default_gtf` is a read-only property.** It was a settable attribute the registration
  path reassigned as it adopted a sole annotation; the registry now decides it, so a caller that
  used to assign to it names the annotation at construction —
  `Genome(assembly, default_gtf=<name>)` — instead.
- **What an assembly name means is a value now, in a module named for it.** *Where do these bytes
  come from* was answered inline by the downloader, which is why reading a chimera's record needed
  four deferred imports to dodge an import cycle and why three module-level functions reached into
  the downloader's privates to ask. `genome.io.source` resolves a name into one of the three
  **Source** kinds — a URL pinned or derived, a path or URL the caller seeded, or a component set —
  and the registration dispatches on which came back. The four ordered checks and their precedence
  are unchanged (ADR-0008), and `genome register <name>` is still one command for all three kinds.
- **The downloader is a registration and nothing else.** It used to inherit from `Downloader` as
  well, whose constructor never ran because its answer to *which directory?* was the wrong one; the
  plain `Downloader` is unchanged and still fetches into a cache directory of its own.
- **An external tool is one module, not five.** Locating a binary, asking its version, running it,
  running it only when its output is stale, and saying what installs it were spread across
  `external.py`, `io/utils.py`, `io/completion.py` and every aligner, with two byte-identical version
  detectors between STAR and chromap. They are now one `ExternalTool` — `path`, `version`, `run`,
  `run_to`, `install_instructions` — with two adapters: the one that shells out, and a recording
  stand-in that runs nothing. Errors, the freshness rule and the version cache are decided once, so
  the two cannot drift.
- **An aligner is given its tool instead of making one, and constructing one runs nothing.** The
  binary was resolved and asked for its version *in the constructor*, which meant a `STAR(...)` could
  not exist on a machine without STAR and every test had to patch two names to get one. Both are now
  answered on first use, and a caller may pass the tool to drive.
- **A missing aligner raises its install instructions rather than printing them to stderr.** The text
  is the exception's message, so the caller that catches it has what to do; a library writing to a
  console its caller may not have was never an error message.
- **`genome doctor` checks the tools the package actually shells out to.** `samtools`, `faToTwoBit`
  and `twoBitInfo` — the three that prepare an assembly — where it used to check `bedtools`, which
  nothing here has ever run, and neither of the two UCSC binaries `prepare_fasta` cannot work
  without. A tool that is installed but rejects `--version`, as those two do, is reported present
  rather than left out or raised on. `pixi add` commands now name the conda package rather than the
  binary, so the command in the error is one that works: `ucsc-fatotwobit`, not `faToTwoBit`.
- **The suite is two lanes, and together they are a partition of it.** `-m aligner` selects the three
  tests that build a real STAR or chromap index; `-m 'not aligner'` selects everything else. The
  aligner lane used to be the *whole* suite re-run in an environment that also had the binaries,
  which meant the three tests skipped silently in the other lane and a skip is green. `pixi run
  test-aligner` now refuses to select until both binaries answer `--version`, so that lane cannot
  report green having built nothing. `_needs` in `tests/test_aligner.py` applies the marker and the
  skip under one name, so they cannot come apart.
- **The other lane refuses to run without its tools too.** `pixi run check` and `pixi run test` now
  front `scripts/require_tools.sh`, which proves `samtools`, `faToTwoBit` and `twoBitInfo` answer
  before either selects a test. With those three off `PATH` the unit lane reported 136 skipped, 694
  passed and exit 0 — every test that writes a real FASTA, `.2bit`, chrom.sizes or annotation
  database skipping itself green, so "check is green" covered 694 of 830 tests. The per-test skips
  stay and are still right; what ends is a lane that passes having built nothing. The two UCSC
  binaries reject `--version`, so the probe runs each bare and reads who answered: the shell's 127
  and 126 mean absent and unrunnable, and any other status came from the binary itself.
- **The suite runs on eight workers and the gate runs its steps concurrently.** `pytest -n auto
  --maxprocesses 8`: 2.8 s against 5.9 s serial, `auto` finding fewer cores on a small CI runner
  where the cap does not bind. `pixi run check` moved off a sequential `depends-on` onto
  `scripts/check.sh`, which runs lint, fmt-check, typecheck and test at once and prints each one's
  output whole, in a fixed order — 4 s against 9.4 s. Measurements in
  `docs/research/test-suite-parallelism-2026-08-14.md`.
- **CI builds the docs inside the lint job** rather than on a runner of its own, `default` already
  carrying the `docs` feature. Deploying stays in `docs.yml`, the one workflow granted write access.
- **A broken registration raises instead of being quietly rebuilt or quietly trusted.** Files with no
  record mean an interrupted run; a record that disagrees with disk means something changed behind
  our back. Both raise, naming the file that differs and the command that repairs it. An absent or
  empty directory is not a broken state — that is a fresh registration.
- **Checksums are taken over unpacked content, never the archive it arrived in.** Gzip bytes change
  under recompression while the FASTA inside does not, so a hand-copied or mirrored file checks
  against the same official row.
- **One fetch step for the whole package.** Every download goes through it, and it is the only
  caller of pooch. pooch is used as a downloader; its own cache is deliberately not relied on,
  because the completion record owns that judgment.
- **`Genome` takes one metadata record instead of six per-field overrides**, bringing the
  constructor from eleven arguments to six.
- **`ucsc_name` may be blank**, for a reference UCSC has never carried. The assembly id is a local
  key and UCSC is the default source rather than the namespace; the schema had not caught up.
- **The default annotation comes from the table's flag**, so everyone reaches for the same one. An
  explicit choice still wins, and the previous rule — the sole registered annotation, otherwise none
  — remains the fallback.
- **`ce11` is sourced from WormBase rather than UCSC.** UCSC spells its chromosomes `chrI`/`chrM`
  while every WormBase annotation spells them `I`/`MtDNA`, so pairing the two would force an
  override on every registration; taking both the FASTA and the annotation from PRJNA13758 WS298
  makes them agree by construction.

- **The documentation site is four pages instead of six**, and each is a tutorial rather than an
  account of why the package is built the way it is. `Home`, `Genome`, `Sequences` and `CLI`, plus
  the generated API reference. Design rationale is not repeated on the site: a decision lives in
  `docs/adr/` and is read there.
- **A glossary term the records settled but the code does not have yet says so**, with a *decided,
  not built* marker naming the record, so the vocabulary can run ahead of the implementation without
  reading as an API that exists. No term carries one now: `Chimera`, `Component` and `Merged
  annotation` are all built.
- **An assembly's annotations are one registry, bound once to that assembly.** Whether an annotation
  is registered, broken, offered but not begun, or nothing at all used to be assembled from the same
  three scans in three separate places — as a `Genome` opened, as `genome annotations` reported, and
  as the error a name nobody registered earned. It is now settled once, and a `Genome` holds one and
  delegates to it instead of keeping three dictionaries in step by hand. The registry carries the
  assembly directory it was opened with, so it cannot file an annotation somewhere other than where
  the caller looking for it is looking.
- **`Genome.register_gtf` over a directory nothing vouches for now names a command a shell can run**
  — `genome register-gtf <assembly> <gtf> <name> --force` — rather than the equivalent Python call.
  A genome knows which assembly it is; only the by-directory `register_gtf`, which does not, still
  names the call.
- **`Genome.metadata` is always a record.** It was `AssemblyMetadata | None`, and an assembly the
  curated table does not list got `None` — so every reader guarded a missing record before reading a
  field off it, eight times on `Genome` alone. Unlisted is now a record whose fields are unknown,
  carrying the assembly's own name and nothing else, which is what a blank cell already means
  everywhere else in that table. Read `genome.metadata.species` and it is `None` when nobody knows,
  as before; the guard has nowhere left to live. A record passed to `Genome(metadata=...)` still
  replaces the row wholesale, and passing none is still optional.
- **Two accessors on the metadata table, because there are two questions.** The new
  `assembly_metadata(assembly)` is total and answers *what is known about this assembly*;
  `lookup_assembly(assembly)` still returns `None` and answers *does the curated table list this
  name*. Only the second question has a `None` answer — it is what tells a chimera's derived name
  from a free-form local key on a machine holding neither, so making it total would read `my_ref` as
  a chimera of `my` and `ref` (ADR-0003, ADR-0008). The downloader now works from the total one:
  `UCSCGenomeDownloader.metadata` is an `AssemblyMetadata` rather than `AssemblyMetadata | None`,
  so the three places that guarded it before reading a field no longer do. Which name the table
  lists is asked elsewhere and is untouched.
- **What a registration answers with has a type.** Nine API functions handed back
  `dict[str, object]`, so the command line — a thin client — re-narrowed every value it read and
  knew the completion record's key names by heart. They now return frozen records:
  `RegisteredAssembly`, `VerifiedAssembly`, `RegisteredAnnotation`, `AnnotationStatus` and its
  `AnnotationStatusRow`, each with an `as_json()` for the `--json` path; `assembly_table_row`
  returns the `AssemblyMetadata` it was always describing. **The JSON is unchanged, key for key and
  in the same order**, and so is every key written to `.completion.json` — the types wrap those
  names and never rename them. A caller that indexed a returned dict reads an attribute instead, or
  calls `as_json()` for the mapping it had before.

### Removed

- **Four annotation steps leave the package surface** — `list_annotations`,
  `list_broken_annotations`, `default_annotation` and `fetch_annotation` are no longer re-exported
  from `genome.io`, so `from genome.io import fetch_annotation` breaks. Nothing in the package calls
  any of them by name; what a caller wants from an assembly's annotations is `AnnotationRegistry`,
  which takes their place in `genome.io`. They stay importable from `genome.io.gtf`.
- **The eight metadata pass-throughs on `Genome`** — `assembly_name`, `species`, `ucsc_name`,
  `ncbi_name`, `ncbi_assembly_id`, `ncbi_taxid`, `source_url` and `sha256`. Each was one line
  guarding a record that is now always there. Read them off the record: `genome.metadata.species`.

- **`genome.external.tool_version` and the loose `_resolve` beside it.** Both are `ExternalTool`
  now — `InstalledTool(name).version` and `.path` — and an aligner's `_detect_version` and
  `install_instructions` are gone with them, the tool answering both.
- **Four tests that asserted nothing the rest of the suite did not.** The smoke test (every module
  imports the package, and the CLI's `version` command is tested on its merits); a `bedtools`
  version check strictly subsumed by `doctor`; and an assertion that `download` re-exports four
  names `registration` defines, which held the module wiring rather than any behaviour and would
  have broken on a refactor that changed nothing observable.
- **`Usage`, `Genome files` and `Annotations & indexes` are gone as site pages.** What a user needs
  from them moved into `Genome` and the new `CLI` page; what remained described machinery a caller
  never invokes directly. Links to those three URLs break.
- **The compressed download no longer sits beside the unpacked FASTA forever.** Downloads land in a
  disposable working area and the archive is deleted once the record is written. It is kept for the
  duration of a run, so an interrupted job repairs without re-downloading.
- **`curl` is no longer shelled out to** and is no longer a required external tool.
- **The three old completion markers are gone**: the `.genome_prepared` sentinel, an index's
  `.success` flag, and the separate `<name>.index.json` sidecar, whose contents the shared record
  now carries.
- **`genome.io.utils._run` and `_run_to`**, a name-addressed layer that restated
  `ExternalTool.run`/`run_to` and kept a second copy of the freshness branch. Each FASTA
  preparation step now holds its own tool and calls `run_to` on it.

### Fixed

- **An annotation is no longer counted as registered merely because its database file exists**, so a
  `gffutils` build killed halfway is reported as broken rather than silently queried.
- **Re-registering a valid annotation is a silent no-op** rather than a warning that returned an
  annotation which might not exist.
- **`genome table-row` reports an existing pinned checksum rather than enforcing it**, since it is
  the command a maintainer runs precisely when an upstream file has changed and the pin must be
  regenerated.
- **A blank identifier cell in the assembly table reads back as unknown**, rather than as the string
  `nan` — or, for a blank taxonomy id, an exception. The row `genome table-row` emits for an
  assembly the table does not list yet leaves the species and the NCBI identifiers blank, so pasting
  that line in is now a working route rather than one that breaks the next lookup.
- **`genomeChrBinNbits` is computed and passed to STAR.** It was never passed at all, so STAR's
  default of 18 always held and every sequence was padded up to a whole multiple of 262,144 bases —
  a ×5.18 inflation of the genome file on an 87-scaffold draft like `ecHT115`, and a silent one,
  since STAR neither warns on this parameter nor clamps it. It is now sized from the mean sequence
  length and the read length `sjdb_overhang` implies, and
  passed even when it lands on 18, so the record says what the build asked for rather than what it
  left out. Passing it yourself still wins. Measured in
  `docs/research/aligner-index-params-and-reference-names.md`.
- **An index record pins the digest of the assembly it was built from.** It recorded the FASTA's
  path and never its bytes, so re-registering a reference left every index over it still reading as
  finished — chromap worst of all, whose index file stays byte-identical while the sequence names
  beneath it change. `details["assembly_sha256"]` is copied from the assembly's own record and
  compared record to record, with no sequence bytes read; a disagreement names both digests and the
  rebuild that repairs it. An index built before this change pins nothing, so it reads as *unknown*
  rather than as wrong and stays unguarded until it is next rebuilt.
- **Registering with `--no-check-chromosomes` is no longer told to register the assembly first**,
  which it may well have done already. An annotation record now says *why* the names went
  unchecked — `details["chromosomes_unchecked_because"]` is `"caller-override"` or
  `"no-chrom-sizes"`, and `None` when they were checked — so both registration commands report the
  state they are in, advice included only where there is any to give. A check that ran and passed
  now says so too, rather than being reported by silence. A record written before that field reads
  as *unknown* rather than as either reason.

### Upgrading

There is no migration. A directory prepared by an older version has no completion record, so it
raises and needs one forced re-registration: `genome register <assembly> --force`, and
`genome register-annotation <assembly> <name> --force` for each annotation. A forced re-registration
keeps an unpacked FASTA whose checksum still matches and rebuilds only the derived files.

Interrupting a *first* registration leaves a directory holding files with no record, which raises
next time and needs the same forced re-registration. That is the accepted trade-off, chosen over
silently resuming.
