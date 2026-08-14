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
- **Test fixtures under `tests/data/`** — real subsampled `sacCer3` bytes, replacing inline fixtures
  for the work that needs real FASTA and GTF content. `tests/data/chimera/` adds four tiny component
  assemblies cut from those same bytes, between them carrying what no shipped assembly can
  demonstrate: a chromosome-name collision, a name that is a strict prefix of another, names already
  holding an underscore, and one holding a doubled underscore.

### Changed

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

### Removed

- **The compressed download no longer sits beside the unpacked FASTA forever.** Downloads land in a
  disposable working area and the archive is deleted once the record is written. It is kept for the
  duration of a run, so an interrupted job repairs without re-downloading.
- **`curl` is no longer shelled out to** and is no longer a required external tool.
- **The three old completion markers are gone**: the `.genome_prepared` sentinel, an index's
  `.success` flag, and the separate `<name>.index.json` sidecar, whose contents the shared record
  now carries.

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
