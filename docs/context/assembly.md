# Assembly

What an assembly is on disk, and what it costs to make it so. This context covers `genome.py`,
`metadata.py`, `external.py` and `io/{download,fasta,twobit,utils}.py`: a name becomes a directory of
prepared files, everything in that directory is derivable from one FASTA, and nothing is derived
twice.

Words every context shares — **Assembly**, **Genome**, **Chromosome**, **Region**, **Data dir** and
the rest — are defined once in the repo-root `CONTEXT-MAP.md`.

## Language

### Naming and placing

**Assembly metadata**:
The curated TSV row cross-referencing one **Assembly** across three naming authorities:
`assembly_name`, `species`, `ucsc_name`, `ncbi_name`, `ncbi_assembly_id`, `ncbi_taxid` — and, where
the lab has pinned them, the **Source** that assembly's **FASTA** is fetched from and the sha256 of
the *unpacked* FASTA that source yields, checked after decompression rather than over the archive.
Both may be blank: an unpinned checksum is unverified rather than wrong. A cross-reference and never
an allow-list — an assembly absent from the table is perfectly legal, and a complete record handed to
`Genome(...)` replaces the row wholesale, every field or none.
_Avoid_: registry, catalog, database, manifest — each implies the table decides what exists

**Assembly dir**:
The single directory holding everything tied to one **Assembly**, `<data dir>/genome/<assembly>/` —
the FASTA and its derived files, plus the `gtf/` and `index/` subtrees other contexts own. Addressed
by assembly name alone, so the name is the only thing a caller ever needs to find any of it. It also
holds a hidden, disposable working area, `.work/`, that a registration downloads into: same
filesystem, so placing an unpacked file is a rename; nothing in it is ever claimed by a **Completion
marker**; and it is discarded once one is written.
_Avoid_: cache dir (this is the lab's reference data, not something an eviction policy may delete);
genome dir, download dir

### Getting the bytes

**Downloader**:
What puts a FASTA in the **Assembly dir**, from the **Source** the assembly's **Assembly metadata**
pins or, failing that, from a UCSC golden-path URL derived from the assembly name. Only the derived
URL is checked against UCSC first, so a typo fails on a `HEAD` request rather than three minutes
into a download.
_Avoid_: fetcher, client, mirror, provider

**Source**:
Where an assembly's **FASTA** comes from when no golden-path URL is derived: the URL its **Assembly
metadata** row pins, or the local path or URL handed to `Genome(path_or_url=...)`. Naming a source
means UCSC is never consulted about the assembly name — validation is a property of the source — and
only a pinned source is checked against a recorded checksum; a hand-supplied one is trusted as given
and degrades the assembly name to a label for the directory.
_Avoid_: input, custom genome, reference

**Genome files**:
The complete prepared set for an assembly: the FASTA plus its three derived files. Prepared as a
set, so a directory holding some of them is unprepared rather than partly prepared, and a **Genome**
either gets all four or refetches.
_Avoid_: outputs, artifacts, assets; and never for the `gtf/` or `index/` subtrees, which are the
Annotation and Index contexts' to name

**Freshness**:
The rule deciding whether a derived file must be rebuilt: it exists, is non-empty, and is no older
than every input — `make`'s rule and nothing cleverer. What makes re-running a whole preparation
free.
_Avoid_: cache hit, staleness check, up-to-date; and never a **Completion marker**, which asserts a
pipeline finished rather than one file being current

### The files themselves

**FASTA**:
The soft-masked whole-genome sequence file an assembly is prepared from — the one input every other
file descends from, kept decompressed. **Soft-masking** is data, so its lower case is carried
through and never normalized away.
_Avoid_: reference, sequence file, `.fa`

**`.fai`**:
The random-access sidecar `samtools faidx` writes beside a FASTA: one line per sequence, giving name,
length, and byte offsets into the FASTA. Its first two columns *are* **chrom.sizes** — the same two
facts, written down twice in two shapes.
_Avoid_: index, FASTA index — an **Index** is an aligner index and nothing else

**chrom.sizes**:
Every sequence's name and length, as a two-column file on disk and as the pandas Series a **Genome**
answers `chromosomes` and `chrom_sizes` from. Redundant against the **`.fai`** on purpose: this is
the shape every **External tool** asks for.
_Avoid_: genome file (bedtools' name for it), sizes, lengths, contig list

**`.2bit`**:
The two-bits-per-base encoding of the assembly, and the file every **Fetch sequence** actually reads
— the FASTA is never parsed for bases. A public artifact, not an internal cache: `Genome.twobit_path`
is a promise, made to be handed to **External tools**.
_Avoid_: index; binary FASTA, compressed FASTA

### Using them

**Fetch sequence**:
Turning a **Region** into bases — the one crossing from coordinates to sequence, and the only place
an assembly yields a `DNA`. Bounds are checked against **chrom.sizes**, **soft-masking** survives,
and a `-` **Strand** returns the reverse complement.
_Avoid_: get sequence, extract, slice, query, lookup

**Doctor**:
Reporting which **External tools** are discoverable and at what versions, before any operation needs
them. Diagnostic, never repair: it names the missing tool and the exact command that installs it,
then stops.
_Avoid_: check, health check, validate, preflight
