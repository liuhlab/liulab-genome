# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Calendar Versioning](https://calver.org/) using
`YYYY.MM.MICRO` (e.g. `2026.6.0`).

## [Unreleased]

### Added

- **`genome assembly list`, and the `assembly_status()` behind it.** The CLI could not answer
  the first question a new user has. `genome assembly register` is the first command anyone
  runs and nothing said what to put after it; the only way to find out was to import
  `genome.assembly` and read the table. The new command prints the eight assemblies the
  metadata table offers **set against the ones prepared on this machine**, and names
  `genome assembly register <name>` and `genome assembly verify <name>` as the next actions.
  It downloads, prepares and builds nothing, and exits `0` on a machine holding none of them.
  **The second half is what nothing in the package could ask before**: `is_prepared` answered
  for one name a caller already held, and nothing enumerated the tree — so an assembly
  registered from the UCSC golden path, which is a legitimate registration with no row behind
  it, was invisible. Such a name now reads as `registered, not offered`, and a directory in
  the tree with no record beside it reads as `here, not registered` rather than being reported
  as absent or as an assembly. `assembly_status`, `AssemblyStatus`, `AssemblyStatusRow` and
  `is_prepared` are exported from `genome.assembly`; the row's `state` is derived from its
  booleans and stays out of the `--json` payload, as the annotation report's does. Integrity
  is deliberately not asked here — `genome assembly verify` owns it, and asking it cheaply
  would report *unchecked* in the words of *checked*.

- **Twelve error types the package raises are now importable.** Every exception class this
  package can hand a caller resolves from a public name, so catching one no longer means
  importing from a module the API reference declares free to move between releases.
  `ToolNotFoundError`, `ShippedTableError`, `MetadataRowError`, `AmbiguousDefaultAnnotationError`,
  `CuratedGeneListError` and `GeneListAssemblyMismatchError` come from `genome`; the last two
  also from `genome.annotation`; and `RegistrationError`, `RegistrationMismatchError`,
  `UnfinishedRegistrationError`, `ChecksumMismatchError`, `PreparedSetNotDownloadedError` and
  `PreparedChecksumError` from `genome.store`, which is now rendered on the API reference page.
  **Two of these change what a caller can express**, rather than only where they import from:
  `ShippedTableError` and `PreparedSetNotDownloadedError` are the parents of eight error types
  that were already public, so `except ShippedTableError` to catch any broken shipped table, and
  `except PreparedSetNotDownloadedError` to catch any set that was never downloaded, could not
  be written before and now can. A guard test holds the invariant rather than the list: an
  exception class added anywhere in the package with no public name fails the suite.
  `genome.store` re-exports these six and stays closed to callables, which is what keeps the
  suite's offline guard able to patch the fetch step on the module.

- **The lab's three static writing gates are here and runnable, and they gate nothing yet.** The
  router states the writing rules in three separate places and nothing checked any of them, so a
  document over its cap, a phrase the lab has ruled out, or an agent-facing tree nobody declared
  was caught by a reviewer or by nobody. `pixi run vale`, `pixi run markdownlint` and
  `pixi run conformance` come across from `liulab-repo-template` at `9c4f8ed` — the four
  `styles/Lab/` rules, `.markdownlint-cli2.yaml` and `scripts/conformance.py` byte for byte,
  `.vale.ini` rewritten around this repository's paths. **None of the three is a step of `check`
  or of CI**, which is the whole reason they land first and alone: the tree fails vale with 5
  errors and markdownlint with 966 across 40 files, and a gate that is red on the day it arrives
  teaches nobody anything. Cleaning the prose, and then wiring the three in, are the commits after
  this one. `pyproject.toml` gains `[tool.liulab.agent-docs]`, this repository's agent-facing set
  spelled with its real paths: `AGENTS.md` and `docs/agents/` take the 1000-word document cap,
  `docs/adr/` the tighter 400-word record cap, `docs/research/` none. **`CONTEXT-MAP.md` and
  `docs/context/` take none either**, and that is a decision rather than an omission — a glossary
  is capped per entry, because how long the file runs reports how many terms the domain happens to
  have, which measures nothing about the writing (ADR-0024). Every `.vale.ini` section names every
  rule, including the two it leaves off: vale's `*` matches `/`, so the sections overlap and the
  per-rule settings pile up with the later section winning, and a rule left out of a section keeps
  whatever an earlier one said instead of defaulting on. That is how a record tree ends up checked
  by nothing with a green run to show for it, and conformance now fails if any of the six declared
  keys loses its section. One rule is waived, `agent-docs-unpublished`: `mkdocs.yml` keeps the
  record, agent, context and research trees out of the built site with `exclude_docs`, so the
  `search: exclude: true` front matter that rule wants would sit inert on 48 files the site never
  builds. The waiver names its reason, suppresses that one rule and nothing else, and prints on
  every run including the green ones — so the day the site generator changes and the premise stops
  holding, it is on screen rather than forgotten. `.gitignore` gains `.cache/` and `/.gemini/` and
  anchors `.codex/`, which matters now that markdownlint hands it the whole exclusion list.

### Changed

- **A missing type annotation now fails the gate, and the type checker is pinned exactly.** The
  full-annotations rule was stated as non-negotiable with nothing behind it — pyright ran at
  `basic`, where a missing annotation is not a diagnostic — and the router's own rule text said so.
  **Raising the mode is not by itself the mechanism**, which is the thing worth writing down:
  measured against pyright 1.1.411's own rule tables and confirmed on a probe file, `standard` adds
  exactly five checks over `basic` — possibly-unbound variable, incompatible method override,
  incompatible variable override, overlapping overload, function member access — and leaves
  `reportMissingParameterType` and `reportUnknownParameterType` at `none` exactly as `basic` does.
  Only `strict` turns those on, and it turns on **556** other diagnostics with them, which is a
  different decision than this one. So the two rules are named beside the mode, and that is what
  makes the rule's text true. **All of it reports zero errors across 156 files**, so there was no
  fallout to fix: every function under `src/`, `tests/` and `scripts/` already carried the
  annotations the rule asked for, and what changes is that the next one that does not fails
  `pixi run check` rather than depending on a reviewer noticing. pyright is pinned `==1.1.411` where
  every other dev dependency floats — the language server an editor runs is this same checker, so a
  version that drifts between the two is the editor and CI disagreeing about one file, and a release
  of it would turn CI red on a commit that changed nothing. It is the version already resolved, so
  `pixi.lock` did not move; bumping the pin becomes a commit of its own.

- **Python narrows to 3.13, so an install on 3.12 now fails rather than succeeding untested.**
  `requires-python` floored at `>=3.12` while the lock held one interpreter, so no lane ever ran
  the version the floor advertised. The two settings that looked like they held it — `ruff
  target-version` and `pyright pythonVersion` — read source syntax and not a runtime, and the
  failure a 3.12 user would actually hit is a 3.13-only behaviour in a dependency, which neither
  can see. Both are **deleted rather than retargeted**: ruff derives its target from
  `requires-python` and pyright takes its version from the environment's interpreter, so Python is
  declared in exactly two places, `requires-python` and the pixi pin. The `Programming Language ::
  Python :: 3.12` classifier is gone with them, which is what makes a 3.12 resolver refuse the
  package instead of installing it. A 3.12 test lane was the alternative and stays the answer if a
  consumer pinned to 3.12 turns up; it lost because nothing is known to need one and every bioconda
  pin would have to resolve on a second interpreter (ADR-0025). No solved package changed — the
  environment was already 3.13.

- **The package's docstring examples now run, in the unit lane.** `AGENTS.md` asks for at least
  one runnable example on every public object and the package carries over twelve hundred of them;
  nothing executed a single one, so the bar was held by review alone. `pytest` now collects
  `src/genome` with `--doctest-modules`, `ELLIPSIS` set at config level, which puts every example
  in `pixi run check` and in CI beside the tests it already runs. Three examples had drifted and
  are fixed. `genome.store.prepared.login_node_help` promised in its `Returns` prose and asserted
  in its example a sentence ending in the quoted command; the command moved mid-sentence and
  neither followed. The two `UCSCGenomeDownloader` examples each had a line whose trailing text
  read as a `# doctest: +SKIP` directive and was prose — so the one line that would fetch a genome
  from UCSC was the one line *not* skipped, harmless only because the line above it had been
  skipped into a `NameError`. **The suite's two autouse guards were given reach over `src/` first,
  and proven, before either directive was repaired.** A conftest's fixtures reach the directory it
  sits in and nothing above, so a doctest item collected from `src/` ran behind neither the
  offline guard nor the per-test **Data dir** — and repairing the directive while that was true
  would have turned a latent hazard into a live multi-gigabyte download from CI. Both guards now
  live in `tests/_guards.py`, loaded by a root `conftest.py` sitting above both trees, and a test
  runs a throwaway example in a pytest of its own to prove the reach rather than assert the shape
  of a file. The examples carrying `+SKIP` still execute nothing; whether an example nothing runs
  meets the bar is a separate question and is not answered here. **One suite-wide change came
  with the load**: hypothesis's default 200 ms deadline is a per-test wall-clock budget, and in
  a ten-worker lane where several modules spawn process pools of their own it reads the
  scheduler rather than the code — it fired three times in sixteen runs once the examples were
  collected, and never in eight without them. It is now off for the suite, which is the same
  answer this repository already gave when it declined to assert such a budget itself.

- **The built site carries no architecture-decision-record numbers.** A record number is a
  citation between agents: it resolves for anyone reading the source, because the record tree is
  in the repository beside it, and for nobody reading the site, which excludes that tree on
  purpose. One hundred and eleven numbers were published on the API reference page — seventy-six
  rendered from docstrings, thirty-five inside "Source code in ..." listings. Docstrings keep
  citing records, and `mkdocstrings` now drops the citation as it renders; source listings are no
  longer published, since a listing with the numbers stripped out would disagree with the file it
  names. `genome.__doc__` and `help()` are untouched. Two guards hold it: a record is cited as a
  trailing parenthetical — `(ADR-0006)` — which the test suite checks so the removal is total,
  and `pixi run docs-build` fails if a number reaches the built site by any route.

- **The API reference no longer embeds source listings.** Use the source links on GitHub. The
  page still carries every signature, argument type and attribute exactly as the code declares
  them, which is what it is for.

- **User-facing messages no longer cite architecture-decision-record numbers.** Nine exception
  messages, one CLI command's `--help`, and the `limits` string that rides on a *successful*
  symbol match dropped a trailing `(ADR-00xx)`. The records live in a tree kept out of the built
  site on purpose, so the number resolved nowhere a reader could reach; every message still says
  in words what the record decided, and still names the same next action. **The `limits` change is
  the one a caller may be reading**: it is not an error path — `XrefSet.for_symbols(...)`
  `.match_symbols(...).limits`, the `"limits"` key of `as_json()`, and what `genome xref symbols`
  echoes to stderr on every mouse or worm run all carry it, so code matching that string exactly
  needs updating. A guard test holds the invariant rather than the list: a record number in any
  string this package can print, or in any CLI command's docstring, fails the suite. Comments and
  the docstrings of ordinary objects keep their citations, which are between agents and stay.

- **BREAKING — the CLI grows a sub-app per topic, and fourteen commands are renamed.** `genome
  register` is `genome assembly register`, `genome tf-gene-list` is `genome tf gene-list`, `genome
  xref` is `genome xref ids`, `genome match-symbols` is `genome xref symbols`, `genome homologs` is
  `genome homology links`, `genome motif-scan` is `genome motif scan`, and the eight remaining
  moves follow the same rule — a sub-app named for the package its commands ship from, which is
  why the Orthology context's is `homology`. `version`, `revcomp` and `doctor` belong to no topic
  and are unchanged. Each topic package now ships its own CLI module and its own renderers beside
  the result types they render, so what `gene-list` prints and what a `GeneList` holds change in
  one place; `genome.cli` keeps the three commands, the `add_typer` calls and no rendering helper,
  and the console script still points `genome` at `genome.cli:app`. The seam under the commands did
  not move: no command constructs a `Genome` it did not construct before, every one still takes
  `--json`, and every error message and docstring naming a command names the new spelling.
  **Every old spelling still runs for this one release**, hidden from `genome --help` and printing
  a deprecation notice on **stderr** so `--json` on stdout still parses — one function object per
  command registered under both names, from one table and one loop in `genome.cli`. The single
  exception is `genome xref`, whose old spelling is a sub-app's name as well as a command's: the
  `xref` group answers it by handing an unrecognised first token to `ids`, with the same notice on
  stderr. **The aliases go in the release after this one**, deleted as a unit with that table.
- **`io/` retired: every context owns its own I/O, and what is left is a store.** The one
  directory grouped by *kind of work* is gone, and its modules moved to whatever they are about.
  `genome.assembly` is the whole Assembly context — `genome.py`, `chimera.py` and the seven
  `io/` modules its glossary already claimed (`source`, `components`, `download`, `chimera` as
  `chimera_build`, `registration`, `fasta`, `twobit`). `genome.annotation` is the annotation
  package plus the shipped **Curated gene list**, which takes the glossary's own name as
  `annotation/curated.py` so that a module and the registry's `gene_list()` function no longer
  share one namespace. `genome.store` keeps what belongs to no context and is reached by all of
  them: the **Data dir** root (`data_dir.py`), the one fetch step, the **Completion marker**,
  checksumming (`io/utils.py` → `store/checksum.py`) and the **Prepared set** pipeline — and a
  test reads every file under it to hold that none imports a context back. **`metadata.py`
  split**: the assembly table is `genome.assembly.metadata` and the annotation table is
  `genome.annotation.metadata`, the two having shared their shape and nothing else. **Every
  root under the **Data dir** now sits with the code that reads what lives under it** — `motif/` in
  `tf.motif.jaspar`, `xref/` in `xref.xref`, `homology/` in `homology.compara` — and the pipeline
  they share declares none of the three. No deferred import was added; the three that existed
  are carried over with their reasoning. `tests/` mirrors the new tree, one package per package,
  and every import-edge guard moved with the module it defends. **Nothing the package exports
  changed**: `genome.__all__`, the CLI's commands and its `--json` are untouched, and the
  shipped `data/` tree did not move.
- **One module writes every shipped table, where three build scripts each declared the writer
  again.** `genome.shipped_writer` owns the unquoted TSV rendering with its header, the
  deterministic gzip and the provenance merge that replaces a row by its key and re-sorts the
  file; `scripts/build_tf_census.py`, `scripts/build_tf_cofactor.py` and
  `scripts/build_tf_links.py` keep only their publisher's recipe. **The writer is handed the
  reader's own table declaration**, so a file is held to the header, the required columns, the
  flag spellings and the key it will be read under *before it reaches disk* — the column tuples,
  the file names and the value separator are declared once and imported by both halves, and no
  build script restates a reader's constant or spells a file suffix of its own. Each generator
  supplies its own error class and repair, so a refused build names the recipe to fix rather than
  telling whoever ran it to run it again. The byte-stability promise the gzip call carries is
  stated in that one module, where three write functions and three module docstrings used to
  carry it between them. `tests/test_shipped_writer.py` holds it
  to all of that with no download, and re-renders every shipped table's own rows back to the bytes
  that ship. `VALUE_SEPARATOR` moves to `genome.shipped` beside the flag spellings and is
  re-exported from `genome.tf`; each format's declaration is public as `CENSUS_FORMAT`,
  `COFACTOR_FORMAT`, `LINK_FORMAT` and their provenance peers, and the two species-keyed
  provenance tables now declare that key. **No shipped `.tsv` or `.tsv.gz` byte changed.**
- **`io/gtf.py` split four ways, and gene id stem resolution became findable.** The largest
  module in the package held four clusters that barely touched each other, and is now
  `genome.io.annotation`, a package of four: `registration.py` puts an annotation on disk — the
  fetch, the placement, the **Chromosome** check, the repair-command strings, the **Completion
  marker**, the **Merged annotation** a chimera build derives, and the two registrars addressed
  by assembly name; `registry.py` holds `AnnotationRegistry`, the three scans, the **Default
  annotation** rule and the by-assembly-name questions; `stems.py` holds the **Gene id stem**
  crossing that the Xref, Orthology and TF contexts all make; and `database.py` is the `gffutils`
  adapter, sixty lines, **the only module in the package that imports the library** — held by a
  test that reads every file under `src/` rather than only the four. Resolving stems still walks
  the database's gene rows one at a time, now over a generator that yields off the cursor, so a
  GENCODE-sized annotation is still never held in memory. `AnnotationRegistry` stays one class
  with one interface and calls across the four. **Nothing a caller can see changed.** `genome`
  and `genome.io` export exactly the names they exported before, `genome.__all__` is untouched,
  and every command's `--json` emits the same keys in the same order. `tests/test_gtf.py` splits
  the same four ways into `tests/annotation/`, with the shared registration helpers in one
  `conftest.py`; every test that was there is still there, under its own name.
- **Every result type moved beside whatever returns it, and `genome.io.results` retired.** The
  module held frozen dataclasses whose one shared property was being returned; a third of it
  belonged to contexts no module under `io/` imports. `RegisteredAssembly` and `VerifiedAssembly`
  are now in `genome.io.download`, beside the two functions that build them; `RegisteredAnnotation`,
  `AnnotationStatus`, `AnnotationStatusRow`, `GeneList`, `GeneListSource`, `ResolvedGeneIds`,
  `chromosome_check_summary` and `annotation_register_command` in `genome.io.gtf`; `ResolvedStems`,
  `ResolvedXrefIds`, `ResolvedSymbols` and `SymbolMatch` in `genome.xref.xref`; `HomologyLink` and
  `HomologyAnswer` in `genome.homology.compara`; `ResolvedHomologs` in
  `genome.homology.annotation`. `ResolvedGeneIds` is the one type two contexts share and it lands
  beside `resolve_gene_ids`, because the rule is where a type is returned and not where it is read.
  The module had been opened as a seam both halves of `io` could reach, and that reason expired
  once each type sat beside its producer; the guard that closed it is replaced by its assembly-half
  mirror in `tests/test_download.py`, `io.gtf` already holding the annotation half's. The `as_json`
  convention the module's docstring stated is now ADR-0022. **Nothing a caller can see changed.**
  Every name `genome`, `genome.io` and `genome.xref` exported still imports from the same place,
  `genome.__all__` is untouched, and every command's `--json` emits the same keys in the same order.
- **One module reads every shipped table, where six brought their own loader.** `genome.shipped`
  owns resource lookup, gzip, header validation, cell parsing, the blank-cell rules, duplicate-key
  refusal and the shape of the error a broken file raises; `metadata.py`, `xref/metadata.py`,
  `homology/metadata.py`, `tf/gene/census.py`, `tf/cofactor/table.py` and `tf/link.py` keep a
  **table declaration** — resource path, columns, row type, the noun the table is called by and the
  command that repairs it — and nothing else. Six failures are now checked in one place and reach
  every table: an empty file, a header that is not the declared columns, a row with the wrong cell
  count, a blank required cell, a flag spelled a way no table spells one, and a repeated key. The
  last is declared per table, since a motif link table carries many rows per **Gene id stem** by
  design. **The curated assembly, annotation and xref tables gain the header validation they
  lacked**: those two readers went through pandas and checked no header at all, so a renamed or
  missing column reached the cell parser as a blank cell — or read as `None` where the column was
  optional. Every message still names its own noun and its own repair, pinned word for word by
  `tests/test_shipped.py`, and several gained a repair they never carried. `species_slug` and
  `parse_cell` move to live with the reader; `genome.metadata`, `genome.tf.gene` and
  `genome.tf.cofactor` re-export them, so every existing import path still resolves. No shipped
  `.tsv` or `.tsv.gz` byte changed — this is the reader half only.
- **One module prepares every release-pinned set, where three had each written the pipeline out.**
  A **Prepared set** — a **Motif set**, an **Xref set**, a **Homology set** — is files pinned to a
  **Release**, belonging to no assembly, filed beside the assembly tree under the **Data dir**.
  `genome.io.prepared` now owns the whole of preparing one: the set's directory, the working area,
  the fetch, the digest, the staged rename, the **Completion marker** and the one sentence that
  sends a caller to a login node. A source declares a URL, a checksum and how to slice or parse
  what arrives, and nothing else — a test prepares a fictitious fourth set end to end on exactly
  that. The three roots under the **Data dir** are declared together there (`homology/` included,
  which used to be spelled in `homology/compara.py`), and each context keeps its own reader, its
  own answer types and its own not-downloaded error quoting its own exact prepare command.
- **Where a checksum is enforced now follows what it covers, and decompress-while-hashing has one
  implementation.** A pin over the **unpacked** bytes (ADR-0006) is checked as the file is streamed
  and unpacked; a pin over the archive as served — Ensembl Compara's own md5 of its `.gz` — is
  pooch's `known_hash`, as before. The same streaming step digests the stored slice and every
  re-read of it, so the whole-slice-into-memory read in `xref/xref.py` is gone and the four
  spellings of *unpack while hashing* are one. The working area is now kept exactly when something
  can vouch for what is in it: with an archive pin pooch re-checks a leftover download, so an
  interrupted 110 MB fetch still costs no second download; with no such pin a leftover is
  unverifiable and is swept before fetching rather than adopted.
- **A JASPAR release now writes a Completion marker, and is prepared one set per directory.**
  It used to write none, substituting an atomic rename and a motif count on the grounds that the
  files are under a megabyte. That covered *is this finished* and missed the other half of what a
  record is for — the only answer to *how was this made*: the URL, the package version, when, and
  what the bytes hashed to. JASPAR publishes no checksum to pin, so what is recorded is the digest
  of what was stored and every re-read is held to it; the motif count and the base-id check stay,
  because they say the file is the *right release, whole*, which no digest of ours can. A release
  is therefore prepared in `motif/jaspar/<release>/<tax group>/` rather than flat under
  `motif/jaspar/`; **a file prepared by an earlier version is left where it lies and prepared again
  under the new layout**, one download of under a megabyte. A fetch that fails now raises
  `MotifSetNotDownloadedError` naming the call to make on a login node, where it used to surface
  pooch's own transport error.
- **The TF context moved out of the Annotation module, and the registry kept one seam.**
  `AnnotationRegistry.resolve_gene_ids` is now the only identifier surface `genome.io.gtf` exposes:
  it answers *which gene ids does this **Gene id stem** name here* and knows nothing about what it
  is handed a list of. Roughly 900 lines of TF code moved out from behind it — `TFGene`,
  `TFGeneList`, `NoTFCensusError` and the census plumbing to the new `genome.tf.gene.annotation`;
  `TFCofactor`, `TFCofactorList`, `NoCofactorTableError` and the table plumbing to
  `genome.tf.cofactor.annotation`; and `UnknownSpeciesError`, shared by both halves, to
  `genome.tf.species`. `genome.io.gtf` now imports nothing under
  `genome.tf`, held by a guard: those import lines used to pull in sixteen `genome.tf`
  modules — both shipped-table readers, the motif link table and the whole motif tree down to the
  scan, its worker pool and its Parquet sink — so neither module drags the motif tree in any
  longer. The package still re-exports the five TF names eagerly, so `import genome` loads
  `genome.tf` as it always did — what closed is the edge from `io/`, not the cost of the
  re-export. **Nothing a caller can see changed.** `import genome` still yields `TFGeneList`,
  `TFCofactorList`, `NoCofactorTableError`, `NoTFCensusError` and `UnknownSpeciesError` under those
  names and `genome.__all__` is untouched; `Genome.tf_gene_list()` and `Genome.tf_cofactor_list()`
  answer as before and now delegate into `genome.tf` rather than into the registry; `genome
  tf-gene-list` and `genome tf-cofactor-list` keep their stdout, their stderr attribution lines,
  their `--json` records and their exit codes. Adding a second bio topic that wants an annotation's
  own gene ids is now one directory under `src/genome/` rather than an edit to four modules, three
  of which have no stake in it.
- **Worker resolution left the motif package.** `genome.tf.motif.workers` is now
  `genome.workers`. It answers exactly the same one question — how many worker processes a run
  may use: an explicit count as given, else the Slurm allocation, else this process's affinity,
  else the machine — and it always had no domain coupling: about a hundred lines that import
  nothing else from this package and name no motif, assembly or region. `genome.cli` already
  reached past the motif package to get at it, which was the tell. `genome.tf.motif` goes on
  re-exporting `DEFAULT_WORKERS`, `SLURM_CPU_VARS` and `resolve_workers`, so
  `from genome.tf.motif import resolve_workers` keeps working unchanged and no answer changes.
- **CI no longer recompiles memelite's JIT on every run, and the suite is a third of its size.**
  The `test` lane's largest single item was not a test: `memelite`'s scan and compare engines are
  `@numba.njit(cache=True)`, numba writes that cache inside the pixi env, and `setup-pixi` saves
  its env tarball before any test has run — so every run recompiled from source, 27.5 s, 47% of
  the lane's wall. `NUMBA_CACHE_DIR` now points at a workspace path cached on the `pixi.lock`
  hash. Separately the suite was cut **2425 → 933 tests with coverage unchanged at 98%** (60 → 61
  missed statements of 4149), by merging sibling one-assert tests that shared a fixture and
  collapsing parametrize tables that walked one code path with many rows of data; all 169 test
  classes and all 154 distinct refusal messages survive, checked as whole-tree set differences
  against the previous commit. Serial **56.3 s → 24.0 s**, parallel **17–19 s → 8.4 s**, the
  `check` gate **24 s → 10 s**. Measurements and method:
  `docs/research/test-suite-cost-and-parallelism-2026-08-30.md`.
- **The unit lane distributes by group rather than by load.** Nine tests spawn their own process
  pools; distributed by load, several could fork pools at once and the lane's wall went bimodal —
  13.5 s or 16.1 s on the same commit, decided by scheduling alone. `--dist=loadgroup` pins each
  such module to one worker and is faster on a 4-core runner too (9.4 s against 12.3 s). The
  worker cap moves 8 → 10; on a CI runner `auto` finds four and the cap does not bind.
- **The docs site is rewritten for a reader who writes Python and is not a genomicist.** Three long
  pages become eleven task pages under three nav sections: Get started, then *Working with a genome*
  (Assembly, Sequences and regions, Annotations, Aligner indexes), *Topics* (Transcription factors,
  Motifs, Gene identifiers, Homology) and *Interfaces* (a CLI overview and two command pages, plus
  the API reference). `docs/genome.md`, `docs/sequences.md` and `docs/cli.md` are gone, split into
  the pages that replace them. **`xref` and `homology` get Python coverage for the first time**, and
  so do the JASPAR motif surface and the shipped TF-to-motif link table, which appeared nowhere
  before. Biology and file-format context get a sentence; design rationale and ADR numbers get
  nothing, since `docs/adr/`, `docs/context/` and `docs/research/` are excluded from the built site.
  Every example is read-verified against the source, and every quoted output is one that was
  actually produced. **`docs/reference.md` drops from 59 module directives to eight package
  sections**, rendering only what each package's `__init__.py` re-exports, so the page no longer
  publishes some fifty private modules to a reader who should never reach for them; `genome.store`
  re-exports nothing on purpose and so gets no section. The README is rewritten to the same pitch
  and points at the site. Nothing under `src/` changed.

### Added

- **Which genes in another species a gene is homologous to, on Ensembl Compara's own trees.**
  `genome.homology` is the Orthology context and a peer of `genome.tf`. `HomologySet(species,
  other_species)` pins **Compara release 116**, fetches the per-species dump that holds that pair
  once into `$LIULAB_DATA/homology/` — a sibling of the assembly tree, beside `motif/` — checks the
  publisher's own md5 against the bytes as they arrive, slices the pair out a row at a time, and
  records a **Completion marker** carrying both checksums: the publisher's as provenance and the
  derived slice's own sha256 as the integrity check, re-verified on every read. All three pairings
  among human, mouse and worm answer, in either direction, off one prepared file per pair. The slice
  is a plain gzipped TSV with Compara's own header and Compara's own rows, so a collaborator reads it
  in R without this package. **Nothing here is computed**: every field of every link is a cell of
  that file, and this package publishes no quality score, no ranking and no "best ortholog" of its
  own.
- **The Compara partition guard, which is the trap this exists for.** The per-species dumps are a
  de-duplicated partition *at the pair level* — Compara's own README calls each file "an arbitrary
  subset of orthologies involving the given genome" — and which file holds a pair is **not stable
  across releases**: counted on 116, the human file holds **0** human↔mouse rows and 23,982
  human↔worm; the mouse file holds 23,764 human and 25,006 worm; the worm file holds neither. (On
  110 the human file held 16,242 mouse rows, so the assignment really does move.) Which file holds
  which pair is therefore a **measurement** in the shipped `homology_metadata.tsv`, re-taken on every
  prepare: a slice that comes back empty **raises naming the other file** and writes nothing, rather
  than answering empty. A pair is never partially present, which is what makes zero trustworthy.
- **A Homology type is the publisher's and is never recomputed** (ADR-0020). `HomologyLink` carries
  Compara's `homology_type` verbatim, its high-confidence flag and both quality scores; resolving an
  answer into a registered **Annotation** through `resolve_gene_ids` — used unchanged, with no
  `Genome` mixin and no `Genome`-level convenience — leaves the label alone and reports the
  **Dropped partner**s separately, so a view that *looks* one-to-one is never mistaken for one.
  Orthologs are the default and paralogs come back only on request. **Measured, and worth stating:
  release 116 publishes no cross-species paralogy at all** — zero `between_species_paralog` in the
  whole 4.0 M-row human dump, and every `other_paralog` (128,020), `within_species_paralog` (13,144)
  and `gene_split` (9) row relates two genes of *one* species — so on this release the switch
  changes nothing for these three pairs. It is kept because it is where such a row would land and
  because *not an ortholog* must stay distinguishable from *absent* (ADR-0013).
- **A quality score that is null for a whole species pair says so before a filter empties.**
  `goc_score` and `wga_coverage` are `NULL` on **100%** of the rows of *either* worm pairing — all
  23,982 human↔worm and all 25,006 mouse↔worm — where the human↔mouse pair carries real values.
  Which columns a set holds nothing in is measured over the prepared slice, not listed against a
  pair, and rides on the set and on every answer. Compara 116 is the release pinned because it is the
  newest that publishes both an `MD5SUM` and the compressed naming: only releases 90 and 116 publish
  an md5 at all, and 113 ships these dumps uncompressed — so each row carries a URL read off the
  release's own listing rather than one built from a template, and a row with no checksum is refused.
  **Orthology is served and never consumed** (ADR-0019): no table this package publishes is derived
  through homology, held structurally — nothing outside `genome.homology` may import it, which a test
  asserts by reading every module's source.

- **A gene symbol reaches the package's answers, and a retired spelling no longer drops silently.**
  `XrefSet.match_symbols(symbols)` matches **approved, previous and alias** spellings, answers with
  every **Gene id stem** any of them names, and says on each hit which kind it matched — ambiguity is
  the return type and not an edge case, so `ADCY3`, HGNC's approved symbol for one gene and a symbol
  it retired from another, answers with both. **The measured failure it prevents:** of EpiFactors
  v2.0's 801 human rows all 801 carry an HGNC id and all 801 resolve, but **31 still spell the gene
  by a symbol HGNC has retired** — `ARNTL` for `BMAL1`, `C11orf30` for `EMSY` — and an approved-only
  join mis-keys or drops exactly those. Matching is **exact by default**, so `Brca1` asked of a human
  set matches nothing rather than half-working; `case_insensitive=True` is opt-in and still returns
  **every** gene matched rather than picking one. **From a shell it is `genome match-symbols
  <species> <symbols>...`** — a command of its own, because a symbol is a verb of its own and not a
  third direction of `genome xref`: it parses arguments, makes that one call and renders, so
  `import genome` and the shell hit one code path and `--json` is `as_json()` verbatim. The matches
  go to stdout tab-separated in the four columns `asked, symbol, gene_id_stem, kind` — the last three
  being the keys `SymbolMatch.as_json()` writes, in its order, so the two renderings cannot drift —
  with the heading, the URL, the counts and what this source could not have matched on stderr.
  **Every symbol passed leaves with at least one row**, the ones that matched nothing getting one
  with the other columns empty. `--case-insensitive` is the opt-in fold and still prints every match;
  `--source` names an **Xref source**, and one named that carries no symbols exits `1` naming the one
  that does.
- **A symbol lookup needs no source named, because a default is per species *and per question***
  (ADR-0021). The curated table flags a **Default xref source** twice — once for identifiers, once
  for symbols — because the two are usually different rows: human's identifiers come from `alliance`,
  whose 25 MB cross-reference file publishes **no human symbol at all**, and its symbols from `hgnc`;
  mouse's and worm's from `alliance_bgi`, a third source again. So `genome match-symbols "Homo
  sapiens" ARNTL` and `genome xref "Homo sapiens" --from-stems symbol ENSG00000133794` answer on the
  first try, where both previously exited `1` telling you what to pass. **The fill-in is the API's
  and the shell decides none of it**: `XrefSet.for_symbols(species)` is the ordinary constructor with
  the question named, and `XrefSet.for_namespace(species, namespace)` is that same fill-in for a
  caller holding a namespace rather than a verb — which is what a caller who read one off a flag
  holds, so `genome xref` hands its `--from-stems` namespace straight over. One code path, so the
  shell and `import genome` cannot resolve the default differently. **A named source is never swapped
  and a built set never borrows**: `XrefSet("Homo sapiens")` still matches no symbol and raises,
  because it holds one publisher's bytes and answering from another's is the merge that one query
  reading exactly one set forbids (ADR-0017). **Both** symbol questions — matching a symbol, and
  labelling a stem with one — miss on such a set with that same message, and it names *this species'*
  source and `for_symbols`, rather than every source that publishes symbols for anybody, or, as the
  labelling direction did, nowhere at all.
- **The two directions are deliberately not mirror images.** Away from the hub,
  `from_stems(stems, "symbol")` gives the authority's **single current approved symbol** — this is
  labelling a figure axis, and it is one-to-one by the authority's own construction. Toward it,
  `to_stems(ids, "symbol")` **raises** naming `match_symbols` rather than answering on approved
  spellings alone, which is the failure above. The shell says the same thing in its own words:
  `genome xref --from-stems symbol` labels, `genome xref --to-stems symbol` exits `2` naming
  `genome match-symbols` — a shell user cannot act on the name of a Python call — and the
  `--to-stems` help no longer lists the namespace the command refuses, so help and behaviour agree.
  The kind of each spelling is stored in the set's own plain gzipped TSV — `symbol`,
  `previous_symbol`, `alias_symbol` in the namespace column — so a collaborator reads it in a shell
  with `awk`.
- **Two further Xref sources, and both pin.** `hgnc` is human's, from HGNC's **quarterly archive file
  of 2026-07-07**, and is the only one that publishes previous and alias spellings *typed*. Its
  reader **parses by header name and never by position**, because the schema has drifted from **52**
  columns in 2020 to **54** now, and a named test reverses the whole header to prove the answer does
  not move. Its pin **names a file read out of the bucket listing** rather than one built from a
  date: the archive's dates are irregular — `2024-07-02`, `2025-01-06`, and both `2026-07-03` and
  `2026-07-07` in one quarter. HGNC's remembered EBI FTP path is a live 404; the live one is a Google
  Cloud Storage bucket with a **doubled path segment**, `…/hgnc/archive/archive/quarterly/tsv/…`.
- **Mouse and worm get current symbols, and are told why they get no others.** `alliance_bgi` is the
  Alliance's per-species gene submission — MGI's for mouse, WormBase's own **WS298** for worm — read
  a record at a time out of 72 MB and 74 MB of JSON. It publishes the current approved symbol alone,
  because each record's `synonyms` list is **undifferentiated**: WormBase files the genuine former
  name `daf-17` beside the sequence names `R13H8.1` and `CELE_R13H8.1`, and typing them would put a
  kind on a claim no publisher made. So every answer carries `kinds` and `limits` saying which kinds
  that source could match and why the rest are missing — *this gene is absent* and *this source
  cannot match that spelling* must not both be silence. It is the Alliance's copy because neither
  authority can be reached directly: MGI keeps no dated archive (ADR-0018), and
  `downloads.wormbase.org` answers **403 to an automated client**, measured from two networks.
- **Ensembl is a second Xref source, selectable and pinned to a release of its own — and it is not
  the equal of the first.** `XrefSet("Homo sapiens", "ensembl", "116")` answers off Ensembl's
  per-species TSV dumps for human and mouse, pinned to release **116** — the release the lab's
  registered `gencode_v50` corresponds to — independently of Alliance's `9.0.0`, which stays the
  **Default xref source**. Adding it was what the design promised: a reader module, one entry in the
  reader table and two rows in the shipped metadata table. **Two publishers are two answers and
  nothing merges them** (ADR-0017): Entrez GeneID `79166` names **two** stems in Alliance 9.0.0 and
  **seventy-two** in Ensembl 116, each answer carrying the source and release that produced it, and
  the narrower one is never widened nor the wider one voted down. **The fan-out is stated where the
  source is chosen** rather than in a note — in the constructor's `source` parameter, in the
  reader's own module and in the shipped attribution — because it is the reason to choose
  deliberately: the two publishers agree on only **57.6%** of human gene-level (GeneID, ENSG) pairs,
  NCBI's mapping being a sequence match at a published overlap threshold and near-one-to-one where
  Ensembl's reaches 72 stems for one GeneID and **208 GeneIDs for one stem**. **The intuitive
  quality filter raises rather than answering nothing**: every human `EntrezGene` row release 116
  publishes carries `info_type=DEPENDENT` and **not one** carries `DIRECT` — 552,633 rows, zero
  direct, and mouse the same at 358,853 — so `evidence="DIRECT"` empties the set rather than
  narrowing it, and is met with an error naming what the release actually carries. The filter is a
  real capability and not only a guard: it selects which of the publisher's rows are read, a
  filtered set is prepared beside the unfiltered one rather than over it, and a source whose file
  grades nothing refuses a filter instead of ignoring it. Two conventions that cannot be assumed
  from one another are now recorded side by side: Alliance publishes an md5 of the **unpacked** TSV,
  Ensembl a BSD `sum` of the **served** `.gz` — 16 bits, and no integrity check for a 6 MB file — so
  the rows pin a digest of the unpacked bytes and the attribution records Ensembl's own value and
  what it covers. No worm row, and not by oversight: Ensembl files *C. elegans* under Ensembl
  Genomes' numbering, so release-116's worm directory holds a file stamped **63**, and worm is
  answered by the Alliance where the hop is the identity.
- **An `XrefSet` carries the curated row it actually resolved to**, as `provenance`, so what is
  cited is what answered rather than what a second lookup would resolve the defaults to today.
- **Which genes are transcription factors, answered by a published census that ships in the
  wheel.** `genome.tf.gene` is the gene half of the TF context, keyed by gene where the motif half
  is keyed by motif. Name an assembly and `Genome.tf_gene_list()` answers with the genes one census
  judges transcription factors, **in that annotation's own gene ids**, so the answer joins to a
  counts matrix with nothing left to normalise. Human is **Lambert et al. 2018** v_1.01 (2,765
  genes, 1,639 judged transcription factors, PMID 29425488) and mouse is **AnimalTFDB 4.0** (1,611
  genes, all positive, PMID 36268869) — not human TFs mapped through orthologs, which reaches at
  best 88.6% and loses 236 C2H2 zinc fingers (ADR-0014). **Nothing here decides what a TF is**: every
  verdict travels with the census that reached it, and the answer names the publisher, version and
  PubMed id to cite. **Absence is not emptiness** throughout — a gene the census never assessed and
  one it assessed and rejected are different answers, the 1,126 rejected genes ship so that "is this
  a TF?" can be answered *no*, and an assembly whose species has no census **raises and names the
  species that do** rather than returning an empty list. The species is read off the assembly's own
  metadata row and never passed in, so asking for human transcription factors while holding a mouse
  assembly is not expressible (ADR-0003). Thirteen of Lambert's twenty-nine published columns ship,
  gzipped, at 52 KB against the publisher's 326 KB; the four uniform across every census — the
  **Gene id stem**, the symbol, the TF flag and the **DBD family** — lead every table, and every
  other column keeps a snake_case spelling of its own publisher's name and is **never compared
  across tables**. The two family vocabularies are deliberately un-crosswalked: 75 values under
  Lambert's `DBD`, 72 under AnimalTFDB's `Family`, and `ARID/BRIGHT` and `ARID` are not asserted
  equivalent.
- **An annotation resolves a Gene id stem into the gene ids it actually spells.**
  `AnnotationRegistry.resolve_gene_ids()` is general and not TF-specific, and is the first thing in
  this package to open the annotation database every registration has always built — one indexed
  query, a row at a time, never the whole annotation in memory. **It answers with every gene id a
  stem names and never picks one**: measured against the registered annotations, `gencode_v50` and
  `gencode_vM39` collide for none and `gencode_v50lift37` names two gene ids for nine stems, eight
  of them pseudoautosomal-Y, where taking the first would be wrong nine times and silently. A stem
  carrying no version resolves to itself, so an annotation whose ids are not Ensembl-shaped is
  untouched by an Ensembl-shaped assumption. **Stems that resolve to nothing ride back on the
  answer** rather than being dropped, which is how the three entries Lambert records as UniProt
  entry names rather than Ensembl ids stay visible instead of vanishing.
- **Which JASPAR motifs answer for a TF gene, in a stated order.** `genome.tf.link` is the join
  neither half owns — it imports both and neither imports it — and `motif_links(gene, species)`
  hands back one gene's **Motif link**s, each saying what it is a motif *of*: `monomer` where the
  profile names one gene, `complex` otherwise with its partners recorded, so a heterodimer matrix
  is never read as a monomer's, and AHR, DDIT3, TAL1 and TLX1 are linked rather than reported
  motif-less. Links arrive in **Attribution specificity** order — monomer before complex,
  species-matched before **Cross-species link**, then higher total **Information content**, then
  **Motif id** — which is total and stable, so "the motif for this factor" means the same thing on
  two machines and in two releases. **It is an ordering and not a quality score**, and no quality
  score is computed or shipped anywhere here: JASPAR publishes none, matrix depth is normalised per
  assay so ranking on it ranks the assay, and the canonical AP-1 motif is a complex that describes
  JUN's binding better than any JUN monomer matrix. A caller who disagrees re-sorts on the
  attributes the link already carries. **Cross-species profiles are kept and marked, not excluded**
  (ADR-0013) — JASPAR files an orthologous pair's matrix under whichever species was assayed, so
  excluding them costs human 108 genes and mouse 552 of 689 — and `cross_species=False` drops them
  for a question that demands species-matched profiles. **Provenance is captured before any
  filtering**, since what comes out of a filter is no longer the release it came from. A versioned
  gene id is **refused rather than stemmed**, because a stem may name two gene ids and this package
  never picks a gene the caller did not name.
- **The TF-to-motif mapping ships as a curated table rather than a computed rule** (ADR-0015). One
  gzipped TSV per species per **Release** under `data/tf_link/` — human and mouse × JASPAR 2024 and
  2026, CORE vertebrates, 70 KB in all — readable in R or a shell **without importing this
  package**, because a TF-to-motif mapping is wanted outside here as often as inside. Bulk ships
  gzipped and small metadata ships plain, as everywhere else in this package's data, so the
  three-row alias table beside them stays a readable diff. Twelve columns carrying the release and
  species on every row so tables concatenate. The link is made by upper-casing the **Motif name**,
  splitting on `::` and matching each part against the census's own symbol column, with an alias
  table keyed on **gene id** and not on symbol for the profiles JASPAR renamed after Lambert was
  published — `TBXT` for `T`, `SCAND3` for `ZBED9`, `ZFTA` for `C11orf95`. The `EWSR1-FLI1` fusion
  names no gene and stays unlinked by design. Only assessed-positive genes receive links: human 876
  genes / 1,085 links on 2026 and 745 / 946 on 2024, mouse 693 / 896 and 653 / 851, with 732 of
  mouse's 896 cross-species — which is the coverage argument made concrete. **No test proves a
  shipped table still matches JASPAR and none can**, since regenerating one needs a download CI
  cannot make; the pinned counts convert a silent drift into a loud failure, which is the most that
  is available, and it is written down rather than papered over.
- **`genome tf-gene-list <assembly>` prints an assembly's TF genes**, shaped exactly like
  `gene-list` because a caller who learned one has learned the other. Gene ids to stdout one per
  line, the heading and the census attribution to stderr so the output pipes; `--annotation` names
  which registered annotation to ask; `--json` emits the whole record — the genes with their **TF
  assessment** and **DBD family**, the census's provenance, and the unresolved **Gene id stem**s.
  Non-zero exit with a message naming the next action for each of three distinct failures: the
  annotation is not registered, no census ships for this assembly's species, and nothing says what
  species the assembly is. There is deliberately **no `--include-rejected` flag** — stdout is a bare
  id list with nowhere to carry which verdict an id holds, so the flag would feed assessed-negative
  ids into a pipeline that reads them as TFs; the widened answer is expressible in Python, where
  each id travels with its flag. And **no command wraps the Motif link table**, which is already a
  file anyone can read — the reason for shipping it.
- **The censuses and the link tables are rebuilt by committed generators**, `scripts/build_tf_census.py`
  and `scripts/build_tf_links.py`, so a publisher's re-release or a new JASPAR release is a re-run
  and a reviewable diff rather than a manual edit. They **fail loudly and name the column** when a
  publisher re-spells one rather than dropping it, write **byte-stable** output so an unchanged input
  produces no diff, and live outside the shipped package — the wheel carries the data and not the
  build tooling. The link generator reads JASPAR's SQLite dump for per-profile species, which is
  why **the motif subpackage needed no change at all**: no new field on **Motif**, no change to the
  loader, parser or scan path.
- **The human cofactor list, and it is this package that publishes it.** `cofactor_table("Homo
  sapiens")` answers with **1,466 genes — 354 both publishers list, 670 AnimalTFDB 4.0 alone, 442
  EpiFactors v2.0 alone** — a union neither publisher releases and therefore, uniquely in
  `genome.tf`, **nobody's verdict but ours** (ADR-0016). Everywhere else here a verdict travels with
  the census that reached it, and classification still does: the row carries AnimalTFDB's family and
  category *and* EpiFactors' function, target, modification and complex under namespaced columns,
  filled only by the publisher that actually named the gene, with **nothing crosswalked between the
  two in either direction** (ADR-0014) — so a `source` of `both` is agreement about membership and
  about nothing else. EpiFactors keys on HGNC ids and publishes no Ensembl ids at all, so the stems
  of the 442 genes only it lists come from **a pinned dated HGNC monthly archive** and never the
  rolling current file, which is what makes them reproducible; the join is on the id and **never on
  the symbol**, because 31 of EpiFactors' 801 rows still name their gene by a symbol HGNC has
  retired — `ACINU` for `ACIN1`, `ARNTL` for `BMAL1` — and human's `symbol` column is HGNC's
  approved spelling throughout. Five genes carry two EpiFactors rows each and ship as one with their
  cells unioned and deduplicated, at a cost stated rather than hidden: for those five the pairing
  between a function and its own modification is lost. **The TF list and the cofactor list overlap —
  151 human genes are both a Lambert-positive TF gene and a cofactor**, 57 from the AnimalTFDB side
  and 122 from the EpiFactors side, so a caller who unions the two answers double-counts them; being
  a cofactor never suppresses a motif the census already reached. Human's provenance carries **three**
  source rows, HGNC's among them, because a source that earns 442 stems earns a citation. Multi-valued
  cells split on `;`, the separator the rest of this package already uses, and the build refuses a
  published value that already contains one. Every curation rule lives in
  `scripts/build_tf_cofactor.py` and none of it in the wheel.
- **Which genes a publisher lists as transcription cofactors, answered by a table that ships in the
  wheel.** `genome.tf.cofactor` is the third part of the TF context and a peer of the gene and motif
  halves: the gene half answers whether a gene is a **TF gene** and of what **DBD family**, this one
  answers whether it is a **Transcription cofactor** and of what class. It is keyed the same way, by
  **Gene id stem**. Mouse and worm ship **AnimalTFDB 4.0** (PMID 36268869) — 970 genes across 84 of
  the publisher's own families, *C. elegans* 317 across 57, and the same six categories in each — as
  gzipped TSVs under `data/tf_cofactor/`, found by enumerating that directory so that adding a
  species is dropping in a file. **Nothing here decides what a cofactor is** for those two:
  membership and classification both travel with the publisher that reached them, and the answer
  names the publisher, version and PubMed id to cite. Four uniform columns lead every table —
  the stem, the symbol, the cofactor flag and a closed-vocabulary `source` validated on read — and
  everything after them is one publisher's own column under a namespaced name, never compared with
  another's (ADR-0014). `is_cofactor` reads `yes` on every row today and is kept anyway: dropping it
  would make presence in the file the verdict, and a future source could then not record a rejection
  without a format change. **Worm ships although no publisher has censused worm transcription
  factors**, so a worm assembly answers here while the TF gene half has nothing to say — the
  publishers' shape and not a defect, stated beside the data and pinned in the tests. Provenance is
  **two** plain tables rather than one, keyed by species and by species-and-source, because one row
  cannot describe a table built from several publishers and joining them positionally inside a cell
  is the shape that breaks quietly. `parse_cofactor_table` is public and takes the table as text, so
  every way a shipped file can be malformed is reachable without writing a broken one into the
  package, and every refusal names the file and the command that regenerates it.
  `scripts/build_tf_cofactor.py` is the third committed generator: it takes file paths and downloads
  nothing, joins AnimalTFDB's own two files through five hand-written family spellings **whose
  arithmetic it re-runs on every build**, and fails loudly both when a family survives that map with
  no category and when the publisher's own counts stop reconciling — so a release that renames a
  family is a broken build rather than a quietly blanked column.
- **An assembly's transcription cofactors, in its own annotation's gene ids.**
  `Genome.tf_cofactor_list()` is the **TF cofactor list**, the counterpart of `tf_gene_list()` and
  the same three layers: a method on the genome, one on the annotation registry taking a **Registered
  name**, and a module-level `tf_cofactor_list(assembly)` for a caller who has not opened one — one
  code path, so a shell surface over it adds no second. Every **Gene id stem** the **Cofactor
  table** is keyed by is resolved through the annotation's own gene ids, so the answer joins to a
  counts matrix with nothing left to normalise; a stem naming two gene ids answers with **both** and
  never picks one; and the stems this annotation carries no gene for **ride back on the answer**
  rather than being dropped. The species is read off the assembly's own metadata row and never
  passed in, so asking for human cofactors while holding a mouse assembly is not expressible
  (ADR-0003). Each entry carries the four uniform columns as fields of its own — the stem, the
  symbol, the cofactor flag and which publisher listed the gene — with every publisher's own
  vocabulary reachable beside them under that publisher's namespaced column name, exactly as a **TF
  gene** carries its **DBD family** and the census's other judgements. **Absence is not emptiness
  here either:** an assembly whose species has no cofactor table raises the new
  `NoCofactorTableError` **naming the species that do**, rather than answering with none, so that
  *nobody has published one for this species* can never be read as *this species has no cofactors*.
  `UnknownSpeciesError` is now shared by both halves and says which of them was asked for. **The two
  halves do not raise for the same assemblies**: a worm assembly answers here and raises from
  `tf_gene_list()`, because a publisher assessed worm cofactors and none has released a worm TF
  census.
- **`genome tf-cofactor-list <assembly>` prints an assembly's transcription cofactors**, shaped
  exactly like `genome tf-gene-list` because a caller who learned one has learned the other. Gene
  ids to stdout one per line, the heading and the publishers' attribution to stderr so the output
  pipes; `--annotation` names which registered annotation to ask; `--json` emits the whole record —
  every gene with the publisher that listed it and that publisher's own classification, one
  provenance entry per publisher to cite, and the unresolved **Gene id stem**s. Non-zero exit with a
  message naming the next action for each of three distinct failures: the annotation is not
  registered, no cofactor table ships for this assembly's species — **and the message names the
  species that do** — and nothing says what species the assembly is. The command computes nothing:
  `tf_cofactor_list()` is the one code path, so the shell and a notebook cannot drift. **A worm
  assembly is answered here and refused by `genome tf-gene-list`**, pinned in a test of its own so
  that the asymmetry is not mistaken for a bug and quietly "fixed".
- **Asking a transcription cofactor for its motifs now answers that it is one.** `motif_links()`
  met every gene no census assessed with the same `GeneNotAssessedError` — *this census never
  assessed this gene* — which reads as *nothing here knows this gene* when a publisher does list it,
  as a **Transcription cofactor** that recognises no sequence of its own and so has no motif to look
  for. The lookup now has an **order**, and the order is what keeps it correct. The census is asked
  first, so the **151 human genes that are both a TF gene and a cofactor** — TBP, KMT2A and DNMT1
  among them — come back with exactly the links they always did, because a second table must never
  suppress an answer the census already reached; a gene it assessed and turned down still answers
  with that verdict. Only then is that species' **Cofactor table** asked, and a gene it lists raises
  the new `TranscriptionCofactorError`, whose message names the census that did not assess the gene,
  the publisher that lists it as a cofactor and that there is no motif here to look for. A gene
  neither knows raises `GeneNotAssessedError` unchanged, and so does every gene of a species that
  ships no cofactor table. **The new error subclasses the old one**, so an `except` clause written
  before any cofactor table shipped keeps covering every gene it covered — a sibling under
  `LookupError` would read more cleanly and would silently stop covering exactly the genes this
  added knowledge about. The subclass relation is literally true, no TF census having assessed the
  gene, and claims nothing about biology: a cofactor is not a kind of transcription factor.
- **The words for transcription cofactors, settled before the code that will use them.** The TF gene
  context glossary now covers the whole TF context at `docs/context/tf.md`, and defines
  **Transcription cofactor**, **Cofactor table** and **TF cofactor list** — so that a bare
  "cofactor", which names NAD+ and heme to most of biology, never stands on its own in an issue, a
  test name or an error message, and so that a table saying two publishers listed a gene is read as
  agreement on membership only and never on classification. **`genome.tf` no longer says cofactors
  are out of scope**: they are a carve-out with a subpackage named for them, `genome.tf.cofactor`,
  and the reason a cofactor has no motif is that it recognises no sequence — not that nothing here
  knows about it. ADR-0016 records that **this package publishes the human cofactor list**, the
  first place it decides anything rather than relaying a publisher's verdict: 1,466 genes unioned
  from AnimalTFDB 4.0 and EpiFactors v2.0, whose stems come through a pinned dated HGNC archive,
  with both costs stated rather than smoothed over — 151 human genes are both a **TF gene** and a
  cofactor, and the two classification vocabularies are deliberately not crosswalked (ADR-0014).
  Vocabulary and records only; no behaviour changes.
- **The words for cross-database gene identifiers and for cross-species homology, settled before the
  code that will use them.** Two new bounded contexts get glossaries — `docs/context/xref.md` and
  `docs/context/orthology.md` — so that an **Xref set** is never called a crosswalk, an id map or a
  mapping table in an issue, a test name or an error message, and so that a **Homology link**'s
  **Homology type** is read as its publisher's claim about evolution rather than as a count of rows
  in whatever file came back. Xref coins **Xref set**, **Xref source**, **Default xref source**,
  **Namespace** and **Symbol match**; orthology coins **Homology set**, **Homology link**,
  **Homology type**, **Dropped partner** and **Paralogy link**. The context map lists both and draws
  the seven edges joining them to Annotation, TF and Motif — two of which are **prohibitions rather
  than calls**: nothing shipped here is keyed by a foreign **Namespace**, and no list this package
  publishes is derived through homology. **Four records settle what the design decided.** ADR-0017:
  the hub is the **Gene id stem**, a query reads exactly one **Xref set**, and nothing composes two
  hops or merges two publishers — NCBI and Ensembl agree on only 57.6% of human gene-level (GeneID,
  ENSG) pairs, so a merged table would decide nearly half its rows by a rule nobody published.
  ADR-0018: only a publisher keeping dated releases at stable URLs is eligible, and its bytes are
  downloaded rather than shipped, at two costs stated rather than smoothed over — neither set
  answers offline on a fresh install, and mouse gets current symbols only, previous and alias
  spellings being MGI's and MGI publishing no dated archive. ADR-0019: **orthology is served and
  never consumed**, so worm still has no TF census and mouse still has no assessed-negative genes
  even with a homology table on disk. ADR-0020: a **Homology type** is the publisher's tree-derived
  label and is never recomputed after filtering. ADR-0014 gains a `**Status.**` line, its cost line
  having said this package builds no ortholog or homology support at all, and the rule against
  assuming an assembly, a coordinate system or a strand is extended in place — a gene id never
  crosses species implicitly. **The shared kernel gains an eleventh word**: **Gene id stem** moves
  out of the TF glossary, four contexts now keying on the identifier every **Xref set** hangs off,
  and its definition is unchanged. Vocabulary and records only; no behaviour changes.
- **A column of Entrez, UniProt, HGNC, MGI or WormBase ids reaches this package's answers, with no
  assembly and no genome open.** `genome.xref` is a peer of `genome.tf`, and an **Xref set** is one
  species, one **Xref source** and one pinned **Release**: `XrefSet("Homo sapiens")` fetches the
  publisher's file once into `$LIULAB_DATA/xref/`, a sibling of the assembly tree beside `motif/`,
  slices it to that species and re-reads it thereafter — the shape `JasparDatabase` established,
  with a `cache_dir` override naming the directory itself. **Two verbs and only two** (ADR-0017):
  `to_stems` toward the hub and `from_stems` away from it, so a caller wanting Entrez → HGNC makes
  both calls and owns the join, and a query reads exactly one set, which makes merging two
  publishers inexpressible rather than merely discouraged. Both answer in the shape
  `resolve_gene_ids` already established — ask order kept, **every** id a key names and never a
  chosen one, no resolved value ever empty, what named nothing riding back in `unresolved`, a
  flattener that documents what flattening loses, and `as_json()` — and every answer names the
  species, source and release that produced it, because an answer that did not would be
  unreproducible a year later. **Alliance of Genome Resources 9.0.0 is the first source and the
  default for all three species**, gene level only, across Ensembl, Entrez, UniProt and each
  species' own authority; naming a release is enough to fetch it, the curated row in
  `data/xref/xref_metadata.tsv` knowing the publisher, the version, the URL and the checksum. **The
  three species' hops turn out to have three different shapes**, which is the argument for this
  being an object a caller opens: for worm it is the identity — all 46,926 `WB:WBGene…` genes carry
  `ENSEMBL:WBGene…`, the same string, with zero differing — for mouse a real join onto `ENSMUSG…`,
  and for human a join through HGNC in which **2,535 of 40,665 genes carry two or more Ensembl
  cross-references**, so 6.2% of HGNC ids name two stems and nothing picks one. **Every incoming id
  is reduced to its Gene id stem on ingest**, a publisher's and a caller's alike, because joining a
  versioned id to a bare one returns zero matches and says nothing — the most error-prone detail in
  this landscape — and each **Namespace**'s CURIE prefix is accepted whether or not it is written,
  so `HGNC:11998` and `11998` are one identifier. Alliance's duplication is deduplicated **on the
  key and never on the row**: 2,659,704 rows reduce to 1,811,267 distinct
  `(GeneID, GlobalCrossReferenceID, TaxonID)`, 31.9% redundant, and a whole-row `uniq` removes none
  of it. **The stored form is a plain gzipped TSV** of `namespace`, `xref_id`, `gene_id_stem`,
  sorted, unique and written with no gzip timestamp, so two machines slicing one release produce
  identical bytes and a collaborator who does not use Python reads it in R. The **Completion
  marker** beside it carries **both** checksums — the publisher's own md5 as provenance and the
  slice's own sha256 as the integrity check, since what is stored is a derived slice rather than the
  publisher's bytes — and either one disagreeing with what is on disk means unfinished rather than
  present. The publisher's md5 is over its **unpacked** bytes (ADR-0006), which is Alliance's own
  convention and a trap: hashing the `.tsv.gz` as it arrives mismatches every time. Four errors name
  their next action: a set that is not downloaded names the call to make on a login node, an
  unsupported species names the three that have a set, a **Namespace** the source does not carry
  names the ones it does, and a file that does not match its pin refuses rather than answering with
  silently fewer genes.
- **`genome xref` — a species, a set of ids and a direction, and the answer comes back with the
  misses still on it.** The shell surface over an **Xref set**, and a thin one: it parses arguments,
  makes one API call and renders, so `import genome` and the shell hit one code path and `--json`
  is `as_json()` verbatim. **The direction is named and never inferred** — `--to-stems NAMESPACE`
  reads the ids as that namespace and answers in **Gene id stem**s, `--from-stems NAMESPACE` reads
  them as stems and answers in that namespace — and each flag carries the namespace, so a direction
  without one is not expressible and naming neither or both exits `2`. Inferring it from the id
  strings would be a judgement the API does not make: `HGNC:11998` asked the wrong way answers
  *nothing found* rather than quietly turning around. **The pairs go to stdout, tab-separated, so
  the output pipes** — `cut -f2` is the answer and `cut -f1` says what asked for it — with the
  heading, the publisher's URL and the counts on stderr; an id naming two genes prints two rows
  rather than whichever came first, and **an id that resolved to nothing gets a row of its own with
  an empty second column**, which is exactly what a hand-rolled join drops without saying so. Every
  id passed therefore leaves with at least one row, in both renderings. `--source` names an **Xref
  source** and omitting it answers from the species' **Default xref source**; either way the answer
  names the source and release that produced it. **Every failure exits non-zero naming the next
  action**: a species no set exists for names the three that do, a set that is not here and cannot
  be fetched names the call to make on a login node, a **Namespace** the source does not carry names
  the ones it does, and a directory an interrupted download left unfinished names the repair.
- **`genome homologs` — a species pair, a set of stems, and the publisher's own label on every
  row.** The shell surface over a **Homology set**, and as thin as its sibling: it parses
  arguments, makes one API call and renders, so `import genome` and the shell hit one code path.
  Any pairing among human, mouse and worm answers, either way round. **The links go to stdout,
  tab-separated, so the output pipes** — the seven columns are the keys `as_json()` writes, in its
  order, so the text rendering and `--json` cannot drift and every cell printed is a value the API
  put in the answer. **Every stem passed leaves with at least one row**: one with three homologs
  prints three, and one this release names no homolog for gets a row with the other columns empty
   — which is not `NULL`, Compara's own word for a cell it recorded nothing in on a link that does
  exist. The heading, the attribution and two qualifications go to stderr: the **Dropped
  partner**s, counted *and* named so a link that merely looks one-to-one stays distinguishable
  from one the publisher called one-to-one, and whichever quality columns the set holds no value
  in anywhere — `goc_score` and `wga_coverage` are null on every link of *either* worm pairing, so
  a shell user is told before `awk` empties their filter rather than after. `--paralogs` returns
  every link the publisher wrote and a **Paralogy link** is marked by its own `homology_type`
  rather than excluded (ADR-0013); release 116 publishes none cross-species, so on it the flag
  changes nothing and the heading says which question was asked either way. **Every failure exits
  non-zero naming the next action**, the wrong-file case most of all: a pair fetched from the
  Compara dump that no longer holds it raises naming the other file rather than answering empty.
  A set that cannot be fetched now names the call to make on a login node, as an **Xref set**
  already did.

- **`genome motif-scan` — a FASTA in, a Parquet file out, a summary on standard output.** The
  batch case, and the one motif operation that belongs in a shell script and a scheduler job;
  listing, plotting and comparing motifs get no command, because they are notebook work. It takes
  the FASTA and the output path, plus `--release`, `--tax-group`, `--threshold`, `--background`
  and `--workers`, and prints what the run was — release, tax group, motifs scanned, motifs
  skipped, the background actually used, the threshold, sequences scanned, hits written, workers,
  and where the hits went. **The hits go to the named file and the summary to standard output**,
  so `--json` is never corrupted by table data, and any progress display is suppressed under it.
  **It defaults to every core the allocation granted** — the Slurm allocation first, then process
  affinity, then the machine — where the library defaults to one worker: a console script is a
  proper entry point, so the process-pool hazard behind that default does not apply here.
  `--background` takes the three modes; four frequencies of your own stay a Python call, since a
  mistyped one on a command line would change every cutoff and look like a scan that simply found
  fewer hits. The summary is read off the Parquet footer and never by reading the hits back, which
  is what `provenance_of(path)` and `hit_count(path)` are for — what a written scan was, and how
  many hits it holds, at the same cost on 550 million rows as on none.
- **Scanning regions takes the same arguments every other scan does.**
  `Genome.scan_regions` now forwards `background=` and `workers=` to the scan underneath it,
  so the case an HPC user hits first — a background derived from the peak set's own composition,
  scanned across the whole allocation — is reachable from a region scan and not only from a
  sequence one. Two workers produce the **identical table** a serial scan does after the lift as
  well as before it. `output=` is the one argument that is *not* forwarded and is now refused by
  name: a scan that streams to Parquet hands back a path, and a path holds no coordinates to lift
  into the assembly's frame, so the refusal names both ways to get one written instead.
- **A scan can use more than one core, and says the same thing when it does.** `workers=` on
  `MotifSet.scan`, `scan_sequences` and `scan_fasta` shards the work across processes — MOODS holds
  the GIL, so parallelism here means processes and not threads, each worker constructing its own
  engine once and keeping it. **Serial and parallel produce the identical table**, row for row,
  dtypes and provenance included: the shards of one sequence are put back into the order a serial
  scan would have emitted them before the batch is handed on, so choosing two workers is a choice
  about wall time and about nothing else. **The library default is one worker**, so importing this
  package and scanning never starts a process unasked — under the spawn start method a pool
  re-imports the caller's script, and an unguarded script would re-execute itself; the command-line
  entry point will pass `None` instead, which resolves the count with `resolve_workers()`: **the
  Slurm allocation first** (`SLURM_CPUS_PER_TASK`, then `SLURM_CPUS_ON_NODE`), then process
  affinity, then the machine's cores — never the last alone, which would put fourteen workers into
  a two-CPU job. A sequence long enough to be worth cutting is split with an **overlap of one less
  than the longest matrix**, and each shard keeps only hits *starting* inside the region it owns, so
  a hit lying across a boundary is reported exactly once. Shards are submitted a bounded number
  ahead rather than all at once, so a genome FASTA is still streamed rather than cut up in advance.
  The engine is also now built once per scan rather than once per sequence, which is the same
  automaton and the same answer, built a thousand times less often on a thousand-record FASTA.
- **A scan too large to hold streams to Parquet and hands back the path.** Passing `output=` to
  `MotifSet.scan`, `scan_sequences` or `scan_fasta` writes the hits to Parquet and returns the path
  rather than a table. Batches are written **as they are produced** — one row group per named
  sequence — so the whole result is never materialised: hg38 against a full vertebrate release is
  about 550 million rows, which at the 19 bytes a row the fixed dtypes cost is 10.5 GB and is not a
  DataFrame on any machine in the lab. **There is no row-count guard and no refusal**; a
  genome-scale scan is the caller's decision. `read_hits(path)` reads one back and is the reader to
  use: what comes off the disk equals the in-memory table for the same scan, **dtypes included**,
  down to the order of a categorical column's categories — the writer pins `int32` dictionary
  indices so batches of four sequence names and of four hundred agree on one schema, and the reader
  sorts the categories back into the order `astype("category")` produces. **The provenance travels
  in the file**, under its own key in the Parquet metadata, because `frame.attrs` does not survive
  pandas' round trip: `pandas.read_parquet` gives the rows and drops what the scan was, and
  `read_hits` puts the background, threshold, release, tax group and motif lists back on `attrs`.
  This adds `pyarrow` to the core dependency table.
- **A scan derives its background from what it was handed, and stops converting the same
  thresholds twice.** The background is now **automatic**: derived from the input when the input
  holds at least 10 000 unambiguous bases, uniform below that, since a composition estimated from
  fewer would distort the very cutoffs it sets — at that floor the standard error on each base
  frequency is about 0.004, which moves a per-position log-odds term by under 0.02 nats. This
  matters more than any other scan parameter: switching from uniform to a real chromosome's
  composition changed the hit count by 2.5% but turned over **26%** of the hit set. `background=`
  still takes four frequencies and now also takes `"uniform"` to pin it, `"derive"` to derive
  whatever the input holds, and `"auto"`, which is what omitting it means; **whichever it was, the
  background actually used is recorded on the result**, so handing the recorded value back
  reproduces the scan exactly. Deriving reads a **bounded prefix** of the input and hands those
  records back in front of the rest, so a FASTA is still read once and a generator source is still
  drained once — the estimate is a head sample rather than the whole input, which is exactly the
  accuracy the floor was chosen for. Every background, however it arrived, is **rounded onto a
  0.001 grid** and that one value builds the matrices, sets the cutoffs, keys the cache and is what
  gets recorded. Converting the threshold into one cutoff per motif is the engine's one slow step —
  a few seconds for a full vertebrate release — and a pure function of `(matrices, background, p)`,
  so it is **cached on disk** under `<LIULAB_DATA>/motif/thresholds/`, shared by every project on
  the machine exactly as the JASPAR files are; the rounding is what lets two peak sets from one
  genome share an entry. A cache is a speed-up and never a dependency: a corrupt, truncated or
  older entry is a miss, and a data root that cannot be written to makes a scan slow rather than
  broken.
- **Scan regions of a genome and get hits in that assembly's own coordinates.**
  `Genome.scan_regions(motifs, regions)` fetches each region's bases, scans them with a motif set
  and lifts every hit into the assembly's frame, so the translation that used to be redone in every
  notebook — and is where the off-by-ones live — exists in exactly one place. `sequence_name`
  carries the chromosome as *this* assembly spells it, and a region naming one the assembly does not
  carry raises rather than being reconciled. **The arithmetic is written out in the docstring so a
  reader can check it without running it**: for a `+` region (and a `.` one with it) a hit found at
  local `[s, e)` is at `[S + s, S + e)` and keeps its strand; for a `-` region, whose bases are
  fetched reverse-complemented, it is at `[E - e, E - s)` with the hit strand flipped — the two ends
  swap, and the `- 1`s cancel exactly because coordinates are 0-based half-open, which is why the
  same flip written for a 1-based inclusive interval is wrong. An unknown strand is **not** promoted
  to `+`: the fetch hands back forward bases for it, so there is nothing to flip. Regions may
  overlap and several may sit on one chromosome — each is scanned in its own right, and a hit seen
  from two regions is two rows, since deduplicating would be deciding for the caller which peak a
  site belongs to. A locus *string* is refused, because it carries no strand and the strand is the
  whole question. It arrives as a mixin on `Genome`, following the aligner's, and **the dependency
  runs Genome to motif and never back** — the motif modules name `Genome` under type checking alone,
  since a motif belongs to no assembly and a motif set is usable with no genome open. **The raw form
  is untouched**: scanning a mapping of sequences still answers in region-local coordinates and
  names no assembly. What comes back is the same hit table with the same columns and dtypes, and its
  provenance is **extended rather than replaced** — the assembly joins the background, threshold,
  release, tax group and the two motif lists on `frame.attrs`, because a chromosome name and an
  interval mean nothing until something says which reference they are in (ADR-0003).
- **Scan DNA with a motif set and get one hit table back, whatever you scanned.**
  `MotifSet.scan(sequence, name="sequence")`, `MotifSet.scan_sequences({name: bases})` and
  `MotifSet.scan_fasta(path)` answer with the same table — the same columns, in the same order,
  with the same compact dtypes (`motif_id`, `motif_name`, `sequence_name` and `strand` categorical,
  `start` and `end` `int32`, `score` `float16`) — so nothing downstream branches on how the scan was
  called. **The dtypes are the contract and not an optimisation**: 19 bytes a row against about 100
  for what pandas would infer. Coordinates are **0-based half-open and always in the forward frame**,
  both strands are scanned, and the strand is `+` or `-` and never unknown, because a scan knows
  which of the two it scored. The engine is MOODS, which has no strand concept: the adapter doubles
  the matrix list with reverse complements, computes a cutoff per entry, and splits the results back
  by index — and a position reported for a reverse-complement matrix is already a forward-frame
  start, which the suite asserts against the engine rather than trusting. **The threshold is a
  per-position p-value**, 1e-4 by default, converted per motif against the background; the `score`
  column carries log-odds **in bits** (MOODS scores in nats and the one conversion lives in the
  adapter) and no per-hit p-value is computed. **Motifs shorter than 7 positions are not scanned and
  are named on the result** rather than called at a looser cutoff than was asked for: a 6-mer has
  4096 possible words, so its best match has p = 2.44e-4, and an engine asked for 1e-4 anyway clamps
  and over-calls in silence on roughly 98 of 879 vertebrate motifs. **Input is upper-cased, so a
  soft-masked sequence yields exactly the hits its upper-case equivalent does**, with no option to
  honour the masking — the one place a scan contradicts a shared-kernel term, recorded as ADR-0012.
  A FASTA record's name is its header **up to the first whitespace**, matching what STAR and chromap
  write into an alignment made from the same file, so a hit table joins against that alignment with
  nobody renaming anything; plain or gzipped, and records are read one at a time. The table carries
  its provenance on `frame.attrs` — the background, the threshold, the release, the tax group, and
  which motifs were scanned and which skipped — read off the set, so a `JasparDatabase` names its
  release and a filtered one answers `None` for both, since a filtered release is no longer that
  release. This release is serial with a uniform background; the scan produces batches internally,
  one per named sequence, so a Parquet sink and a parallel source attach without restructuring it.
- **Name a JASPAR release and query it like a dictionary that filters itself.**
  `JasparDatabase("2024", "vertebrates")` fetches that release's transfac file once into
  `<LIULAB_DATA>/motif/jaspar/`, flat and with the release and tax group in the name, and every
  construction after re-reads it and fetches nothing. Both releases (2024, 2026) and all eight tax
  groups are supported, vertebrates by default, and `"all"` selects the union file — whose
  published name drops the tax group, which is why the cached name is built rather than copied.
  Motif data is the **first thing filed beside the assembly tree rather than inside it**, since a
  motif belongs to no assembly. The transfac serialization is read rather than `.jaspar`: it
  carries the count matrix and all six annotations in one file and does not round counts to
  integers. **Annotation values are separated by a semicolon and never by a comma** — commas live
  inside single values (`C3H(C),C2HC zinc-fingers like factors`, `PBM, CSA and/or DIP-chip`), so
  splitting on one would silently corrupt about fifty records per release. A `MotifSet` holds the
  motifs — built from *any* motifs, so a model's de novo matrices get the same API — and indexing
  it resolves a matrix id, a bare base id or a unique factor name to **exactly one** motif, never a
  union type. An ambiguous name raises `AmbiguousMotifNameError` naming every matching id and the
  call that returns them all, since 66 names collide in 2024 and 71 in 2026; `by_name` always hands
  back a tuple, of one where the name is unique, and absence raises rather than answering with
  nothing. Filtering takes annotation keywords (prose matched as a case-insensitive substring, ids
  matched exactly) or an arbitrary predicate and returns a plain `MotifSet`, because **a filtered
  release is no longer that release**. There is **no completion record here and that is
  deliberate** — the files are under 1 MB, so integrity is a download to a temporary name renamed
  into place only on success, plus a parsed-motif count checked against a constant on every read,
  which turns a truncated file into an error rather than half a release. A non-redundant release
  shipping one version of each matrix is asserted rather than assumed, because it is what makes a
  bare base id address one motif.
- **`Motif.tf_class` and `Motif.tf_family` are tuples, not strings.** They are genuinely
  multi-valued — a dimer carries one of each per half — so they join `uniprot_ids` and `pubmed_ids`
  as plural annotations, and a bare string handed to any of the four is refused rather than stored
  letter by letter. `tax_group` and `data_type` stay single strings.
- **A motif, built from counts and able to answer for itself.** `genome.tf.motif.Motif` is a frozen
  4 × L count matrix plus the id and name it is addressed by and the six annotations JASPAR
  publishes — tax group, class, family, UniProt accessions, PubMed ids, data type — and it needs no
  file, no download and no network to be built or asked anything. **The counts are the single
  source of truth**: probabilities are a column normalisation, and log-odds take a background and a
  pseudocount as *arguments and never fields*, so one motif scored against two backgrounds stays
  one motif and two motifs sharing an id can never differ. Counts are float rather than int,
  because the source's own records carry fractional values. It reports its information content per
  position in bits, its consensus as a typed `DNA` rather than a `str`, and its sequence logo
  through `plot()`, which takes an optional axes and returns the axes so a grid of motifs is one
  figure — with the y-axis in bits, the same quantity trimming thresholds on. **Trimming acts only
  on the ends**, so a degenerate spacer in the middle of a dimer can never split it in two; it
  honours an optional maximum length, never goes below the minimum scannable length of 7 — a 6-mer
  has 4096 possible words, so its best match has p = 2.44e-4 and cannot reach the default 1e-4
  threshold — and returns a motif with the same id, the same name and an offset such that a
  position in the trimmed frame plus the offset is the position in the full one. Trimming a trimmed
  motif composes. The matrix is copied and marked read-only at construction, so a frozen motif is
  frozen all the way down.
- **The motif libraries, in every environment.** MOODS, logomaker, matplotlib and xarray join the
  core conda dependency table and `memelite` the PyPI one, ahead of the motif subpackage that will
  import them. None gets its own feature or environment: what earns a feature is a binary the
  package does not ship and checks for at call time, which a library is not, and the GTF library
  already sits in the core table on the same reasoning — so **the suite stays two lanes**. It is
  `matplotlib-base` rather than `matplotlib`, which is what logomaker itself depends on and drops a
  GUI stack nothing here opens. None of them is added to the PEP 621 `dependencies`, following
  `gffutils`: that list is what `pip install liulab-genome` can actually satisfy, and MOODS
  publishes no wheel — only an sdist needing SWIG and a C++ toolchain — so declaring it there would
  turn a working install into a failing build. **The scan engine was chosen by measurement**, not
  preference: `memelite`'s `fimo` silently returns zero hits for 97 of the 879 JASPAR 2024
  vertebrate motifs, never scans the final window of a sequence, and cannot be driven without
  holding the whole file; MOODS is also 9.7× faster per core once thresholds are prepared. The
  numbers and the method are in `docs/research/motif-scan-engine-2026-08-28.md`.
- **An annotation can name the genes in a category.** `Genome.gene_list("rRNA")` returns the gene
  ids, and `Genome.gene_lists()` every category that annotation declares — from the CLI, `genome
  gene-list <assembly> <category>` and `genome gene-categories <assembly>`, both with `--json`.
  The answer comes from a curated gene list shipped inside the package, one per annotation, and not
  from the GTF's own biotype attribute, which four publishers spell two ways over three taxonomies
  that disagree and which `sacCer3/ensgene_v101` omits altogether (ADR-0011). Every annotation
  declares one category today, `rRNA`, holding **everything rRNA-derived it carries** — mature
  genes, pseudogene copies, mitochondrial rRNAs, and yeast's 35S precursor. It is drawn for
  counting rRNA-derived reads as a QC metric rather than for describing rRNA biology, so it is
  inclusive and its genes may overlap; each list says in its own prose what it holds and what it
  under-reports. A merged annotation answers with one source per contributing component, so a
  chimera's genes stay attributable. **An annotation that cannot answer raises rather than
  returning nothing**: `NoGeneCategoriesError` for one no list ships for and
  `GeneCategoryNotDeclaredError` for one whose list declares other categories, both `LookupError`s
  so they can be caught together and still told apart. No declared category is ever empty and no
  call ever hands back an empty collection. Lists ship for all seven registered annotations.
- **An intron bound on the assembly table.** A row may now carry `intron_length_cap` — the longest
  gap a spliced aligner should take for an intron on that assembly — and
  `intron_length_cap_rationale`, which says why that number and not another. Registered for `ce11`
  (50,000), `ecHT115` (1), `hg38` and `mm39` (1,000,000 each), and blank everywhere else, including
  `sacCer3`: a blank cap reads back as `None` and means nobody has chosen a bound, which changes no
  alignment. Every value is set by hand and none is derived from an annotation, whose longest intron
  is a floor on what the organism does rather than a ceiling on it (ADR-0010). Nothing in this
  package reads either column; they are curated here for whoever configures the aligner.
- **Ask what a motif looks like: `MotifSet.compare`.** Hand it one `Motif`, several, or a whole
  `MotifSet` and it compares them against the motifs of the set it is called on — read
  `release.compare(de_novo)` as *compare these against this release*. The use case is naming: a
  chromBPNet or TF-MoDISco run hands back matrices with no names on them, and this says which
  published motif each one most resembles. The comparison is `memelite`'s `tomtom`, handed the same
  4 × L probability matrices `Motif.probabilities` already produces, so **nothing transposes,
  permutes or rescales on the way in**; this is the only place `memelite` is used, its scanner
  having been measured against MOODS and rejected. What comes back is a `MotifComparison` wrapping
  an `xarray.Dataset` **indexed by motif id on both axes**, so `data.sel(query="pattern_0",
  target="MA0139.2")` asks about one pair and `data["neg_log10_p"]` is a similarity matrix ready to
  cluster; it carries negative log10 p, score, offset, overlap and strand, the strand always `+` or
  `-` because a comparison knows which of the two it aligned. **Negative log10 p is stored and raw
  p never is**: the array holds half precision, whose smallest normal value is 6.1e-5, so a p of
  1e-20 stored raw would flush to zero and take every motif's best match down with it — stored as
  20.0 it is an ordinary small number, and an underflowed p reads back as infinite rather than as a
  warning. `to_frame()` flattens to one row per pair, defaulting to the single best target per
  query and taking a larger limit or none at all, ranked by p with score and then the target set's
  own order breaking ties — which is also the order the engine's fast path returns, so the two
  agree on which target is best. **Passing `top=n` is not a convenience over the complete answer**:
  it takes `tomtom`'s nearest-neighbour path, which never scores the targets that lose, so the
  result is *ragged* — dimensions `(query, rank)` with the target ids riding along as a variable,
  because rank 0 names a different motif for each query. Such a result **cannot be widened without
  recomputing, and that is accepted rather than a defect**; `RaggedComparisonError` says so and
  names the call that would. A `top` above the number of targets is refused up front, because the
  engine answers that one with a `SystemError` out of numba that names nothing a caller can act on.
  **A motif compared against itself aligns to itself perfectly** — offset 0, the whole length
  overlapping, on `+`, and no target scores higher. It is usually ranked first too, and the
  exception is documented rather than papered over: TOMTOM's p-value rewards a short dense
  alignment, so a long motif that embeds a shorter one can rank the shorter one above itself, which
  both of the 31- and 33-column CTCF matrices do with the 15-column `MA0139.2` inside them.

### Removed

- **pre-commit, a second toolchain the environment manager did not control.**
  `.pre-commit-config.yaml` pinned `ruff-pre-commit` at `v0.15.16` and let pre-commit build that
  copy in an ephemeral environment of its own, while `pixi.lock` resolved the ruff the gate runs;
  two versions formatting one file can disagree, so a commit the hook passed could fail
  `pixi run check` — and the gate is what CI runs. The hook file also carried its own rules —
  trailing whitespace, end-of-file, YAML and TOML parses, a 500 KB file cap — enforced at commit
  time on one machine and by nothing in CI, which is a check that reads as green because it never
  ran there. Nothing in the repository told a reader to install the hooks except the removed file's
  own header. `pre-commit` leaves the `dev` feature and the lock loses it with `cfgv`, `identify`
  and their build dependencies; the entry point is `pixi run check`, unchanged.

### Fixed

- **CI's numba cache is reused on every run, not only when the runner's CPU model happens to
  match.** numba keys each compiled overload on a tuple whose middle element is
  `(target_triple, cpu_name, cpu_features)`, so the host CPU **model** was part of the key. GitHub's
  runners are heterogeneous and a job does not choose the one it lands on, which made reuse luck:
  the same commit took the test lane to 32 s on one run and 72 s on another, the slow one having
  logged `Cache restored from key: …` and then recompiled four tests at 38.96 s, 26.06 s, 23.63 s
  and 23.60 s. Nothing reported it, because the two caches disagree about what a hit means —
  `actions/cache` succeeds on the *files*, numba then rejects the *contents* per overload — so the
  log said hit and only `--durations` said otherwise. The `test` job now pins
  `NUMBA_CPU_NAME=generic`, collapsing that element to `('…', 'generic', '')`, which is stable
  across machines of one architecture. Giving up CPU-specific vectorisation was measured rather
  than assumed, and on the architecture that matters: on a Sapphire Rapids Xeon with the full
  AVX-512 set — the widest vectors x86-64 offers, so the largest gap the pin can cost — memelite's
  scan engine was a dead heat with host codegen and stayed one as the scanned sequence grew 333×,
  and its compare engine landed within ±4% on either side of host. A genuinely genome-scale scan
  is still unmeasured and deserves its own dated measurement, but nothing found so far makes host
  codegen worth the cache bug. The pin is on the CI job and not in the pixi task, so a developer's
  build keeps host codegen. A new `numba cache key` step prints the target description and the
  source stamps numba also keys on, so a recurrence of either is readable in the log instead of
  inferred from wall-clock.
- **A release named without a source is honoured instead of quietly ignored.** `lookup_xref` — and
  so `XrefSet` — returned the default source's newest release as soon as no source was named, never
  consulting the release asked for. Harmless while every species listed exactly one release, and
  wrong the moment a second source arrives with a numbering of its own: a caller pinning a release
  would have been handed another one, under a release string saying they had not been, which is the
  whole of what pinning is for. It now answers with that release or raises naming the ones the
  default source actually has.
- **`normalise_id` is idempotent, including where whitespace hides behind a version separator.**
  Stripping ran before the version was dropped, so `"7157\r."` stemmed to `"7157\r"` on the first
  pass and only reached `"7157"` on the second — two spellings of one id settling on different
  strings, which joins to nothing and says nothing about it. The suite's own hypothesis property
  found it, and it was a flaky-CI landmine besides: it passed until a machine's example database
  happened to find the case. Whitespace now goes on both sides of the version drop.
- **Neither cross-species subpackage imports the TF half any more, and a test holds it there.**
  `genome.homology` and `genome.xref` both reached into `genome.tf.gene.census` for `species_slug`,
  which the context map forbids in both directions — *Orthology → TF gene* is "a prohibition rather
  than a call" and "the xref half reads no census". The helper names files after a species and is no
  TF concept, so it now lives in `genome.metadata`, the module that owns the curated tables a species
  is spelled in and that every context may read. Every previous import path still answers, and a new
  guard reads both subpackages' source — deferred and `TYPE_CHECKING` imports included — for anything
  naming `genome.tf`.
- **A `cache_dir` handed to a `HomologySet` names the directory it prepares in**, as it already did
  for an `XrefSet` and a `JasparDatabase`. It was being treated as a homology *root* with
  `ensembl_compara/<release>/<pair>/` re-applied beneath it, which is the one exception a caller
  reading the other two would not expect.
- **A malformed quality cell raises instead of reading as "the publisher recorded nothing."**
  `goc_score` and `wga_coverage` fell back to `None` for anything `int()` or `float()` refused,
  putting *Compara scored nothing here* and *this package could not read the score* under one value —
  in the very columns a whole species pair is measured null in. A cell that is neither a number nor
  Compara's own `NULL` now raises `ComparaFileError` naming the file, the column and the value.
- **A crossing into an annotation stops dropping two things it was handed.** `resolve_homologs`
  built its answer with no `null_quality_scores`, though the measurement rides on every other answer,
  and it *replaced* the answer's `dropped_partners` with what the annotation had dropped — losing
  every partner a **Homology type** filter removed before the crossing. A **Dropped partner** is one
  the answer no longer names whichever step removed it, so both causes are now counted together.
- **The repair for an unfinished homology set names the call that rebuilds it**, not only the `rm
  -rf` that empties the directory — the shape the xref half already had.
- **The Orthology glossary no longer says the code does not exist.** `docs/context/orthology.md` and
  the context map's Orthology row both carried the *(decided, not built)* marker this branch made
  false.
- **The project URLs name the account that actually holds the repository.** `Homepage`, `Issues` and
  `Changelog` all said `github.com/lhqing/liulab-genome`, which is not where this project lives —
  and those three are the only copy of that address a package index shows, so the Homepage and
  Issues links on the index page went somewhere that is not this project. All three now name
  `liuhlab`, agreeing with the site config, which had it right. Two manifest tidies ride along,
  neither user-visible: the version fallback is gone, so a checkout with no tags in history fails
  the build loudly instead of stamping a wheel `0.0.0+dev` — both workflows that build already
  fetch full history, so it never fired — and the `tests/**` lint ignore drops `S101`, a rule
  belonging to a rule set `select` does not name and which therefore fired never.

## [2026.8.0] - 2026-08-17

The first tagged release, and it is the whole package: everything below is what `genome` is, not what
changed since a predecessor. Cut because a consumer needs to pin a release rather than a branch.

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
- **Registering an annotation by name.** `genome.annotations.register("gencode_v50")` fetches,
  verifies, checks chromosome names, builds the database and records it. The path-based
  `genome.annotations.register_path(gtf, name)` remains for an unlisted GTF.
- **A chromosome-name check that runs before the database build**, so a GTF spelled for the wrong
  assembly fails in seconds rather than after many minutes. Every sequence name in the GTF must
  exist in the assembly's `chrom.sizes`; the reverse is not required, since an assembly may carry
  scaffolds an annotation never mentions. `check_chromosomes=False`, or
  `--no-check-chromosomes`, overrides it, and the record says which annotations were checked.
- **`genome.annotations.offered`**, the table's rows for this assembly, answering a different
  question from `genome.annotations.registered`, which is "registered on this machine".
- **`genome.annotations.broken`**, the complement of `.registered`: between the two, every
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
- **`Genome.separator` and `Genome.component_annotations`**, both read off the completion record and
  both `None` for an assembly that is not a chimera, exactly as `components` is. The first is the
  underscore run this chimera's names were actually suffixed with, so a caller reading a suffixed
  name back never assumes `__` and never splits an escalated name in the wrong place; the second is
  the registered annotation each component contributed to the merge, or `None` for one that
  contributed none — which the merged name cannot say, since it names only the contributors, and
  which decides the annotation a per-component count is taken against.
- **A chimera's sequence order is a published contract**, not an artefact of the concatenation: one
  contiguous block per component, components in the sorted order the derived name spells, and each
  component's own declared order inside its block. It was already verified after every build and is
  now stated in the Chimera glossary entry, because a consumer that filters one component's sequences
  back out of an alignment header recovers a single-assembly header only while it holds.
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
- **`AnnotationRegistry` is a package export** — `from genome import AnnotationRegistry`. It was an
  `io` internal, and `Genome.annotations` now hands one back, so callers hold one.

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
  well, whose constructor never ran because its answer to *which directory?* was the wrong one. That
  base has since gone entirely — see Removed.
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
- **Registering by path over a directory nothing vouches for now names a command a shell can run**
  — `genome register-gtf <assembly> <gtf> <name> --force` — rather than the equivalent Python call.
  Every way in now knows which assembly it is registering for, so there is no route left that has to
  name a call instead.
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
- **An index record's `parameters` means one thing: every tuning knob that determined the build**,
  caller-supplied and package-computed alike, whether or not the command line also spells it. STAR
  already recorded that superset; chromap recorded only its extra flags and rendered its command
  line from them, so the same key held two contracts. Flags are now rendered from the caller's
  keywords alone. Records already on disk are not migrated. `Genome.build_star_index`,
  `build_chromap_index` and `get_index` also take `tool=` now, forwarding it to the aligner.
- **A chimera's completion-record shape is `genome.io.components`, and the one fetch step is
  `genome.io.fetch`** — split out of `genome.io.source` and `genome.io.download`, which closes the
  `download → chimera → gtf → download` import cycle. `ChimeraDetails` and friends are no longer
  importable from `genome.io.chimera` or `genome.io.source`, and `download.fetch_url` is now
  `fetch.fetch_url` — the patch target for taking a test offline is `genome.io.fetch.fetch_url`.
- **`ChimeraDetails.from_details` takes the assembly name and refuses a broken record.** It returned
  `None` both for *not a chimera* and for *this record is malformed*, so a half-written chimera
  record read back as an ordinary assembly nothing ever checked against its components. Malformed
  now raises `RegistrationError` naming `genome register <assembly> --force` (ADR-0007); *not a
  chimera* is still `None`.
- **The five records a registration answers with are `genome.io.results`** — moved whole out of
  `io/gtf.py` and `io/download.py` with `chromosome_check_summary` and the `EXPECTED_FROM_*`
  constants, and re-exported from `genome.io` as before. The JSON is unchanged, key for key and in
  order; `AnnotationStatusRow.state` and `AnnotationStatus.default_summary` are new, being the two
  sentences the command line used to derive for itself.
- **A binary is asked its version once per process, not once per build step**, remembered against
  the path it was located at — so a build stops re-probing its tools for provenance that cannot
  change. Recorded versions are unchanged; `genome.external.clear_version_cache()` forgets them.
- **`Genome.annotations` is the `AnnotationRegistry`, not `list[str]`.** Everything about the
  *collection* of annotations is now one object — `.registered`, `.broken`, `.offered`,
  `.path(name)`, `.register(name)`, `.register_path(gtf, name)` — while the default annotation, the
  everyday read, stays on `Genome` as `default_gtf` and `default_gtf_path`. It is deliberately not a
  list: a registry settles a four-way state (registered, broken, offered, nothing) and a `__len__`
  or `__iter__` over any one of them would hide which set is being walked, so there is no `len()`,
  no `in` and no iteration. **Every stale use therefore fails loudly — except truthiness.** A
  registry object is always truthy, so `if genome.annotations:` silently flips from *has annotations
  registered* to *always*. That is this change's one silent break, and the fix is
  `if genome.annotations.registered:`. Grep for it.
- **`gtf.register_annotation_by_path` is `gtf.register_gtf`.** Same signature, same
  `RegisteredAnnotation` return; the mouthful existed only because the name was taken. The module
  pair now matches the two CLI commands exactly — `gtf.register_annotation(assembly, name)` behind
  `genome register-annotation`, `gtf.register_gtf(assembly, gtf, name)` behind `genome register-gtf`.
- **`Genome.build_star_index` defaults its annotation.** `gtf` is optional and falls back to
  `default_gtf`, so `chimera.build_star_index(threads=8)` is the everyday call rather than
  `build_star_index(gtf=chimera.default_gtf, threads=8)`. Naming a key explicitly still works and
  still names `index/star_<gtf>/`. With no default at all it raises `ValueError` naming both fixes.
  `build_chromap_index` takes no annotation and is unchanged.

### Removed

- **Four annotation steps leave the package surface** — `list_annotations`,
  `list_broken_annotations`, `default_annotation` and `fetch_annotation` are no longer re-exported
  from `genome.io`, so `from genome.io import fetch_annotation` breaks. Nothing in the package calls
  any of them by name; what a caller wants from an assembly's annotations is `AnnotationRegistry`,
  which takes their place in `genome.io`. The first three stay importable from `genome.io.gtf`;
  `fetch_annotation` is gone altogether, having become one line of delegation to
  `AnnotationRegistry.register`.
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
- **`Downloader`**, the pooch-cache wrapper no caller in the package used.
  `UCSCGenomeDownloader` deliberately did not subclass it, so nothing moves anywhere: reach
  `genome.io.fetch.fetch_url` with a directory of your own instead.
- **Five of `Genome`'s eight annotation members.** The registry `Genome.annotations` returns already
  answered every one of them, so each was a pass-through restating its docstring:

  | Gone from `Genome` | Reached as |
  |---|---|
  | `broken_annotations` | `genome.annotations.broken` |
  | `offered_annotations` | `genome.annotations.offered` |
  | `register_annotation(name)` | `genome.annotations.register(name)` |
  | `register_gtf(gtf, name)` | `genome.annotations.register_path(gtf, name)` |
  | `get_gtf_path(name)` | `genome.annotations.path(name)` |

  `default_gtf`, `default_gtf_path` and `annotations` stay, because the default annotation is the
  usage and the collection is the rare case.
- **The directory-addressed `gtf.register_gtf(assembly_dir, gtf, name, …)`**, which had no caller
  anywhere in the package and was the one registration form that knew no assembly name — so it had
  to be handed a `chrom_sizes` or check nothing, and could only name a Python call as its repair.
  Address the assembly instead: `gtf.register_gtf(assembly, gtf, name, cache_dir=…)`, which now
  carries that name, or `AnnotationRegistry(assembly_dir).register_path(gtf, name)` for a directory
  you hold.

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

**One code change will not announce itself.** `Genome.annotations` returns a registry rather than a
list, and a registry object is always truthy, so `if genome.annotations:` keeps running and stops
meaning anything — it is now true even for an assembly with nothing registered. Every other stale
use raises. Search your code for `.annotations` in a boolean position and write
`.annotations.registered` instead.
