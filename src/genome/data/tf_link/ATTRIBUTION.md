# Attribution — the shipped Motif link tables

Every file in this directory says which JASPAR motifs answer for which **TF gene**. The
profiles are JASPAR's and the verdicts on which genes are transcription factors are the
censuses' — **none of either is this package's**. The tables are plain TSV so that the
mapping is readable in R or a shell by collaborators who never import this package, which
is what it is most often wanted for (ADR-0015).

## What ships

`<species slug>.jaspar<release>.motif_link_table.tsv` — one table per species per
**Release**, CORE `vertebrates` only. Two keys name one table, so both are in the name,
and a new release or a new species is a file dropped in rather than a line of code.

| File | Genes | Links | Cross-species |
|---|---|---|---|
| `homo_sapiens.jaspar2024.motif_link_table.tsv` | 745 | 946 | 161 |
| `homo_sapiens.jaspar2026.motif_link_table.tsv` | 876 | 1,085 | 162 |
| `mus_musculus.jaspar2024.motif_link_table.tsv` | 653 | 851 | 690 |
| `mus_musculus.jaspar2026.motif_link_table.tsv` | 693 | 896 | 732 |

Twelve columns, the same twelve in every table, so two tables concatenate into one frame
that still says what each row came from: `release`, `species`, `gene_id_stem`, `symbol`,
`motif_id`, `motif_name`, `role`, `partners`, `motif_tax_ids`, `is_cross_species`,
`total_information_content`, `rank`. `partners` and `motif_tax_ids` join their values
with `;`; a flag is spelled `yes` or `no`, as a census spells its own.

`motif_name_alias.tsv` beside them is the alias table — the **Motif name** parts no
census spells that way.

## Attribution

The profiles are JASPAR's, <https://jaspar.elixir.no/>, read from the SQLite dump of one
release, `https://jaspar.elixir.no/download/database/JASPAR<year>.sqlite`. **Cite the
JASPAR release you used.** The genes are the censuses', which carry their own attribution
in `../tf_gene/ATTRIBUTION.md` — Lambert *et al.* 2018 for human, AnimalTFDB 4.0 for
mouse. What is redistributed here is JASPAR's identifiers, its published names and its
per-profile species, plus one number computed from its matrices; no matrix is copied.

## How a link is made

The **Motif name** is upper-cased and split on `::`, and each part is matched against the
census's own symbol column. **Only assessed-positive genes receive links**, so a gene a
census assessed and turned down is unlinked however a profile is named — human `SMAD2`
is the case to know: `MA1964.2` is named for it, Lambert assessed it `No`, and it is
correctly unlinked rather than aliased.

The parts no census spells that way are handled by `motif_name_alias.tsv`, which is
**keyed on the Gene id stem and not on a symbol**: a symbol is the thing that moved, so
keying on one keys on the moving part. All three rows are human profiles JASPAR renamed
after Lambert's census was published — `TBXT` for `T`, `SCAND3` for `ZBED9`, `ZFTA` for
`C11orf95` — and the generator checks each one against the census it names, so a stale
row is an error rather than a comment nobody reads. `SCAND3` in particular was got wrong
twice by guessing at symbol history, which is why the key is the gene id.

A profile that names no gene at all stays unlinked **by design**: `MA0149.1`
`EWSR1-FLI1` is an oncogenic fusion, and asserting a gene for it would be inventing one.

Two mouse-tagged profiles resolve to no mouse gene and are shipped unlinked:
`MA2503.1 Banp`, which AnimalTFDB does not list and Lambert assessed and turned down; and
`MA0611.3 Dux`, which AnimalTFDB spells `Duxf3`. The second is a rename of the kind the
alias table exists for and has not been curated into it — the three rows above are what
was verified, and adding a fourth is a change to the counts pinned below and so a
reviewable one.

## Role, ordering, and the species flag

`role` is `monomer` where the profile names one gene and `complex` otherwise, with the
other named parts in `partners`, so a heterodimer matrix is never read as a monomer's:
`FOS::JUN` links to both genes as a complex and to neither as a monomer. A gene whose
only motifs are complexes is linked rather than reported motif-less — AHR, DDIT3, TAL1
and TLX1 are exactly those genes on the 2026 release.

`rank` is dense from one within a gene and encodes **Attribution specificity**: `monomer`
before `complex`, species-matched before a **Cross-species link**, then higher
`total_information_content`, then `motif_id`, so the order is total and stable across
machines and releases. It says what a matrix is attributable to and explicitly **not**
which motif is better — the canonical AP-1 matrix is a complex and describes JUN's
binding better than any JUN monomer does. **No quality score is computed or shipped**:
JASPAR publishes none, and matrix depth is normalised per assay, so ranking on depth
ranks the assay.

`total_information_content` is the matrix's **Information content** summed over its
columns, in bits — the package's own `Motif.information_content`, computed from the
dump's count matrix rather than reimplemented. The dump's counts and the transfac file
the package's loader reads agree to within 1e-6 bits on every profile of both releases.

`is_cross_species` marks a profile measured on a vertebrate other than the gene's own
species. Such links are **kept and marked, never excluded** (ADR-0013): JASPAR files an
orthologous pair's matrix under whichever species was assayed, so the species field is an
artefact of the experiment rather than a claim about the factor, and excluding them costs
mouse 732 of its 896 links on 2026. A profile the dump records no species for at all —
`MA0108`, TBP — has a blank `motif_tax_ids` and is marked cross-species, because the row
cannot claim a species match it has no evidence for.

## Rebuilding

`python scripts/build_tf_links.py <species> <release> <that release's SQLite dump>`. The
generator lives outside the wheel, fails loudly when JASPAR re-spells a table or a
column, checks the profile count against what the package's own loader pins for that
release, and writes byte-stable output — so a re-release is a re-run and a reviewable
diff. It reads the dump rather than the transfac file because the transfac file carries
no per-profile species, which is why **the JASPAR loader needs no change** to make this
join possible.

## What no test proves

**No test proves a shipped table still matches JASPAR, and none can.** Regenerating one
needs a download and CI has no network — the same limitation ADR-0011 already accepts for
the curated gene lists, recorded here rather than papered over. `tests/test_tf_link.py`
pins the counts in the table above, which converts a silent drift into a loud failure.
That is the most that is available.

## Licensing

JASPAR is an open-access database published for reuse and asks to be cited; the censuses
state no licence and are redistributed under the attribution above. The tables here are a
derived join carrying identifiers, names and one computed number, and no count matrix.
