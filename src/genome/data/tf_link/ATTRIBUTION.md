# Attribution — the shipped Motif link tables

Which JASPAR motifs answer for which **TF gene**: the profiles are JASPAR's, the verdicts on which
genes are transcription factors are the censuses', **and neither is this package's**.

| | |
|---|---|
| Motifs | JASPAR, <https://jaspar.elixir.no/>, read from one release's SQLite dump, `https://jaspar.elixir.no/download/database/JASPAR<year>.sqlite`. **Cite the JASPAR release you used.** |
| Genes | The censuses, attributed in `../tf_gene/ATTRIBUTION.md` — Lambert *et al.* 2018 v_1.01 (PMID 29425488) for human, AnimalTFDB 4.0 (PMID 36268869) for mouse. **Cite the one you used.** |
| Licensing | Redistributed here: JASPAR's identifiers, published names and per-profile species, plus one number computed from its matrices — no count matrix. JASPAR is open-access, published for reuse, and asks to be cited; the censuses state no licence and are redistributed under the attribution above. |

`<species slug>.jaspar<release>.motif_link_table.tsv.gz` ships, one per species per **Release**, CORE
`vertebrates` only, twelve columns identical across tables so two concatenate and still say what each
row came from. No quoting; multi-value cells join on `;`. Gzipped because they are bulk, where the
three-row `motif_name_alias.tsv` beside them is plain, as every small metadata table here is.

| Table | Genes | Links | Cross-species |
|---|---|---|---|
| `homo_sapiens.jaspar2024` | 745 | 946 | 161 |
| `homo_sapiens.jaspar2026` | 876 | 1,085 | 162 |
| `mus_musculus.jaspar2024` | 653 | 851 | 690 |
| `mus_musculus.jaspar2026` | 693 | 896 | 732 |

## How a link is made, and what stays unlinked

Upper-case the **Motif name**, split it on `::`, match each part against the census's own symbol
column. **Only assessed-positive genes receive links**: `MA1964.2` is named for human `SMAD2`, which
Lambert turned down, so it is unlinked rather than aliased. Parts no census spells that way go
through `motif_name_alias.tsv`, **keyed on the Gene id stem, not on a symbol** — three rows, the
human profiles JASPAR renamed after Lambert (`TBXT` for `T`, `SCAND3` for `ZBED9`, `ZFTA` for
`C11orf95`), each checked against the census as the table is built. `role` and `partners` say whether
the matrix is a monomer's or a complex's, so `FOS::JUN` links to both genes as a complex and to
neither as a monomer; `rank` encodes **Attribution specificity** and is **not** a quality score;
`is_cross_species` marks a profile measured on another vertebrate, kept and marked (ADR-0013).

Three profiles stay unlinked, each correctly: `MA0149.1 EWSR1-FLI1` is an oncogenic fusion naming no
gene; `MA2503.1 Banp` is not in AnimalTFDB; and **`MA0611.3 Dux` is not a rename** — JASPAR's `Dux`
is UniProt A1JVI8, **MGI:3703875**, AnimalTFDB's `Duxf3` is ENSMUSG00000075046, **MGI:1921649**, and
two MGI accessions are two genes, so an alias row would assert an identity MGI denies. That question
is closed, and `tests/test_tf_link.py` pins it.

## Rebuilding, and what no test proves

`python scripts/build_tf_links.py <species> <release> <that release's SQLite dump>` — the generator
lives outside the wheel, fails loudly when JASPAR re-spells a table or a column, and writes
byte-stable output, so a re-release is a re-run and a reviewable diff. **No test proves a shipped
table still matches JASPAR, and none can**: regenerating one needs a download and CI has no network,
the limitation ADR-0011 already accepts. The counts above are pinned in `tests/test_tf_link.py`,
turning silent drift into a loud failure — the most that is available.
