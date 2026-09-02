# Gene identifiers

One gene carries many names: an Ensembl stable id, an NCBI Entrez GeneID, a UniProt
accession, an HGNC id, and the symbol a paper prints. An **xref set** is one publisher's
table of which of those names go together, for one species, at one pinned release. Name a
species and the set is prepared on this machine, then ask it in either direction.

```python
from genome.xref import XrefSet

human = XrefSet("Homo sapiens")
human.to_stems(["7157"], "entrez").gene_id_stems    # ['ENSG00000141510']
```

No genome is opened and no assembly is named. What a gene is called in another species is a
separate question, answered on [Homology](homology.md).

## Preparing a set

The first `XrefSet("Homo sapiens")` downloads the publisher's file, slices it to that one
species and writes the slice as a plain gzipped TSV. Every construction after reads the
slice and fetches nothing.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

The first use downloads. The lab's CPU compute nodes have no internet, so run `genome xref ids` once on a login node before submitting a job that needs it.

A set prints which species, source and release it holds, and how many genes it covers.
`namespaces` is read off the slice, so it lists what this set can actually answer in:

```python
human
# XrefSet(species='Homo sapiens', source='alliance', release='9.0.0', stems=43867)

human.namespaces    # ('ensembl', 'entrez', 'uniprot', 'hgnc')
```

Each species has its own authority, so a human set carries `hgnc`, a mouse set `mgi` and a
worm set `wormbase`. Asking for one the set does not carry raises
`NamespaceNotCarriedError`, and the message lists the ones it does.

## Sources and releases

Three species have a set, and more than one publisher offers each of them:

```python
from genome.xref import xref_releases, xref_sources, xref_species

xref_species()
# ('Homo sapiens', 'Mus musculus', 'Caenorhabditis elegans')

xref_sources("Homo sapiens")                 # ('alliance', 'ensembl', 'hgnc')
xref_releases("Homo sapiens", "alliance")    # ('9.0.0',)
```

Leave the source out and the species' default answers, so everyone in the lab reaches for
the same one. **The default is per species and per question**, because the publisher that
carries a species' ids is usually not the one that carries its symbols. `lookup_xref` says
which row will answer, without downloading anything:

```python
from genome.xref import lookup_xref

lookup_xref("Homo sapiens").source                      # 'alliance'
lookup_xref("Homo sapiens", for_symbols=True).source    # 'hgnc'
lookup_xref("Mus musculus", for_symbols=True).source    # 'alliance_bgi'
```

Name a source, as `XrefSet("Homo sapiens", "ensembl")`, and it is honoured rather than
swapped for a default. Two publishers can disagree about which Ensembl gene an Entrez id
names, so pick one deliberately when it matters; every answer records the source and
release that produced it. Naming a species, source or release nothing is prepared for
raises `NoXrefSetError`.

There is more here than the common path uses: an `evidence=` filter for publishers that
grade their rows, the `normalise_id` and `normalise_symbol` helpers, and a reader module per
publisher. All of it is in the [API reference](../reference.md).

## Converting ids

Everything is joined through one namespace. A **gene id stem** is an Ensembl gene id with
any version suffix removed, and it is what both directions have in common. `to_stems` reads
ids in a namespace and answers in stems; `from_stems` reads stems and answers in a
namespace. Which one you want is named rather than guessed from the string:

```python
human.to_stems(["7157", "672"], "entrez").resolved
# {'7157': ('ENSG00000141510',), '672': ('ENSG00000012048',)}

human.from_stems(["ENSG00000141510"], "uniprot").resolved
# {'ENSG00000141510': ('P04637',)}
```

The keys are your own spelling of what you asked about, in the order you asked, so the
answer lines up against your table row for row. **A value is every gene the publisher
asserts, never one picked for you.** An HGNC id that names two genes answers with both:

```python
human.to_stems(["HGNC:11998", "HGNC:13666"], "hgnc").resolved
# {'HGNC:11998': ('ENSG00000141510',),
#  'HGNC:13666': ('ENSG00000094914', 'ENSG00000291836')}
```

`gene_id_stems` flattens those values into one list. It drops which id named which gene, so
read `resolved` whenever that matters. There is no direct hop from one namespace to
another: Entrez to UniProt is two calls with the stem in the middle, and you join them.

### What did not resolve

Ids this release names nothing for ride back on the answer instead of being dropped:

```python
answer = human.to_stems(["7157", "999999999"], "entrez")
answer.resolved      # {'7157': ('ENSG00000141510',)}
answer.unresolved    # ('999999999',)
```

**Read `unresolved` on every answer**, because a hand-written join drops these rows without
saying so. It is in ask order, and a stem that carries no id in the namespace you asked for
lands there too. Nothing here separates an id a publisher retired from one it never
carried.

## Versioned ids and stems

GENCODE and Ensembl's own GTFs spell a gene id with a version suffix,
`ENSG00000141510.18`. Most cross-reference files spell it bare. **Joining a versioned column
to a bare one matches nothing and raises no error**, which is the usual cause of an empty
lookup with nothing to explain it.

The version is dropped on both sides here, so the two spellings are one identifier:

```python
human.from_stems(["ENSG00000141510.18", "ENSG00000141510"], "uniprot").resolved
# {'ENSG00000141510.18': ('P04637',), 'ENSG00000141510': ('P04637',)}
```

Two keys, because the keys are your spelling, and the same value under both. For a column
you are cleaning up before a join of your own, `gene_id_stem` does the reduction on its own:

```python
from genome.xref import gene_id_stem

gene_id_stem("ENSG00000141510.18")    # 'ENSG00000141510'
```

Worm ids are published unversioned, so `WBGene00003001` is already a stem and the reduction
changes nothing.

## Matching symbols

The two directions are different questions. Starting from a stem, a gene has one current
approved symbol, which is what a figure axis wants, and `from_stems` gives it:

```python
symbols = XrefSet.for_symbols("Homo sapiens")
symbols.from_stems(["ENSG00000133794"], "symbol").resolved    # {'ENSG00000133794': ('BMAL1',)}
```

Starting from a symbol, the answer is genes, and the match is made against spellings the
authority has since retired as well as the current one. That is `match_symbols`, and every
hit says which kind of spelling matched:

```python
answer = symbols.match_symbols(["ARNTL", "ADCY3", "Brca1"])
answer.resolved
# {'ARNTL': (SymbolMatch(symbol='ARNTL', gene_id_stem='ENSG00000133794', kind='previous'),),
#  'ADCY3': (SymbolMatch(symbol='ADCY3', gene_id_stem='ENSG00000138031', kind='approved'),
#            SymbolMatch(symbol='ADCY3', gene_id_stem='ENSG00000155897', kind='previous'))}
answer.unresolved    # ('Brca1',)
```

`ARNTL` is a spelling HGNC retired and it still reaches `BMAL1`. `ADCY3` is HGNC's approved
symbol for one gene and a symbol it took away from another, so it answers with both.
**Matching is exact by default**, so `Brca1`, a mouse spelling asked of the human
authority, matches nothing. Pass `case_insensitive=True` to fold both sides.

`XrefSet.for_symbols` fills in the source that carries symbols; the plain constructor fills
in the identifier one, which for human carries no symbol at all, so
`XrefSet("Homo sapiens").match_symbols(...)` raises `NamespaceNotCarriedError`. Asking
`to_stems` for the symbol namespace raises `SymbolDirectionError`, whose message names
`match_symbols`.

### What a source could not match

Only HGNC labels its previous and alias spellings as such. The mouse and worm authorities
each publish one approved symbol per gene beside an undifferentiated synonyms list, so
those sets match approved spellings and nothing else. The pair that works for human does
not work for mouse:

```python
mouse = XrefSet.for_symbols("Mus musculus")
answer = mouse.match_symbols(["Bmal1", "Arntl"])
answer.resolved
# {'Bmal1': (SymbolMatch(symbol='Bmal1', gene_id_stem='ENSMUSG00000055116', kind='approved'),)}
answer.unresolved    # ('Arntl',)
answer.kinds         # ('approved',)
```

**A symbol this release does not have and a spelling this source could never have matched
both land in `unresolved`**, and the entry alone does not say which happened. `kinds` says
which kinds the set can match, and `limits` says why the rest are missing. `limits` is a
long string, trimmed here after its first clause:

```python
print(answer.limits)
# this source publishes one current approved symbol per gene beside an undifferentiated
# synonyms list, and a spelling in that list does not say whether the authority retired
# it or merely records it ...
```

For a human symbol set, `kinds` is `('approved', 'previous', 'alias')` and `limits` is
`None`.

## Citing what answered

A set carries the curated row it resolved to, so what you cite is what answered rather than
what you looked up afterwards:

```python
print(symbols.provenance.attribution())
# HGNC 2026-07-07 (PMID 41287213) —
# https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-07-07.txt
```

Publisher, their own version string, the paper to cite where there is one, and the file the
bytes came from. `provenance` is the whole record, including the species, the NCBI taxid
and the publisher's checksum. Each command below prints its own source URL on every answer.

## From the command line

Neither question needs a genome, and both have a command. `genome xref ids` converts in the
direction the flag names:

```console
$ genome xref ids "Homo sapiens" --to-stems entrez 7157 672 999999999
entrez ids -> gene id stems for Homo sapiens (alliance 9.0.0)
  source  https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz
  2 resolved, 2 gene id stems, 1 this release names none for
7157 ENSG00000141510
672 ENSG00000012048
999999999
```

The pairs go to stdout, tab-separated, so `cut -f2` is the answer and the rest pipes. The
heading and the counts go to stderr. Every id you passed gets a row, and one that resolved
to nothing gets an empty second column.

`genome xref symbols` matches, and prints the kind of spelling in the fourth column:

```console
$ genome xref symbols "Homo sapiens" ARNTL ADCY3 Brca1
gene symbols -> gene id stems for Homo sapiens (hgnc 2026-07-07)
  source   https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-07-07.txt
  columns  asked, symbol, gene_id_stem, kind
  matching exact, on approved, previous, alias spellings
  2 resolved, 3 matches, 1 this release matched nothing for
ARNTL ARNTL ENSG00000133794 previous
ADCY3 ADCY3 ENSG00000138031 approved
ADCY3 ADCY3 ENSG00000155897 previous
Brca1
```

Both take `--source` and `--json`. Both exit non-zero when a species, source or namespace is
not one this package has, rather than printing an answer shorter than you asked for. The
full option list is on [CLI: lookup commands](../cli/lookups.md).
