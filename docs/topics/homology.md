# Homology

Two genes in two species are **homologous** when they descend from one gene in a common
ancestor. Ensembl Compara builds a gene tree for each protein family across the species it
covers and reads the relationships off that tree. Name a species pair here and you get
those published links for whatever gene ids you ask about.

```python
from genome.homology import HomologySet

worm_human = HomologySet("Caenorhabditis elegans", "Homo sapiens")
worm_human.homologs(["WBGene00020462"]).homolog_gene_id_stems
# ['ENSG00000177479']
```

This page answers what a gene is in another species. What the same gene is called
elsewhere in its own species is [Gene identifiers](gene-identifiers.md).

## Preparing a set

Constructing a `HomologySet` prepares it. The first construction downloads the Compara
dump that holds the pair, checks the bytes against the publisher's own md5 as they arrive,
slices them to the two species and records that it finished. Every construction after that
reads the slice and fetches nothing.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

The first use downloads. The lab's CPU compute nodes have no internet, so run `genome homology links` once on a login node before submitting a job that needs it.

Three species are pinned, at Compara release 116, and any pair of them answers. `len()`
says how many links the prepared set holds:

```python
from genome.homology import homology_species

homology_species()
# ('Caenorhabditis elegans', 'Homo sapiens', 'Mus musculus')
len(worm_human)    # 23982
```

The pair is unordered for the download and ordered for the question.
`HomologySet("Homo sapiens", "Caenorhabditis elegans")` reads the same prepared file as its
reverse and answers about human genes rather than worm ones. Asking about a species this
package does not prepare raises `UnknownHomologySpeciesError`, and asking for a pair the
shipped table does not pin at that release raises `NoHomologyPairError`.

## Asking for homologs

`homologs()` takes gene id stems and returns a `HomologyAnswer`. `resolved` maps every stem
that named at least one link to all of them, in the order you asked. A stem this release
names no homolog for goes into `unresolved` instead of vanishing:

```python
answer = worm_human.homologs(
    ["WBGene00008317", "WBGene00012385", "WBGene00008352"]
)
{stem: len(links) for stem, links in answer.resolved.items()}
# {'WBGene00008317': 1, 'WBGene00012385': 29}
answer.unresolved
# ('WBGene00008352',)
```

**Ids go in bare, with no version suffix.** Compara writes its gene ids without one, so
`WBGene00020462.1` would match nothing and is refused with `VersionedGeneIdError`, whose
message names the stem to pass.

Each entry of `resolved` is a `HomologyLink` holding the two gene ids and what the
publisher said about them:

```python
link = answer.resolved["WBGene00008317"][0]
link.homolog_gene_id_stem    # 'ENSG00000164074'
link.homology_type           # 'ortholog_one2one'
```

Orthologs come back by default. `paralogs=True` also returns the links whose label calls a
duplication rather than a speciation, of which release 116 publishes none for these pairs.
The rest of the answer, `as_json()` included, is in the [API reference](../reference.md).

## One gene, many genes, or none

A gene can have several counterparts in the other species, and it can have none. The answer
above holds both cases, and `WBGene00012385` is the crowded one:

```python
len(answer.resolved["WBGene00012385"])         # 29
answer.resolved["WBGene00012385"][0].homolog_gene_id_stem
# 'ENSG00000083838'
len(answer.homolog_gene_id_stems)              # 30
```

`homolog_gene_id_stems` flattens the whole answer to the other species' ids and loses which
asked stem reached which partner. **Taking the first entry as "the" human gene is wrong for
most worm genes**, because most links in this set are many-to-many. Read `resolved` when
you need to know which gene an id came back for.

Absence is as common. This set links 6,663 worm genes to human ones, and `ce11`'s WormBase
annotation carries 46,926 genes, so most worm genes name no human homolog at all.

## Homology type

`homology_type` is Compara's own tree-derived label, carried through verbatim. This package
never recomputes it and never derives one of its own.

| Label | What the tree says |
| --- | --- |
| `ortholog_one2one` | One gene on each side, with no duplication since the species split. |
| `ortholog_one2many` | One gene on this side, several on the other, duplicated after the split. |
| `ortholog_many2many` | Both sides duplicated, so neither gene has a single counterpart. |

Of the 23,982 links in the worm-human set, 2,764 are `ortholog_one2one`, 3,073 are
`ortholog_one2many` and 18,145 are `ortholog_many2many`.

Because the label is the publisher's, **a view that has been filtered down to one partner
still reads `ortholog_one2many`**. Partners that fell away are counted in
`dropped_partners` rather than folded into a corrected label.

## Quality scores

Compara publishes three confidence fields on a link, and all three come through as it wrote
them. `is_high_confidence` is its own flag, true on 2,979 of this set's 23,982 links.
`goc_score` is gene order conservation and `wga_coverage` is whole-genome-alignment
coverage.

**Compara records neither `goc_score` nor `wga_coverage` on any link of this set**, so a
filter on either empties the answer rather than narrowing it. Both the set and every answer
it gives name the fields they hold no value in, so you can check before writing the filter:

```python
worm_human.null_quality_scores    # ('goc_score', 'wga_coverage')
answer.null_quality_scores        # ('goc_score', 'wga_coverage')
link.goc_score                    # None
link.is_high_confidence           # True
```

Which fields carry values is a property of the pair. Compara scored human against mouse, so
that set lists nothing as null and both columns are there to filter on:

```python
human_mouse = HomologySet("Homo sapiens", "Mus musculus")
human_mouse.null_quality_scores    # ()
human_mouse.homologs(["ENSG00000141510"]).resolved["ENSG00000141510"][0].goc_score
# 100
```

## Gene ids in your own annotation

A `HomologyAnswer` speaks in stems. `resolve_homologs` puts the other species' stems into
one registered annotation's own gene ids, which is what a counts matrix is indexed by:

```python
from genome import Genome
from genome.homology import resolve_homologs

human_worm = HomologySet("Homo sapiens", "Caenorhabditis elegans")
answer = human_worm.homologs(["ENSG00000152670", "ENSG00000177479"])
crossed = resolve_homologs(answer, Genome("ce11").annotations, "wormbase_ws298")
crossed.homolog_gene_ids
# ['WBGene00001598', 'WBGene00001599', 'WBGene00001600', 'WBGene00020462']
```

WormBase writes its gene ids bare, so these come back spelled the same. An annotation that
versions its ids gives you the versioned spellings.

The registry has to belong to an assembly of the answer's `other_species`. **Nothing checks
that**, so handing a human answer a worm registry resolves nothing instead of raising. Omit
the annotation name and the assembly's default annotation answers; naming one that is not
registered there raises `AnnotationNotRegisteredError`.

Partners the annotation carries no gene for are added to `dropped_partners`, and asked
stems left naming nothing end up in `unresolved`.

## Where the links come from

Every link, label and confidence value in a set is Ensembl Compara's, at the release the
shipped provenance table pins. `provenance` is that table row and `attribution()` renders
the line to print beside anything the set answered:

```python
print(worm_human.provenance.attribution())
# Ensembl Compara release 116 (PMID 26896847) — https://ftp.ensembl.org/pub/
# release-116/tsv/ensembl-compara/homologies/homo_sapiens/
# Compara.116.protein_default.homologies.tsv.gz
```

The URL is one line and is wrapped here to fit. `source_url` is that address on its own,
and `path` is the sliced TSV on this machine, a gzipped file with the publisher's own
header and rows.

## Nothing is mapped for you

Nothing crosses species unless you call `homologs()` yourself. No table this package
publishes is derived through homology, no answer is silently species-mapped, and a species
with no data of its own never borrows another's.
