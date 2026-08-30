# Transcription factors

A **transcription factor** is a protein that binds a particular DNA sequence and changes how
nearby genes are transcribed. A published census says which genes those are, one census
per species. Name an assembly and you get that census back in the annotation's own gene
ids, so the answer joins straight onto a counts matrix.

```python
from genome import Genome

tfs = Genome("hg38").tf_gene_list()
```

Nothing is downloaded and nothing is written to answer that: the censuses ship inside the
package. Loading a JASPAR release and scanning sequence for where its matrices occur is on
[Motifs](motifs.md).

## Which genes are transcription factors

`tf_gene_list()` answers against one registered annotation, the assembly's default unless
you name another with `annotation=`. Three things come back. `genes` holds one entry per
gene the census judged a transcription factor and this annotation carries. `gene_ids`
flattens those into the ids the annotation itself spells, version suffix included where it
carries one. `unresolved` holds the census gene id stems this annotation has no gene for,
so what the crossing drops is countable rather than silent.

A census is keyed by **gene id stem**, a gene id with its version suffix dropped:
`ENSG00000137203` and not `ENSG00000137203.12`. One stem can name two gene ids in one
annotation, a pseudoautosomal gene and its `_PAR_Y` copy among them, and both come back
rather than one being picked.

**The species is not a parameter.** It is read from the assembly's own metadata row, so
you cannot ask a mouse assembly for human transcription factors. Registering the
annotation the ids come back in is on [Annotations](../genome/annotations.md).

## What the census recorded about each gene

Four fields are the same in every census: the gene id stem, the symbol, the TF flag, and the
**DBD family**, which is the DNA-binding domain the publisher classifies the factor under.
Everything else a publisher recorded stays under that publisher's own column name in
`judgements`, and how much there is varies a lot:

```python
from genome.tf.gene import tf_gene_table

tf_gene_table("Homo sapiens").columns[4:]
# ('tf_assessment', 'binding_mode', 'motif_status', 'interpro_ids',
#  'vaquerizas_2009_classification', 'cisbp_considers_it_a_tf',
#  'tfcat_classification', 'is_a_go_tf', 'is_c2h2_zf_krab')
tf_gene_table("Mus musculus").columns[4:]    # ()
```

Those are the keys of `judgements` on every gene in the answer, so tightening or loosening
the criteria is a filter over what you already hold rather than a second call. Lambert's
`tf_assessment` takes the values `Known motif`, `Likely to be sequence specific TF` and
`Inferred motif` across the genes it accepted:

```python
known = [gene for gene in tfs.genes if gene.judgements["tf_assessment"] == "Known motif"]
```

**Group by `dbd_family` within one species and never across two.** The publishers did not
harmonise their vocabularies, so Lambert's `ARID/BRIGHT` and AnimalTFDB's `ARID` are two
spellings rather than one family counted twice.

## Assessed, rejected and never assessed

Three states, and two of them are easy to confuse. Lambert assessed 2,765 human genes
and judged 1,639 of them transcription factors. The other 1,126 are a verdict: assessed and
turned down. A gene in neither number was never looked at.

```python
census = tf_gene_table("Homo sapiens")
len(census), len(census.assessed_positive)    # (2765, 1639)
```

`tf_gene_list()` answers with the assessed-positive genes, which is what most callers
want. `include_rejected=True` widens it to every gene the census assessed, each carrying
the census's verdict in `is_tf`:

```python
everything = Genome("hg38").tf_gene_list(include_rejected=True)
```

**A gene absent from either answer was never assessed, which is not the same as rejected.**
Whether a census records rejections at all is its publisher's choice: AnimalTFDB lists only
genes it accepts, so mouse's 1,611 rows are all positive and `include_rejected=True` changes
nothing for a mouse assembly.

## Which genes are cofactors

A **transcription cofactor** acts on transcription without binding a DNA sequence of its
own, so it has no motif. Chromatin remodellers, histone-modifying enzymes and Mediator
subunits are cofactors. Which genes those are is published too, and `tf_cofactor_list()`
answers in the same shape as the census half:

```python
from genome import Genome

worm = Genome("ce11")
cofactors = worm.tf_cofactor_list()
len(cofactors.cofactors), len(cofactors.unresolved)    # (317, 0)
cofactors.gene_ids[:3]
# ['WBGene00000064', 'WBGene00000099', 'WBGene00000105']
```

Each entry says which publisher listed the gene, and keeps that publisher's own vocabulary
under that publisher's own column name:

```python
entry = cofactors.cofactors[0]
entry.symbol, entry.source    # ('act-2', 'animaltfdb')
entry.classifications["animaltfdb_category"]    # 'Chromatin Remodeling Factors'
```

`source` is `animaltfdb`, `epifactors`, or `both` where two publishers listed the same gene.
**`both` means the two agree the gene is a cofactor and says nothing about how either
classified it**, so group within one publisher's columns and never across two. Human is the
only species here built from more than one publisher; mouse and worm come from one each.
There is no `include_rejected` on this half, because no publisher shipping today releases a
rejected set.

## When no table answers for a species

The two halves do not cover the same species. Nobody has published a worm TF census, but a
publisher has listed worm cofactors, so `ce11` raises from one call and answers from the
other:

```python
from genome import Genome, NoCofactorTableError, NoTFCensusError

worm = Genome("ce11")

try:
    tfs = worm.tf_gene_list()
except NoTFCensusError:
    tfs = None

try:
    cofactors = worm.tf_cofactor_list()
except NoCofactorTableError:
    cofactors = None

tfs, len(cofactors.cofactors)    # (None, 317)
```

Nothing on a `Genome` says in advance which way either call will go, so the `try` is how you
ask. **Neither error means the species has no transcription factors or no cofactors.** It
means nobody has published for it. Both messages name the species that are covered. A
chimera, and an assembly the curated metadata table does not list, raise
`UnknownSpeciesError` instead, because nothing on them says what species they are.

## Citing the publisher

Every verdict on this page is a publisher's. `provenance.attribution()` renders the line to
print beside anything they answered, and it is on the answers `tf_gene_list()` and
`tf_cofactor_list()` return as well as on the shipped tables:

```python
print(tf_gene_table("Homo sapiens").provenance.attribution())
# Lambert et al. 2018 v_1.01 (PMID 29425488) — https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
```

A cofactor table built from more than one publisher renders one line per publisher, joined,
and `provenance.sources` is what carries them apart:

```python
from genome.tf.cofactor import cofactor_table

[source.publisher for source in cofactor_table("Homo sapiens").provenance.sources]
# ['AnimalTFDB', 'EpiFactors', 'HGNC']
```

## Which motifs answer for a factor

`motif_links` maps one transcription factor to the JASPAR profiles attributed to it. It
reads a curated table shipped inside the package. No assembly is involved and nothing is
downloaded:

```python
from genome.tf import motif_links

jun = motif_links("JUN", "Homo sapiens")
len(jun), jun.release    # (12, '2026')
jun.motif_ids[:3]
# ('MA0488.2', 'MA0489.3', 'MA1131.2')
```

Name the gene by its census symbol or by its gene id stem, and pass the species, which is
never guessed from the id. A versioned gene id raises `VersionedGeneIdError`, and the
message names the stem to pass instead.

Each `MotifLink` says what the matrix is a motif *of*. `role` is `monomer` where the profile
names this gene alone and `complex` where it names a dimer, with `partners` holding the other
genes in the name:

```python
link = jun[2]
link.motif_id, link.motif_name    # ('MA1131.2', 'FOSL2::JUN')
link.role, link.partners          # ('complex', ('FOSL2',))
```

Links come back monomers first, then species-matched profiles before ones measured on
another vertebrate, then higher information content. **That order says what the matrix is
attributable to, not which matrix is better**, and every field it sorts on is on the link,
so you can re-sort. `cross_species=False` drops the profiles measured on another vertebrate,
which is how a question that needs species-matched matrices asks for them:

```python
motif_links("Jun", "Mus musculus", cross_species=False).motif_ids    # ('MA0489.3',)
```

An empty answer has two causes. A gene the census turned down receives no links at all, and
a gene it accepted may still have no profile in this release. `is_tf` tells the two apart.
The rest of the link's fields are in the [API reference](../reference.md).

```python
smad2 = motif_links("SMAD2", "Homo sapiens")
len(smad2), smad2.is_tf    # (0, False)
```

## A gene with no motif to look for

Two absences raise here, and `TranscriptionCofactorError` subclasses `GeneNotAssessedError`,
so **the narrower clause has to come first** or it never runs:

```python
from genome.tf import GeneNotAssessedError, TranscriptionCofactorError, motif_links

for symbol in ("JUN", "WDR5", "GAPDH"):
    try:
        links = motif_links(symbol, "Homo sapiens")
    except TranscriptionCofactorError:
        print(symbol, "is a cofactor, so no motif answers for it")
    except GeneNotAssessedError:
        print(symbol, "was never assessed")
    else:
        print(symbol, "has", len(links), "motifs")
# JUN has 12 motifs
# WDR5 is a cofactor, so no motif answers for it
# GAPDH was never assessed
```

The subclassing follows the censuses, not the biology: a cofactor is not a kind of
transcription factor, but no census assessed one, so it lands in the never-assessed set.
`motif_links` asks the census first, so a gene that is both, TBP among them, is answered by
the census rather than raised over. The cofactor error is reached only where no census
assessed the gene at all.
