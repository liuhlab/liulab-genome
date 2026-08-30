# Annotations

An **annotation** says where the genes, transcripts and exons sit on an assembly's
coordinates. It is a GTF file, one line per feature, and registering it builds a SQLite
database next to that text so later queries do not have to re-read it. A genome can hold
several annotations, each under a short name.

```python
from genome import Genome

worm = Genome("ce11")
worm.annotations.registered    # ['wormbase_ws298']
worm.default_gtf               # 'wormbase_ws298'
```

Matching gene ids from a published table onto an annotation's own is covered on
[Gene identifiers](../topics/gene-identifiers.md).

## Registering an annotation

Name one the annotation table lists for this assembly:

```python
worm.annotations.register("wormbase_ws298")
```

That fetches the GTF from the URL the table's row pins, checks the unpacked file against the
sha256 in the same row, checks its chromosome names against the assembly's, and builds the
database. The record saying all of it finished is written last.

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

Each annotation gets its own directory, `<assembly dir>/gtf/<name>/`, holding the GTF, the
`.db` and that record.

How long it takes depends on the size of the GTF. Yeast's `ensgene_v101` is a 5 MB download
and the whole registration finishes in about a second. Worm's `wormbase_ws298` is 183 MB and
builds a 452 MB database, which is 21 seconds of work once the file is on disk, plus the
download. Registering a name that is already registered returns straight away: nothing is
fetched, nothing is rebuilt.

The same work runs from a shell:

```console
$ genome annotation register ce11 wormbase_ws298
registered wormbase_ws298 for ce11 in /Users/hanqing/liulab_data/genome/ce11/gtf/wormbase_ws298
  source  https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf.gz
  sha256  f6e4569453a80908e1182039e9af5fc0184d2607db582e57788e979b40b8ca81
  files   wormbase_ws298.db, wormbase_ws298.gtf
  chromosomes checked — every name the GTF uses is one the assembly carries
```

The first use downloads. The lab's CPU compute nodes have no internet, so run `genome annotation register hg38 gencode_v50` once on a login node before submitting a job that needs it.

## What is registered, offered or broken

`worm.annotations` holds every annotation this assembly has. Three lists say where each one
stands:

- `registered` — built on this machine and trustworthy
- `offered` — listed in the annotation table for this assembly, whether or not anyone here
  has built it
- `broken` — a directory that is here and cannot be trusted; each entry says what is wrong
  with it and names the command that repairs it

```python
worm.annotations.registered                       # ['wormbase_ws298']
worm.annotations.broken                           # []
[row.name for row in worm.annotations.offered]    # ['wormbase_ws298']
```

A name in none of the three is the fourth state: not here, and not offered either. There is
no `len()`, no `in` and no iterating over `worm.annotations`, so name the list you mean.

`genome annotation list` prints what is here against what the table offers:

```console
$ genome annotation list ce11
annotations for ce11 in /Users/hanqing/liulab_data/genome/ce11
  wormbase_ws298  registered  WormBase WS298
default: wormbase_ws298
```

`path` gives the GTF's location, ready to hand to another tool:

```python
worm.annotations.path("wormbase_ws298")
# PosixPath('/Users/hanqing/liulab_data/genome/ce11/gtf/wormbase_ws298/wormbase_ws298.gtf')
```

Asking for an annotation that is not registered raises `AnnotationNotRegisteredError`. The
message names the command that repairs it.

## Registering your own GTF

For a GTF the table does not list, hand over the path and a name to address it by:

```python
worm.annotations.register_path("/data/my_genes.gtf.gz", "my_genes")
```

A `.gz` source is decompressed on the way in and a plain GTF is copied as it is. No checksum
is compared, because an unlisted GTF has none pinned for it. Everything else matches a
listed annotation: the same directory, the same chromosome check, the same database, the
same record. `my_genes` is addressed the same way from then on, and `genome annotation list`
shows it as registered but not offered.

The shell spelling takes the assembly, the file and the name, in that order:

```console
$ genome annotation register-gtf ce11 /data/my_genes.gtf my_genes
```

## The default annotation

`default_gtf` is what a call gets when it names no annotation: the `default_gtf=` passed at
construction, else the one the table flags for this assembly, else the sole registered
annotation, else `None`.

```python
worm.default_gtf    # 'wormbase_ws298'
worm.default_gtf_path
# PosixPath('/Users/hanqing/liulab_data/genome/ce11/gtf/wormbase_ws298/wormbase_ws298.gtf')
```

Naming a default does not register it. Opening a genome starts no download, so on a machine
where nothing has been registered yet the default points at an annotation that is not there.
That is the usual state of a fresh machine and not an error. `genome annotation list` reports
the gap without raising:

```console
$ genome annotation list hg38
annotations for hg38 in /Users/hanqing/liulab_data/genome/hg38
  gencode_v50  offered, not registered  GENCODE v50
default: gencode_v50 — not registered here; register it with `genome annotation register hg38 gencode_v50`
```

`default_gtf_path` resolves that name to a file, so asking for a default that is not
registered raises `AnnotationNotRegisteredError`. To pick a different default, pass it at
construction:

```python
Genome("ce11", default_gtf="my_genes")
```

## When a registration is refused

**Every chromosome name the GTF uses has to be one the assembly carries.** The check is one
streaming pass over the GTF and runs before the database build, so a mismatch fails in
seconds. An Ensembl-spelled GTF against a UCSC-spelled assembly is the common case:

```python
sacCer3 = Genome("sacCer3")
sacCer3.annotations.register_path("ensembl.gtf", "ensembl_spelled")
# ChromosomeMismatchError: the GTF for 'ensembl_spelled' names 2 chromosomes the
# assembly does not carry: 1, MT. An annotation and its assembly must spell
# chromosomes the same way, and the usual cause is a UCSC-versus-Ensembl mismatch
# ('chr1' against '1', 'chrM' against 'MtDNA'). The assembly carries: chrI, chrII,
# chrIII, chrIV, chrIX, chrM, chrV, chrVI, chrVII, chrVIII (and 7 more). Register
# the annotation built for this assembly, or pass check_chromosomes=False —
# --no-check-chromosomes from a shell — to register this one anyway.
```

Nothing was created, and the annotation directory is left as it was found. The other
direction is fine: an assembly may carry scaffolds the annotation never mentions.

Gene and transcript inference is off by default. GENCODE, Ensembl, RefSeq and WormBase GTFs
declare `gene` and `transcript` features already, and rebuilding them from exon lines is
gffutils' slow path. Pass `disable_infer_genes=False`, or `--infer-genes` from a shell, for a
bare exon-level GTF, which otherwise registers as a database of exons and nothing else.

## Gene categories

A **gene category** is a named set of an annotation's own gene ids. Ask for one by name and
the ids come back in the spelling that annotation uses:

```python
rrna = worm.gene_list("rRNA")
rrna.gene_ids[:3]     # ['WBGene00004512', 'WBGene00004513', 'WBGene00004567']
len(rrna.gene_ids)    # 23
```

`gene_lists()` says what may be asked for. Nothing in the package holds a vocabulary of
categories; they are whatever the curated list for that annotation declares, and today every
one of them declares `rRNA` and nothing else:

```python
[answer.category for answer in worm.gene_lists()]    # ['rRNA']
```

`rRNA` holds everything rRNA-derived that annotation carries: the mature genes, the
pseudogene copies, the mitochondrial rRNAs. It is built for counting rRNA-derived reads as a
QC metric rather than for describing rRNA biology, so a pseudogene copy that captures reads
counts the same as a functional gene. Each answer carries a description of what it holds and
where the ids came from. The wording differs per annotation, so read it before using the ids
for a measurement. The description below is trimmed after the first two clauses:

```python
rrna.sources[0].description
# 'Every rRNA-derived gene WormBase WS298 annotates on ce11: the chrI 45S unit as
# rrn-1.1 and rrn-1.2 (18S), rrn-2.1 (5.8S) and rrn-3.1 (26S), the sixteen rrn-4.*
# 5S genes, ...'
```

**Features in a category may overlap**, so count over the union of intervals rather than
summing per-gene counts. In `ce11`, `rrn-3.56` is 851 of 852 bases identical to the 3' end of
`rrn-3.1`, and the chrI rDNA array is collapsed to about one and a half of the roughly 100
copies the real genome carries, so reads from every copy land on these few features.

The ids are curated and ship inside the package. They are not read from the GTF's own biotype
attribute, which yeast's `ensgene_v101` does not carry at all, so deriving the list yourself
from that attribute returns no rRNA genes for yeast and no warning either.

A `GeneList` carries more than the ids: the assembly and annotation it answers for, and one
`GeneListSource` per contributing curated list, so a chimera's genes stay attributable to the
component they came from. Both classes are in the [API reference](../reference.md).

### When an annotation cannot answer

There is no empty answer here. An annotation nothing ships a curated list for, and one whose
list does not declare the category asked for, are different facts and each raises its own
error. Neither can be checked ahead of the call, because `gene_lists()` raises the first of
the two itself:

```python
from genome import GeneCategoryNotDeclaredError, NoGeneCategoriesError

try:
    genes = worm.gene_list("tRNA").gene_ids
except NoGeneCategoriesError:
    print("no curated list ships for this annotation, so nothing can be asked of it")
except GeneCategoryNotDeclaredError as error:
    print("declared here:", ", ".join(error.declared))
# declared here: rRNA
```

`ce11` takes the second branch: a curated list ships for `wormbase_ws298` and it declares
`rRNA` only. Both errors are `LookupError`s, so `except LookupError` catches the pair where
the two need no telling apart. No declared category is ever empty and neither call ever
returns an empty collection.
