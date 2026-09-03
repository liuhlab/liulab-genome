# Sequences and regions

A prepared genome answers by locus: name a chromosome and a range of coordinates, get the
bases. The result is a `DNA`, a `str` subclass carrying the biological transforms, so
complementing it or measuring its GC content is a method call away.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.fetch_sequence("chrIV:0-10")    # DNA('ACACCACACC')
```

Preparing an assembly so there is something to read is on [Assembly](assembly.md).

## Coordinates are 0-based and half-open

Every coordinate in this package is 0-based and half-open. **The first base of a
chromosome is at position `0`, and `end` is one past the last base you want.** So
`chrIV:0-10` is the first ten bases:

```python
len(sacCer3["chrIV:0-10"])    # 10
sacCer3["chrIV:0-1"]          # DNA('A')
```

This is the convention BED files use. Genome browsers, GTF and GFF files, VCF and SAM all
count from 1 and include both ends, so a feature shown as `chrIV:1-10` in IGV is
`chrIV:0-10` here. The conversion happens where those files are read and written. Nothing
you pass to `fetch_sequence`, and nothing you build a `Region` from, is ever 1-based.
Thousands separators are tolerated, so `chrIV:1,000-1,010` and `chrIV:1000-1010` are the
same request.

## Fetching sequence

`fetch_sequence` takes a locus string, a bare chromosome name, or a `Region`, and returns
a `DNA`. Indexing the genome is sugar for the same call:

```python
sacCer3.fetch_sequence("chrIV:0-10")    # DNA('ACACCACACC')
sacCer3["chrIV:0-10"]                   # DNA('ACACCACACC')
```

A name with no coordinates means the whole sequence:

```python
len(sacCer3["chrM"])    # 85779
```

Because the result is a `DNA`, the transforms are already on it:

```python
sacCer3["chrIV:0-10"].reverse_complement()    # DNA('GGTGTGGTGT')
sacCer3["chrIV:0-1000"].gc_content            # 0.435
```

**Lower-case bases come back as they are.** Lower case is soft-masking, the repeat
annotation the FASTA carries, and nothing upper-cases it for you.

An `end` past the chromosome length raises rather than truncating quietly. An `end` equal
to the length is valid, because the end is exclusive:

```python
sacCer3["chrIV:0-999999999"]
# ValueError: region chrIV:0-999999999: end (999999999) exceeds chrIV length (1531933).
# Coordinates are 0-based half-open, so the maximum valid end is 1531933.
```

A chromosome name this reference does not carry, and a locus that is not
`chrom:start-end`, both raise `ValueError` too.

## Strand

A locus string carries no strand, so it is read on the forward strand. Pass a `Region`
with strand `"-"` and you get the reverse complement of that interval:

```python
from genome import Region

sacCer3[Region("chrIV", 1000, 1010, "+")]    # DNA('TCTATAGTCA')
sacCer3[Region("chrIV", 1000, 1010, "-")]    # DNA('TGACTATAGA')
```

**Strand `"."` means unknown, and nothing promotes it to `"+"`.** Fetching with `"."`
returns the forward-strand bases, which is what the reference stores, but it makes no
claim about which strand the feature is on. Set the strand when it matters.

## Regions

`Region` is the coordinate type the rest of the package passes around: frozen, validated,
0-based half-open, with an explicit strand.

```python
r = Region("chrIV", 0, 10)    # Region(chrom='chrIV', start=0, end=10, strand='.')
len(r)                        # 10
str(r)                        # 'chrIV:0-10'
```

`len(r)` and `r.length` are the same number, `end - start`. `str(r)` gives back the locus
string in the coordinates it stores, so it goes straight back into `fetch_sequence`.

Construction enforces `start >= 0`, `end >= start`, and a `strand` of `"+"`, `"-"` or
`"."`; anything else raises `ValueError`. The class is frozen, so a `Region` you hand to a
function comes back unchanged.

`Region.from_string` builds one from a locus string. The strand is a separate argument,
since the string carries none:

```python
Region.from_string("chrIV:100-200", strand="-")
# Region(chrom='chrIV', start=100, end=200, strand='-')
```

To take a locus string apart without building a `Region`, `genome.region.parse_region`
returns the three pieces. A bare chromosome name gives `None` for both coordinates, which
is how a caller knows to resolve it against a table of chromosome lengths:

```python
from genome.region import parse_region

parse_region("chrIV:1,000-2,000")    # ('chrIV', 1000, 2000)
parse_region("chrM")                 # ('chrM', None, None)
```

## Typed sequences

`DNA` and `RNA` subclass `str`. A `DNA` is a string anywhere a string is expected, and it
carries the transforms besides.

| Class | Alphabet | Transforms |
| --- | --- | --- |
| `DNA` | `A C G T` | `complement`, `reverse_complement`, `transcribe`, `gc_content` |
| `RNA` | `A C G U` | `complement`, `reverse_complement`, `back_transcribe`, `gc_content` |

Build one from any string. The value is stored verbatim, case included:

```python
from genome import DNA, RNA

DNA("ATCG")         # DNA('ATCG')
DNA("aTcG")         # DNA('aTcG')
RNA("AUCG")         # RNA('AUCG')
```

**The alphabet is not checked at construction.** Scanning every character costs too much
on a whole chromosome, so `DNA("ATCX")` builds and gives you `DNA('ATCX')`. IUPAC
ambiguity codes are out of scope, and an `N` run from a reference FASTA arrives as
ordinary characters.

### Checking the alphabet

Reject bad input where it enters, and ask the class rather than spelling out `ACGT`
yourself. `outside_alphabet` returns the offending characters, distinct and sorted, and an
empty list when nothing offends:

```python
DNA.outside_alphabet("ATCX")    # ['X']
DNA.outside_alphabet("aTcG")    # []
RNA.outside_alphabet("ATCG")    # ['T']
```

Case is not an offence, since lower case is soft-masking. `DNA.ALPHABET` is the frozenset
itself, for an error message that has to name it.

### Slicing stays typed

Indexing and slicing return the same subclass:

```python
DNA("ATCGATCG")[2:5]     # DNA('CGA')
RNA("AUCGAUCG")[1:3]      # RNA('UC')
```

Nothing else does. `upper`, `lower`, `replace`, `+` and every other inherited `str` method
return a plain `str`. Wrap the result to get the type back:

```python
DNA(DNA("aTcG").upper())    # DNA('ATCG')
```

## Transforms

Each transform returns a new typed value and preserves case:

```python
DNA("aTcG").complement()            # DNA('tAgC')
DNA("aTcG").reverse_complement()    # DNA('CgAt')
DNA("ATCG").transcribe()            # RNA('AUCG')
RNA("AUCG").back_transcribe()       # DNA('ATCG')
```

`gc_content` is the fraction of bases that are `G` or `C`, counted case-insensitively, and
`0.0` for an empty sequence:

```python
DNA("aTcG").gc_content    # 0.5
```

There is no translation to protein: protein sequences live in `liulab-protein`, not here.
The full surface of both types is in the [API reference](../reference.md).

## From the command line

`genome revcomp` reverse-complements its argument, with the same case handling:

```console
$ genome revcomp aTcG
CgAt
```

It rejects a character outside `DNA.ALPHABET` and names it. The rest of the command line
is on the [CLI overview](../cli/index.md).
