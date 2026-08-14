# Sequences

`genome.seq` gives you three typed wrappers around `str`:

| Class | Alphabet | Transforms |
|-------|----------|------------|
| `DNA` | `A C G T` | `complement`, `reverse_complement`, `transcribe`, `gc_content` |
| `RNA` | `A C G U` | `complement`, `reverse_complement`, `back_transcribe`, `gc_content` |
| `Protein` | the 20 standard amino acids | _(none yet)_ |

They subclass `str`, so a `DNA` **is** a string anywhere a string is expected — and
`Genome.fetch_sequence` hands you one, so the transforms are always in reach.

## Construction

The value is stored verbatim, case included. Lower case is soft-masking and means
something, so it is never normalised away.

```python
from genome import DNA, RNA, Protein

DNA("ATCG")     # DNA('ATCG')
DNA("aTcG")     # DNA('aTcG')  — case preserved
DNA("")         # DNA('')      — empty is valid
```

The alphabet in the table is documentation, not a runtime check: scanning every character
costs too much on a whole chromosome, so any string is accepted. Reject bad input where it
enters — the CLI does this for `genome revcomp`. IUPAC ambiguity codes (`N`, `R`, `Y`, …)
are out of scope.

## Transforms

Each returns a new typed value, case preserved:

```python
>>> DNA("aTcG").complement()
DNA('tAgC')
>>> DNA("aTcG").reverse_complement()
DNA('CgAt')
>>> DNA("ATCG").transcribe()
RNA('AUCG')
>>> RNA("AUCG").back_transcribe()
DNA('ATCG')
```

`gc_content` is a `float` in `[0.0, 1.0]`, counted case-insensitively, and `0.0` for the
empty sequence:

```python
>>> DNA("aTcG").gc_content
0.5
```

## Slicing stays typed

Indexing and slicing return the same subclass:

```python
>>> DNA("ATCGATCG")[2:5]
DNA('CGA')
```

Nothing else does. Every other inherited `str` method (`upper`, `lower`, `replace`,
`split`, `+`, …) returns a plain `str` — wrap it yourself if you want the type back:

```python
>>> DNA(DNA("aTcG").upper())
DNA('ATCG')
```

## From the command line

```console
$ genome revcomp aTcG
CgAt
```

See the [CLI](cli.md) page.
