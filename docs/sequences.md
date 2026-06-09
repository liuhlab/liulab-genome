# Typed biological sequences

`genome.seq` provides three thin, typed wrappers around `str`:

| Class | Alphabet | Transforms |
|-------|----------|------------|
| `DNA` | `A C G T` | `complement`, `reverse_complement`, `transcribe`, `gc_content` |
| `RNA` | `A C G U` | `complement`, `reverse_complement`, `back_transcribe`, `gc_content` |
| `Protein` | 20 standard amino acids (`ACDEFGHIKLMNPQRSTVWY`) | _(none yet)_ |

These classes subclass `str`, so a `DNA` **is** a string anywhere a string is expected — but
construction validates the alphabet, slicing keeps the type, and the biological methods return
the right typed result with case preserved.

> **IUPAC ambiguity codes (`N`, `R`, `Y`, …) are intentionally out of scope.** Add them only
> when the project explicitly needs them.

## Construction and validation

Validation runs in `__new__` (because `str` is immutable), is **case-insensitive**, and
**preserves case** in the stored value. Lower case is meaningful — it is the conventional
soft-masking marker — so `aTcG` is stored verbatim, not normalised.

```python
from genome import DNA, RNA, Protein

DNA("ATCG")        # → DNA('ATCG')
DNA("aTcG")        # → DNA('aTcG')   (case preserved)
DNA("")            # → DNA('')       (empty is valid)
```

An invalid character raises `ValueError` naming the offenders:

```python
>>> DNA("ATCX")
Traceback (most recent call last):
    ...
ValueError: DNA contains characters outside alphabet {ACGT}: ['X']

>>> RNA("ATCG")    # T isn't part of RNA
Traceback (most recent call last):
    ...
ValueError: RNA contains characters outside alphabet {ACGU}: ['T']
```

The private base class `_Seq` is abstract and can't be instantiated:

```python
>>> from genome.seq import _Seq
>>> _Seq("ATCG")
Traceback (most recent call last):
    ...
TypeError: _Seq is abstract; instantiate DNA, RNA, or Protein instead.
```

## Slicing stays typed

`__getitem__` is overridden so indexing and slicing return the same subclass — not a bare `str`.

```python
>>> s = DNA("ATCGATCG")
>>> s[2:5]
DNA('CGA')
>>> s[0]
DNA('A')
>>> type(s[2:5]).__name__
'DNA'
```

## Biological transforms

Each transform returns a new typed instance of the appropriate class with case preserved.

```python
>>> DNA("ATCG").complement()
DNA('TAGC')
>>> DNA("aTcG").complement()
DNA('tAgC')

>>> DNA("ATCG").reverse_complement()
DNA('CGAT')
>>> DNA("aTcG").reverse_complement()
DNA('CgAt')

>>> DNA("ATCG").transcribe()
RNA('AUCG')
>>> RNA("AUCG").back_transcribe()
DNA('ATCG')
```

Empty sequences are well-behaved everywhere:

```python
>>> DNA("").reverse_complement()
DNA('')
```

## GC content

`gc_content` is a `float` in `[0.0, 1.0]`, defined as `0.0` for the empty sequence.

```python
>>> DNA("GGCC").gc_content
1.0
>>> DNA("ATAT").gc_content
0.0
>>> DNA("aTcG").gc_content   # case-insensitive counting
0.5
>>> DNA("").gc_content
0.0
```

## What _doesn't_ stay typed (by design)

Only `__getitem__` and the explicit biological methods preserve the subclass. All other
inherited `str` methods (`upper`, `lower`, `replace`, `strip`, `split`, `+`, …) return plain
`str`:

```python
>>> type(DNA("aTcG").upper()).__name__
'str'
```

If you want a typed variant, wrap the result yourself:

```python
>>> DNA(DNA("aTcG").upper())
DNA('ATCG')
```

This is deliberate — silently overriding every inherited method would surprise callers and
require validation on results that `str` makes no promises about. Slicing is the one
exception because the typed return is unambiguously safe (a slice of valid DNA is valid DNA).

## Reverse-complementing from the command line

The CLI exposes a `revcomp` command — a one-line wrapper over `DNA.reverse_complement`:

```bash
$ genome revcomp ATCG
CGAT

$ genome revcomp aTcG --json
{"input": "aTcG", "reverse_complement": "CgAt"}
```

Invalid input exits with code `2` and an error message naming the offending characters.

## Checking the native toolchain

The package shells out to `samtools` and `bedtools` for the bulk of file-format work. Verify
they're on `PATH` (resolved by pixi) with:

```bash
$ genome doctor
samtools: samtools 1.21 ...
bedtools: bedtools v2.31.1
```

Or programmatically:

```python
>>> from genome.external import doctor
>>> doctor()
{'samtools': 'samtools 1.21 ...', 'bedtools': 'bedtools v2.31.1'}
```

A missing tool raises `ToolNotFoundError` with an actionable message pointing at `pixi add`.
