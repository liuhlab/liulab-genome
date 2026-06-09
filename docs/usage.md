# Usage

## Installation

The package depends on native tools from bioconda (`samtools`, `bedtools`), so
the recommended path is [pixi](https://pixi.sh):

```bash
git clone https://github.com/lhqing/liulab-genome.git
cd liulab-genome
pixi install --locked        # solve from the committed pixi.lock
pixi shell                   # activate the env in your shell
```

If you only need the Python API and already have `samtools` / `bedtools` on
your `PATH` from somewhere else, a pip install also works:

```bash
pip install liulab-genome
```

You can verify the native toolchain is reachable at any time:

```bash
$ genome doctor
samtools: samtools 1.21 ...
bedtools: bedtools v2.31.1
```

`genome doctor` exits with code `1` and an actionable message if either tool
is missing.

## Command-line interface

The `genome` command is a thin Typer wrapper over the Python API — every
command has a corresponding function in the package, so scripts and notebooks
can call the same code paths without shelling out.

```bash
$ genome --help
$ genome version
$ genome revcomp ATCG
$ genome revcomp aTcG --json
$ genome doctor [--json]
```

`--json` toggles machine-readable output. Errors go to stderr and use non-zero
exit codes (`2` for invalid input, `1` for missing native tools).

## Python API quickstart

```python
from genome import DNA, RNA, Protein

# Construction preserves the value verbatim (including case). The alphabet is
# NOT enforced — scanning every char is too costly on large sequences.
seq = DNA("aTcGN")  # accepted and stored as-is; no validation
seq = DNA("aTcG")   # OK

# Biological transforms return the right typed result.
seq.complement()                # DNA('tAgC')
seq.reverse_complement()        # DNA('CgAt')
seq.transcribe()                # RNA('aUcG')
seq.gc_content                  # 0.5

# Slicing stays typed.
seq[1:3]                        # DNA('Tc')

# Other str methods (upper, lower, replace, +) return plain str by design;
# wrap explicitly if you want a typed result back.
DNA(seq.upper())                # DNA('ATCG')
```

See [Sequences](sequences.md) for the full contract, including what happens
on invalid input and which `str` methods do or do not preserve the type.

## Running the gates locally

The repository ships with a single `pixi run check` that runs lint, format
check, type check, and the full test suite — the same commands CI runs:

```bash
pixi run check
pixi run test-cov              # tests + coverage report
pixi run docs                  # live MkDocs server at http://127.0.0.1:8000
```

To enable the pre-commit hooks in your local clone (one-off):

```bash
pixi run -- pre-commit install
```
