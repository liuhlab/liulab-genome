# CLI: genome commands

Eight commands in two groups: one prepares a reference on this machine, the other registers
annotations against it and reads what they carry. Both take an assembly name as their first
argument.

```console
$ genome assembly register sacCer3
registered sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3
  source  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
  sha256  6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3
  files   sacCer3.2bit, sacCer3.chrom.sizes, sacCer3.fa, sacCer3.fa.fai
```

The files land in the shared [data directory](../index.md#the-data-directory), where every project on the machine reads the same copy.

Each synopsis below is its group's own `--help`, so it leaves the leading `genome` off:
`assembly register sacCer3` is run as `genome assembly register sacCer3`. `--json`, exit
codes and the split between stdout and stderr are on the [CLI overview](index.md).

## genome assembly

`register` prepares an assembly and prints where it landed, as above. `verify` re-reads the
FASTA and recomputes its sha256, which is the full-file check; registering and reopening go
by size and are instant. `table-row` downloads an assembly and prints the line to paste into
the shipped metadata table, which is how a new assembly gets a pinned checksum.
[Assembly](../genome/assembly.md) says what the four files hold and what a chimera is.

```console
$ genome assembly verify sacCer3
/Users/hanqing/liulab_data/genome/sacCer3/sacCer3.fa: sha256 6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3 matches the digest the metadata table pins for it
```

Registering an assembly that is already here downloads nothing and prints the same summary
off its record. A directory whose files changed after they were registered exits `1`
instead. **`--force` is the repair**, and it keeps whatever is provably good: a FASTA whose
checksum still matches is reused and only the derived files are rebuilt.

The three commands, with every argument and option:

::: mkdocs-typer2
    :module: genome.assembly.cli
    :name: assembly
    :engine: native

## genome annotation

`register` fetches one of the annotations the table lists for an assembly and builds its
gffutils database; `register-gtf` does the same for a GTF of your own, which nothing is
pinned for. `list` puts what is offered beside what is registered here, and downloads
nothing to answer, so it works for an assembly you have never prepared. `gene-list` and
`gene-categories` read the curated gene categories that ship with the package. What a
category is and where it comes from is on [Annotations](../genome/annotations.md).

```console
$ genome annotation list ce11
annotations for ce11 in /Users/hanqing/liulab_data/genome/ce11
  wormbase_ws298  registered  WormBase WS298
default: wormbase_ws298
```

The last line of a registration always says whether the GTF's chromosome names were checked
against the assembly's, and when they were not, why. **Registering an annotation before its
assembly is allowed**, and reports that there was no `chrom.sizes` to check against;
registering the assembly is what fixes that.

`gene-list` prints ids and nothing else on stdout, so the output pipes or redirects. The
heading and the per-source counts go to stderr:

```console
$ genome annotation gene-list ce11 rRNA > rrna.txt
rRNA for ce11 / wormbase_ws298
  wormbase_ws298  23
```

Asking for a category the annotation does not declare exits `1`, and so does asking an
annotation that is not registered here. Neither prints an empty list. `genome annotation
gene-categories ce11` says which categories may be asked for.

Every argument and option of the five:

::: mkdocs-typer2
    :module: genome.annotation.cli
    :name: annotation
    :engine: native
