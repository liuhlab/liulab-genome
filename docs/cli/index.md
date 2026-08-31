# CLI overview

Installing the package puts a `genome` command on `PATH`. Every command calls the same
function the Python API calls, so a shell script and a notebook take one code path.

```console
$ genome annotation list sacCer3
annotations for sacCer3 in /Users/hanqing/liulab_data/genome/sacCer3
  ensgene_v101  registered  UCSC ensGene.v101
default: ensgene_v101
```

[CLI: genome commands](genome.md) covers everything that names an assembly, and
[CLI: lookup commands](lookups.md) covers the rest, which name a species or a file instead.

## The command tree

Three commands belong to no topic. Everything else hangs off a group named for the part of
the package it ships from:

| Group | Commands |
|---|---|
| `genome assembly` | `register`, `list`, `verify`, `table-row` |
| `genome annotation` | `register`, `register-gtf`, `list`, `gene-list`, `gene-categories` |
| `genome tf` | `gene-list`, `cofactor-list` |
| `genome xref` | `ids`, `symbols` |
| `genome homology` | `links` |
| `genome motif` | `scan` |

`genome --help` prints the same tree, and `--help` on a group lists that group's commands.

Building a STAR or chromap index has no command here. It is a Python call, and
[Aligner indexes](../genome/aligner.md) is the page for it.

## JSON output

Every command takes `--json`, which replaces the text output with one JSON object:

```console
$ genome version --json
{"version": "2026.8.1.dev13+g476734570.d20260830"}
```

Six commands answer with a list: `annotation gene-list`, `tf gene-list`, `tf
cofactor-list`, `xref ids`, `xref symbols` and `homology links`. Those send the ids or the
pairs to stdout and the heading, the attribution and the counts to stderr, so `> file.tsv`
captures the answer and leaves the commentary on screen. `--json` puts both in one object.

## Exit codes

`0` is success. `2` is invalid input, refused before anything is downloaded or read. `1` is
everything else: a missing native tool, a failed download, a checksum mismatch, a lookup
that found nothing.

Errors go to stderr with an `error: ` prefix, and the message names the next action rather
than only the problem:

```console
$ genome annotation register-gtf ce11 /nonexistent.gtf wormbase_ws298
error: GTF file not found: /nonexistent.gtf. Pass the path of an existing .gtf or .gtf.gz, or register a listed annotation by name instead.
```

## genome version

Print the installed version.

```console
$ genome version
2026.8.1.dev13+g476734570.d20260830
```

## genome revcomp

Reverse-complement a DNA sequence typed on the command line. Case is preserved, so a
lower-case base comes back lower-case.

```console
$ genome revcomp ATCG
CGAT
$ genome revcomp aTcG --json
{"input": "aTcG", "reverse_complement": "CgAt"}
```

A character outside `ACGT` exits `2` and names the offending character. It is for a primer
you are holding. A file of sequence is Python work, on
[Sequences and regions](../genome/sequences.md).

## genome doctor

Report which native tools are on `PATH`, and at what version. It exits `1` naming the
install command when one of them is missing.

```console
$ genome doctor
samtools: samtools 1.22.1
faToTwoBit: installed; reports no version
twoBitInfo: installed; reports no version
```

A tool that is installed but refuses `--version`, as the UCSC binaries do, is listed anyway.
Presence is the question. **STAR and chromap are not checked here.** They are optional, and
each one checks for itself when you ask it to build an index.

## The old flat spellings

Fourteen commands used to hang off the root rather than off a group: `genome verify`,
`genome tf-gene-list`, `genome motif-scan` and the rest. Each still runs and calls the same
function its group command calls, hidden from `genome --help` and printing a deprecation
notice on stderr, so `--json` on stdout still parses for a script that has not moved yet.
**They are removed in the next release.**
