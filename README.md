# liulab-genome

Reference genomes on disk, so a script does not have to keep them. Name an assembly and the
package fetches its FASTA, derives the companion files other tools expect (`.fai`, `.2bit`,
`chrom.sizes`), and hands back an object that reads sequence out of it. It also registers
GTF annotations against an assembly, builds STAR and chromap indexes, and carries published
tables of transcription factors, motifs, gene identifiers and orthologs.

The import name is `genome`.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.fetch_sequence("chrIV:0-10")    # DNA('ACACCACACC')
sacCer3.chrom_sizes["chrIV"]            # 1531933
```

Coordinates in the Python API are 0-based and half-open, so `chrIV:0-10` is the first ten
bases. The same work runs from a shell, starting with `genome assembly register sacCer3`.

Docs: <https://liuhlab.github.io/liulab-genome/>

## Development

The work is done by native binaries, so [pixi](https://pixi.sh) is the supported path. One
lock file brings the Python package, `samtools`, `bedtools` and the rest of the tooling
together from `conda-forge` and `bioconda`.

```bash
pixi install --locked              # install the default env from pixi.lock
pixi shell                         # activate it
pixi run check                     # the gate: lint, fmt-check, typecheck, test, run concurrently
pixi run -e aligners test-aligner  # the other lane: the tests that build a real STAR/chromap index
```

STAR and chromap are large and most work never touches them, so they live in a second pixi
environment named `aligners`. `pixi run check` skips those tests, so run both lanes.

See [`AGENTS.md`](./AGENTS.md) (`CLAUDE.md` symlinks to it) for the full contributor and
agent working agreement.

## License

MIT. See [`LICENSE`](./LICENSE).
