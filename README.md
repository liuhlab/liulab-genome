# liulab-genome

Reference genomes on disk, ready to query. Name an assembly and it is fetched, prepared
(`.fai`, `.2bit`, `chrom.sizes`) and answering sequence queries — plus GTF annotations
registered against it and STAR/chromap indexes built from it.

Import name: `genome`.

```python
from genome import Genome

sacCer3 = Genome("sacCer3")
sacCer3.fetch_sequence("chrIV:0-10")   # DNA('ACACCACACC') — 0-based, half-open
```

Docs: <https://liuhlab.github.io/liulab-genome/>

## Development

This project uses [pixi](https://pixi.sh) with `conda-forge` + `bioconda` channels. Native
deps (`samtools`, `bedtools`) and Python tooling are all managed by pixi.

```bash
pixi install                     # solve & install the default env (resolves from pixi.lock if present)
pixi shell                       # activate the env
pixi run check                   # lint + fmt-check + typecheck + test, run concurrently (the gate)
pixi run -e aligners test-aligner # the other lane: the tests that build a real STAR/chromap index
```

See [`AGENTS.md`](./AGENTS.md) (`CLAUDE.md` symlinks to it) for the full contributor/agent working agreement.

## License

MIT — see [`LICENSE`](./LICENSE).
