# liulab-genome

Handling genomic files — metadata management, file processing, feature extraction.

Import name: `genome`.

## Status

Early scaffolding. The first real feature (typed biological sequences in `genome.seq`)
lands in Phase 3.

## Development

This project uses [pixi](https://pixi.sh) with `conda-forge` + `bioconda` channels. Native
deps (`samtools`, `bedtools`) and Python tooling are all managed by pixi.

```bash
pixi install            # solve & install the default env (resolves from pixi.lock if present)
pixi shell              # activate the env
pixi run check          # lint + fmt-check + typecheck + test (the CI gate)
```

See [`CLAUDE.md`](./CLAUDE.md) for the full contributor/agent working agreement.

## License

MIT — see [`LICENSE`](./LICENSE).
