# 1. The environment manager is pixi, not uv or pip

`samtools` and `bedtools` are runtime dependencies and are not Python packages, so no PyPI-only
resolver — uv included — can install them; that is the whole reason the modern default loses here.
Channels are `conda-forge` then `bioconda` and the order is priority, the manifest is `[tool.pixi.*]`
in `pyproject.toml`, and `pixi.lock` is committed. What it costs: a smaller ecosystem, and a slower
and less familiar tool than uv.
