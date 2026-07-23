"""Chromap aligner — genome index construction.

`Chromap <https://github.com/haowenz/chromap>`_ is a fast aligner and
preprocessor for chromatin profiles (ATAC-seq/scATAC-seq, ChIP-seq, Hi-C). Its
genome index is built with ``chromap --build-index``, which writes a single
index file from the reference FASTA alone — unlike a splice-aware RNA index,
chromap needs no gene annotation, so one index serves every use of an assembly.

:class:`Chromap` exposes only the two minimizer knobs tuned in practice (k-mer
length and window size); every other ``--build-index`` option is reachable
through ``**kwargs``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

from genome.aligner.aligner import Aligner


class Chromap(Aligner):
    """Chromap aligner index builder.

    A chromap index is a single minimizer table built from the reference FASTA
    alone; it carries no annotation, so — unlike
    :class:`~genome.aligner.star.STAR` — there is exactly one index per assembly.
    The index is written to ``.../index/chromap/chromap.index`` and
    :attr:`index_path` returns that file (chromap consumes it via ``-x/--index``).

    Parameters
    ----------
    genome : genome.genome.Genome
        The genome whose reference FASTA will be indexed.
    """

    name = "chromap"
    binary = "chromap"

    def install_instructions(self) -> str:
        """Return how to install chromap (bioconda)."""
        return (
            "chromap is not installed. Install it from bioconda, e.g.:\n"
            "    pixi add chromap         # into the project environment\n"
            "See https://github.com/haowenz/chromap for details."
        )

    def _detect_version(self) -> str:
        """Return the version reported by ``chromap --version`` (e.g. ``0.3.2-r518``)."""
        result = subprocess.run(
            [self._executable, "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        return (result.stdout or result.stderr).strip()

    @property
    def _artifact(self) -> Path:
        """The single index file chromap loads via ``-x/--index``."""
        return self.index_dir / f"{self.name}.index"

    def index(
        self,
        *,
        kmer: int | None = None,
        window: int | None = None,
        overwrite: bool = False,
        **kwargs: Any,
    ) -> Path:
        """Build the chromap index for the bound assembly and return :attr:`index_path`.

        Output goes to ``<LIULAB_DATA>/genome/<assembly>/index/chromap/chromap.index``.
        When a successful index already exists it is reused unless ``overwrite=True``.
        chromap needs only the reference FASTA — no gene annotation — so one index
        serves every use of the assembly.

        Only the two minimizer knobs are named below. Any other ``--build-index``
        option may be passed as a keyword argument using chromap's flag name with
        underscores for hyphens (e.g. ``min_frag_length=30`` -> ``--min-frag-length
        30``); for their meaning see ``chromap --help``.

        Parameters
        ----------
        kmer : int, optional
            ``-k/--kmer``: minimizer k-mer length. chromap's own default is used
            when omitted.
        window : int, optional
            ``-w/--window``: minimizer window size. chromap's own default is used
            when omitted.
        overwrite : bool, default False
            Rebuild even if a successful index already exists.
        **kwargs : Any
            Extra ``--build-index`` options forwarded verbatim as chromap flags.

        Returns
        -------
        pathlib.Path
            The built index file (also available as :attr:`index_path`).
        """
        if self._flag_path.is_file() and not overwrite:
            return self.index_path

        self.index_dir.mkdir(parents=True, exist_ok=True)
        self._flag_path.unlink(missing_ok=True)  # invalidate any stale marker

        fasta = self._genome.files.fasta

        parameters: dict[str, Any] = {}
        if kmer is not None:
            parameters["kmer"] = kmer
        if window is not None:
            parameters["window"] = window
        parameters.update(kwargs)

        args: list[str] = [
            "--build-index",
            "--ref",
            str(fasta),
            "--output",
            str(self._artifact),
        ]
        args += _kwargs_to_flags(parameters)

        self._run(args)
        self._write_metadata(command=[self.binary, *args], parameters=parameters)
        self._mark_success()
        return self.index_path


def _kwargs_to_flags(kwargs: dict[str, Any]) -> list[str]:
    """Turn ``{"min_frag_length": 30}`` into ``["--min-frag-length", "30"]``.

    chromap's long options are hyphenated, so underscores in keyword names become
    hyphens. List/tuple values become multiple space-separated arguments after the
    flag.
    """
    flags: list[str] = []
    for key, value in kwargs.items():
        flag = f"--{key.replace('_', '-')}"
        if isinstance(value, (list, tuple)):
            flags += [flag, *(str(item) for item in value)]
        else:
            flags += [flag, str(value)]
    return flags
