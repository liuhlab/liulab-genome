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

from functools import partial
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
    genome : genome.assembly.genome.Genome
        The genome whose reference FASTA will be indexed.
    tool : genome.external.ExternalTool, optional
        As :class:`~genome.aligner.aligner.Aligner`.
    """

    name = "chromap"
    binary = "chromap"
    # chromap's long options are hyphenated, so a Python keyword's underscores become
    # hyphens on the way to the command line.
    _flag_separator = "-"

    @property
    def _artifact(self) -> Path:
        """The single index file chromap loads via ``-x/--index``."""
        return self.index_dir / f"{self.name}.index"

    @property
    def _build_arguments(self) -> str:
        """None — one chromap index serves the whole assembly, so nothing selects it."""
        return ""

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
        An index whose completion record says it finished is reused unless
        ``overwrite=True``; a directory holding index files that no record vouches for
        raises rather than being silently rebuilt. chromap needs only the reference
        FASTA — no gene annotation — so one index serves every use of the assembly.

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
            Rebuild even if a finished index already exists, and rebuild over a
            directory that cannot be trusted rather than raising on it.
        **kwargs : Any
            Extra ``--build-index`` options forwarded verbatim as chromap flags.

        Returns
        -------
        pathlib.Path
            The built index file (also available as :attr:`index_path`).

        Raises
        ------
        genome.store.completion.RegistrationError
            If the index directory holds files without a record, a record that
            disagrees with them, or a record pinning a different assembly digest
            than the one registered now. Pass ``overwrite=True`` to rebuild.
        RuntimeError
            If chromap exits non-zero.
        """
        compose = partial(self._compose, kmer=kmer, window=window, **kwargs)
        return self._build(compose, overwrite=overwrite)

    def _compose(
        self, *, kmer: int | None, window: int | None, **kwargs: Any
    ) -> tuple[list[str], dict[str, Any]]:
        """Return the ``--build-index`` command line, and the knobs that determined it.

        Reached only when a build is going to run — see
        :meth:`~genome.aligner.aligner.Aligner._build`.
        """
        fasta = self._genome.files.fasta

        args: list[str] = [
            "--build-index",
            "--ref",
            str(fasta),
            "--output",
            str(self._artifact),
        ]

        # The two named knobs are spelled here and recorded below, rather than rendered
        # out of what is recorded: chromap has no default of this package's choosing, so
        # each reaches the command line only when the caller asked for it.
        parameters: dict[str, Any] = {}
        if kmer is not None:
            parameters["kmer"] = kmer
            args += ["--kmer", str(kmer)]
        if window is not None:
            parameters["window"] = window
            args += ["--window", str(window)]
        parameters.update(kwargs)

        args += self._flags(kwargs)

        return args, parameters
