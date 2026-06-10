"""STAR aligner — genome index construction.

`STAR <https://github.com/alexdobin/STAR>`_ is a splice-aware RNA-seq aligner.
Its genome index is built with ``STAR --runMode genomeGenerate``, which writes a
directory of binary files (the *genomeDir*) that STAR later loads for mapping.

:class:`STAR` exposes only the handful of ``genomeGenerate`` options that are
tuned in practice; every other STAR flag is reachable through ``**kwargs``.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genome.aligner.aligner import Aligner

if TYPE_CHECKING:
    from genome.genome import Genome


class STAR(Aligner):
    """STAR aligner index builder.

    A STAR index is splice-junction-aware: it is built against one gene
    annotation, so each annotation gets its own *genomeDir*. The bound GTF key
    selects that annotation (its path resolved via
    :meth:`~genome.genome.Genome.get_gtf_path`) and names the index directory
    ``star_<gtf_key>``. STAR's index is the *genomeDir* directory itself, so
    :attr:`index_path` returns :attr:`index_dir`.

    Parameters
    ----------
    genome : genome.genome.Genome
        The genome whose reference FASTA will be indexed.
    gtf : str
        Name of a GTF annotation registered on ``genome`` (see
        :meth:`~genome.genome.Genome.register_gtf`).
    """

    name = "star"
    binary = "STAR"

    def __init__(self, genome: Genome, *, gtf: str) -> None:
        self._gtf_key = gtf
        super().__init__(genome)

    @property
    def index_dir(self) -> Path:
        """Per-annotation genome directory ``.../index/star_<gtf_key>/``."""
        base = super().index_dir
        return base.with_name(f"{base.name}_{self._gtf_key}")

    def install_instructions(self) -> str:
        """Return how to install STAR (bioconda)."""
        return (
            "STAR is not installed. Install it from bioconda, e.g.:\n"
            "    pixi add star            # into the project environment\n"
            "See https://github.com/alexdobin/STAR for details."
        )

    def _detect_version(self) -> str:
        """Return the version reported by ``STAR --version`` (e.g. ``2.7.11b``)."""
        result = subprocess.run(
            [self._executable, "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
        return (result.stdout or result.stderr).strip()

    @property
    def _artifact(self) -> Path:
        """STAR loads the genome directory directly."""
        return self.index_dir

    def index(
        self,
        *,
        sjdb_overhang: int = 100,
        threads: int = 1,
        overwrite: bool = False,
        **kwargs: Any,
    ) -> Path:
        """Build the STAR genome index for the bound assembly and annotation.

        Output goes to ``<LIULAB_DATA>/genome/<assembly>/index/star_<gtf_key>/``.
        When a successful index already exists it is reused unless
        ``overwrite=True``. The annotation GTF is resolved from the bound ``gtf``
        key via :meth:`~genome.genome.Genome.get_gtf_path` and passed to STAR as
        ``--sjdbGTFfile`` for splice-junction-aware indexing.

        Only the most commonly tuned options are named below. Any other STAR
        ``genomeGenerate`` option may be passed as a keyword argument using its
        STAR name without the leading ``--`` (e.g. ``genomeSAindexNbases=11``);
        for the meaning of those options see the STAR manual / ``STAR --help``.

        Parameters
        ----------
        sjdb_overhang : int, default 100
            ``--sjdbOverhang``: ideally ``read_length - 1``.
        threads : int, default 1
            ``--runThreadN``: number of threads to build with.
        overwrite : bool, default False
            Rebuild even if a successful index already exists.
        **kwargs : Any
            Extra ``genomeGenerate`` options forwarded verbatim as STAR flags.

        Returns
        -------
        pathlib.Path
            The genome directory (also available as :attr:`index_path`).
        """
        if self._flag_path.is_file() and not overwrite:
            return self.index_path

        self.index_dir.mkdir(parents=True, exist_ok=True)
        self._flag_path.unlink(missing_ok=True)  # invalidate any stale marker

        fasta = self._genome.files.fasta
        # STAR runs with the index dir as CWD, so resolve the annotation path.
        gtf_file = self._genome.get_gtf_path(self._gtf_key).resolve()

        # STAR requires a reduced suffix-array index size for small genomes:
        #   min(14, log2(genomeLength) / 2 - 1). Default it unless overridden.
        if "genomeSAindexNbases" not in kwargs:
            genome_length = int(self._genome.chrom_sizes.sum())
            kwargs["genomeSAindexNbases"] = min(14, max(2, int(math.log2(genome_length) / 2 - 1)))

        parameters: dict[str, Any] = {
            "threads": threads,
            "gtf": self._gtf_key,
            "sjdb_gtf_file": str(gtf_file),
            "sjdb_overhang": sjdb_overhang,
            **kwargs,
        }

        args: list[str] = [
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            str(self.index_dir),
            "--genomeFastaFiles",
            str(fasta),
            "--runThreadN",
            str(threads),
            "--sjdbGTFfile",
            str(gtf_file),
            "--sjdbOverhang",
            str(sjdb_overhang),
        ]
        args += _kwargs_to_flags(kwargs)

        self._run(args)
        self._write_metadata(command=[self.binary, *args], parameters=parameters)
        self._mark_success()
        return self.index_path


def _kwargs_to_flags(kwargs: dict[str, Any]) -> list[str]:
    """Turn ``{"genomeSAindexNbases": 11}`` into ``["--genomeSAindexNbases", "11"]``.

    List/tuple values become multiple space-separated arguments after the flag.
    """
    flags: list[str] = []
    for key, value in kwargs.items():
        flag = f"--{key}"
        if isinstance(value, (list, tuple)):
            flags += [flag, *(str(item) for item in value)]
        else:
            flags += [flag, str(value)]
    return flags
