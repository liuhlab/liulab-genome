"""STAR aligner — genome index construction.

`STAR <https://github.com/alexdobin/STAR>`_ is a splice-aware RNA-seq aligner.
Its genome index is built with ``STAR --runMode genomeGenerate``, which writes a
directory of binary files (the *genomeDir*) that STAR later loads for mapping.

:class:`STAR` exposes only the handful of ``genomeGenerate`` options that are
tuned in practice; every other STAR flag is reachable through ``**kwargs``.
"""

from __future__ import annotations

import math
from functools import partial
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genome.aligner.aligner import Aligner

if TYPE_CHECKING:
    from genome.assembly.genome import Genome
    from genome.external import ExternalTool


class STAR(Aligner):
    """STAR aligner index builder.

    A STAR index is splice-junction-aware: it is built against one gene
    annotation, so each annotation gets its own *genomeDir*. The bound GTF key
    selects that annotation (its path resolved via
    :meth:`~genome.annotation.registry.AnnotationRegistry.path`) and names the index directory
    ``star_<gtf_key>``. STAR's index is the *genomeDir* directory itself, so
    :attr:`index_path` returns :attr:`index_dir`.

    Parameters
    ----------
    genome : genome.assembly.genome.Genome
        The genome whose reference FASTA will be indexed.
    gtf : str
        Name of a GTF annotation registered on ``genome`` (see
        :meth:`~genome.annotation.registry.AnnotationRegistry.register_path`).
    tool : genome.external.ExternalTool, optional
        As :class:`~genome.aligner.aligner.Aligner`.
    """

    name = "star"
    binary = "STAR"
    # STAR's long options are camel-case words with no separator, so a Python keyword
    # reaches the command line unchanged.
    _flag_separator = "_"

    def __init__(self, genome: Genome, *, gtf: str, tool: ExternalTool | None = None) -> None:
        self._gtf_key = gtf
        super().__init__(genome, tool=tool)

    @property
    def index_dir(self) -> Path:
        """Per-annotation genome directory ``.../index/star_<gtf_key>/``."""
        base = super().index_dir
        return base.with_name(f"{base.name}_{self._gtf_key}")

    @property
    def _artifact(self) -> Path:
        """STAR loads the genome directory directly."""
        return self.index_dir

    @property
    def _build_arguments(self) -> str:
        """The annotation key, which is what picks one STAR index out of several."""
        return f"gtf={self._gtf_key!r}"

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
        An index whose completion record says it finished is reused unless
        ``overwrite=True``; a directory holding index files that no record
        vouches for raises rather than being silently rebuilt. The annotation GTF
        is resolved from the bound ``gtf`` key via
        :meth:`~genome.annotation.registry.AnnotationRegistry.path` and passed to STAR as
        ``--sjdbGTFfile`` for splice-junction-aware indexing.

        Only the most commonly tuned options are named below. Any other STAR
        ``genomeGenerate`` option may be passed as a keyword argument using its
        STAR name without the leading ``--`` (e.g. ``genomeSAindexNbases=11``);
        for the meaning of those options see the STAR manual / ``STAR --help``.
        Two of them are sized from the assembly rather than left to STAR's
        defaults — the suffix-array index size ``genomeSAindexNbases``, from the
        total sequence length, and the genome-storage bin size
        ``genomeChrBinNbits``, from the mean sequence length and the read length
        ``sjdb_overhang`` implies. Passing either one yourself wins.

        Parameters
        ----------
        sjdb_overhang : int, default 100
            ``--sjdbOverhang``: ideally ``read_length - 1``, which is also where
            the computed ``genomeChrBinNbits`` reads the read length back out of.
        threads : int, default 1
            ``--runThreadN``: number of threads to build with.
        overwrite : bool, default False
            Rebuild even if a finished index already exists, and rebuild over a
            directory that cannot be trusted rather than raising on it.
        **kwargs : Any
            Extra ``genomeGenerate`` options forwarded verbatim as STAR flags.

        Returns
        -------
        pathlib.Path
            The genome directory (also available as :attr:`index_path`).

        Raises
        ------
        genome.store.completion.RegistrationError
            If the index directory holds files without a record, a record that
            disagrees with them, or a record pinning a different assembly digest
            than the one registered now. Pass ``overwrite=True`` to rebuild.
        RuntimeError
            If STAR exits non-zero.
        """
        compose = partial(self._compose, sjdb_overhang=sjdb_overhang, threads=threads, **kwargs)
        return self._build(compose, overwrite=overwrite)

    def _compose(
        self, *, sjdb_overhang: int, threads: int, **kwargs: Any
    ) -> tuple[list[str], dict[str, Any]]:
        """Return the ``genomeGenerate`` command line, and the knobs that determined it.

        Reached only when a build is going to run, which is what keeps resolving the
        annotation out of the path that reuses a finished index — see
        :meth:`~genome.aligner.aligner.Aligner._build`.
        """
        fasta = self._genome.files.fasta
        # STAR runs with the index dir as CWD, so resolve the annotation path.
        gtf_file = self._genome.annotations.path(self._gtf_key).resolve()

        # STAR requires a reduced suffix-array index size for small genomes:
        #   min(14, log2(genomeLength) / 2 - 1). Default it unless overridden.
        if "genomeSAindexNbases" not in kwargs:
            genome_length = int(self._genome.chrom_sizes.sum())
            kwargs["genomeSAindexNbases"] = min(14, max(2, int(math.log2(genome_length) / 2 - 1)))

        # Each sequence is padded up to a whole number of 2^genomeChrBinNbits bases, so a
        # many-sequence reference left at STAR's default of 18 wastes up to 262,144 bases
        # per sequence — silently, because STAR neither warns on this one nor clamps it.
        # Scale it as STAR's manual recommends,
        #   min(18, log2(max(genomeLength / references, readLength))),
        # taking the read length back out of sjdb_overhang, which is read_length - 1 by
        # definition. That inner max is already the floor — a bin shorter than one read
        # helps nobody — so no separate guard. Passed even when it lands on 18, so the
        # completion record says what the build asked for rather than what it left out.
        # Truncated only to agree with the line above: there truncation reproduces STAR's
        # own C++ where its manual rounds, but for this knob STAR computes nothing at all,
        # so there is no source to match and the neighbour is the whole argument.
        if "genomeChrBinNbits" not in kwargs:
            chrom_sizes = self._genome.chrom_sizes
            mean_length = float(chrom_sizes.sum()) / len(chrom_sizes)
            bin_nbits = int(math.log2(max(mean_length, sjdb_overhang + 1)))
            kwargs["genomeChrBinNbits"] = min(18, bin_nbits)

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
        args += self._flags(kwargs)

        return args, parameters
