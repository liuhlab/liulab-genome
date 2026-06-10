"""Mixin adding aligner-index construction to :class:`~genome.genome.Genome`."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from genome.genome import Genome


class AlignerMixin:
    """Build aligner genome indexes for a :class:`~genome.genome.Genome`.

    Each ``build_<aligner>_index`` method instantiates the corresponding
    :class:`~genome.aligner.aligner.Aligner`, which checks that the aligner is
    installed before doing anything. Index files land under
    ``<LIULAB_DATA>/genome/<assembly>/index/<aligner>/``.
    """

    def build_star_index(self, gtf: str, **kwargs: Any) -> Path:
        """Build a STAR genome index for this assembly against annotation ``gtf``.

        A thin entry point onto :meth:`genome.aligner.star.STAR.index`; the
        remaining keyword arguments are forwarded there (see its docstring for
        the exposed options and how to pass arbitrary STAR flags).

        Parameters
        ----------
        gtf : str
            Name of a GTF annotation registered on this genome (see
            :meth:`~genome.genome.Genome.register_gtf`). Its path is resolved via
            :meth:`~genome.genome.Genome.get_gtf_path` and passed to STAR; the
            index is written to a per-annotation directory ``index/star_<gtf>/``,
            so different annotations build independent indexes.

        Returns
        -------
        pathlib.Path
            The built STAR genome directory.
        """
        from genome.aligner.star import STAR

        return STAR(cast("Genome", self), gtf=gtf).index(**kwargs)
