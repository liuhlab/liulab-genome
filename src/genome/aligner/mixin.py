"""Mixin adding aligner-index construction to :class:`~genome.genome.Genome`."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from genome.aligner.aligner import Aligner
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

    def build_chromap_index(self, **kwargs: Any) -> Path:
        """Build a chromap genome index for this assembly.

        A thin entry point onto :meth:`genome.aligner.chromap.Chromap.index`; the
        keyword arguments are forwarded there (see its docstring for the exposed
        options and how to pass arbitrary chromap flags). Unlike
        :meth:`build_star_index`, chromap needs no gene annotation, so there is no
        ``gtf`` argument and one index serves the whole assembly. The index is
        written to ``index/chromap/chromap.index``.

        Returns
        -------
        pathlib.Path
            The built chromap index file.
        """
        from genome.aligner.chromap import Chromap

        return Chromap(cast("Genome", self)).index(**kwargs)

    def get_index(self, aligner: str, **kwargs: Any) -> Path:
        """Return the path of an already-built index, for use in aligner commands.

        Locates the index a prior ``build_<aligner>_index`` produced and returns
        the file or prefix to hand to the aligner's own command line (e.g. STAR's
        ``--genomeDir``). Nothing is built here; if no successful index exists yet
        this raises, directing you to build it first.

        The aligner-specific selectors that identify *which* index are passed as
        keyword arguments, mirroring the matching ``build_*`` method — for STAR,
        the annotation ``gtf`` key that named the index.

        Parameters
        ----------
        aligner : str
            Aligner identifier, case-insensitive (e.g. ``"star"``).
        **kwargs : Any
            Aligner-specific selectors forwarded to the aligner constructor to
            pin down the index. STAR requires ``gtf``.

        Returns
        -------
        pathlib.Path
            The built index file or prefix, ready to drop into the aligner command.

        Raises
        ------
        ValueError
            If ``aligner`` is not a known aligner.
        RuntimeError
            If no successful index exists yet — build it first with the
            corresponding ``build_<aligner>_index`` method.
        """
        aligner_cls = _resolve_aligner(aligner)
        return aligner_cls(cast("Genome", self), **kwargs).index_path

    def get_star_index(self, gtf: str) -> Path:
        """Return the path of the STAR *genomeDir* built for annotation ``gtf``.

        Convenience wrapper over :meth:`get_index` for STAR — the returned
        directory is what STAR's ``--genomeDir`` expects at mapping time.

        Parameters
        ----------
        gtf : str
            Name of the GTF annotation the index was built against (the same key
            passed to :meth:`build_star_index`).

        Returns
        -------
        pathlib.Path
            The STAR genome directory.

        Raises
        ------
        RuntimeError
            If no successful STAR index exists yet for ``gtf`` — build it first
            with :meth:`build_star_index`.
        """
        return self.get_index("star", gtf=gtf)

    def get_chromap_index(self) -> Path:
        """Return the path of the chromap index built for this assembly.

        Convenience wrapper over :meth:`get_index` for chromap — the returned file
        is what chromap's ``-x/--index`` expects at mapping time. A chromap index
        carries no annotation, so no selector is needed.

        Returns
        -------
        pathlib.Path
            The chromap index file.

        Raises
        ------
        RuntimeError
            If no successful chromap index exists yet for this assembly — build it
            first with :meth:`build_chromap_index`.
        """
        return self.get_index("chromap")


def _resolve_aligner(name: str) -> type[Aligner]:
    """Look up an :class:`~genome.aligner.aligner.Aligner` subclass by its name.

    Raises :class:`ValueError` with the known names if ``name`` is unregistered.
    """
    from genome.aligner.chromap import Chromap
    from genome.aligner.star import STAR

    registry: dict[str, type[Aligner]] = {STAR.name: STAR, Chromap.name: Chromap}
    try:
        return registry[name.lower()]
    except KeyError:
        raise ValueError(f"Unknown aligner {name!r}; known aligners: {sorted(registry)}.") from None
