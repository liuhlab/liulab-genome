"""Download and cache large genomic files with `pooch <https://www.fatiando.org/pooch/>`_.

This is an I/O boundary module: it reaches out to the network. :class:`Downloader`
is a thin, reusable wrapper over :func:`pooch.retrieve` that fetches a file once
and caches it on disk, so repeated requests for the same URL are served locally.
:class:`UCSCGenomeDownloader` specializes it for reference-genome FASTA files
served from the UCSC golden path.

Notes
-----
- Downloads are cached under :func:`pooch.os_cache` (per-user cache directory)
  by default; pass ``cache_dir`` to override (e.g. a shared lab scratch path).
- ``known_hash`` is optional. When omitted, the download is **not** verified and
  pooch logs the computed hash — pin that value on subsequent calls to detect a
  corrupted or upstream-changed file. UCSC publishes ``md5sum.txt`` next to its
  ``bigZips`` downloads if you want to supply a hash.

Examples
--------
>>> from genome.io.download import UCSCGenomeDownloader
>>> dl = UCSCGenomeDownloader("hg38")            # doctest: +SKIP
>>> fasta = dl.fetch_fasta()                     # downloads + decompresses, cached
>>> fasta.name                                   # doctest: +SKIP
'hg38.fa'
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pooch

# A pooch post-processor: called with (fname, action, pooch_instance) and
# returns the path (or paths) to use as the result of the download.
Processor = Callable[..., object]


class Downloader:
    """Download and cache large files via pooch.

    Parameters
    ----------
    cache_dir : str or pathlib.Path, optional
        Directory under which downloads are stored. Defaults to the per-user
        cache location for the ``genome`` application
        (``pooch.os_cache("genome")``). The directory is created on first use.

    Attributes
    ----------
    cache_dir : pathlib.Path
        Resolved cache directory used for all downloads.

    Examples
    --------
    >>> dl = Downloader()                         # doctest: +SKIP
    >>> path = dl.fetch("https://example.org/big.bed.gz")   # doctest: +SKIP
    """

    def __init__(self, cache_dir: str | Path | None = None) -> None:
        self.cache_dir: Path = Path(pooch.os_cache("genome") if cache_dir is None else cache_dir)

    def fetch(
        self,
        url: str,
        *,
        known_hash: str | None = None,
        fname: str | None = None,
        processor: Processor | None = None,
        progressbar: bool = True,
    ) -> Path:
        """Download ``url`` into the cache and return the local path.

        If a previously downloaded copy is already present (and matches
        ``known_hash`` when one is given), it is reused without hitting the
        network.

        Parameters
        ----------
        url : str
            The file URL to download.
        known_hash : str, optional
            Expected hash as ``"<algorithm>:<hexdigest>"`` (e.g.
            ``"md5:8f3c..."``). If ``None``, verification is skipped and pooch
            logs the computed hash so you can pin it next time.
        fname : str, optional
            Local file name to save as. Defaults to the basename of ``url``.
        processor : callable, optional
            A pooch post-processor applied after download, such as
            :class:`pooch.Decompress` or :class:`pooch.Untar`. Its return value
            becomes the path returned by this method.
        progressbar : bool, default True
            Show a textual download progress bar (requires ``tqdm``).

        Returns
        -------
        pathlib.Path
            Absolute path to the cached (and, if ``processor`` was given,
            processed) file.

        Raises
        ------
        requests.exceptions.HTTPError
            If the download fails (e.g. the URL 404s).
        ValueError
            If ``known_hash`` is given and the downloaded file does not match.
        """
        result = pooch.retrieve(
            url=url,
            known_hash=known_hash,
            fname=fname,
            path=self.cache_dir,
            processor=processor,
            progressbar=progressbar,
        )
        return Path(result)


class UCSCGenomeDownloader(Downloader):
    """Download reference-genome FASTA files from the UCSC golden path.

    UCSC serves per-assembly downloads under
    ``https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/bigZips/``. This
    fetches the soft-masked, gzipped whole-genome FASTA
    (``<assembly>.fa.gz``) and, by default, decompresses it to
    ``<assembly>.fa`` inside the cache.

    Parameters
    ----------
    assembly : str
        UCSC assembly name, e.g. ``"hg38"``, ``"hg19"``, ``"mm39"``.
    cache_dir : str or pathlib.Path, optional
        See :class:`Downloader`.

    Attributes
    ----------
    assembly : str
        The assembly name passed at construction.

    Examples
    --------
    >>> dl = UCSCGenomeDownloader("hg38")
    >>> dl.fasta_url
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    >>> fasta = dl.fetch_fasta()                  # doctest: +SKIP
    """

    BASE_URL: str = "https://hgdownload.soe.ucsc.edu/goldenPath"

    def __init__(self, assembly: str, cache_dir: str | Path | None = None) -> None:
        super().__init__(cache_dir)
        self.assembly = assembly

    @property
    def fasta_url(self) -> str:
        """URL of the gzipped whole-genome FASTA for this assembly."""
        return f"{self.BASE_URL}/{self.assembly}/bigZips/{self.assembly}.fa.gz"

    def fetch_fasta(
        self,
        *,
        known_hash: str | None = None,
        decompress: bool = True,
        progressbar: bool = True,
    ) -> Path:
        """Download (and optionally decompress) the genome FASTA.

        Parameters
        ----------
        known_hash : str, optional
            Expected hash of the **downloaded ``.fa.gz``** (computed before
            decompression); see :meth:`Downloader.fetch`.
        decompress : bool, default True
            If ``True``, gunzip the download to ``<assembly>.fa`` and return
            that path. If ``False``, keep and return the ``.fa.gz``.
        progressbar : bool, default True
            Show a download progress bar (requires ``tqdm``).

        Returns
        -------
        pathlib.Path
            Path to the decompressed ``<assembly>.fa`` (or the ``.fa.gz`` when
            ``decompress=False``).
        """
        processor: Processor | None = (
            pooch.Decompress(method="gzip", name=f"{self.assembly}.fa") if decompress else None
        )
        return self.fetch(
            self.fasta_url,
            known_hash=known_hash,
            processor=processor,
            progressbar=progressbar,
        )
