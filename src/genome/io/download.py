"""Download and cache large genomic files with `pooch <https://www.fatiando.org/pooch/>`_.

I/O boundary module: it reaches out to the network. :class:`Downloader` is a thin,
reusable wrapper over :func:`pooch.retrieve` that fetches a URL once and caches it;
:class:`UCSCGenomeDownloader` specializes it for reference-genome FASTA files from
the UCSC golden path. See each class for caching, storage layout, and hashing.

Examples
--------
>>> from genome.io.download import UCSCGenomeDownloader
>>> dl = UCSCGenomeDownloader("hg38")            # doctest: +SKIP
>>> files = dl.fetch_genome()                    # download + decompress + prepare
>>> files.chrom_sizes.name                       # doctest: +SKIP
'hg38.chrom.sizes'
"""

from __future__ import annotations

import os
import shutil
import subprocess
from collections.abc import Callable
from pathlib import Path
from urllib.parse import urlparse

import pooch
import requests

from genome.external import _resolve
from genome.io.fasta import GenomeFiles, prepare_fasta
from genome.io.utils import _gunzip

# A pooch post-processor: called with (fname, action, pooch_instance) and
# returns the path (or paths) to use as the result of the download.
_Processor = Callable[..., object]

#: Environment variable naming the lab data root directory.
LIULAB_DATA_ENV = "LIULAB_DATA"

#: Well-known lab data roots, tried in order when ``LIULAB_DATA`` is unset.
DEFAULT_LIULAB_DATA_PATHS = [
    "/share/lhqlab/liulab_data",
    "/large_storage/zhoulab/hanliu/liulab_data",
]

#: Sentinel file written in an assembly's cache dir once its FASTA pipeline
#: (download/seed → faidx → 2bit → chrom.sizes) has completed successfully.
#: Its presence lets repeat constructions skip fetching entirely.
PREPARED_MARKER = ".genome_prepared"


def liulab_data_dir() -> Path:
    """Return the root directory for lab reference data.

    The location is read from the ``LIULAB_DATA`` environment variable. When that
    is unset (or empty), each entry in :data:`DEFAULT_LIULAB_DATA_PATHS` is checked
    in order and the first that exists is used as the root. If none exist, it falls
    back to ``~/liulab_data``. The path is expanded (``~`` resolved) but **not**
    created here — callers create the specific subdirectory they need on first write.

    Returns
    -------
    pathlib.Path
        The resolved lab data root.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> liulab_data_dir()
    PosixPath('/scratch/liulab')
    >>> del os.environ["LIULAB_DATA"]
    """
    env = os.environ.get(LIULAB_DATA_ENV)
    if env:
        return Path(env).expanduser()
    for candidate in DEFAULT_LIULAB_DATA_PATHS:
        path = Path(candidate).expanduser()
        if path.exists():
            return path
    return (Path.home() / "liulab_data").expanduser()


def assembly_data_dir(assembly: str) -> Path:
    """Return the directory holding all reference files for ``assembly``.

    Every file tied to a reference assembly (FASTA, indexes, annotations, …)
    lives under ``<liulab_data>/genome/<assembly>/`` so they stay co-located.

    Parameters
    ----------
    assembly : str
        Assembly name, e.g. ``"hg38"``.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/genome/<assembly>``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> assembly_data_dir("hg38")
    PosixPath('/scratch/liulab/genome/hg38')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / "genome" / assembly


def _looks_like_url(source: str) -> bool:
    """Return whether ``source`` is an http(s)/ftp URL rather than a local path."""
    return urlparse(source).scheme.lower() in {"http", "https", "ftp", "ftps"}


def _curl_download(url: str, dest: Path, *, progressbar: bool = True) -> Path:
    """Download ``url`` to ``dest`` with ``curl``, writing atomically via a ``.part`` file."""
    curl = _resolve("curl")
    part = dest.with_name(dest.name + ".part")
    verbosity = "-#" if progressbar else "-sS"
    args = [curl, "-fL", verbosity, "--retry", "3", "-o", str(part), url]
    try:
        subprocess.run(args, check=True)
    except subprocess.CalledProcessError as err:
        part.unlink(missing_ok=True)
        raise RuntimeError(
            f"curl failed (exit {err.returncode}) downloading {url!r}. Check the URL "
            f"and your network/proxy settings (HTTP(S)_PROXY)."
        ) from err
    part.replace(dest)
    return dest


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
        processor: _Processor | None = None,
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
    ``<assembly>.fa``.

    Unless ``cache_dir`` is given, files are stored under the per-assembly
    reference directory ``<LIULAB_DATA>/genome/<assembly>/`` (see
    :func:`assembly_data_dir`), keeping all reference files for an assembly
    together.

    Parameters
    ----------
    assembly : str
        UCSC assembly name, e.g. ``"hg38"``, ``"hg19"``, ``"mm39"``.
    cache_dir : str or pathlib.Path, optional
        Override the storage directory. Defaults to
        :func:`assembly_data_dir(assembly) <assembly_data_dir>`.

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
        if cache_dir is None:
            cache_dir = assembly_data_dir(assembly)
        super().__init__(cache_dir)
        self.assembly = assembly

    @property
    def assembly_url(self) -> str:
        """URL of the UCSC golden-path directory for this assembly."""
        return f"{self.BASE_URL}/{self.assembly}/"

    @property
    def fasta_url(self) -> str:
        """URL of the gzipped whole-genome FASTA for this assembly."""
        return f"{self.BASE_URL}/{self.assembly}/bigZips/{self.assembly}.fa.gz"

    @property
    def prepared_marker(self) -> Path:
        """Path to the FASTA-pipeline-done sentinel for this assembly.

        See :data:`PREPARED_MARKER`. Lives directly in :attr:`cache_dir`.
        """
        return self.cache_dir / PREPARED_MARKER

    def _expected_genome_files(self) -> GenomeFiles:
        """Paths the FASTA pipeline produces for this assembly (whether or not they exist).

        Both :meth:`fetch_genome` and :meth:`fetch_genome_from` materialize the
        FASTA as ``<assembly>.fa`` and derive identically named companions, so a
        single layout describes either entry point.
        """
        fasta = self.cache_dir / f"{self.assembly}.fa"
        return GenomeFiles(
            fasta=fasta,
            fai=fasta.with_name(fasta.name + ".fai"),
            twobit=self.cache_dir / f"{self.assembly}.2bit",
            chrom_sizes=self.cache_dir / f"{self.assembly}.chrom.sizes",
        )

    def _completed_genome(self, *, overwrite: bool) -> GenomeFiles | None:
        """Return the cached GenomeFiles when the pipeline is already done, else ``None``.

        Skips short-circuiting when ``overwrite`` is set. Treats the
        :attr:`prepared_marker` as authoritative, but still confirms every
        derived file is present so a partially deleted cache falls through to a
        fresh fetch instead of returning dangling paths.
        """
        if overwrite or not self.prepared_marker.is_file():
            return None
        files = self._expected_genome_files()
        if all(p.is_file() for p in (files.fasta, files.fai, files.twobit, files.chrom_sizes)):
            return files
        return None

    def _mark_prepared(self) -> None:
        """Write the :attr:`prepared_marker` flag recording a completed FASTA pipeline."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.prepared_marker.write_text(
            f"FASTA pipeline complete for {self.assembly!r}.\n"
            "Presence of this file makes Genome construction skip re-fetching.\n"
        )

    def validate_assembly(self, *, timeout: float = 30.0) -> None:
        """Check that ``assembly`` is a real golden-path directory at UCSC.

        Sends an HTTP ``HEAD`` to :attr:`assembly_url` so a typo in the assembly
        name fails fast with a clear message, rather than surfacing later as an
        opaque 404 on the FASTA file itself.

        Parameters
        ----------
        timeout : float, default 30.0
            Seconds to wait for the server before giving up.

        Raises
        ------
        ValueError
            If no directory exists for ``assembly`` (HTTP 404) — the assembly
            name is almost certainly wrong.
        requests.exceptions.RequestException
            If the request fails for any other reason (network error, timeout,
            or an unexpected non-success status).

        Examples
        --------
        >>> UCSCGenomeDownloader("hg38").validate_assembly()    # doctest: +SKIP
        >>> UCSCGenomeDownloader("nope99").validate_assembly()  # doctest: +SKIP
        Traceback (most recent call last):
        ValueError: Unknown UCSC assembly 'nope99': no directory at ...
        """
        response = requests.head(self.assembly_url, timeout=timeout, allow_redirects=True)
        if response.status_code == 404:
            raise ValueError(
                f"Unknown UCSC assembly {self.assembly!r}: no directory at "
                f"{self.assembly_url}. Check the name against {self.BASE_URL}/ "
                f"(e.g. 'hg38', 'mm39', 'sacCer3')."
            )
        response.raise_for_status()

    def fetch_fasta(
        self,
        *,
        known_hash: str | None = None,
        decompress: bool = True,
        progressbar: bool = True,
    ) -> Path:
        """Download (and optionally decompress) the genome FASTA.

        The assembly is always confirmed to exist at UCSC via
        :meth:`validate_assembly` before any download is attempted, so a bad
        name fails fast with a clear message.

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

        Raises
        ------
        ValueError
            If the assembly is unknown to UCSC.
        """
        self.validate_assembly()
        processor: _Processor | None = (
            pooch.Decompress(method="gzip", name=f"{self.assembly}.fa") if decompress else None
        )
        return self.fetch(
            self.fasta_url,
            known_hash=known_hash,
            processor=processor,
            progressbar=progressbar,
        )

    def fetch_genome(
        self,
        *,
        known_hash: str | None = None,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> GenomeFiles:
        r"""Download and fully prepare the reference genome in one call.

        Chains :meth:`fetch_fasta` and :func:`genome.io.fasta.prepare_fasta`:
        download ``<assembly>.fa.gz`` from the UCSC golden path, decompress it, then
        build the ``.fai`` index, ``.2bit`` encoding, and ``.chrom.sizes``. All
        outputs land in :attr:`cache_dir` (``<LIULAB_DATA>/genome/<assembly>/`` by
        default), co-located with the kept ``.fa.gz`` download. Every step is cached;
        pass ``overwrite=True`` to force the preparation steps to rerun.

        Parameters
        ----------
        known_hash : str, optional
            Expected hash of the **downloaded ``.fa.gz``** (before decompression);
            see :meth:`Downloader.fetch`. When ``None``, verification is skipped.
        progressbar : bool, default True
            Show a download progress bar (requires ``tqdm``).
        overwrite : bool, default False
            Force the preparation steps (faidx, 2bit, chrom.sizes) to rerun even
            when their outputs look fresh. The pooch download/decompression cache is
            unaffected.

        Returns
        -------
        genome.io.fasta.GenomeFiles
            Paths to the decompressed FASTA and its three derived files.

        Raises
        ------
        requests.exceptions.HTTPError
            If the download fails (e.g. a wrong assembly name 404s).
        ValueError
            If the assembly is unknown to UCSC, or if ``known_hash`` is given
            and the download does not match.
        genome.external.ToolNotFoundError
            If ``samtools``, ``faToTwoBit``, or ``twoBitInfo`` are not on ``PATH``.
        RuntimeError
            If any native preparation tool exits non-zero.

        Examples
        --------
        >>> dl = UCSCGenomeDownloader("hg38")         # doctest: +SKIP
        >>> files = dl.fetch_genome()                 # download + decompress + prepare
        >>> files.fai.name, files.twobit.name, files.chrom_sizes.name   # doctest: +SKIP
        ('hg38.fa.fai', 'hg38.2bit', 'hg38.chrom.sizes')
        """
        cached = self._completed_genome(overwrite=overwrite)
        if cached is not None:
            return cached
        fasta = self.fetch_fasta(
            known_hash=known_hash,
            decompress=True,
            progressbar=progressbar,
        )
        files = prepare_fasta(fasta, overwrite=overwrite)
        self._mark_prepared()
        return files

    def fetch_genome_from(
        self,
        source: str | Path,
        *,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> GenomeFiles:
        """Prepare the genome from a user-provided FASTA instead of downloading from UCSC.

        Use this to seed an assembly from a file you already have or a non-UCSC
        URL — handy when the UCSC golden path is unreachable (firewall/proxy) or
        for a custom reference. ``source`` is either a **local filesystem path**
        (copied into :attr:`cache_dir`) or an **http(s)/ftp URL** (downloaded with
        ``curl``). Gzipped sources (``.gz``) are decompressed. The resulting
        ``<assembly>.fa`` is then indexed/2bit/chrom.sizes-prepared exactly as
        :meth:`fetch_genome` does. UCSC is never contacted.

        Parameters
        ----------
        source : str or pathlib.Path
            Local FASTA path or http(s)/ftp URL. ``.gz`` is decompressed.
        progressbar : bool, default True
            Show curl's progress meter while downloading a URL (ignored for a
            local copy).
        overwrite : bool, default False
            Re-fetch the source and rerun preparation even when a fresh
            ``<assembly>.fa`` is already cached.

        Returns
        -------
        genome.io.fasta.GenomeFiles
            Paths to the prepared FASTA and its three derived files.

        Raises
        ------
        FileNotFoundError
            If ``source`` is a local path that does not exist.
        genome.external.ToolNotFoundError
            If ``curl`` (for a URL) or a preparation tool is not on ``PATH``.
        RuntimeError
            If the curl download or any native preparation tool fails.
        """
        cached = self._completed_genome(overwrite=overwrite)
        if cached is not None:
            return cached
        fasta = self._materialize_fasta(source, progressbar=progressbar, overwrite=overwrite)
        files = prepare_fasta(fasta, overwrite=overwrite)
        self._mark_prepared()
        return files

    def _materialize_fasta(
        self,
        source: str | Path,
        *,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> Path:
        """Place ``source`` into the cache as ``<assembly>.fa`` and return that path.

        Copies a local file or downloads a URL with curl, decompressing a ``.gz``
        source. A fresh existing ``<assembly>.fa`` is reused unless ``overwrite``.
        """
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        fasta = self.cache_dir / f"{self.assembly}.fa"
        if fasta.is_file() and not overwrite:
            return fasta

        src = str(source)
        gzipped = src.endswith(".gz")
        download = self.cache_dir / (f"{self.assembly}.fa.gz" if gzipped else f"{self.assembly}.fa")

        if _looks_like_url(src):
            _curl_download(src, download, progressbar=progressbar)
        else:
            local_source = Path(src).expanduser()
            if not local_source.is_file():
                raise FileNotFoundError(
                    f"local FASTA source not found: {local_source}. Pass an existing "
                    f"file path or an http(s)/ftp URL."
                )
            if local_source.resolve() != download.resolve():
                shutil.copy2(local_source, download)

        if gzipped:
            _gunzip(download, fasta)
        return fasta


if __name__ == "__main__":
    downloader = UCSCGenomeDownloader("sacCer3")
    files = downloader.fetch_genome()
    print(files)
