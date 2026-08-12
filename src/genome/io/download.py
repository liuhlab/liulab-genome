"""Download and cache large genomic files with `pooch <https://www.fatiando.org/pooch/>`_.

I/O boundary module: it reaches out to the network. :func:`fetch_url` is the package's
single fetch step — every download goes through it and it is the only call site of
:func:`pooch.retrieve`. :class:`Downloader` binds that step to a cache directory, and
:class:`UCSCGenomeDownloader` specializes it for reference-genome FASTA files from the
UCSC golden path. See each for caching, storage layout, and hashing.

An assembly listed in the curated metadata table brings its own source URL and, where
the lab has pinned one, the sha256 of its **unpacked** FASTA. That digest is checked
after decompression rather than by pooch: pooch hashes the bytes it downloaded, which
are the ``.fa.gz``, and gzip bytes change under recompression while the FASTA inside
does not (ADR-0006).

pooch is used as a downloader and nothing more — its own cache is deliberately not
relied on. Downloads land in the assembly's working area, the unpacked FASTA is moved
out of it, and the whole area is discarded once
:class:`~genome.io.completion.CompletionRecord` says the registration finished. That
record owns the judgment pooch's cache would otherwise make: what is already here, and
whether it is usable.

Examples
--------
>>> from genome.io.download import UCSCGenomeDownloader
>>> dl = UCSCGenomeDownloader("hg38")            # doctest: +SKIP
>>> files = dl.fetch_genome()                    # download + decompress + verify + prepare
>>> files.chrom_sizes.name                       # doctest: +SKIP
'hg38.chrom.sizes'
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable
from pathlib import Path
from urllib.parse import urlparse

import pooch
import requests

from genome.io.completion import (
    build_record,
    clear_work_dir,
    disagreements,
    read_record,
    work_dir,
    write_record,
)
from genome.io.fasta import PREPARATION_TOOLS, GenomeFiles, prepare_fasta
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import METADATA_FIELDS, AssemblyMetadata, lookup_assembly

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

#: URL schemes a source may use, matching the downloaders pooch ships. Anything
#: else is a local path (or an error, if it carries a scheme we cannot fetch).
_URL_SCHEMES = frozenset({"http", "https", "ftp", "sftp"})


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
    """Return whether ``source`` is a fetchable URL rather than a local path."""
    return urlparse(source).scheme.lower() in _URL_SCHEMES


def fetch_url(
    url: str,
    dest_dir: Path,
    *,
    known_hash: str | None = None,
    fname: str | None = None,
    processor: _Processor | None = None,
    progressbar: bool = True,
) -> Path:
    """Download ``url`` into ``dest_dir`` and return the local path.

    The package's one fetch step: every download it performs goes through here, and
    this is the only call site of :func:`pooch.retrieve`. pooch picks its transport
    from the URL scheme, so ``http``, ``https``, ``ftp`` and ``sftp`` all work — the
    last additionally needs ``paramiko`` installed, which pooch will say for itself.
    A file already sitting at the destination is reused (verified against
    ``known_hash`` when one is given) and no network call is made.

    Reach this function through the module, never by importing the name: write
    ``from genome.io import download`` and call ``download.fetch_url(...)``, so that a
    single ``monkeypatch.setattr(download, "fetch_url", ...)`` takes every download in
    the package offline. Callers inside this module call it as a module global for the
    same reason.

    Parameters
    ----------
    url : str
        The file URL to download, including its scheme.
    dest_dir : pathlib.Path
        Directory the file is written into. Created when a download actually happens.
    known_hash : str, optional
        Expected hash as ``"<algorithm>:<hexdigest>"`` (e.g. ``"md5:8f3c..."``), or a
        bare hex digest for sha256. If ``None``, verification is skipped and pooch logs
        the computed hash so you can pin it next time.
    fname : str, optional
        Local file name to save as. Defaults to a hash-prefixed unique name pooch
        derives from ``url``.
    processor : callable, optional
        A pooch post-processor applied after the download, such as
        :class:`pooch.Decompress` or :class:`pooch.Untar`. Its return value becomes the
        path returned here.
    progressbar : bool, default True
        Show a textual download progress bar (requires ``tqdm``).

    Returns
    -------
    pathlib.Path
        Absolute path to the downloaded (and, if ``processor`` was given, processed)
        file.

    Raises
    ------
    requests.exceptions.HTTPError
        If an http(s) download fails (e.g. the URL 404s).
    ValueError
        If ``known_hash`` is given and the file does not match it, or if no downloader
        exists for the URL's scheme.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io import download
    >>> download.fetch_url(                                  # doctest: +SKIP
    ...     "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
    ...     Path("/scratch/liulab/genome/sacCer3"),
    ...     fname="sacCer3.fa.gz",
    ... )
    PosixPath('/scratch/liulab/genome/sacCer3/sacCer3.fa.gz')
    """
    result = pooch.retrieve(
        url=url,
        known_hash=known_hash,
        fname=fname,
        path=dest_dir,
        processor=processor,
        progressbar=progressbar,
    )
    return Path(result)


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
        """Download ``url`` into :attr:`cache_dir` and return the local path.

        :func:`fetch_url` bound to this downloader's cache directory — see it for the
        arguments, the reuse-without-network behaviour, the hashing and what raises.

        Returns
        -------
        pathlib.Path
            Absolute path to the cached (and, if ``processor`` was given, processed)
            file.
        """
        return fetch_url(
            url,
            self.cache_dir,
            known_hash=known_hash,
            fname=fname,
            processor=processor,
            progressbar=progressbar,
        )


class UCSCGenomeDownloader(Downloader):
    """Download reference-genome FASTA files from the UCSC golden path.

    UCSC serves per-assembly downloads under
    ``https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/bigZips/``. This
    fetches the soft-masked, gzipped whole-genome FASTA
    (``<assembly>.fa.gz``) and, by default, decompresses it to
    ``<assembly>.fa``.

    The assembly's metadata row is consulted first: when it pins a source URL that URL
    is fetched instead of the derived golden-path one, and when it pins a sha256 the
    unpacked FASTA is checked against it. An assembly the table does not list keeps
    working exactly as before — the table is a cross-reference, not an allow-list.

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
    metadata : genome.metadata.AssemblyMetadata, optional
        A complete metadata record to use *instead of* the curated table's row for
        ``assembly``. Omit it and the row is looked up here, so a downloader used on
        its own still gets the pinned source and checksum.

    Attributes
    ----------
    assembly : str
        The assembly name passed at construction.
    metadata : genome.metadata.AssemblyMetadata or None
        The record this downloader works from — the one passed in, else the curated
        table's row, else ``None`` for an assembly the table does not list.

    Examples
    --------
    >>> dl = UCSCGenomeDownloader("hg38")
    >>> dl.fasta_url
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    >>> fasta = dl.fetch_fasta()                  # doctest: +SKIP
    """

    BASE_URL: str = "https://hgdownload.soe.ucsc.edu/goldenPath"

    def __init__(
        self,
        assembly: str,
        cache_dir: str | Path | None = None,
        *,
        metadata: AssemblyMetadata | None = None,
    ) -> None:
        if cache_dir is None:
            cache_dir = assembly_data_dir(assembly)
        super().__init__(cache_dir)
        self.assembly = assembly
        self.metadata: AssemblyMetadata | None = (
            metadata if metadata is not None else lookup_assembly(assembly)
        )

    @property
    def assembly_url(self) -> str:
        """URL of the UCSC golden-path directory for this assembly."""
        return f"{self.BASE_URL}/{self.assembly}/"

    @property
    def _pinned_source_url(self) -> str | None:
        """The source URL this assembly's metadata pins, or ``None`` when it pins none."""
        return self.metadata.source_url if self.metadata else None

    @property
    def _expected_sha256(self) -> str | None:
        """The sha256 this assembly's metadata pins for the unpacked FASTA, or ``None``."""
        return self.metadata.sha256 if self.metadata else None

    @property
    def fasta_url(self) -> str:
        """URL of the gzipped whole-genome FASTA — the pinned source, else the golden path.

        A metadata row that pins a source URL answers this outright; otherwise the URL
        is derived from the assembly name and UCSC's golden-path layout, as it always
        was.
        """
        pinned = self._pinned_source_url
        if pinned:
            return pinned
        return f"{self.BASE_URL}/{self.assembly}/bigZips/{self.assembly}.fa.gz"

    @property
    def _work_dir(self) -> Path:
        """The disposable working area this assembly downloads into.

        Inside :attr:`cache_dir` on purpose: same filesystem, so placing the unpacked
        FASTA is a rename rather than a copy, and it goes when the assembly does. See
        :func:`~genome.io.completion.work_dir`.
        """
        return work_dir(self.cache_dir)

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
        """Return the prepared GenomeFiles when the record says so, else ``None``.

        The completion record is the only thing consulted: it must be there, and every
        file it claims must be present at the size it claims. That is one ``stat`` per
        file and no file contents, so reopening a prepared genome is instant. Anything
        else — no record, or a record disagreeing with disk — reads as unfinished here
        and falls through to a fresh registration, exactly as the absence of the old
        marker did. ``overwrite`` skips the question entirely.
        """
        if overwrite:
            return None
        record = read_record(self.cache_dir)
        if record is None or disagreements(self.cache_dir, record):
            return None
        return self._expected_genome_files()

    def _record_completion(
        self, files: GenomeFiles, *, source_url: str | None, sha256: str | None
    ) -> None:
        """Write this assembly's completion record, then discard the working area.

        Called last, once every derived file exists — writing the record is what makes
        the registration finished, and the archive is only disposable after that. An
        interrupted run therefore leaves its download in place and repairs without
        fetching a whole genome again.
        """
        record = build_record(
            self.cache_dir,
            kind="genome",
            name=self.assembly,
            files=[files.fasta, files.fai, files.twobit, files.chrom_sizes],
            source_url=source_url,
            sha256=sha256,
            tools=PREPARATION_TOOLS,
        )
        write_record(self.cache_dir, record)
        clear_work_dir(self.cache_dir)

    def _place_fasta(self, unpacked: Path) -> Path:
        """Move ``unpacked`` out of the working area to ``<assembly>.fa`` and return it.

        A rename within one filesystem, since the working area sits inside
        :attr:`cache_dir` — no second copy of a whole genome is ever made.
        """
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        fasta = self._expected_genome_files().fasta
        if unpacked.resolve() != fasta.resolve():
            unpacked.replace(fasta)
        return fasta

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

    def verify_fasta(self, fasta: Path | None = None) -> str:
        """Return the sha256 of the unpacked FASTA, raising when it is not the pinned one.

        The digest is taken over the **unpacked** ``<assembly>.fa``, never over the
        ``.fa.gz`` it arrived in, which is why pooch's own ``known_hash`` cannot do this
        job: pooch hashes what it downloaded. Gzip bytes change under recompression
        while the FASTA inside does not, so a content digest also matches a copy taken
        from a mirror or handed over by hand (ADR-0006). The file is streamed, so a
        whole-genome FASTA is never held in memory.

        Parameters
        ----------
        fasta : pathlib.Path, optional
            The file to check. Defaults to this assembly's ``<assembly>.fa`` in
            :attr:`cache_dir`.

        Returns
        -------
        str
            The computed hex digest. An assembly whose metadata pins no sha256 — or
            that the table does not list — has nothing to disagree with, so the value
            is simply reported back for recording or for pinning later.

        Raises
        ------
        genome.io.utils.ChecksumMismatchError
            If the metadata pins a sha256 and the file's digest is a different one;
            the message names both values.
        FileNotFoundError
            If the file does not exist.

        Examples
        --------
        >>> UCSCGenomeDownloader("sacCer3").verify_fasta()      # doctest: +SKIP
        '6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3'
        """
        target = self._expected_genome_files().fasta if fasta is None else fasta
        actual = sha256_file(target)
        expected = self._expected_sha256
        if expected is not None and actual != expected:
            raise ChecksumMismatchError(target, expected, actual)
        return actual

    def fetch_fasta(
        self,
        *,
        known_hash: str | None = None,
        decompress: bool = True,
        progressbar: bool = True,
    ) -> Path:
        """Download (and optionally decompress) the genome FASTA into the working area.

        Both files land in the assembly's working area rather than beside its prepared
        files: nothing there is claimed by a completion record, and the whole area is
        discarded once one is written. :meth:`fetch_genome` is what moves the unpacked
        FASTA out of it.

        When the URL had to be derived from the assembly name, the assembly is first
        confirmed to exist at UCSC via :meth:`validate_assembly`, so a typo fails fast
        with a clear message. When the metadata row pins a source URL that check is
        skipped: validation is a property of the source (ADR-0003), and a pinned URL
        *is* the source, so there is nothing left to guess.

        Parameters
        ----------
        known_hash : str, optional
            Expected hash of the **downloaded ``.fa.gz``** (checked by pooch, before
            decompression); see :meth:`Downloader.fetch`. Unrelated to the metadata
            row's ``sha256``, which covers the unpacked FASTA — see
            :meth:`verify_fasta`.
        decompress : bool, default True
            If ``True``, gunzip the download to ``<assembly>.fa`` and return
            that path. If ``False``, keep and return the ``.fa.gz``.
        progressbar : bool, default True
            Show a download progress bar (requires ``tqdm``).

        Returns
        -------
        pathlib.Path
            Path **inside the working area** to the decompressed ``<assembly>.fa`` (or
            to the ``<assembly>.fa.gz`` when ``decompress=False``). The archive is kept
            there for the duration of the run, so an interrupted registration repairs
            without downloading a whole genome again.

        Raises
        ------
        ValueError
            If the URL was derived and the assembly is unknown to UCSC.
        """
        if self._pinned_source_url is None:
            self.validate_assembly()
        processor: _Processor | None = (
            pooch.Decompress(method="gzip", name=f"{self.assembly}.fa") if decompress else None
        )
        # The archive is named after the assembly, not after the URL: a pinned source
        # need not be UCSC, so its file name says nothing about which assembly this is.
        return fetch_url(
            self.fasta_url,
            self._work_dir,
            known_hash=known_hash,
            fname=f"{self.assembly}.fa.gz",
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

        Chains :meth:`fetch_fasta`, :meth:`verify_fasta` and
        :func:`genome.io.fasta.prepare_fasta`: download ``<assembly>.fa.gz`` from the
        assembly's source into the working area, decompress it, check the unpacked FASTA
        against the sha256 its metadata pins, move it to :attr:`cache_dir`, then build
        the ``.fai`` index, ``.2bit`` encoding, and ``.chrom.sizes`` beside it
        (``<LIULAB_DATA>/genome/<assembly>/`` by default).

        A completion record is written last, once all four files exist, holding the URL
        fetched, the digest computed, every file with its size, the tool versions, the
        package version and the time. That record is what makes a later call cheap: it
        is read, its claims are checked against disk by size alone, and nothing is
        fetched. The archive is deleted with the rest of the working area at that point
        — and only then, so an interrupted run still repairs from it. Pass
        ``overwrite=True`` to register again from scratch.

        Parameters
        ----------
        known_hash : str, optional
            Expected hash of the **downloaded ``.fa.gz``** (before decompression);
            see :meth:`Downloader.fetch`. When ``None``, pooch verifies nothing — which
            is independent of the metadata row's ``sha256`` over the unpacked FASTA,
            always checked here.
        progressbar : bool, default True
            Show a download progress bar (requires ``tqdm``).
        overwrite : bool, default False
            Register again from scratch: the assembly's completion record is not
            consulted, and the preparation steps (faidx, 2bit, chrom.sizes) rerun even
            when their outputs look fresh. An archive still sitting in the working area
            is reused rather than downloaded again.

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
        genome.io.utils.ChecksumMismatchError
            If the metadata pins a sha256 and the unpacked FASTA is not it.
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
        downloaded = self.fetch_fasta(
            known_hash=known_hash,
            decompress=True,
            progressbar=progressbar,
        )
        # Checked before it is placed, so a FASTA that is not the pinned one never
        # reaches the assembly dir; recorded from here, so a whole genome is not
        # hashed a second time just to write the record.
        digest = self.verify_fasta(downloaded)
        fasta = self._place_fasta(downloaded)
        files = prepare_fasta(fasta, overwrite=overwrite)
        self._record_completion(files, source_url=self.fasta_url, sha256=digest)
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
        (copied into the working area) or a **URL** (fetched with
        :func:`fetch_url`, so http(s), ftp and sftp all work). Gzipped sources
        (``.gz``) are decompressed. The resulting ``<assembly>.fa`` is then
        indexed/2bit/chrom.sizes-prepared exactly as :meth:`fetch_genome` does, and a
        completion record is written last, recording the source it was given and the
        digest of what arrived. UCSC is never contacted. Neither is the metadata row's
        pinned source or checksum: a seeded FASTA is whatever the caller handed over —
        it is recorded, never compared — and the assembly name is only a label for the
        directory it lands in.

        Parameters
        ----------
        source : str or pathlib.Path
            Local FASTA path or http(s)/ftp/sftp URL. ``.gz`` is decompressed.
        progressbar : bool, default True
            Show a download progress bar while fetching a URL (ignored for a
            local copy).
        overwrite : bool, default False
            Re-read the source and rerun preparation even when this assembly's
            completion record already says it is registered.

        Returns
        -------
        genome.io.fasta.GenomeFiles
            Paths to the prepared FASTA and its three derived files.

        Raises
        ------
        FileNotFoundError
            If ``source`` is a local path that does not exist.
        ValueError
            If ``source`` carries a URL scheme no downloader handles.
        genome.external.ToolNotFoundError
            If a preparation tool is not on ``PATH``.
        RuntimeError
            If any native preparation tool fails.
        """
        cached = self._completed_genome(overwrite=overwrite)
        if cached is not None:
            return cached
        fasta = self._materialize_fasta(source, progressbar=progressbar, overwrite=overwrite)
        files = prepare_fasta(fasta, overwrite=overwrite)
        self._record_completion(files, source_url=str(source), sha256=sha256_file(fasta))
        return files

    def _materialize_fasta(
        self,
        source: str | Path,
        *,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> Path:
        """Place ``source`` at ``<assembly>.fa`` in :attr:`cache_dir` and return that path.

        Copies a local file or fetches a URL through :func:`fetch_url` into the working
        area, decompresses a ``.gz`` source there, then moves the FASTA into place. A
        fresh existing ``<assembly>.fa`` is reused unless ``overwrite``, which also
        discards any kept download so the source is read again rather than reused.
        """
        fasta = self._expected_genome_files().fasta
        if fasta.is_file() and not overwrite:
            return fasta

        work = self._work_dir
        work.mkdir(parents=True, exist_ok=True)
        src = str(source)
        gzipped = src.endswith(".gz")
        name = f"{self.assembly}.fa.gz" if gzipped else f"{self.assembly}.fa"
        downloaded = work / name

        if _looks_like_url(src):
            if overwrite:
                downloaded.unlink(missing_ok=True)
            fetch_url(src, work, fname=name, progressbar=progressbar)
        else:
            local_source = Path(src).expanduser()
            if not local_source.is_file():
                scheme = urlparse(src).scheme.lower()
                if scheme:
                    raise ValueError(
                        f"cannot fetch {src!r}: no downloader for the {scheme!r} scheme. "
                        f"Pass a local file path, or a URL using one of: "
                        f"{', '.join(sorted(_URL_SCHEMES))}."
                    )
                raise FileNotFoundError(
                    f"local FASTA source not found: {local_source}. Pass an existing "
                    f"file path or an http(s)/ftp/sftp URL."
                )
            shutil.copy2(local_source, downloaded)

        if gzipped:
            downloaded = _gunzip(downloaded, work / f"{self.assembly}.fa")
        return self._place_fasta(downloaded)


def assembly_table_row(
    assembly: str,
    *,
    cache_dir: str | Path | None = None,
    progressbar: bool = True,
) -> dict[str, object]:
    r"""Fetch ``assembly``'s FASTA and return the metadata table row describing it.

    What makes filling in the table's checksum column a copy-paste rather than a manual
    hashing chore: the FASTA is downloaded and unpacked, its sha256 is computed over the
    **unpacked** file, and the assembly's row comes back with ``source_url`` set to the
    URL that was actually fetched and ``sha256`` to that digest. Every other field is
    the curated table's own — or blank for an assembly the table does not list, since
    those identifiers are ones only a person can supply.

    Only the FASTA is fetched: no ``.fai``, ``.2bit`` or ``chrom.sizes`` is built, so
    this needs no native tools. Nothing is registered either — both the download and its
    unpacked form stay in the assembly's working area, so hashing a genome for the table
    never leaves an unregistered FASTA among that assembly's own files, and re-running it
    reuses what is already there.

    Parameters
    ----------
    assembly : str
        The assembly to fetch, e.g. ``"sacCer3"``.
    cache_dir : str or pathlib.Path, optional
        Override which assembly directory the download works in. Defaults to
        :func:`assembly_data_dir(assembly) <assembly_data_dir>`.
    progressbar : bool, default True
        Show a download progress bar (requires ``tqdm``).

    Returns
    -------
    dict
        Field name to value over :data:`~genome.metadata.METADATA_FIELDS`, with
        ``None`` for anything still unknown. Hand it to
        :func:`~genome.metadata.format_table_row` for the line to paste into
        ``data/assembly_metadata.tsv``.

    Raises
    ------
    ValueError
        If the URL had to be derived and the assembly is unknown to UCSC.
    genome.io.utils.ChecksumMismatchError
        If the row already pins a sha256 and the FASTA that arrived is not it.

    Examples
    --------
    >>> from genome.metadata import format_table_row
    >>> format_table_row(assembly_table_row("sacCer3"))       # doctest: +SKIP
    'sacCer3\tSaccharomyces cerevisiae\t...\t6ff72f07...'
    """
    downloader = UCSCGenomeDownloader(assembly, cache_dir)
    fasta = downloader.fetch_fasta(progressbar=progressbar)
    record = downloader.metadata
    # Every column the record knows; all blank when the table lists no row at all, in
    # which case the name is the only thing that does not need a person to supply it.
    row: dict[str, object] = {name: getattr(record, name, None) for name in METADATA_FIELDS}
    if record is None:
        row["assembly_name"] = assembly
    row["source_url"] = downloader.fasta_url
    row["sha256"] = downloader.verify_fasta(fasta)
    return row


if __name__ == "__main__":
    downloader = UCSCGenomeDownloader("sacCer3")
    files = downloader.fetch_genome()
    print(files)
