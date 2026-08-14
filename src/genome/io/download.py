"""Download and cache large genomic files with `pooch <https://www.fatiando.org/pooch/>`_.

I/O boundary module: it reaches out to the network. :func:`fetch_url` is the package's
single fetch step — every download goes through it and it is the only call site of
:func:`pooch.retrieve`. :class:`Downloader` binds that step to a cache directory, and
:class:`UCSCGenomeDownloader` specializes it for reference-genome FASTA files from the
UCSC golden path. See each for caching, storage layout, and hashing.

Only the fetching is here. The rest of registering an assembly — the **Assembly dir**
layout, the working area, placing the FASTA and writing the record — is
:class:`~genome.io.registration.AssemblyRegistration`, which the downloader extends and
which knows nothing about where the bytes came from. The layout names it now owns stay
importable from this module, which is where they used to live.

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

A directory that cannot be trusted stops the work rather than being quietly rebuilt:
:func:`~genome.io.completion.check_registration` turns files-without-a-record and a
record-that-disagrees into errors naming ``genome register <assembly> --force``, and
that command is what repairs them (ADR-0007).

Examples
--------
>>> from genome.io.download import UCSCGenomeDownloader
>>> dl = UCSCGenomeDownloader("hg38")            # doctest: +SKIP
>>> files = dl.fetch_genome()                    # download + decompress + verify + prepare
>>> files.chrom_sizes.name                       # doctest: +SKIP
'hg38.chrom.sizes'
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from urllib.parse import urlparse

import pooch
import requests

from genome.io.completion import RECORD_NAME, RegistrationError, read_record
from genome.io.fasta import GenomeFiles, prepare_fasta

# Re-exported, not used here: the **Assembly dir** layout moved to `registration`
# because the shared registration steps are written in terms of it, and these are the
# spellings the rest of the package and its docs already import from this module.
from genome.io.registration import ANNOTATIONS_SUBDIR as ANNOTATIONS_SUBDIR
from genome.io.registration import INDEXES_SUBDIR as INDEXES_SUBDIR
from genome.io.registration import AssemblyRegistration
from genome.io.registration import assembly_data_dir as assembly_data_dir
from genome.io.registration import liulab_data_dir as liulab_data_dir
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import METADATA_FIELDS, AssemblyMetadata, lookup_assembly

# A pooch post-processor: called with (fname, action, pooch_instance) and
# returns the path (or paths) to use as the result of the download.
_Processor = Callable[..., object]

#: URL schemes a source may use, matching the downloaders pooch ships. Anything
#: else is a local path (or an error, if it carries a scheme we cannot fetch).
_URL_SCHEMES = frozenset({"http", "https", "ftp", "sftp"})

#: What answered *which digest should this FASTA have?* — the assembly's curated
#: metadata row, or the completion record its own registration wrote. Reported by
#: :func:`verify_assembly` so that being held to a pin and being held only to what this
#: machine last produced are never read as the same result.
_EXPECTED_FROM_TABLE = "table"
_EXPECTED_FROM_RECORD = "record"


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


class UCSCGenomeDownloader(AssemblyRegistration, Downloader):
    """Download reference-genome FASTA files from the UCSC golden path.

    UCSC serves per-assembly downloads under
    ``https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/bigZips/``. This
    fetches the soft-masked, gzipped whole-genome FASTA
    (``<assembly>.fa.gz``) and, by default, decompresses it to
    ``<assembly>.fa``.

    An :class:`~genome.io.registration.AssemblyRegistration` whose FASTA arrives over
    the network: the base owns the assembly's directory and the steps that finish a
    registration in it, and everything added here is about *where the bytes come from*
    — the URL, the name check, the pinned digest. :class:`Downloader` is the second base
    only for :meth:`Downloader.fetch`; ``cache_dir`` is the registration's.

    The assembly's metadata row is consulted first: when it pins a source URL that URL
    is fetched instead of the derived golden-path one, and when it pins a sha256 the
    unpacked FASTA is checked against it. An assembly the table does not list keeps
    working exactly as before — the table is a cross-reference, not an allow-list.

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
        # Resolves cache_dir to the assembly's own directory, which is the one thing
        # Downloader.__init__ would otherwise do — and it would default it to pooch's
        # cache instead. The bases are not cooperative; the assembly-aware one wins.
        super().__init__(assembly, cache_dir)
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

    def _expected_digest(self) -> tuple[str | None, str | None]:
        """Return the digest this assembly's FASTA is held to, and what supplied it.

        A **fallback**, and deliberately not a question about what kind of assembly this
        is: the curated row's pin answers whenever there is one, and otherwise the digest
        the assembly's own completion record already holds does. The second is the only
        thing that can answer for an assembly nothing was downloaded for — a chimera pins
        nothing in the table by design (ADR-0008) — but nothing here asks whether it is
        one, so any row that pins no digest gets the same treatment.

        For verification alone. Fetching still consults :attr:`_expected_sha256` by
        itself, because a record's digest is what this machine last produced, and holding
        a fresh download to it would refuse exactly when an upstream file has legitimately
        changed — which is the moment re-registering is what a caller wants.

        Returns
        -------
        tuple of (str or None, str or None)
            The digest to expect and where it came from — ``"table"`` for the curated
            row, ``"record"`` for the completion record — or ``(None, None)`` when
            neither pins one and there is nothing to check against.
        """
        pinned = self._expected_sha256
        if pinned is not None:
            return pinned, _EXPECTED_FROM_TABLE
        record = read_record(self.cache_dir)
        if record is not None and record.sha256 is not None:
            return record.sha256, _EXPECTED_FROM_RECORD
        return None, None

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

    def _proven_fasta(self) -> tuple[Path, str] | None:
        """Return the FASTA already on disk with its digest, when it is provably right.

        What makes repairing cheap: a re-registration that can prove the unpacked FASTA
        is the pinned one keeps it and rebuilds only the derived files, rather than
        pulling a whole genome down again. ``None`` — fetch the source again — in all
        three of the cases where it cannot be proven: the FASTA is missing, its digest
        is a different one, or **this assembly pins no digest at all**, since with
        nothing to compare against there is no way to show what is there is right.
        """
        expected = self._expected_sha256
        fasta = self._expected_genome_files().fasta
        if expected is None or not fasta.is_file():
            return None
        actual = sha256_file(fasta)
        return (fasta, actual) if actual == expected else None

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
        — and only then, so an interrupted run still repairs from it.

        A directory that cannot be trusted — files with no record, or a record that
        disagrees with what is on disk — **raises** rather than being rebuilt or
        trusted, naming ``genome register <assembly> --force`` (ADR-0007). That is what
        ``overwrite=True`` is: it skips the question, keeps the unpacked FASTA when its
        digest can be shown to be the pinned one, and fetches the source again when it
        cannot (see :meth:`_proven_fasta`). An absent or empty directory is not a broken
        state — it is a fresh registration and proceeds normally.

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
        genome.io.completion.UnfinishedRegistrationError
            If the assembly directory holds files but no record.
        genome.io.completion.RegistrationMismatchError
            If its record disagrees with what is on disk.
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
        registered = self._completed_genome(overwrite=overwrite, repair=self._repair_command())
        if registered is not None:
            return registered
        kept = self._proven_fasta()
        if kept is not None:
            fasta, digest = kept
        else:
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
        genome.io.completion.RegistrationError
            If the assembly directory holds files but no record, or a record that
            disagrees with what is on disk. The message names this same call as the
            repair — ``genome register <assembly> --force --source <source>`` — rather
            than the plain one, which would fetch from somewhere this assembly never
            came from.
        FileNotFoundError
            If ``source`` is a local path that does not exist.
        ValueError
            If ``source`` carries a URL scheme no downloader handles.
        genome.external.ToolNotFoundError
            If a preparation tool is not on ``PATH``.
        RuntimeError
            If any native preparation tool fails.
        """
        registered = self._completed_genome(
            overwrite=overwrite, repair=self._repair_command(source)
        )
        if registered is not None:
            return registered
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


def register_assembly(
    assembly: str,
    *,
    source: str | Path | None = None,
    force: bool = False,
    cache_dir: str | Path | None = None,
    progressbar: bool = True,
    metadata: AssemblyMetadata | None = None,
) -> dict[str, object]:
    """Prepare ``assembly`` on disk and return the record of what that did.

    Naming an assembly is enough: where its FASTA comes from and which digest it must
    match are the metadata table's to know. The whole pipeline runs — fetch, unpack,
    verify, index, derive — and the completion record lands last. An assembly that is
    already registered is returned from its record without fetching anything.

    A directory that cannot be trusted raises instead (see :meth:`
    UCSCGenomeDownloader.fetch_genome`); ``force=True`` is what repairs one, and it
    keeps an unpacked FASTA it can prove is the pinned one rather than downloading a
    whole genome again.

    Parameters
    ----------
    assembly : str
        The assembly to register, e.g. ``"hg38"``.
    source : str or pathlib.Path, optional
        Seed the assembly from this FASTA — a local path or an http(s)/ftp/sftp URL —
        instead of fetching the source its metadata pins. See
        :meth:`UCSCGenomeDownloader.fetch_genome_from`.
    force : bool, default False
        Register again from scratch, repairing a directory that raises.
    cache_dir : str or pathlib.Path, optional
        Override which directory the assembly is registered in. Defaults to
        :func:`assembly_data_dir(assembly) <assembly_data_dir>`.
    progressbar : bool, default True
        Show a download progress bar (requires ``tqdm``).
    metadata : genome.metadata.AssemblyMetadata, optional
        A complete metadata record to use instead of the curated table's row.

    Returns
    -------
    dict
        The completion record's own fields — ``files``, ``source_url``, ``sha256``,
        ``tool_versions``, ``completed_at`` and the rest — plus ``assembly`` and the
        ``directory`` they live in. Ready to serialize as it is.

    Raises
    ------
    genome.io.completion.RegistrationError
        If the directory holds a build that cannot be trusted as finished, or (with
        ``force``) if the run somehow left no record behind.
    genome.io.utils.ChecksumMismatchError
        If the metadata pins a sha256 and the unpacked FASTA is not it.
    genome.external.ToolNotFoundError
        If ``samtools``, ``faToTwoBit`` or ``twoBitInfo`` are not on ``PATH``.

    Examples
    --------
    >>> register_assembly("sacCer3")                       # doctest: +SKIP
    {'kind': 'genome', 'name': 'sacCer3', 'files': {...}, ...}
    """
    downloader = UCSCGenomeDownloader(assembly, cache_dir, metadata=metadata)
    if source is None:
        downloader.fetch_genome(progressbar=progressbar, overwrite=force)
    else:
        downloader.fetch_genome_from(source, progressbar=progressbar, overwrite=force)
    record = read_record(downloader.cache_dir)
    if record is None:
        raise RegistrationError(
            f"{assembly} was prepared in {downloader.cache_dir} but no {RECORD_NAME} is "
            f"there, so nothing can vouch for it. Register it again with "
            f"`{downloader._repair_command(source)}`."
        )
    payload: dict[str, object] = dict(asdict(record))
    payload["assembly"] = assembly
    payload["directory"] = str(downloader.cache_dir)
    return payload


def verify_assembly(
    assembly: str,
    *,
    fasta: str | Path | None = None,
    cache_dir: str | Path | None = None,
    metadata: AssemblyMetadata | None = None,
) -> dict[str, object]:
    """Re-read a FASTA and check its sha256 against the digest expected of it.

    The one operation that reads bytes rather than sizes. Registering an assembly and
    reopening it both go by presence and size, which is what makes them instant; this
    is the deliberate re-verification for when integrity is actually in doubt, and it
    costs a full pass over the file.

    What is expected of it comes from the assembly's curated row, and **failing that,
    from the completion record its own registration wrote** — a fallback rather than a
    question about what kind of assembly this is (see
    :meth:`UCSCGenomeDownloader._expected_digest`). ``expected_from`` says which
    answered, because being held to a pinned digest and being held only to what this
    machine last produced are different results and a caller must be able to tell them
    apart.

    With no ``fasta`` it verifies the assembly's own registered FASTA, and the
    registration must be intact — a directory that cannot be trusted raises here as it
    does anywhere else. Two things are then checked beside the digest, both by reading
    records rather than bytes: that each component this assembly was built from is still
    the one it was built from, and that each annotation merged into its own is still the
    one that was merged. Neither costs an assembly without components anything. Point
    ``fasta`` at any file to check that instead: a copy taken from a mirror or handed
    over by hand is checkable before anything is built on it, and nothing needs to be
    registered first — nothing is then asked about components, since the assembly's own
    registration is not what is being verified.

    Parameters
    ----------
    assembly : str
        The assembly whose row supplies the digest to check against, e.g. ``"hg38"``.
    fasta : str or pathlib.Path, optional
        A FASTA to check instead of the assembly's registered one.
    cache_dir : str or pathlib.Path, optional
        Override which directory the assembly is registered in.
    metadata : genome.metadata.AssemblyMetadata, optional
        A complete metadata record to use instead of the curated table's row.

    Returns
    -------
    dict
        ``assembly``, the ``fasta`` that was read, its computed ``sha256``, the
        ``expected`` digest (``None`` when nothing pins one), ``expected_from`` —
        ``"table"``, ``"record"``, or ``None`` for the same case — and ``verified``,
        whether there was anything to check against at all. A digest that disagrees
        raises rather than reporting ``False``.

    Raises
    ------
    genome.io.utils.ChecksumMismatchError
        If a digest is expected and the file's is a different one.
    genome.io.completion.RegistrationError
        If the assembly's directory holds a build that cannot be trusted as finished —
        including one whose components were registered again underneath it.
    FileNotFoundError
        If there is no file to read — nothing registered for ``assembly``, or no file
        at an explicit ``fasta``.

    Examples
    --------
    >>> verify_assembly("sacCer3")                          # doctest: +SKIP
    {'assembly': 'sacCer3', 'fasta': '...', 'sha256': '6ff72f07...', ...}
    >>> verify_assembly("sacCer3", fasta="/tmp/copied.fa")  # doctest: +SKIP
    {'assembly': 'sacCer3', 'fasta': '/tmp/copied.fa', ...}
    """
    # Deferred: `genome.io.chimera` reaches this module through `genome.io.gtf`, so
    # importing it at the top would be a cycle. Nothing else here needs it.
    from genome.io.chimera import check_components_unchanged

    downloader = UCSCGenomeDownloader(assembly, cache_dir, metadata=metadata)
    if fasta is not None:
        target = Path(fasta).expanduser()
        if not target.is_file():
            raise FileNotFoundError(
                f"no FASTA at {target}: pass the path of a file to check against "
                f"{assembly}'s row, or omit it to check {assembly}'s own registered FASTA."
            )
    else:
        target = downloader._expected_genome_files().fasta
        registered = downloader._completed_genome(
            overwrite=False, repair=downloader._repair_command()
        )
        if registered is None or not target.is_file():
            raise FileNotFoundError(
                f"{assembly} is not registered in {downloader.cache_dir}, so there is "
                f"nothing to verify. Register it with `genome register {assembly}`, or "
                f"pass the FASTA to check with --fasta."
            )
        # Beside the digest, never instead of it: this one reads records rather than
        # bytes, and answers what a digest of this assembly's own bytes cannot.
        check_components_unchanged(downloader.cache_dir, assembly)
    expected, expected_from = downloader._expected_digest()
    actual = sha256_file(target)
    if expected is not None and actual != expected:
        raise ChecksumMismatchError(target, expected, actual)
    return {
        "assembly": assembly,
        "fasta": str(target),
        "sha256": actual,
        "expected": expected,
        "expected_from": expected_from,
        "verified": expected is not None,
    }


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

    Notes
    -----
    The digest is **reported, never enforced**. A row that already pins one is not
    consulted, because this is the command to reach for when an upstream file has
    legitimately changed and the pin has to be regenerated — refusing on a mismatch
    would refuse exactly when the command is needed. Checking a FASTA you already
    hold against the official row is
    :meth:`UCSCGenomeDownloader.verify_fasta`'s job, not this one's.

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
    # Reported, never enforced: this is the command a maintainer runs *because* the
    # pinned digest needs regenerating, so comparing against the stale one would
    # refuse exactly when it is needed. Checking a FASTA against the official row is
    # what verifying an assembly is for.
    row["sha256"] = sha256_file(fasta)
    return row


if __name__ == "__main__":
    downloader = UCSCGenomeDownloader("sacCer3")
    files = downloader.fetch_genome()
    print(files)
