"""Download and cache large genomic files with `pooch <https://www.fatiando.org/pooch/>`_.

I/O boundary module: it reaches out to the network. The package's single fetch step is
:func:`genome.io.fetch.fetch_url`, reached through the module so that one rebinding takes
every download offline; :class:`UCSCGenomeDownloader` drives it for reference-genome
FASTA files from the UCSC golden path. See it for caching, storage layout, and hashing.

Only the fetching is here. The rest of registering an assembly — the **Assembly dir**
layout, the working area, placing the FASTA and writing the record — is
:class:`~genome.io.registration.AssemblyRegistration`, which the downloader extends and
which knows nothing about where the bytes came from. The layout names it now owns stay
importable from this module, which is where they used to live.

**What a name means is settled before any of it**, by
:func:`~genome.io.source.resolve_source`: one name becomes a **Source** of one of three
kinds, and :meth:`UCSCGenomeDownloader.fetch_genome` dispatches on which came back rather
than asking again. A URL is fetched here; a component set is handed to
:func:`~genome.io.chimera.build_chimera`. That is what makes ``genome register <name>``
one command for every kind of **Source** — and what keeps the command line a thin client,
since neither the resolution nor its refusals are written there.

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
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import pooch
import requests

from genome.io import fetch
from genome.io.completion import RECORD_NAME, CompletionRecord, RegistrationError, read_record
from genome.io.components import ChimeraDetails, components_status
from genome.io.fasta import GenomeFiles, prepare_fasta
from genome.io.fetch import _Processor

# Re-exported, not used here: the **Assembly dir** layout moved to `registration`
# because the shared registration steps are written in terms of it, and these are the
# spellings the rest of the package and its docs already import from this module.
from genome.io.registration import ANNOTATIONS_SUBDIR as ANNOTATIONS_SUBDIR
from genome.io.registration import INDEXES_SUBDIR as INDEXES_SUBDIR
from genome.io.registration import AssemblyDir, AssemblyRegistration, assembly_repair_command
from genome.io.registration import assembly_data_dir as assembly_data_dir
from genome.io.registration import liulab_data_dir as liulab_data_dir
from genome.io.source import (
    ComponentSource,
    FetchedSource,
    SeededSource,
    fetched_source,
    resolve_source,
)
from genome.io.utils import ChecksumMismatchError, _gunzip, sha256_file
from genome.metadata import AssemblyMetadata, assembly_metadata

#: URL schemes a source may use, matching the downloaders pooch ships. Anything
#: else is a local path (or an error, if it carries a scheme we cannot fetch).
_URL_SCHEMES = frozenset({"http", "https", "ftp", "sftp"})

#: What answered *which digest should this FASTA have?* — the assembly's curated
#: metadata row, or the completion record its own registration wrote. Reported by
#: :func:`verify_assembly` so that being held to a pin and being held only to what this
#: machine last produced are never read as the same result. Public because the CLI keys
#: its two sentences on them: a surface that spelled the strings again would print the
#: raw status the day one of these was renamed, rather than failing.
EXPECTED_FROM_TABLE = "table"
EXPECTED_FROM_RECORD = "record"


@dataclass(frozen=True)
class RegisteredAssembly:
    """What preparing an assembly on disk produced: its record, and where that landed.

    :func:`register_assembly`'s answer — what ``genome register`` prints, and what its
    ``--json`` serializes. The **Completion marker** the run wrote *is* the answer, so it
    is carried whole rather than copied out field by field, and the two questions a
    surface then asks — which files are claimed, and is this a **Chimera** — are answered
    from that one record instead of by reading the directory again.

    Attributes
    ----------
    assembly : str
        The **Assembly** that was registered, under the name the caller asked for.
    directory : pathlib.Path
        Its **Assembly dir** — where those files and that record are.
    record : genome.io.completion.CompletionRecord
        The record the registration wrote, read back.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io.completion import CompletionRecord
    >>> registered = RegisteredAssembly(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38"),
    ...     record=CompletionRecord(
    ...         kind="genome",
    ...         name="hg38",
    ...         files={"hg38.fa.fai": 21, "hg38.fa": 12},
    ...         source_url="https://example.org/hg38.fa.gz",
    ...         sha256="1a2b3c",
    ...         tool_versions={},
    ...         package_version="2026.8.0",
    ...         completed_at="2026-08-12T09:00:00+00:00",
    ...         details={},
    ...     ),
    ... )
    >>> registered.file_names
    ['hg38.fa', 'hg38.fa.fai']
    >>> registered.chimera is None
    True
    >>> registered.as_json()["directory"]
    '/data/genome/hg38'
    """

    assembly: str
    directory: Path
    record: CompletionRecord

    @property
    def source_url(self) -> str | None:
        """Where the bytes were fetched from, or ``None`` when nothing was — a chimera's."""
        return self.record.source_url

    @property
    def sha256(self) -> str | None:
        """Digest of the unpacked FASTA, or ``None`` when none was computed."""
        return self.record.sha256

    @property
    def file_names(self) -> list[str]:
        """Every file the record claims, sorted — a fresh list each call."""
        return sorted(self.record.files)

    @property
    def chimera(self) -> ChimeraDetails | None:
        """What the build recorded about its components, or ``None`` for anything else.

        The record is what says an assembly is a **Chimera**, here as everywhere else —
        and the record is already in hand, so a surface reporting the registration that
        just happened never reads the same file a second time to find out.
        """
        return ChimeraDetails.from_record(self.record)

    def as_json(self) -> dict[str, Any]:
        """Return this registration as the payload ``--json`` serializes.

        The record's own fields under the record's own names, then the ``assembly`` asked
        for and the ``directory`` it landed in — the two facts a record does not hold
        about itself. The names are the ones written on disk and are never respelled here.

        Returns
        -------
        dict
            The record's fields, followed by ``assembly`` and ``directory``.
        """
        return {**asdict(self.record), "assembly": self.assembly, "directory": str(self.directory)}


@dataclass(frozen=True)
class VerifiedAssembly:
    """What re-reading a FASTA proved: its digest, what that was held to, and the components.

    :func:`verify_assembly`'s answer, and three results a caller must be able to tell
    apart, so each is a field of its own: the digest computed, *what supplied* the digest
    it was held to — being held to the lab's pin and being held to what this machine last
    produced are different answers — and, for a **Chimera**, what comparing its components
    settled. A digest that disagreed raises rather than arriving here, so this is what
    nothing refused.

    Attributes
    ----------
    assembly : str
        The **Assembly** whose row supplied the digest to check against.
    fasta : pathlib.Path
        The file that was read.
    sha256 : str
        The digest computed over it.
    expected : str or None
        The digest it was held to, or ``None`` when nothing pinned one.
    expected_from : str or None
        What answered with ``expected`` — :data:`EXPECTED_FROM_TABLE`,
        :data:`EXPECTED_FROM_RECORD`, or ``None`` when nothing did.
    components : str or None
        :data:`~genome.io.components.COMPONENTS_UNCHANGED` or
        :data:`~genome.io.components.COMPONENTS_UNKNOWN` for a chimera, and ``None`` for
        anything else — including every ``fasta`` checked on its own.

    Examples
    --------
    >>> from pathlib import Path
    >>> checked = VerifiedAssembly(
    ...     assembly="sacCer3",
    ...     fasta=Path("/data/genome/sacCer3/sacCer3.fa"),
    ...     sha256="6ff72f07",
    ...     expected="6ff72f07",
    ...     expected_from=EXPECTED_FROM_TABLE,
    ...     components=None,
    ... )
    >>> checked.verified
    True
    >>> checked.as_json()["expected_from"]
    'table'
    """

    assembly: str
    fasta: Path
    sha256: str
    expected: str | None
    expected_from: str | None
    components: str | None

    @property
    def verified(self) -> bool:
        """Whether there was a digest to check against at all, rather than merely one computed."""
        return self.expected is not None

    def as_json(self) -> dict[str, Any]:
        """Return this verification as the payload ``--json`` serializes.

        Returns
        -------
        dict
            Every attribute above, with ``fasta`` rendered as text and :attr:`verified`
            written out beside the fields it is read from.
        """
        return {
            "assembly": self.assembly,
            "fasta": str(self.fasta),
            "sha256": self.sha256,
            "expected": self.expected,
            "expected_from": self.expected_from,
            "verified": self.verified,
            "components": self.components,
        }


def _looks_like_url(source: str) -> bool:
    """Return whether ``source`` is a fetchable URL rather than a local path."""
    return urlparse(source).scheme.lower() in _URL_SCHEMES


def _expected_digest(
    assembly_dir: AssemblyDir, metadata: AssemblyMetadata
) -> tuple[str | None, str | None]:
    """Return the digest this assembly's FASTA is held to, and what supplied it.

    A **fallback**, and deliberately not a question about what kind of assembly this is:
    the curated row's pin answers whenever there is one, and otherwise the digest the
    assembly's own completion record already holds does. The second is the only thing that
    can answer for an assembly nothing was downloaded for — a chimera pins nothing in the
    table by design (ADR-0008) — but nothing here asks whether it is one, so any row that
    pins no digest gets the same treatment.

    For verification alone. Fetching consults the row by itself, because a record's digest
    is what this machine last produced, and holding a fresh download to it would refuse
    exactly when an upstream file has legitimately changed — which is the moment
    re-registering is what a caller wants.

    Returns
    -------
    tuple of (str or None, str or None)
        The digest to expect and where it came from — ``"table"`` for the curated row,
        ``"record"`` for the completion record — or ``(None, None)`` when neither pins one
        and there is nothing to check against.
    """
    pinned = metadata.sha256
    if pinned is not None:
        return pinned, EXPECTED_FROM_TABLE
    record = assembly_dir.read_record()
    if record is not None and record.sha256 is not None:
        return record.sha256, EXPECTED_FROM_RECORD
    return None, None


class UCSCGenomeDownloader(AssemblyRegistration):
    """Download reference-genome FASTA files from the UCSC golden path.

    UCSC serves per-assembly downloads under
    ``https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/bigZips/``. This
    fetches the soft-masked, gzipped whole-genome FASTA
    (``<assembly>.fa.gz``) and, by default, decompresses it to
    ``<assembly>.fa``.

    An :class:`~genome.io.registration.AssemblyRegistration` whose FASTA arrives over
    the network: the base owns the assembly's directory and the steps that finish a
    registration in it, and everything added here is about *where the bytes come from*
    — the URL, the name check, the pinned digest. The fetch itself is
    :func:`~genome.io.fetch.fetch_url`, given this assembly's working area: where a
    download lands is a decision the assembly's directory already made, so nothing here
    binds it to a cache of its own.

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
    metadata : genome.metadata.AssemblyMetadata
        The record this downloader works from — the one passed in, else what
        :func:`~genome.metadata.assembly_metadata` knows about ``assembly``. **Total**:
        an assembly the table does not list has a record whose every identifier is
        unknown, never no record, so nothing here asks whether there is one before
        reading a field off it. Whether the table *lists* a name is the other question
        and is not asked here at all (ADR-0003).

    Examples
    --------
    >>> dl = UCSCGenomeDownloader("hg38")
    >>> dl.fasta_url
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz'
    >>> UCSCGenomeDownloader("no_such_assembly").metadata.sha256 is None
    True
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
        super().__init__(assembly, cache_dir)
        self.metadata: AssemblyMetadata = (
            metadata if metadata is not None else assembly_metadata(assembly)
        )

    @property
    def assembly_url(self) -> str:
        """URL of the UCSC golden-path directory for this assembly."""
        return f"{self.BASE_URL}/{self.assembly}/"

    @property
    def _golden_path_fasta_url(self) -> str:
        """The FASTA URL derived from the assembly name and UCSC's golden-path layout."""
        return f"{self.BASE_URL}/{self.assembly}/bigZips/{self.assembly}.fa.gz"

    @property
    def _expected_sha256(self) -> str | None:
        """The sha256 this assembly's metadata pins for the unpacked FASTA, or ``None``."""
        return self.metadata.sha256

    def _fetched_source(self) -> FetchedSource:
        """Where this assembly's FASTA is downloaded from, and whether that URL was derived.

        Asked without reading the disk, so a caller that only wants the URL pays nothing:
        this is the fourth of :func:`~genome.io.source.resolve_source`'s checks on its own.
        """
        return fetched_source(self.metadata, self._golden_path_fasta_url)

    def _source(self) -> FetchedSource | ComponentSource:
        """Resolve this assembly's **Source** — see :func:`~genome.io.source.resolve_source`."""
        return resolve_source(
            self.dir, metadata=self.metadata, golden_path_url=self._golden_path_fasta_url
        )

    @property
    def fasta_url(self) -> str:
        """URL of the gzipped whole-genome FASTA — the pinned source, else the golden path.

        A metadata row that pins a source URL answers this outright; otherwise the URL
        is derived from the assembly name and UCSC's golden-path layout, as it always
        was.
        """
        return self._fetched_source().url

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
            decompression); see :func:`~genome.io.fetch.fetch_url`. Unrelated to the
            metadata row's ``sha256``, which covers the unpacked FASTA — see
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
        source = self._fetched_source()
        if source.derived:
            self.validate_assembly()
        processor: _Processor | None = (
            pooch.Decompress(method="gzip", name=f"{self.assembly}.fa") if decompress else None
        )
        # The archive is named after the assembly, not after the URL: a pinned source
        # need not be UCSC, so its file name says nothing about which assembly this is.
        return fetch.fetch_url(
            source.url,
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
        r"""Prepare the reference genome in one call — by download, or by concatenation.

        The choke point every way in funnels through, and therefore where a name is
        resolved: :func:`~genome.io.source.resolve_source` says whether this assembly's
        FASTA is fetched or concatenated from components already on this disk, and only the
        first of those is what the rest of this method does. A
        :class:`~genome.io.source.ComponentSource` is handed to
        :func:`~genome.io.chimera.build_chimera` instead, which is why ``genome register
        <name> --force`` is one command for all three kinds of **Source** rather than a
        download that fails on two of them.

        A finished registration is returned from its record either way, and a chimera's
        components are checked against their own records as it is handed back — the one
        failure a digest of this assembly's own bytes cannot show, since a component
        registered again underneath leaves those bytes untouched and no longer a copy of
        anything that exists. An assembly with no components pays nothing for the question.

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
            see :func:`~genome.io.fetch.fetch_url`. When ``None``, pooch verifies nothing —
            which
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
            If its record disagrees with what is on disk, or if a component of this
            chimera was registered again since it was built.
        FileNotFoundError
            If the name resolves to a chimera whose components are not all prepared here,
            or spells them in an order that is not the canonical one. Both messages name
            the command to run instead.
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
        source = self._source()
        registered = self._completed_genome(overwrite=overwrite, repair=self._repair_command())
        if registered is not None:
            # Asked of every assembly and answered instantly for one with no components,
            # so that opening a chimera by name is held to what building one is held to —
            # the check used to be reachable only through the builder and the verifier.
            # The refusal is the point here; the answer is for a surface that prints one.
            components_status(self.dir)
            return registered
        if isinstance(source, ComponentSource):
            # The one deferred import left here, and it is the layering and not a dodge:
            # building a chimera opens each component as a whole `Genome`, which is the top
            # of the stack. `genome.io.source` answers what a name *is* without any of that,
            # which is why the resolution above needs no such apology.
            from genome.io.chimera import build_chimera

            return build_chimera(self.dir, source.components, overwrite=overwrite)
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
        self._record_completion(files, source_url=source.url, sha256=digest)
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
        :func:`~genome.io.fetch.fetch_url`, so http(s), ftp and sftp all work). Gzipped sources
        (``.gz``) are decompressed. The resulting ``<assembly>.fa`` is then
        indexed/2bit/chrom.sizes-prepared exactly as :meth:`fetch_genome` does, and a
        completion record is written last, recording the source it was given and the
        digest of what arrived. UCSC is never contacted. Neither is the metadata row's
        pinned source or checksum: a seeded FASTA is whatever the caller handed over —
        it is recorded, never compared — and the assembly name is only a label for the
        directory it lands in.

        The :class:`~genome.io.source.SeededSource` kind, and the one of the three nothing
        resolves: the caller answered before anything was asked, so no name is read and no
        record is consulted for what this assembly *is*. Its record is still what says
        whether the work is already done — that check runs first, exactly as it does for a
        name that had to be resolved.

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
        seeded = SeededSource(source)
        registered = self._completed_genome(
            overwrite=overwrite, repair=self._repair_command(seeded.location)
        )
        if registered is not None:
            return registered
        fasta = self._materialize_fasta(
            seeded.location, progressbar=progressbar, overwrite=overwrite
        )
        files = prepare_fasta(fasta, overwrite=overwrite)
        self._record_completion(files, source_url=str(seeded.location), sha256=sha256_file(fasta))
        return files

    def _materialize_fasta(
        self,
        source: str | Path,
        *,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> Path:
        """Place ``source`` at ``<assembly>.fa`` in :attr:`cache_dir` and return that path.

        Copies a local file or fetches a URL through :func:`~genome.io.fetch.fetch_url` into the working
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
            fetch.fetch_url(src, work, fname=name, progressbar=progressbar)
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
) -> RegisteredAssembly:
    """Prepare ``assembly`` on disk and return the record of what that did.

    Naming an assembly is enough: where its FASTA comes from and which digest it must
    match are the metadata table's to know. The whole pipeline runs — fetch, unpack,
    verify, index, derive — and the completion record lands last. An assembly that is
    already registered is returned from its record without fetching anything.

    A directory that cannot be trusted raises instead (see :meth:`
    UCSCGenomeDownloader.fetch_genome`); ``force=True`` is what repairs one, and it
    keeps an unpacked FASTA it can prove is the pinned one rather than downloading a
    whole genome again.

    **The name is the whole interface**, chimeras included: a name whose parts are a
    prepared or listed assembly each is concatenated from those components instead of
    being fetched, so this one call prepares all three kinds of **Source** and ``force``
    repairs all three. See :func:`~genome.io.source.resolve_source` for the order the
    checks run in and what each one settles.

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
    RegisteredAssembly
        The completion record the run wrote — ``files``, ``source_url``, ``sha256``,
        ``tool_versions``, ``completed_at`` and the rest — with the ``assembly`` and the
        ``directory`` it lives in. :meth:`RegisteredAssembly.as_json` serializes it.

    Raises
    ------
    genome.io.completion.RegistrationError
        If the directory holds a build that cannot be trusted as finished, or (with
        ``force``) if the run somehow left no record behind.
    FileNotFoundError
        If the name resolves to a chimera this machine cannot build — a component that is
        not prepared here, or the components in an order that is not the canonical one.
    genome.io.utils.ChecksumMismatchError
        If the metadata pins a sha256 and the unpacked FASTA is not it.
    genome.external.ToolNotFoundError
        If ``samtools``, ``faToTwoBit`` or ``twoBitInfo`` are not on ``PATH``.

    Examples
    --------
    >>> register_assembly("sacCer3")                       # doctest: +SKIP
    RegisteredAssembly(assembly='sacCer3', directory=PosixPath('...'), record=...)
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
            f"`{assembly_repair_command(assembly, source)}`."
        )
    return RegisteredAssembly(assembly=assembly, directory=downloader.cache_dir, record=record)


def verify_assembly(
    assembly: str,
    *,
    fasta: str | Path | None = None,
    cache_dir: str | Path | None = None,
    metadata: AssemblyMetadata | None = None,
) -> VerifiedAssembly:
    """Re-read a FASTA and check its sha256 against the digest expected of it.

    The one operation that reads bytes rather than sizes. Registering an assembly and
    reopening it both go by presence and size, which is what makes them instant; this
    is the deliberate re-verification for when integrity is actually in doubt, and it
    costs a full pass over the file.

    What is expected of it comes from the assembly's curated row, and **failing that,
    from the completion record its own registration wrote** — a fallback rather than a
    question about what kind of assembly this is (see :func:`_expected_digest`).
    ``expected_from`` says which answered, because being held to a pinned digest and being
    held only to what this machine last produced are different results and a caller must be
    able to tell them apart.

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
    VerifiedAssembly
        The digest computed, what it was held to and what supplied that, and — for a
        chimera — what comparing the components settled. A digest that disagrees raises
        rather than reporting :attr:`~VerifiedAssembly.verified` ``False``; components
        that disagree likewise.

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
    >>> verify_assembly("sacCer3").verified                  # doctest: +SKIP
    True
    >>> verify_assembly("sacCer3", fasta="/tmp/copied.fa")  # doctest: +SKIP
    VerifiedAssembly(assembly='sacCer3', fasta=PosixPath('/tmp/copied.fa'), ...)
    """
    components: str | None = None
    assembly_dir = AssemblyDir.locate(assembly, cache_dir)
    row = metadata if metadata is not None else assembly_metadata(assembly)
    if fasta is not None:
        target = Path(fasta).expanduser()
        if not target.is_file():
            raise FileNotFoundError(
                f"no FASTA at {target}: pass the path of a file to check against "
                f"{assembly}'s row, or omit it to check {assembly}'s own registered FASTA."
            )
    else:
        target = assembly_dir.genome_files.fasta
        registered = assembly_dir.completed_files(repair=assembly_repair_command(assembly))
        if registered is None or not target.is_file():
            raise FileNotFoundError(
                f"{assembly} is not registered in {assembly_dir.path}, so there is "
                f"nothing to verify. Register it with `genome register {assembly}`, or "
                f"pass the FASTA to check with --fasta."
            )
        # Beside the digest, never instead of it: this one reads records rather than
        # bytes, and answers what a digest of this assembly's own bytes cannot. Its
        # answer is reported as well as enforced, so that "nothing was comparable" is
        # never handed back looking like "everything agreed".
        components = components_status(assembly_dir)
    expected, expected_from = _expected_digest(assembly_dir, row)
    actual = sha256_file(target)
    if expected is not None and actual != expected:
        raise ChecksumMismatchError(target, expected, actual)
    return VerifiedAssembly(
        assembly=assembly,
        fasta=target,
        sha256=actual,
        expected=expected,
        expected_from=expected_from,
        components=components,
    )


def assembly_table_row(
    assembly: str,
    *,
    cache_dir: str | Path | None = None,
    progressbar: bool = True,
) -> AssemblyMetadata:
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
    genome.metadata.AssemblyMetadata
        The row itself, with ``None`` for anything still unknown — the same type the
        curated table parses into, since that is what this computes. Hand
        :func:`dataclasses.asdict` of it to
        :func:`~genome.metadata.format_table_row` for the line to paste into
        ``data/assembly_metadata.tsv``.

    A **Chimera** is **refused before anything is fetched**, because it has no work here to
    do: its row pins nothing, so there is no digest to compute and no source to record, and
    the refusal describes that row rather than printing it (ADR-0008).

    Raises
    ------
    ValueError
        If the URL had to be derived and the assembly is unknown to UCSC, or if the name
        resolves to a chimera, which has no row to compute.
    FileNotFoundError
        If the name spells a chimera's components in an order that is not the canonical
        one; the message names the spelling to use.

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
    >>> from dataclasses import asdict
    >>> from genome.metadata import format_table_row
    >>> format_table_row(asdict(assembly_table_row("sacCer3")))   # doctest: +SKIP
    'sacCer3\tSaccharomyces cerevisiae\t...\t6ff72f07...'
    """
    downloader = UCSCGenomeDownloader(assembly, cache_dir)
    source = downloader._source()
    if isinstance(source, ComponentSource):
        raise ValueError(
            f"{assembly} is a chimera, of {', '.join(source.components)}, so there is no row for "
            f"this command to compute and nothing was downloaded. Its row is one line "
            f"carrying the name and nothing else: no source URL, because nothing is "
            f"fetched, and no sha256, because a chimera's bytes are derived by this package "
            f"from components that are themselves pinned, and pinning them again would turn "
            f"our own concatenation into a contract that fails on every disk (ADR-0008). "
            f"Build it with `genome register {assembly}` and check it with `genome verify "
            f"{assembly}`, which compares the components rather than a pin."
        )
    fasta = downloader.fetch_fasta(progressbar=progressbar)
    # Every identifier the row knows, and all of them blank when the table lists no row at
    # all, in which case the name is the only one that does not need a person to supply it.
    # The downloader's record is total, so there is no absent one to stand in for here.
    # The digest is reported, never enforced: this is the command a maintainer runs
    # *because* the pinned one needs regenerating, so comparing against the stale one
    # would refuse exactly when it is needed. Checking a FASTA against the official row
    # is what verifying an assembly is for.
    return replace(downloader.metadata, source_url=downloader.fasta_url, sha256=sha256_file(fasta))


if __name__ == "__main__":
    downloader = UCSCGenomeDownloader("sacCer3")
    files = downloader.fetch_genome()
    print(files)
