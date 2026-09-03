"""One **Prepared set**: the fetch, the working area, the digest and the record.

I/O boundary module, and the one implementation of a pipeline three contexts each own an
instance of — a **Motif set**, an **Xref set** and a **Homology set**. Each is files pinned
to a **Release**, belonging to no **Assembly**, filed beside the assembly tree under the
**Data dir** and prepared by construction; each was written out again where it was needed,
over the same skeleton, three times.

**A source declares three things and the rest is here.** A URL, a checksum, and how to
slice or parse what arrives — :class:`PreparedSource` names them beside the directory the
set lives in and the command that prepares it. What this module owns is everything around
them: the working area, the fetch, the digest, the staged rename, the **Completion marker**
and the one sentence that sends a caller to a login node. Adding a fourth pinned set is
therefore a URL, a checksum, a directory name and a reader — no fetch, no digest, no marker
and no pipeline of its own.

**Where each set is filed is its own context's.** The root a set lives under is one line —
:func:`~genome.store.data_dir.prepared_data_dir` of a name — and it is written where the
code that reads the set is, so ``motif/`` is spelled in :mod:`genome.tf.motif.jaspar`,
``xref/`` in :mod:`genome.xref.xref` and ``homology/`` in :mod:`genome.homology.compara`.
Nothing here knows the three names, which is what stops this module from acquiring a
fourth when a fourth set arrives.

**Every prepared set writes a marker, with no exception.** It is the only thing that says a
build finished, and it is also the only answer to *how was this made* — the URL, the digest,
the package version and whatever the reader measured. A set whose files are small enough to
arrive in one step still gets one, because a directory that answers *is this finished* two
different ways is what this module exists to end.

**Where a checksum is enforced follows what it covers** (ADR-0006). A pin over the
**unpacked** bytes is checked here, as the file is streamed and unpacked, which costs one
pass rather than two — pooch would compare it against the compressed bytes and reject every
download. A pin over the **archive** as served is pooch's own ``known_hash``, checked as the
bytes arrive; that is where Ensembl Compara's published md5 belongs, since its ``MD5SUM``
covers the ``.gz``. A source that pins nothing says so, and the digest of what it stored is
recorded rather than compared.

**The working area is kept only when something can vouch for what is in it.** With an
archive pin, pooch re-checks a leftover download and fetches again when it does not match,
so an interrupted run costs no second download of a 110 MB dump. With no pin over the bytes
as served, a leftover is unverifiable and would be adopted as a finished download — so it is
cleared before the fetch rather than picked up.

**Reach the fetch step through its module**, never by importing the name: this is the
package's third caller of :func:`genome.store.fetch.fetch_url` and it spells it
``fetch.fetch_url(...)`` for the same reason the other two do — one
``monkeypatch.setattr`` on that module takes the whole package offline, and a module
holding an imported reference would keep downloading while the suite believed itself
offline.

Examples
--------
>>> from genome.store import prepared
>>> prepared.login_node_help("genome motif scan --help")[:24]
'Nothing else in this pac'
"""

from __future__ import annotations

import gzip
import hashlib
import shlex
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, field
from io import BufferedIOBase
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol
from urllib.parse import urlparse

from genome.store import fetch
from genome.store.completion import (
    CompletionRecord,
    build_record,
    check_registration,
    clear_work_dir,
    work_dir,
    write_record,
)

if TYPE_CHECKING:  # pragma: no cover - typeshed's name for a live hash object
    from hashlib import _Hash as Hash

#: A checksum taken over the file's **unpacked** content, which is what a pin normally
#: covers here (ADR-0006) and what is verified as the bytes are streamed.
UNPACKED = "unpacked"

#: A checksum taken over the archive **as served**, which is what pooch's ``known_hash``
#: compares and where a publisher's own md5 of a ``.gz`` belongs.
ARCHIVE = "archive"

#: What a checksum may cover. Two, because two is what publishers actually pin.
CHECKSUM_SCOPES: tuple[str, ...] = (UNPACKED, ARCHIVE)

#: What the stored form is called while it is still being written, inside the working area.
#: A suffix rather than the stored name itself, because a set that stores what it fetched
#: would otherwise write over the file it is reading.
PART_SUFFIX = ".part"


class PreparedSetNotDownloadedError(RuntimeError):
    """A **Prepared set** is not on disk here and its bytes could not be fetched.

    A :class:`RuntimeError`, because nothing about the call was wrong: the bytes are simply
    not here and this machine could not go and get them. Each context raises a subclass of
    its own naming its own set and quoting its own prepare command — three specific
    messages rather than one vague one — and the sentence that sends the caller to a login
    node is :func:`login_node_help`, written once.

    Examples
    --------
    >>> issubclass(PreparedSetNotDownloadedError, RuntimeError)
    True
    """


class PreparedChecksumError(ValueError):
    """What arrived is not what the set pins, so nothing is prepared from it.

    A :class:`ValueError`: the bytes are a bad value, and a truncated download is not a
    smaller release. The message names the file, both digests, which bytes the pin covers
    and the command that prepares the set again.

    Examples
    --------
    >>> issubclass(PreparedChecksumError, ValueError)
    True
    """


class PreparedDecodeError(ValueError):
    """The publisher's file is not text, and what reads it here yields text.

    A :class:`ValueError`, for the same reason the checksum error is one: the bytes are a
    bad value for what is being asked of them. Every reader in this package is handed
    decoded lines and all three publishers ship text, so a publisher that ships bytes is a
    design question — whether a source should declare its form, and whether the reader
    protocol should carry bytes — rather than something a caller repairs on disk. The
    message names the file so that the wall is a wall and not a decoder's own error
    surfacing from mid-stream naming nothing. What a *reader* stores is not this: a stored
    form is digested as bytes and never decoded.

    Examples
    --------
    >>> issubclass(PreparedDecodeError, ValueError)
    True
    """


def login_node_help(command: str) -> str:
    """Return the message that sends a caller to a login node, quoting ``command``.

    Written once and quoted by every **Prepared set**'s own not-downloaded error, because
    the fact is one fact: fetching is the only step in this package that needs the network
    and the lab's CPU cluster compute nodes have none.

    Parameters
    ----------
    command : str
        The call that prepares this set, spelled by the set itself so the message names the
        exact repair rather than a general one.

    Returns
    -------
    str
        Two sentences: what needs the network, then where to run ``command`` — quoted —
        and where the set is read from once it has.

    Examples
    --------
    >>> message = login_node_help("genome xref --help")
    >>> "`genome xref --help`" in message
    True
    >>> message.endswith("shared by every project on the machine.")
    True
    """
    return (
        f"Nothing else in this package needs the network, so this is the one step that does. "
        f"Prepare it on a machine with internet — a login node, since the lab's compute nodes "
        f"have none — with `{command}`, after which it is read from the Data dir and shared by "
        f"every project on the machine."
    )


@dataclass(frozen=True)
class Checksum:
    """A digest a fetched file is held to, and which of its bytes the digest covers.

    Attributes
    ----------
    algorithm : str
        The hash algorithm the publisher used, as :func:`hashlib.new` names it — ``"md5"``
        for the two publishers pinned here, whose own files say so.
    digest : str
        The expected hex digest.
    covers : str, default ``"unpacked"``
        :data:`UNPACKED` or :data:`ARCHIVE` — which bytes the digest was taken over, and so
        where it is enforced. Unpacked is the norm (ADR-0006); archive is what a publisher's
        md5 of its own ``.gz`` covers.

    Examples
    --------
    >>> Checksum.parse("md5:d41d8cd98f00b204e9800998ecf8427e")
    Checksum(algorithm='md5', digest='d41d8cd98f00b204e9800998ecf8427e', covers='unpacked')
    >>> str(Checksum("md5", "d41d8cd98f00b204e9800998ecf8427e", covers=ARCHIVE))
    'md5:d41d8cd98f00b204e9800998ecf8427e'
    """

    algorithm: str
    digest: str
    covers: str = UNPACKED

    @classmethod
    def parse(cls, pinned: str, *, covers: str = UNPACKED) -> Checksum:
        """Read a ``"<algorithm>:<hexdigest>"`` pin, algorithm and all.

        Parameters
        ----------
        pinned : str
            The pin as a curated table spells it, algorithm first so that it travels with
            the digest.
        covers : str, default ``"unpacked"``
            Which bytes it was taken over.

        Returns
        -------
        Checksum
            The parsed pin.

        Raises
        ------
        ValueError
            If the algorithm is missing, or ``covers`` is not one of
            :data:`CHECKSUM_SCOPES`.

        Examples
        --------
        >>> Checksum.parse("sha256:" + "0" * 64).algorithm
        'sha256'
        """
        algorithm, separator, digest = pinned.partition(":")
        if not separator:
            raise ValueError(
                f"the checksum {pinned!r} names no algorithm: a pin is spelled "
                f"'<algorithm>:<hexdigest>' so that the algorithm travels with the digest. "
                f"Fix the cell that declares it."
            )
        if covers not in CHECKSUM_SCOPES:
            raise ValueError(
                f"a checksum covers {' or '.join(CHECKSUM_SCOPES)} bytes, not {covers!r}. "
                f"Which of the two decides where the pin is enforced, so it is named rather "
                f"than guessed."
            )
        return cls(algorithm=algorithm, digest=digest, covers=covers)

    def __str__(self) -> str:
        """Return the pin as ``"<algorithm>:<hexdigest>"``, which is what pooch takes."""
        return f"{self.algorithm}:{self.digest}"


class SourceReader(Protocol):
    """Turns the publisher's unpacked lines into the file a **Prepared set** is read from.

    What a source *is*, beside its URL and its checksum: a function from the publisher's
    lines to the stored form, which is what makes adding a set data plus a reader rather
    than a refactor. It fetches nothing and verifies nothing — by the time one is called the
    bytes are on disk and are being streamed past it.
    """

    def __call__(self, lines: Iterator[str], staged: Path, *, origin: str) -> Mapping[str, Any]:
        """Write the stored form at ``staged``, and return what the marker should record.

        ``staged`` is inside the working area, so the set's final name stays unoccupied
        until the whole file exists. ``origin`` is what the lines came out of, for a
        message that has to name the file it refused. Everything returned lands in the
        **Completion marker**'s details, which is where a measurement taken while slicing
        belongs — how many rows were kept, which columns held nothing.
        """
        ...


@dataclass(frozen=True)
class PreparedSource:
    """What one **Prepared set** declares, and the whole of what it declares.

    Attributes
    ----------
    url : str
        Where the publisher's file is fetched from.
    directory : pathlib.Path
        The set's own directory. One set per directory, because each carries a **Completion
        marker** of its own.
    stored_name : str
        The file the set is read from, inside :attr:`directory`.
    kind : str
        What the marker calls what it recorded — ``"motif"``, ``"xref"``, ``"homology"``.
    name : str
        The name the set is addressed by in its own tree.
    prepare_command : str
        The call that prepares this set, in the one spelling the context keeps, quoted
        verbatim into every error a caller repairs by running it.
    description : str
        How an error names this set — ``"the alliance 9.0.0 xref set for 'Homo sapiens'"``.
    read : SourceReader
        How to slice or parse what arrives.
    not_downloaded : type of PreparedSetNotDownloadedError
        The context's own class for *the bytes are not here and could not be fetched*.
    checksum : Checksum, optional
        What the fetched file must hash to, and which of its bytes that covers. ``None``
        for a publisher that pins nothing, whose stored digest is then recorded rather than
        compared.
    download_name : str, optional
        What the fetched file is called inside the working area. The URL's own last segment
        when omitted; named when two sets of one release would otherwise collide on it.
    details : mapping, optional
        What the marker records about the set beyond the pipeline's own fields — the
        publisher, the release, the species. The reader's measurements are merged over it.

    Examples
    --------
    >>> from pathlib import Path
    >>> source = PreparedSource(
    ...     url="https://example.org/beasts.tsv.gz",
    ...     directory=Path("/scratch/liulab/beast/1.0"),
    ...     stored_name="beasts.tsv.gz",
    ...     kind="beast",
    ...     name="beasts/1.0",
    ...     prepare_command="python -c 'from beasts import Beasts; Beasts()'",
    ...     description="the beast 1.0 set",
    ...     read=write_through,
    ...     not_downloaded=PreparedSetNotDownloadedError,
    ... )
    >>> source.path
    PosixPath('/scratch/liulab/beast/1.0/beasts.tsv.gz')
    >>> source.repair.startswith("rm -rf /scratch/liulab/beast/1.0 && ")
    True
    """

    url: str
    directory: Path
    stored_name: str
    kind: str
    name: str
    prepare_command: str
    description: str
    read: SourceReader
    not_downloaded: type[PreparedSetNotDownloadedError]
    checksum: Checksum | None = None
    download_name: str | None = None
    details: Mapping[str, Any] = field(default_factory=dict)

    @property
    def path(self) -> Path:
        """Where the prepared set is read from, whether or not it exists yet."""
        return self.directory / self.stored_name

    @property
    def repair(self) -> str:
        """The command that prepares this set again from nothing.

        Both halves, because deleting the directory on its own leaves the caller with
        nothing and no way back — and the second half is the call that needs a login node,
        which is worth reading before the first half has run.
        """
        return f"rm -rf {shlex.quote(str(self.directory))} && {self.prepare_command}"


@dataclass(frozen=True)
class Prepared:
    """One **Prepared set** on disk: the file it is read from, and what its marker says.

    Attributes
    ----------
    path : pathlib.Path
        The stored file, which exists.
    record : genome.store.completion.CompletionRecord
        The marker beside it — the one already there when nothing was prepared, or the one
        just written.

    Examples
    --------
    >>> prepare(source)                                  # doctest: +SKIP
    Prepared(path=PosixPath('...'), record=CompletionRecord(kind='xref', ...))
    """

    path: Path
    record: CompletionRecord


def prepare(source: PreparedSource, *, progressbar: bool = True) -> Prepared:
    """Return the finished set, fetching and slicing it once if it is not there.

    The whole pipeline, in the order every build in this package takes it: look for a
    marker, fetch into the working area, slice or parse into a staged file, hold the bytes
    to what the source pins, digest the stored form while it is still staged, place it
    under its final name, write the marker, and only then empty the working area. The
    marker is written last, after the file it claims exists, so a run killed anywhere in
    between leaves a directory that reads as unfinished rather than as present.

    **Everything that can fail happens while the directory is still fresh.** What follows
    the move is one ``stat`` and one small record — the digest, which is a full pass over a
    file that may be very large, is taken before it. An interrupted run leaving a claimed
    file with no marker is a state this package refuses to guess at and asks to be rebuilt
    (ADR-0007); a step of its own falling over after the move would manufacture that same
    state from a foreseeable failure, which is a different thing and not one to accept.

    Parameters
    ----------
    source : PreparedSource
        What this set declares.
    progressbar : bool, default True
        Show the download's progress bar. Nothing is drawn when the set is already there,
        since nothing is fetched.

    Returns
    -------
    Prepared
        The stored file and its **Completion marker**.

    Raises
    ------
    genome.store.completion.UnfinishedRegistrationError
        If the directory holds files but no marker — an interrupted run. The message names
        :attr:`PreparedSource.repair`.
    genome.store.completion.RegistrationMismatchError
        If the marker disagrees with what is on disk.
    PreparedSetNotDownloadedError
        As the subclass the source named, if the bytes are not here and could not be
        fetched.
    PreparedChecksumError
        If a pin over the unpacked bytes does not match what arrived. A pin over the archive
        is pooch's to enforce and raises there instead.
    PreparedDecodeError
        If the publisher's file is not valid UTF-8. Reading it yields text, which is what
        every reader here takes; what a reader *stores* is digested as bytes and is not
        this.

    Examples
    --------
    >>> prepare(source, progressbar=False)               # doctest: +SKIP
    Prepared(path=PosixPath('/scratch/liulab/xref/alliance/9.0.0/homo_sapiens/...'), ...)
    """
    directory = source.directory
    existing = check_registration(directory, repair=source.repair)
    if existing is not None:
        return Prepared(path=source.path, record=existing)

    work = work_dir(directory)
    if source.checksum is None or source.checksum.covers != ARCHIVE:
        # Nothing here could vouch for what an interrupted run left, and pooch serves a
        # file already sitting at the destination — so it is swept rather than adopted.
        clear_work_dir(directory)
    fetched = _fetch(source, work=work, progressbar=progressbar)

    # A name of its own, never the stored one: a set that stores what it fetched under the
    # publisher's own name would otherwise have its reader truncate the file it is reading.
    staged = work / f"{source.stored_name}{PART_SUFFIX}"
    digest = _live_digest(source.checksum)
    measured = source.read(unpacked_lines(fetched, digest), staged, origin=str(fetched))
    _check_checksum(source, digest, path=fetched)
    # Digested here, where a failure costs the caller a re-run and not a wedged directory.
    # Packed-ness comes off the STORED name and never off the staged path: the staged file
    # wears the working suffix, so its own name says nothing, and reading it there would
    # start hashing gzip bytes for a set that stores a `.gz` — which is neither what its
    # marker records nor what a pin covers.
    stored_digest = unpacked_digest(staged, packed=_is_packed(source.stored_name))

    directory.mkdir(parents=True, exist_ok=True)
    staged.replace(source.path)
    record = build_record(
        directory,
        kind=source.kind,
        name=source.name,
        files=[source.path],
        source_url=source.url,
        sha256=stored_digest,
        details={**source.details, **measured},
    )
    write_record(directory, record)
    clear_work_dir(directory)
    return Prepared(path=source.path, record=record)


def unpacked_lines(path: Path, digest: Hash | None = None) -> Iterator[str]:
    r"""Yield ``path``'s unpacked lines, digesting the unpacked bytes as they go.

    The one implementation of decompress-while-hashing: stream, unpack, digest what came
    out, never hold the file. Gzip is undone here rather than in a reader, so a source that
    ships plain text needs no branch of its own, and the digest is fed the bytes exactly as
    the publisher wrote them — which is what makes it comparable with a pin taken over the
    unpacked content (ADR-0006).

    Parameters
    ----------
    path : pathlib.Path
        The file to read. Unpacked when its name ends in ``.gz``, read as it is otherwise.
    digest : hashlib hash object, optional
        Updated with every byte yielded. ``None`` to read without hashing.

    Yields
    ------
    str
        One line, its line ending kept, decoded as UTF-8.

    Raises
    ------
    PreparedDecodeError
        If the bytes are not valid UTF-8. This path yields text by decision: every reader
        here takes ``str``, so a publisher shipping bytes is refused by name rather than
        by a decoder's own error from mid-stream.

    Examples
    --------
    >>> import hashlib, gzip
    >>> from pathlib import Path
    >>> path = Path("/tmp/beasts.tsv.gz")               # doctest: +SKIP
    >>> digest = hashlib.sha256()                        # doctest: +SKIP
    >>> [line for line in unpacked_lines(path, digest)]  # doctest: +SKIP
    ['name\tcount\n', 'aardvark\t1\n']
    """
    with _unpacked(path, packed=_is_packed(path.name)) as stream:
        yield from _decoded(stream, digest, origin=path)


def _unpacked(path: Path, *, packed: bool) -> BufferedIOBase:
    """Open ``path``'s unpacked bytes, undoing gzip when ``packed``. The caller closes it."""
    return gzip.open(path, "rb") if packed else path.open("rb")


def _decoded(lines: Iterable[bytes], digest: Hash | None, *, origin: Path) -> Iterator[str]:
    """Decode each line, feeding the bytes to ``digest`` exactly as they were read."""
    for raw in lines:
        if digest is not None:
            digest.update(raw)
        try:
            line = raw.decode("utf-8")
        except UnicodeDecodeError as error:
            raise PreparedDecodeError(
                f"{origin} is not valid UTF-8 text: {error.reason}. Every reader here is "
                f"handed decoded lines, so a publisher that ships bytes cannot be read by "
                f"this package today. Check that the URL the source declares still serves "
                f"the text file it documents — a binary form served under a text name "
                f"reads exactly like this. If the publisher has genuinely moved to one, "
                f"reading it needs a reader that takes bytes, which is a change to this "
                f"package rather than anything to repair on disk."
            ) from error
        yield line


def _is_packed(name: str) -> bool:
    """Say whether ``name`` names a gzip, asked of a **name** because a path may not be one.

    The stored form is digested while it still wears the working suffix, and what its
    digest covers is decided by the name it will be read under, never by the staged path
    (ADR-0006).
    """
    return Path(name).suffix == ".gz"


def unpacked_digest(path: Path, algorithm: str = "sha256", *, packed: bool | None = None) -> str:
    """Return the digest of ``path``'s unpacked content, streaming it once.

    Bytes throughout, decoded nowhere: what a **Prepared set** stores is its reader's
    business, and a stored form that is not text is digested like any other. The file is
    never held, so a stored slice costs one pass over it.

    Parameters
    ----------
    path : pathlib.Path
        The file to digest.
    algorithm : str, default ``"sha256"``
        The hash algorithm, as :func:`hashlib.new` names it.
    packed : bool, optional
        Whether to unpack before digesting. ``None`` reads it off ``path``'s own name,
        which is what every caller holding the file under its final name wants. Name it
        for a file whose path does not carry it — a staged file wearing the working
        suffix, whose stored name is the authority on what its digest covers (ADR-0006).

    Returns
    -------
    str
        The lower-case hex digest of the unpacked bytes.

    Examples
    --------
    >>> from pathlib import Path
    >>> unpacked_digest(Path("tests/data/tiny.fa.gz"))   # doctest: +SKIP
    '9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981'

    A staged file, whose own name says nothing about what it holds:

    >>> import gzip, hashlib, tempfile
    >>> with tempfile.TemporaryDirectory() as scratch:
    ...     staged = Path(scratch) / "beasts.tsv.gz.part"
    ...     _ = staged.write_bytes(gzip.compress(b"aardvark"))
    ...     unpacked_digest(staged, packed=True) == hashlib.sha256(b"aardvark").hexdigest()
    True
    """
    unpack = _is_packed(path.name) if packed is None else packed
    with _unpacked(path, packed=unpack) as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def write_through(lines: Iterator[str], staged: Path, *, origin: str = "") -> Mapping[str, Any]:
    r"""Write the publisher's own bytes through unchanged — the reader that slices nothing.

    What a set whose stored form *is* the published file declares, so that storing the
    bytes as they arrived is still one path through the pipeline rather than a second one
    beside it.

    Parameters
    ----------
    lines : iterator of str
        The publisher's unpacked lines, line endings and all.
    staged : pathlib.Path
        Where to write them, inside the working area.
    origin : str, optional
        What the lines came out of. Unused: nothing here can refuse a line, so there is
        nothing to name.

    Returns
    -------
    mapping of str to object
        How many lines were written, for the marker to record.

    Examples
    --------
    >>> from pathlib import Path
    >>> write_through(iter(["one\\n", "two\\n"]), Path("/tmp/two-lines.txt"))
    {'lines': 2}
    """
    staged.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with staged.open("w", encoding="utf-8", newline="") as out:
        for line in lines:
            out.write(line)
            written += 1
    return {"lines": written}


def _fetch(source: PreparedSource, *, work: Path, progressbar: bool) -> Path:
    """Download the publisher's file into the working area, or say what to do instead.

    ``known_hash`` is handed to pooch for an archive pin and withheld for an unpacked one:
    the latter would be compared against the compressed bytes and would reject every
    download (ADR-0006), and it is checked as the file is read instead.
    """
    checksum = source.checksum
    known_hash = str(checksum) if checksum is not None and checksum.covers == ARCHIVE else None
    name = source.download_name or Path(urlparse(source.url).path).name or source.stored_name
    try:
        return fetch.fetch_url(
            source.url, work, known_hash=known_hash, fname=name, progressbar=progressbar
        )
    except OSError as error:
        # `OSError` alone: a pin that does not match names bytes that arrived, which is a
        # different fact needing a different next action, and re-running the fetch would
        # not repair it.
        raise source.not_downloaded(
            f"{source.description} is not prepared here and {source.url} could not be "
            f"fetched: {error}. {login_node_help(source.prepare_command)}"
        ) from error


def _live_digest(checksum: Checksum | None) -> Hash | None:
    """Return the hash a pin over the unpacked bytes is accumulated in, or ``None``."""
    if checksum is None or checksum.covers != UNPACKED:
        return None
    return hashlib.new(checksum.algorithm)


def _check_checksum(source: PreparedSource, digest: Hash | None, *, path: Path) -> None:
    """Hold what arrived to the digest the source pins over its unpacked bytes."""
    checksum = source.checksum
    if checksum is None or digest is None:
        return
    found = digest.hexdigest()
    if found == checksum.digest:
        return
    raise PreparedChecksumError(
        f"{path} hashes to {found} where {source.description} pins {checksum.digest} over its "
        f"{checksum.covers} bytes. A truncated download is not a smaller release, and preparing "
        f"from it would answer with silently less than the release holds. Prepare it again with "
        f"`{source.repair}`. If {source.url} has genuinely been re-published under the same "
        f"name, that source cannot be pinned and does not belong in a curated table."
    )
