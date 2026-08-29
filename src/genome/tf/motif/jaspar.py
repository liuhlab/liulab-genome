"""One JASPAR **Release** on disk, parsed into a **Motif set** you can query.

I/O boundary module, and the only one in this context: name a release and a **Tax group**
and the matching file is fetched once into the ``motif/`` subtree of the **Data dir**,
shared by every project on the machine, then read back as a :class:`JasparDatabase`. A
second construction re-reads what is there and fetches nothing. The types it hands back
are :mod:`genome.tf.motif.motif`'s, which reads no file and knows nothing about any of
this.

The file fetched is JASPAR's **transfac** serialization rather than the ``.jaspar`` one.
It carries the **Count matrix** and all six annotations in a single file and it does not
round counts to integers — nine of the 1019 vertebrate matrices in the 2026 release carry
fractional counts that ``.jaspar`` discards.

**There is no Completion marker here, and that is deliberate.** Every other multi-step
build in this package writes one, and a reader who knows that will otherwise read its
absence as an omission. These files are all under 1 MB and arrive in a single step, so
the two things a record would buy are bought more cheaply: the file is downloaded under a
temporary name in the working area and renamed into place only on success, so an
interrupted download never occupies the final name; and every construction — a fresh
download and a cache re-read alike — counts the motifs it parsed against
:data:`MOTIF_COUNTS`, so a short file raises rather than quietly scanning with a partial
release. There is nothing a record could add that those two do not already say.

Examples
--------
>>> from genome.tf.motif.jaspar import JasparDatabase, jaspar_url
>>> jaspar_url("2024", "nematodes")                  # doctest: +ELLIPSIS
'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_nematodes_...txt'
>>> database = JasparDatabase("2024", "vertebrates")             # doctest: +SKIP
>>> database["CTCF"]                                             # doctest: +SKIP
Motif(motif_id='MA0139.2', motif_name='CTCF', length=15, offset=0)
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from types import MappingProxyType

import numpy as np

from genome.io import fetch
from genome.io.completion import clear_work_dir, work_dir
from genome.io.registration import motif_data_dir
from genome.tf.motif.motif import BASES, Motif, MotifSet, base_id

#: Where JASPAR publishes its versioned flat files. The REST API is not used: it silently
#: ignores its release parameter and answers a historical query with current-release data,
#: so only these files can be pinned to a **Release**.
JASPAR_BASE_URL = "https://jaspar.elixir.no/download/data"

#: Subdirectory of the **Data dir**'s ``motif/`` tree holding JASPAR's files, flat. One
#: source, one directory: a second motif source would get its own beside this.
JASPAR_SUBDIR = "jaspar"

#: The **Release**s this package prepares, oldest first.
JASPAR_RELEASES: tuple[str, ...] = ("2024", "2026")

#: The **Release** used when none is named — the newest, so a new analysis starts on
#: current data and reproducing an old one is the call that says which.
DEFAULT_RELEASE = "2026"

#: JASPAR's **Tax group**s. ``all`` is exactly the union of the other seven, which
#: :data:`MOTIF_COUNTS` is checked against by the test suite.
JASPAR_TAX_GROUPS: tuple[str, ...] = (
    "vertebrates",
    "plants",
    "insects",
    "nematodes",
    "fungi",
    "urochordates",
    "diatoms",
    "all",
)

#: The **Tax group** used when none is named: the lab's common case is human and mouse.
DEFAULT_TAX_GROUP = "vertebrates"

#: How many motifs each ``(release, tax group)`` file holds, counted from the published
#: files themselves. **This is the integrity check**, standing where a **Completion
#: marker** stands elsewhere: a parse that yields any other number raises, so a truncated
#: download or a half-written cache is an error rather than a quiet partial release. Note
#: ``diatoms``, which really does hold one motif in both releases.
MOTIF_COUNTS: Mapping[tuple[str, str], int] = MappingProxyType(
    {
        ("2024", "vertebrates"): 879,
        ("2024", "plants"): 805,
        ("2024", "insects"): 286,
        ("2024", "nematodes"): 103,
        ("2024", "fungi"): 178,
        ("2024", "urochordates"): 94,
        ("2024", "diatoms"): 1,
        ("2024", "all"): 2346,
        ("2026", "vertebrates"): 1019,
        ("2026", "plants"): 927,
        ("2026", "insects"): 296,
        ("2026", "nematodes"): 103,
        ("2026", "fungi"): 193,
        ("2026", "urochordates"): 94,
        ("2026", "diatoms"): 1,
        ("2026", "all"): 2633,
    }
)

#: What separates one annotation value from the next in a transfac ``CC`` line — a
#: semicolon and never a comma. Commas live *inside* single values, in a class
#: (``"C3H(C),C2HC zinc-fingers like factors"``), a family (``"Zinc finger, BED-type"``)
#: and a data type (``"PBM, CSA and/or DIP-chip"``), so splitting on one corrupts the
#: annotation of roughly fifty records per release without failing.
_VALUE_SEPARATOR = ";"

#: The transfac line prefixes this parser reads. Anything else — ``DE``, which repeats
#: the accession and the name, and the ``XX`` spacers — is skipped.
_ACCESSION, _NAME, _ANNOTATION, _HEADER, _END = "AC ", "ID ", "CC ", "PO", "//"


class TransfacError(ValueError):
    r"""A transfac record cannot be read, so no motif is made from it.

    A bad *file*, not a bad call: the message names the record it stopped on — by its
    accession, or by its position when it has none — and what was wrong with it.

    Examples
    --------
    >>> try:
    ...     parse_transfac("AC MA0001.1\nXX\nID x\nXX\n//\n")
    ... except TransfacError as error:
    ...     print("MA0001.1" in str(error))
    True
    """


class JasparReleaseError(ValueError):
    """A file read as a **Release** is not the release it should be.

    Not a parse failure — every record read cleanly — but the file as a whole is wrong:
    the wrong number of motifs, or two versions of one matrix where a non-redundant
    release ships exactly one. Both mean the same thing to a caller, which is why they are
    one class: what is on disk is not what was asked for, and the repair is to delete it
    and construct again. The message names the file, so the repair is one ``rm``.

    Examples
    --------
    >>> from pathlib import Path
    >>> try:
    ...     _check_count((), release="2024", tax_group="diatoms", path=Path("/tmp/x.txt"))
    ... except JasparReleaseError as error:
    ...     print("holds 0 motifs" in str(error))
    True
    """


def jaspar_url(release: str, tax_group: str = DEFAULT_TAX_GROUP) -> str:
    """Return the URL of one **Release**'s transfac file for one **Tax group**.

    Parameters
    ----------
    release : str
        One of :data:`JASPAR_RELEASES`.
    tax_group : str, default ``"vertebrates"``
        One of :data:`JASPAR_TAX_GROUPS`.

    Returns
    -------
    str
        The published file's URL.

    Raises
    ------
    ValueError
        If the release or the tax group is not one this package prepares.

    Notes
    -----
    The union file's own name carries **no** tax group segment where the seven named
    groups do, which is the one irregularity in JASPAR's naming and the reason this is a
    function rather than a format string spelled at each call site.

    Examples
    --------
    >>> jaspar_url("2026", "insects")
    'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_insects_non-redundant_pfms_transfac.txt'
    >>> jaspar_url("2026", "all")
    'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_transfac.txt'
    """
    release = check_release(release)
    tax_group = check_tax_group(tax_group)
    taxon = "" if tax_group == "all" else f"{tax_group}_"
    return (
        f"{JASPAR_BASE_URL}/{release}/CORE/"
        f"JASPAR{release}_CORE_{taxon}non-redundant_pfms_transfac.txt"
    )


def jaspar_filename(release: str, tax_group: str = DEFAULT_TAX_GROUP) -> str:
    """Return the name one **Release** and **Tax group** is cached under.

    JASPAR's own name for the seven named groups, and *not* for the union: the published
    union file drops the tax group from its name, and a cache directory holding every
    release side by side needs both halves of what a file is on every one of them.

    Parameters
    ----------
    release : str
        One of :data:`JASPAR_RELEASES`.
    tax_group : str, default ``"vertebrates"``
        One of :data:`JASPAR_TAX_GROUPS`.

    Returns
    -------
    str
        The local file name, carrying the release and the tax group.

    Examples
    --------
    >>> jaspar_filename("2024", "nematodes")
    'JASPAR2024_CORE_nematodes_non-redundant_pfms_transfac.txt'
    >>> jaspar_filename("2024", "all")
    'JASPAR2024_CORE_all_non-redundant_pfms_transfac.txt'
    """
    return (
        f"JASPAR{check_release(release)}_CORE_{check_tax_group(tax_group)}"
        f"_non-redundant_pfms_transfac.txt"
    )


def jaspar_data_dir() -> Path:
    """Return the directory JASPAR's files are cached in, flat and shared.

    ``<liulab_data>/motif/jaspar/``. Nothing is created by asking.

    Returns
    -------
    pathlib.Path
        The JASPAR cache directory.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> jaspar_data_dir()
    PosixPath('/scratch/liulab/motif/jaspar')
    >>> del os.environ["LIULAB_DATA"]
    """
    return motif_data_dir() / JASPAR_SUBDIR


def check_release(release: str) -> str:
    """Return ``release`` if this package prepares it, else say which it does.

    Examples
    --------
    >>> check_release("2024")
    '2024'
    """
    if release not in JASPAR_RELEASES:
        raise ValueError(
            f"no JASPAR release {release!r}: this package prepares "
            f"{', '.join(JASPAR_RELEASES)}. Pass one of those as a string."
        )
    return release


def check_tax_group(tax_group: str) -> str:
    """Return ``tax_group`` if JASPAR publishes it, else say which it does.

    Examples
    --------
    >>> check_tax_group("all")
    'all'
    """
    if tax_group not in JASPAR_TAX_GROUPS:
        raise ValueError(
            f"no JASPAR tax group {tax_group!r}: the groups are "
            f"{', '.join(JASPAR_TAX_GROUPS)}. 'all' is the union of the other seven, and "
            f"it selects a file of its own rather than filtering one after it is read."
        )
    return tax_group


def parse_transfac(text: str) -> tuple[Motif, ...]:
    r"""Read transfac text into **Motif**s, in the order the file spells them.

    A pure function from text to motifs: it opens nothing, downloads nothing, and holds no
    opinion about which **Release** the text came from. Records are separated by ``//``,
    the **Count matrix** is the tab-separated block under the ``PO`` header, and the six
    annotations are the ``CC key:value`` lines.

    **Annotation values are separated by a semicolon and never by a comma.** Commas occur
    inside single values — ``"C3H(C),C2HC zinc-fingers like factors"`` is one class and
    ``"PBM, CSA and/or DIP-chip"`` is one data type — so splitting on one would corrupt
    roughly fifty records per release silently. An empty value means the source stated
    nothing and becomes an empty tuple or an empty string, which is common and not an
    error.

    Parameters
    ----------
    text : str
        The whole transfac file. These are under 1 MB, so they are read whole.

    Returns
    -------
    tuple of Motif
        One motif per record, in file order. Empty for empty text.

    Raises
    ------
    TransfacError
        If a record has no accession, no count matrix, a ``PO`` header that is not
        ``A C G T``, or a count row that is not four numbers.

    Examples
    --------
    >>> record = '''AC MA0260.1
    ... XX
    ... ID che-1
    ... XX
    ... PO\tA\tC\tG\tT
    ... 01\t0.0\t0.0\t37.0\t0.0
    ... 02\t37.0\t0.0\t0.0\t0.0
    ... XX
    ... CC tax_group:nematodes
    ... CC tf_class:C2H2 zinc finger factors; Homeo domain factors
    ... XX
    ... //
    ... '''
    >>> (motif,) = parse_transfac(record)
    >>> motif
    Motif(motif_id='MA0260.1', motif_name='che-1', length=2, offset=0)
    >>> motif.tf_class
    ('C2H2 zinc finger factors', 'Homeo domain factors')
    >>> motif.consensus
    DNA('GA')
    """
    motifs: list[Motif] = []
    record: _Record = _Record(position=1)
    for line in text.splitlines():
        if line.startswith(_END):
            if not record.is_empty:
                motifs.append(record.build())
                record = _Record(position=len(motifs) + 1)
        elif line.startswith(_ACCESSION):
            record.accession = line[len(_ACCESSION) :].strip()
        elif line.startswith(_NAME):
            record.name = line[len(_NAME) :].strip()
        elif line.startswith(_ANNOTATION):
            record.annotate(line[len(_ANNOTATION) :])
        elif line.startswith(_HEADER):
            record.start_matrix(line)
        elif record.reading_matrix:
            record.add_row(line)
    if not record.is_empty:
        # A final record whose `//` is missing is still a record, and losing it silently
        # is exactly the truncation the count check exists to catch — so keep it and let
        # the count speak.
        motifs.append(record.build())
    return tuple(motifs)


class _Record:
    """One transfac record being read, and everything a line can add to it."""

    def __init__(self, position: int) -> None:
        self.position = position
        self.accession: str | None = None
        self.name = ""
        self.annotations: dict[str, str] = {}
        self.rows: list[list[float]] = []
        self.reading_matrix = False

    @property
    def is_empty(self) -> bool:
        """Whether nothing has been read into this record yet."""
        return self.accession is None and not self.rows and not self.annotations

    @property
    def _where(self) -> str:
        """How an error names this record: its accession, or where it sits in the file."""
        return repr(self.accession) if self.accession else f"record {self.position}"

    def annotate(self, text: str) -> None:
        """Read one ``CC key:value`` line, keeping the value exactly as it was written."""
        key, _, value = text.partition(":")
        self.annotations[key.strip()] = value.strip()

    def start_matrix(self, header: str) -> None:
        """Read the ``PO`` header, proving its columns are the bases in :data:`BASES` order."""
        columns = tuple(header.split()[1:])
        if columns != BASES:
            raise TransfacError(
                f"{self._where}: the count matrix header names the columns {columns}, and "
                f"this parser reads them as {BASES}. Row order is what every matrix in the "
                f"package is indexed by, so a file in another order is refused rather than "
                f"read transposed."
            )
        self.reading_matrix = True

    def add_row(self, line: str) -> None:
        """Read one count row, or end the matrix at the first line that is not one."""
        fields = line.split("\t")
        if len(fields) != len(BASES) + 1:
            self.reading_matrix = False
            return
        try:
            self.rows.append([float(field) for field in fields[1:]])
        except ValueError as error:
            raise TransfacError(
                f"{self._where}: the count row {line!r} does not hold four numbers. A "
                f"count row is a position label and one count per base, tab separated."
            ) from error

    def build(self) -> Motif:
        """Turn what was read into a :class:`~genome.tf.motif.motif.Motif`."""
        if not self.accession:
            raise TransfacError(
                f"{self._where} has no 'AC' line, so the motif it describes cannot be "
                f"addressed. Every transfac record opens with its accession."
            )
        if not self.rows:
            raise TransfacError(
                f"{self._where} has no count matrix: a record carries a 'PO' header and "
                f"one tab-separated count row per position."
            )
        # The four plural annotations are split and the two singular ones are not, which
        # is the whole of the separator rule as it reaches the Motif type.
        return Motif(
            self.accession,
            self.name,
            np.array(self.rows, dtype=np.float64).T,
            tax_group=self.annotations.get("tax_group", ""),
            tf_class=self._values("tf_class"),
            tf_family=self._values("tf_family"),
            uniprot_ids=self._values("uniprot_ids"),
            pubmed_ids=self._values("pubmed_ids"),
            data_type=self.annotations.get("data_type", ""),
        )

    def _values(self, field: str) -> tuple[str, ...]:
        """Read one plural annotation, split on its semicolons — empty when it says nothing."""
        return _split_values(self.annotations.get(field, ""))


def _split_values(value: str) -> tuple[str, ...]:
    """Split one plural annotation on its semicolons, dropping what says nothing."""
    return tuple(part.strip() for part in value.split(_VALUE_SEPARATOR) if part.strip())


class JasparDatabase(MotifSet):
    """One JASPAR **Release** and **Tax group**, prepared on disk and read into a set.

    A :class:`~genome.tf.motif.motif.MotifSet` that also knows *which* release it is —
    :attr:`release`, :attr:`tax_group`, :attr:`source_url` and the :attr:`path` its bytes
    are cached at — so a **Hit table** produced from it can say what it was scanned with
    months later. Everything a motif set does, it does; and :meth:`~MotifSet.filter` hands
    back a plain motif set rather than another database, because a filtered release is no
    longer that release.

    **Constructing one prepares it**, as opening a :class:`~genome.genome.Genome` does:
    the file is fetched into :func:`jaspar_data_dir` on the first construction and re-read
    on every one after, with nothing fetched. The lab's CPU cluster compute nodes have no
    internet, so the first construction of a release must happen on a login node.

    Parameters
    ----------
    release : str, default ``"2026"``
        One of :data:`JASPAR_RELEASES`.
    tax_group : str, default ``"vertebrates"``
        One of :data:`JASPAR_TAX_GROUPS`. It chooses which file is downloaded rather than
        filtering one afterwards, so a worm scan never pays for a thousand plant matrices.
    cache_dir : str or pathlib.Path, optional
        The directory to cache in, overriding :func:`jaspar_data_dir`. The directory
        itself, not a root to file under.
    progressbar : bool, default True
        Show the download's progress bar. Nothing is drawn when the file is already there.

    Attributes
    ----------
    release : str
        The **Release** this is.
    tax_group : str
        The **Tax group** this is.
    path : pathlib.Path
        The cached file these motifs were read from.
    source_url : str
        Where those bytes came from.

    Raises
    ------
    ValueError
        If the release or tax group is not one this package prepares.
    TransfacError
        If a record in the file cannot be read.
    JasparReleaseError
        If the file holds the wrong number of motifs, or two versions of one matrix.

    Examples
    --------
    >>> from genome.tf.motif import JasparDatabase
    >>> worms = JasparDatabase("2024", "nematodes")               # doctest: +SKIP
    >>> len(worms)                                                # doctest: +SKIP
    103
    >>> worms["MA0260"].motif_name                                # doctest: +SKIP
    'che-1'
    >>> worms.filter(tf_class="zinc finger")                      # doctest: +SKIP
    MotifSet(motifs=27)
    """

    def __init__(
        self,
        release: str = DEFAULT_RELEASE,
        tax_group: str = DEFAULT_TAX_GROUP,
        *,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> None:
        self.release = check_release(release)
        self.tax_group = check_tax_group(tax_group)
        self.source_url = jaspar_url(self.release, self.tax_group)
        directory = Path(cache_dir).expanduser() if cache_dir is not None else jaspar_data_dir()
        self.path = _prepare(
            self.release, self.tax_group, directory=directory, progressbar=progressbar
        )
        motifs = parse_transfac(self.path.read_text(encoding="utf-8"))
        _check_count(motifs, release=self.release, tax_group=self.tax_group, path=self.path)
        _check_base_ids(motifs, path=self.path)
        super().__init__(motifs)

    def __repr__(self) -> str:
        """Return which release this is and how many motifs it holds."""
        return (
            f"JasparDatabase(release={self.release!r}, tax_group={self.tax_group!r}, "
            f"motifs={len(self)})"
        )


def _prepare(release: str, tax_group: str, *, directory: Path, progressbar: bool) -> Path:
    """Return the cached transfac file, downloading it once if it is not there yet.

    The download lands in the directory's working area under a name of its own and is
    renamed into place only once the fetch returned — the rename is atomic within one
    filesystem, and the working area is inside the cache directory so that it is one. So
    the final name is occupied by a whole file or by nothing, which is what stands here in
    place of a **Completion marker**.

    The working area is emptied *before* the fetch as well as after. Every other build in
    this package keeps what an interrupted run downloaded so the repair does not fetch it
    again; that trade is worth nothing under a megabyte, and taking it would mean the fetch
    step could adopt an abandoned part-file as a finished download — the one way bytes from
    an interrupted run could still reach the final name.
    """
    path = directory / jaspar_filename(release, tax_group)
    if path.is_file():
        return path
    clear_work_dir(directory)
    fetched = fetch.fetch_url(
        jaspar_url(release, tax_group),
        work_dir(directory),
        fname=f"{path.name}.part",
        progressbar=progressbar,
    )
    directory.mkdir(parents=True, exist_ok=True)
    fetched.replace(path)
    clear_work_dir(directory)
    return path


def _check_count(motifs: Sequence[Motif], *, release: str, tax_group: str, path: Path) -> None:
    """Hold a parsed release to the motif count :data:`MOTIF_COUNTS` says it has."""
    expected = MOTIF_COUNTS[release, tax_group]
    if len(motifs) != expected:
        raise JasparReleaseError(
            f"{path} holds {len(motifs)} motifs where JASPAR {release} {tax_group} has "
            f"{expected}. A short file is a truncated download and not a smaller release, "
            f"and scanning with it would silently miss motifs — delete {path} and "
            f"construct the database again to fetch it afresh."
        )


def _check_base_ids(motifs: Sequence[Motif], *, path: Path) -> None:
    """Prove the release is non-redundant: one version of each matrix, so base ids resolve.

    Asserted rather than assumed, because it is what lets ``db['MA0139']`` mean the one
    version this release ships without the caller having to remember which.
    """
    seen: dict[str, str] = {}
    for motif in motifs:
        base = base_id(motif)
        first = seen.get(base)
        if first is not None:
            raise JasparReleaseError(
                f"{path} ships two versions of the matrix {base}: {first} and "
                f"{motif.motif_id}. A non-redundant release ships exactly one of each, "
                f"which is what makes a bare base id address one motif — so this file is "
                f"not the non-redundant release it is cached as. Delete it and construct "
                f"the database again."
            )
        seen[base] = motif.motif_id
