"""The shipped table that says where one species pair's homologies are published.

One small hand-maintained TSV ships inside the package,
``data/homology/homology_metadata.tsv``. Keyed by a species pair and a **Release**, it
records the publisher, which per-species file of that release holds the pair, the URL to
fetch and the publisher's own md5 for those bytes — so naming a pair is enough to prepare
it, exactly as **Annotation metadata** makes naming an annotation enough to register one.
No homology data ships here: this table is provenance, and a **Homology set** is
downloaded (ADR-0018's stance on pinnable sources is what makes a shipped checksum still
right a year later).

**The holding species is a measurement, not an assumption.** Compara's per-species files
are a de-duplicated partition at the pair level, and which file a pair lands in is
arbitrary — the human file of release 116 holds 23,982 human↔worm rows and **zero**
human↔mouse rows, which live only in the mouse file. Each row here was written by counting
rows in the published files, and the count is taken again every time a set is prepared: a
pair that comes back empty raises and names the other file rather than answering nothing
(see :class:`~genome.homology.compara.ComparaPartitionError`). Nothing here is trusted to
stay true across releases, which is why a new release is new rows rather than an edit.

**Every row carries a real URL rather than one this package formats.** Compara's naming is
not uniform across releases — release 113 ships these dumps uncompressed, so its file has
no ``.gz`` — and only releases 90 and 116 publish an ``MD5SUM`` at all, with 91 to 112
publishing no checksum of any kind. So the URL is read off the release's own listing and
written down here, and a release that publishes no checksum cannot be added: a row with no
md5 is refused, which is what stops an unpinnable release from being pinned by accident.

The reader is pure — it reads one shipped package resource and nothing else, never the
**Data dir** and never the network. A row that cannot be read raises
:class:`HomologyMetadataError` naming the column, because a table that ships broken is a
defect in this package rather than anything a caller did.

Examples
--------
>>> from genome.homology.metadata import homology_metadata, homology_table
>>> {row.release for row in homology_table()}
{'116'}
>>> row = homology_metadata("Homo sapiens", "Mus musculus", "116")
>>> row.holding_species
'Mus musculus'
>>> print(row.attribution())
Ensembl Compara release 116 (PMID 26896847) — https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/mus_musculus/Compara.116.protein_default.homologies.tsv.gz
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any

#: Directory inside the package holding the homology provenance table. No data file
#: ships beside it: a **Homology set** is downloaded, and this says from where.
HOMOLOGY_SUBDIR = "data/homology"

#: The provenance table keyed by species pair and **Release**. **Plain**, not gzipped —
#: bulk gzipped, small metadata plain, the convention every shipped-data directory here
#: follows.
HOMOLOGY_METADATA_RESOURCE = f"{HOMOLOGY_SUBDIR}/homology_metadata.tsv"

#: The provenance table's columns, in table order.
METADATA_COLUMNS: tuple[str, ...] = (
    "release",
    "species",
    "other_species",
    "holding_species",
    "publisher",
    "pubmed_id",
    "source_url",
    "md5",
)

#: Provenance columns holding a number rather than text.
_NUMERIC_METADATA_COLUMNS = frozenset({"pubmed_id"})


class HomologyMetadataError(ValueError):
    r"""A row of the shipped provenance table cannot be read as a record.

    A defect in this package rather than anything a caller did, so the message names the
    file, the row and the column that refused, and the repair is to fix that cell.

    Examples
    --------
    >>> try:
    ...     read_metadata("release\n", origin="homology_metadata.tsv")
    ... except HomologyMetadataError as error:
    ...     print("homology_metadata.tsv" in str(error))
    True
    """


@dataclass(frozen=True)
class HomologyMetadata:
    """Where one species pair's homologies are published (one row of the table).

    The single declaration of what a **Homology set**'s provenance consists of: the table
    is parsed through these fields, column by column, exactly as
    :class:`~genome.metadata.AssemblyMetadata` parses its own. Every column is required —
    a set nobody can cite is one this package may not point anyone at, and a set with no
    checksum is one a truncated fetch would answer from.

    The pair is unordered: a row is looked up by its two species in either order, since
    which of them a caller asks *about* is a property of the question and not of the file.

    Attributes
    ----------
    release : str
        The Ensembl Compara **Release**, as the publisher numbers it — ``"116"``.
    species : str
        One species of the pair, as the assembly metadata table spells it.
    other_species : str
        The other, likewise.
    holding_species : str
        Whose per-species file actually holds this pair's rows in this release —
        **measured by counting**, never assumed, and re-checked every time a set is
        prepared.
    publisher : str
        Who published it, and who is to be cited.
    pubmed_id : int
        PubMed id of the paper to cite.
    source_url : str
        The published file to fetch: the holding species' own per-species dump.
    md5 : str
        The publisher's own md5 for those bytes, read from the ``MD5SUM`` file beside them
        and checked against the bytes **as they are fetched**. That check is load-bearing
        and not a formality: a resumed download of one of these gzips has been seen to pass
        ``gzip -t`` with the wrong md5, so opening cleanly is no evidence at all. Only
        Compara releases 90 and 116 publish an ``MD5SUM``, which is part of why 116 is what
        is pinned.

    Examples
    --------
    >>> row = homology_metadata("Caenorhabditis elegans", "Homo sapiens", "116")
    >>> row.holding_species, row.md5[:8]
    ('Homo sapiens', '59857f48')
    >>> row.pair
    ('Caenorhabditis elegans', 'Homo sapiens')
    """

    release: str
    species: str
    other_species: str
    holding_species: str
    publisher: str
    pubmed_id: int
    source_url: str
    md5: str

    @property
    def pair(self) -> tuple[str, str]:
        """The two species, sorted — the key a pair is looked up and filed under."""
        first, second = sorted((self.species, self.other_species))
        return first, second

    @property
    def other(self) -> str:
        """The species of the pair whose file does *not* hold it, in this release.

        What a :class:`~genome.homology.compara.ComparaPartitionError` names when the
        recorded file comes back empty: the partition moved, and this is where the pair
        went.
        """
        first, second = self.pair
        return second if first == self.holding_species else first

    @classmethod
    def from_row(cls, row: Mapping[str, str], *, origin: str) -> HomologyMetadata:
        """Read one row of the provenance table into a record.

        Parameters
        ----------
        row : mapping of str to str
            One row, keyed by column name.
        origin : str
            Where the row came from, named in any error.

        Returns
        -------
        HomologyMetadata
            The record.

        Raises
        ------
        HomologyMetadataError
            If a required cell is blank or a numeric one is not a number.

        Examples
        --------
        >>> row = {
        ...     "release": "116",
        ...     "species": "Homo sapiens",
        ...     "other_species": "Mus musculus",
        ...     "holding_species": "Mus musculus",
        ...     "publisher": "Ensembl Compara",
        ...     "pubmed_id": "26896847",
        ...     "source_url": "https://example.invalid/x.tsv.gz",
        ...     "md5": "0" * 32,
        ... }
        >>> HomologyMetadata.from_row(row, origin="test").pubmed_id
        26896847
        """
        return cls(**{name: _parse_cell(name, row, origin=origin) for name in METADATA_COLUMNS})

    def attribution(self) -> str:
        """Return the one line to print beside anything this set answered.

        What a caller owes the publisher, rendered once here so the CLI, a notebook and an
        error message all say it the same way.

        Returns
        -------
        str
            Publisher, release, PubMed id and source URL.

        Examples
        --------
        >>> print(homology_metadata("Homo sapiens", "Mus musculus", "116").attribution())
        Ensembl Compara release 116 (PMID 26896847) — https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/mus_musculus/Compara.116.protein_default.homologies.tsv.gz
        """
        return (
            f"{self.publisher} release {self.release} (PMID {self.pubmed_id}) "
            f"\N{EM DASH} {self.source_url}"
        )


@cache
def homology_table() -> tuple[HomologyMetadata, ...]:
    """Return every row of the shipped provenance table, in table order.

    Which species pairs this package can prepare a **Homology set** for, from which
    release and from whose file. Read once and cached; the records are frozen, so the
    tuple is safe to hold on to.

    Returns
    -------
    tuple of HomologyMetadata
        One record per row.

    Raises
    ------
    HomologyMetadataError
        If a row cannot be read; the message names the column.

    Examples
    --------
    >>> {row.publisher for row in homology_table()}
    {'Ensembl Compara'}
    >>> len(homology_table())
    3
    """
    resource = files("genome").joinpath(HOMOLOGY_METADATA_RESOURCE)
    return read_metadata(resource.read_text(encoding="utf-8"), origin=str(resource))


@cache
def homology_releases() -> tuple[str, ...]:
    """Return every **Release** the shipped table pins, ascending.

    What can be asked for at all, and what an error names when a release cannot be.

    Returns
    -------
    tuple of str
        The release identifiers, sorted.

    Examples
    --------
    >>> homology_releases()
    ('116',)
    """
    return tuple(sorted({row.release for row in homology_table()}))


@cache
def homology_species() -> tuple[str, ...]:
    """Return every species the shipped table names, sorted, as the metadata spells them.

    The species a **Homology set** can be asked about, read off the table rather than
    kept as a second list in code — so adding a species is adding rows.

    Returns
    -------
    tuple of str
        The species names, sorted.

    Examples
    --------
    >>> homology_species()
    ('Caenorhabditis elegans', 'Homo sapiens', 'Mus musculus')
    """
    named: set[str] = set()
    for row in homology_table():
        named.update(row.pair)
    return tuple(sorted(named))


def homology_metadata(species: str, other_species: str, release: str) -> HomologyMetadata | None:
    """Return the table's row for one species pair and **Release**, or ``None``.

    The pair is unordered — ``("Homo sapiens", "Mus musculus")`` and its reverse find the
    same row — because which species a caller asks *about* is a property of the question
    and not of the published file.

    ``None`` is the raw absence, and this is the one place it is how absence is said:
    everything above turns it into an error naming the pairs the table does carry, so
    *nobody pinned this pair* can never be read as *these species share no homologs*.

    Parameters
    ----------
    species : str
        One species of the pair, as the assembly metadata table spells it.
    other_species : str
        The other.
    release : str
        The Compara **Release**.

    Returns
    -------
    HomologyMetadata or None
        The row, or ``None`` when the table pins no such pair in that release.

    Examples
    --------
    >>> homology_metadata("Mus musculus", "Homo sapiens", "116").holding_species
    'Mus musculus'
    >>> homology_metadata("Homo sapiens", "Danio rerio", "116") is None
    True
    """
    wanted = tuple(sorted((species, other_species)))
    for row in homology_table():
        if row.release == release and row.pair == wanted:
            return row
    return None


def read_metadata(text: str, *, origin: str) -> tuple[HomologyMetadata, ...]:
    r"""Read the provenance table from ``text``, holding it to the columns it declares.

    Separate from the resource it came out of, so every way the table can be wrong is
    reachable without writing a broken one into the package.

    Parameters
    ----------
    text : str
        The whole table, tab separated, header first.
    origin : str
        Where the text came from, named in every message.

    Returns
    -------
    tuple of HomologyMetadata
        One record per non-empty row, in table order.

    Raises
    ------
    HomologyMetadataError
        If the header is not :data:`METADATA_COLUMNS`, a row holds the wrong number of
        cells, or a cell cannot be read.

    Examples
    --------
    >>> header = "\t".join(METADATA_COLUMNS)
    >>> read_metadata(header + "\n", origin="test")
    ()
    """
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    header = tuple(lines[0].split("\t")) if lines else ()
    if header != METADATA_COLUMNS:
        raise HomologyMetadataError(
            f"{origin} carries the columns {list(header)} where the provenance table's are "
            f"{list(METADATA_COLUMNS)}. Fix the header, keeping the columns in that order."
        )
    records: list[HomologyMetadata] = []
    for number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        cells = line.split("\t")
        if len(cells) != len(header):
            raise HomologyMetadataError(
                f"{origin} line {number} holds {len(cells)} cells where the header declares "
                f"{len(header)}. The table is a plain TSV with no quoting — fix that line."
            )
        records.append(
            HomologyMetadata.from_row(dict(zip(header, cells, strict=True)), origin=origin)
        )
    return tuple(records)


def _parse_cell(name: str, row: Mapping[str, str], *, origin: str) -> Any:
    """Return one provenance cell, parsed by its column and never blank."""
    text = row.get(name, "").strip()
    if not text:
        raise HomologyMetadataError(
            f"{origin} leaves the {name!r} column blank for the pair "
            f"{row.get('species')!r}/{row.get('other_species')!r}. Every provenance column is "
            f"required: a set nobody can cite is one this package may not point anyone at, and "
            f"one with no checksum is one a truncated fetch would answer from. Fill that cell in."
        )
    if name not in _NUMERIC_METADATA_COLUMNS:
        return text
    try:
        return int(text)
    except ValueError as error:
        raise HomologyMetadataError(
            f"{origin} holds {text!r} in the {name!r} column, which is not a number. Fix that cell."
        ) from error
