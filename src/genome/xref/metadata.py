"""The curated table of **Xref set**s — which exist, and where each one's bytes come from.

One small hand-maintained TSV ships inside the package,
``data/xref/xref_metadata.tsv``, one row per ``(species, Xref source, Release)``. A row
says who published that release, what the publisher calls it, where its file is fetched
from and what the file's unpacked bytes hash to, so **naming a release is enough to fetch
it** — the contract the annotation metadata table already offers one level up. Adding an
**Xref source** is a row here plus a reader, not a branch in code.

Nothing is shipped but the row. The files themselves are 26 MB and up and are downloaded
into the **Data dir** on first construction (ADR-0018), which is also why a publisher is
eligible only if its old releases stay retrievable at stable URLs: a checksum that travels
in the wheel must still match a year later, and a publisher that overwrites its file in
place — or keeps no archive of the release the row names — breaks that.

``source_checksum`` is over the publisher's **unpacked** bytes (ADR-0006), and **which
bytes a publisher's own checksum covers is not a convention any two of them share**:
Alliance's ``md5Sum`` is the digest of the TSV *inside* the gzip, so hashing the
``.tsv.gz`` as it arrives mismatches every time, while Ensembl's ``CHECKSUMS`` is a BSD
``sum`` over the served ``.gz`` — 16 bits, and no integrity check for a 6 MB file. Check
which before pinning a new row; the attribution beside this table records each publisher's
own value and what it covers. The set stored on disk carries a second digest of its own,
since what is stored is a per-species slice rather than these bytes.

``default`` marks the **Default xref source** for a species — the one a caller who names
none is answered by, so everyone in the lab reaches for the same one without discussing
it. It is a default and not a recommendation: naming a source is how the scientific choice
gets made deliberately, and NCBI and Ensembl agree on only 57.6% of human gene-level
(GeneID, ENSG) pairs, so the choice determines nearly half the answer.

``symbol_default`` marks the same thing **for the other question this package answers**
(ADR-0021). A default is per species *and* per question, because the source that carries a
species' ids is often not the one that carries its symbols: human's ids default to the
Alliance, whose cross-reference file publishes no human symbol at all, while HGNC's
quarterly archive does. Two flags rather than one, so ``lookup_xref(species)`` and
``lookup_xref(species, for_symbols=True)`` each read a row somebody wrote down.

**Rows for one ``(species, source)`` are listed oldest release first**, and the last is
what a caller who names no release gets. Ordering is the table's and never parsed out of
the release string: ``"10.0.0"`` sorts before ``"9.0.0"`` as text and after it as a
version, and no reader here is going to be the place that guesses which a publisher meant.

Examples
--------
>>> from genome.xref.metadata import lookup_xref, xref_species, xref_table
>>> "Homo sapiens" in xref_species()
True
>>> row = lookup_xref("Homo sapiens")
>>> row.source, row.publisher
('alliance', 'Alliance of Genome Resources')
>>> row.source_checksum.startswith("md5:")
True
>>> sorted({record.source for record in xref_table()})
['alliance', 'alliance_bgi', 'ensembl', 'hgnc']
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any, get_type_hints

import pandas as pd

from genome.metadata import parse_cell, species_slug

#: Directory inside the package holding the curated **Xref source** table and the
#: attribution beside it. No **Xref set** ships here — only the row that fetches one.
XREF_SUBDIR = "data/xref"

#: The curated table itself. **Plain**, not gzipped: bulk data is compressed and small
#: metadata is not, the convention every shipped-data directory here follows.
XREF_METADATA_RESOURCE = f"{XREF_SUBDIR}/xref_metadata.tsv"


class NoXrefSetError(LookupError):
    """Nothing this package prepares answers for that species, source or release.

    Raised by :func:`lookup_xref` and so by constructing an
    :class:`~genome.xref.xref.XrefSet`. A :class:`LookupError`, because it is a name that
    resolves to nothing rather than a malformed one, and it is the same shape as the
    census and cofactor misses: the message names what *is* available, so a caller who
    guessed a species reads the three that have a set instead of guessing again.

    A species with no Ensembl presence is answered by this permanently rather than
    pending — the registered *E. coli* HT115 assembly has no hub to hang a **Namespace**
    off, so it has no **Xref set** and says so instead of being served a fudged one.

    Examples
    --------
    >>> try:
    ...     lookup_xref("Danio rerio")
    ... except NoXrefSetError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """


@dataclass(frozen=True)
class XrefMetadata:
    """One **Xref set** this package prepares (one row of the curated table).

    The single declaration of what such a row consists of: the table is parsed through
    these fields, in this order, and every one is required except ``pubmed_id`` — a
    publisher with no paper is cited by name and URL, and a set nobody can cite is one
    this package will not fetch.

    Attributes
    ----------
    species : str
        The species, as the assembly metadata table spells it — ``"Homo sapiens"``. Its
        slug names the set's directory, and either spelling is accepted on the way in.
    ncbi_taxid : int
        NCBI taxonomy id, which is how the species' rows are picked out of a publisher's
        multi-species file. Read from the row and never inferred from the species name.
    source : str
        The **Xref source**, lower-cased and stable — ``"alliance"``. It names both the
        reader and the directory the set is filed under.
    release : str
        The pinned **Release**, and the string a caller names this set by.
    publisher : str
        Who published it, and who is to be cited for it.
    version : str
        The publisher's own release identifier, which is often the same string as
        ``release`` and is not required to be: Alliance's file states ``9.0.0`` inside its
        own header, where a quarterly archive states a date.
    pubmed_id : int or None
        PubMed id of the paper to cite, or ``None`` for a publisher with none.
    url : str
        Where the publisher's file is fetched from.
    source_checksum : str
        The publisher's own checksum of that file, as ``"<algorithm>:<hexdigest>"``, taken
        over the **unpacked** bytes (ADR-0006) — which is what Alliance publishes, and
        what this package would have computed anyway. It is provenance rather than the
        integrity check: what is stored on disk is a per-species slice and not the
        publisher's bytes, so the slice carries a digest of its own.
    default : bool
        Whether this is the species' **Default xref source** when the question is
        identifiers.
    symbol_default : bool
        Whether this is it when the question is symbols (ADR-0021). A separate flag because
        the two answers are usually different rows: the source carrying a species' ids
        often publishes none of its symbols.

    Examples
    --------
    >>> record = lookup_xref("Caenorhabditis elegans")
    >>> record.ncbi_taxid, record.default, record.symbol_default
    (6239, True, False)
    >>> record.url.endswith(".tsv.gz")
    True
    >>> lookup_xref("Caenorhabditis elegans", for_symbols=True).symbol_default
    True
    """

    species: str
    ncbi_taxid: int
    source: str
    release: str
    publisher: str
    version: str
    pubmed_id: int | None
    url: str
    source_checksum: str
    default: bool = False
    symbol_default: bool = False

    @classmethod
    def from_row(cls, row: Mapping[str, object]) -> XrefMetadata:
        """Build a record from one row of the curated table.

        Parameters
        ----------
        row : mapping of str to object
            Column name to cell, as the shipped TSV spells one. Keys outside
            :data:`XREF_FIELDS` are ignored.

        Returns
        -------
        XrefMetadata
            The record the row spells.

        Raises
        ------
        genome.metadata.MetadataRowError
            If a cell cannot be read as its column's type, or a column that has no unknown
            is blank. The message names the column.

        Examples
        --------
        >>> XrefMetadata.from_row(
        ...     {
        ...         "species": "Tiny beast",
        ...         "ncbi_taxid": "1",
        ...         "source": "somewhere",
        ...         "release": "1.0",
        ...         "publisher": "Someone et al. 1999",
        ...         "version": "1.0",
        ...         "pubmed_id": "",
        ...         "url": "https://example.org/beast.tsv.gz",
        ...         "source_checksum": "md5:" + "0" * 32,
        ...         "default": "yes",
        ...         "symbol_default": "no",
        ...     }
        ... ).pubmed_id is None
        True
        """
        return cls(**{name: parse_cell(name, row, _XREF_FIELD_TYPES) for name in XREF_FIELDS})

    def attribution(self) -> str:
        """Return the one line to print beside anything this set answered.

        What a caller owes the publisher, rendered once here so a notebook, the CLI and an
        error message all say it the same way.

        Returns
        -------
        str
            Publisher, release, PubMed id where there is one, and the source URL.

        Examples
        --------
        >>> print(lookup_xref("Mus musculus").attribution())
        Alliance of Genome Resources 9.0.0 (PMID 38552170) — https://download.alliancegenome.org/9.0.0/GENECROSSREFERENCE/COMBINED/GENECROSSREFERENCE_COMBINED_11.tsv.gz
        """
        cited = f" (PMID {self.pubmed_id})" if self.pubmed_id is not None else ""
        return f"{self.publisher} {self.version}{cited} \N{EM DASH} {self.url}"


#: Each field's declared type, which parses that field's column of the table.
_XREF_FIELD_TYPES: dict[str, Any] = get_type_hints(XrefMetadata)

#: The field names, in table-column order — the columns every row carries.
XREF_FIELDS: tuple[str, ...] = tuple(_XREF_FIELD_TYPES)


@cache
def xref_table() -> tuple[XrefMetadata, ...]:
    """Return every **Xref set** the shipped table lists, in table order.

    Read once and cached; the records are frozen, so the tuple is safe to hold on to.
    Rows for one ``(species, source)`` are in release order, oldest first.

    Returns
    -------
    tuple of XrefMetadata
        One record per row of ``data/xref/xref_metadata.tsv``.

    Examples
    --------
    >>> sorted({record.release for record in xref_table()})
    ['116', '2026-07-07', '9.0.0']
    """
    resource = files("genome").joinpath(XREF_METADATA_RESOURCE)
    with resource.open("r", encoding="utf-8") as handle:
        frame = pd.read_csv(handle, sep="\t", dtype=str)
    return tuple(XrefMetadata.from_row(dict(row)) for _, row in frame.iterrows())


def xref_species(*, table: Sequence[XrefMetadata] | None = None) -> tuple[str, ...]:
    """Return every species an **Xref set** exists for, as the table spells them.

    What can be asked about at all, and what :class:`NoXrefSetError` names when a species
    cannot be.

    Parameters
    ----------
    table : sequence of XrefMetadata, optional
        The rows to read; the shipped table when omitted.

    Returns
    -------
    tuple of str
        The species names, in first-listed order.

    Examples
    --------
    >>> xref_species()
    ('Homo sapiens', 'Mus musculus', 'Caenorhabditis elegans')
    """
    return tuple(dict.fromkeys(record.species for record in _rows(table)))


def xref_sources(species: str, *, table: Sequence[XrefMetadata] | None = None) -> tuple[str, ...]:
    """Return every **Xref source** that answers for ``species``, in first-listed order.

    Parameters
    ----------
    species : str
        The species, in either the table's spelling or its slug.
    table : sequence of XrefMetadata, optional
        The rows to read; the shipped table when omitted.

    Returns
    -------
    tuple of str
        The source names. Empty for a species no set exists for.

    Examples
    --------
    >>> xref_sources("Mus musculus")
    ('alliance', 'ensembl', 'alliance_bgi')
    >>> xref_sources("Caenorhabditis elegans")
    ('alliance', 'alliance_bgi')
    >>> xref_sources("Danio rerio")
    ()
    """
    slug = species_slug(species)
    return tuple(
        dict.fromkeys(
            record.source for record in _rows(table) if species_slug(record.species) == slug
        )
    )


def xref_releases(
    species: str, source: str, *, table: Sequence[XrefMetadata] | None = None
) -> tuple[str, ...]:
    """Return every **Release** of ``source`` that answers for ``species``, oldest first.

    Parameters
    ----------
    species : str
        The species, in either the table's spelling or its slug.
    source : str
        The **Xref source**.
    table : sequence of XrefMetadata, optional
        The rows to read; the shipped table when omitted.

    Returns
    -------
    tuple of str
        The release strings in table order, which is oldest first — so the last is what a
        caller who names no release gets.

    Examples
    --------
    >>> xref_releases("Homo sapiens", "alliance")
    ('9.0.0',)
    """
    slug = species_slug(species)
    return tuple(
        record.release
        for record in _rows(table)
        if species_slug(record.species) == slug and record.source == source
    )


def lookup_xref(
    species: str,
    source: str | None = None,
    release: str | None = None,
    *,
    for_symbols: bool = False,
    table: Sequence[XrefMetadata] | None = None,
) -> XrefMetadata:
    """Return the curated row for one **Xref set**, filling in the defaults.

    The one lookup, and the one place a miss is turned into an error that says what *is*
    there. It is total in the other direction: it always answers with a row or raises, so
    nothing downstream guards a ``None`` before reading a field.

    **A default is per species and per question** (ADR-0021), which is what ``for_symbols``
    selects. It changes what an *unnamed* source resolves to and nothing else: a named
    source is the caller's deliberate scientific choice and is never swapped for another,
    and a named release is still honoured against whichever source was filled in.

    Parameters
    ----------
    species : str
        The species, in either the table's own spelling (``"Homo sapiens"``) or its slug
        (``"homo_sapiens"``).
    source : str, optional
        The **Xref source**. Omitted, the species' **Default xref source** answers — the
        one flagged for identifiers, or the one flagged for symbols when ``for_symbols``.
    release : str, optional
        The **Release**. Omitted, the newest the table lists for that source answers —
        which is the last row, the table being in release order. Named, it is honoured
        whether or not a source was named too: a release asked for against the default
        source either answers with that release or raises naming the ones that source has,
        and is never quietly swapped for another. Each source numbers its own releases and
        they do not correspond, so ``"116"`` is not a release the default source has.
    for_symbols : bool, default False
        Fill an unnamed source in with the species' symbol-carrying default rather than its
        identifier one. The question a caller is about to ask, named here because a row is
        resolved before anybody knows it — an **Xref set** is built for a species, a source
        and a release, and which of the two verbs it will be asked is not part of that.
    table : sequence of XrefMetadata, optional
        The rows to read; the shipped table when omitted. A caller curating rows of their
        own hands them over here, and nothing is installed by passing them.

    Returns
    -------
    XrefMetadata
        The row, with ``source`` and ``release`` resolved to what actually answered.

    Raises
    ------
    NoXrefSetError
        If no set exists for that species, that source, or that release of it — or, under
        ``for_symbols``, if no row for the species is flagged to answer symbols. The
        message names the species, sources or releases that do exist, whichever missed.

    Examples
    --------
    >>> lookup_xref("homo_sapiens").release
    '9.0.0'
    >>> lookup_xref("Homo sapiens", "alliance", "9.0.0").ncbi_taxid
    9606
    >>> lookup_xref("Homo sapiens", release="9.0.0").source
    'alliance'
    >>> lookup_xref("Homo sapiens", for_symbols=True).source
    'hgnc'
    >>> lookup_xref("Homo sapiens", "alliance", for_symbols=True).source
    'alliance'
    >>> try:
    ...     lookup_xref("Homo sapiens", "alliance", "1.0")
    ... except NoXrefSetError as error:
    ...     print("9.0.0" in str(error))
    True
    """
    rows = _rows(table)
    slug = species_slug(species)
    for_species = [record for record in rows if species_slug(record.species) == slug]
    if not for_species:
        listed = ", ".join(xref_species(table=rows))
        raise NoXrefSetError(
            f"no xref set for {species!r}: this package prepares one for {listed}. Ask for "
            f"one of those, or convert the ids you hold with the publisher's own file — a "
            f"species with no Ensembl presence has no hub to hang a namespace off and is "
            f"unanswerable here by design (ADR-0017)."
        )
    # The default source still has to honour a named release. Resolving the source first
    # and then falling through to the release check keeps one path: naming a release
    # without a source must answer with *that* release or say which ones exist, never
    # quietly with another — two sources numbering their releases differently is exactly
    # where a silently substituted release stops being reproducible. Which of the two
    # defaults is filled in is the only thing `for_symbols` decides.
    if source is None:
        filled = (
            _symbol_default_row(for_species)
            if for_symbols
            else _default_row(for_species, species=species)
        )
        source = filled.source
    for_source = [record for record in for_species if record.source == source]
    if not for_source:
        listed = ", ".join(xref_sources(species, table=rows))
        raise NoXrefSetError(
            f"no xref source {source!r} for {for_species[0].species!r}: this package "
            f"prepares {listed}. Name one of those, or leave the source out to be answered "
            f"by the default one."
        )
    if release is None:
        # Oldest first, so the last row is the newest release — a new analysis starts on
        # current data and reproducing an old one is the call that says which.
        return for_source[-1]
    for record in for_source:
        if record.release == release:
            return record
    listed = ", ".join(record.release for record in for_source)
    raise NoXrefSetError(
        f"no {source} release {release!r} for {for_species[0].species!r}: this package "
        f"prepares {listed}. Name one of those, or leave the release out to be answered by "
        f"the newest."
    )


def _rows(table: Sequence[XrefMetadata] | None) -> Sequence[XrefMetadata]:
    """Return the rows to read — the caller's, else the shipped table's."""
    return xref_table() if table is None else table


def _default_row(for_species: Sequence[XrefMetadata], *, species: str) -> XrefMetadata:
    """Return the species' **Default xref source**'s newest release, or say none is set.

    A species with exactly one source needs no flag to have a default; a species with
    several and no flag has no default, and guessing one would make which publisher
    answered depend on row order rather than on a decision anybody wrote down.
    """
    flagged = [record for record in for_species if record.default]
    chosen = flagged or ([for_species[0]] if len({r.source for r in for_species}) == 1 else [])
    if not chosen:
        listed = ", ".join(dict.fromkeys(record.source for record in for_species))
        raise NoXrefSetError(
            f"no default xref source for {for_species[0].species!r}, which has more than "
            f"one: {listed}. Name the source you want — which publisher answers is a "
            f"scientific choice and is not made for you here (ADR-0017)."
        )
    source = chosen[-1].source
    return [record for record in for_species if record.source == source][-1]


def _symbol_default_row(for_species: Sequence[XrefMetadata]) -> XrefMetadata:
    """Return the newest release of the source flagged to answer symbols for the species.

    **No fallback to a species' only source**, which is where this parts company with
    :func:`_default_row`: a publisher that carries no symbol does not begin carrying them
    by being the only one prepared, and answering from it would hand back a gene list with
    no matches where the honest answer is that nothing here matches symbols for that
    species yet.
    """
    flagged = [record for record in for_species if record.symbol_default]
    if not flagged:
        listed = ", ".join(dict.fromkeys(record.source for record in for_species))
        raise NoXrefSetError(
            f"no xref source answers symbols for {for_species[0].species!r}: the ones "
            f"prepared for it are {listed}, and none of their rows is flagged to answer a "
            f"symbol question. Name the source you want, or flag one in the curated xref "
            f"metadata table's 'symbol_default' column."
        )
    return [record for record in for_species if record.source == flagged[-1].source][-1]
