"""The shipped cofactor tables — which genes a publisher lists as transcription cofactors.

**Attribution.** Three publishers ship here. AnimalTFDB 4.0, Shen *et al.*,
*Nucleic Acids Research* 51(D1):D39-D45, 2023 (PMID 36268869), from
https://guolab.wchscu.cn/, lists cofactors for every species that ships; EpiFactors
v2.0, Marakulina *et al.*, *Nucleic Acids Research* 51(D1):D564-D570, 2023 (PMID
36350659), from https://epifactors.autosome.org/, lists human's beside it; and a
pinned dated HGNC monthly archive (PMID 41287213), from https://www.genenames.org/,
supplies the **Gene id stem** of every gene EpiFactors names. Every classification
is the publisher's who reached it. **Human membership is this package's own** — the
union of the two lists, and so nobody's verdict but ours (ADR-0016) — where mouse
and worm relay one publisher unchanged. Cite the publishers whose table you used;
:meth:`CofactorProvenance.attribution` renders the line to print, and the two
provenance tables beside the data carry the same facts for every table that ships.

A **Transcription cofactor** is a gene that acts on transcription without binding
DNA sequence-specifically — a chromatin remodeller, a histone-modifying enzyme, a
Mediator subunit. It is a different question from the one :mod:`genome.tf.gene`
answers, which is why this is a peer of that subpackage rather than a part of it:
that half says whether a gene is a **TF gene** and of what **DBD family**, this one
says whether it is a cofactor and of what class. It is keyed the same way, by
**Gene id stem**, so a **TF cofactor list** resolves against an **Annotation**
exactly as a TF gene list does.

One **Cofactor table** per species ships inside the package under
``data/tf_cofactor/<species>.cofactor_table.tsv.gz``, and they are found by
enumerating that directory rather than by any list of species in code — so adding a
species is dropping in a file. The name is the species **slug**, the assembly
metadata table's own spelling lower-cased with each run of non-alphanumerics turned
into one underscore, exactly as a census's file is named, and
:func:`cofactor_table` accepts either spelling.

Four columns are uniform across every table and lead every file — the **Gene id
stem**, the symbol, the cofactor flag and the source (:data:`UNIFORM_COLUMNS`).
Everything after them is one publisher's own column under a namespaced snake_case
name — ``animaltfdb_family``, ``animaltfdb_category``, ``epifactors_function`` — so
a table built from two publishers is more columns and one more provenance row,
never a change to the format, and two publishers' vocabularies are never compared
(ADR-0014). A blank cell is that publisher recording nothing and reads back as
``None``, and a cell recording more than one value spells them apart with ``;``,
the separator every multi-valued cell in this package uses.

:data:`SOURCES` is a closed vocabulary, validated as the file is read: it says which
publisher listed the gene, and it asserts agreement on **membership only**, never on
classification. :data:`CITED_SOURCES` is the vocabulary of the *provenance* table's
own ``source`` column and is a different list on purpose — ``both`` is a fact about a
row of a table and describes no publisher, while ``hgnc`` describes a publisher that
contributes identifiers and no membership, so it is cited in its own right and names
no gene's source. ``is_cofactor`` reads ``yes`` on every row that ships today,
because no publisher here releases a rejected set — and it is kept anyway. Dropping
it would make presence in the file the verdict, at which point a future source could
not record a rejection without a format change.

**Worm ships although no publisher has released a worm TF census.** A worm assembly
answers here while the TF gene half raises for it. That asymmetry is the publishers'
shape rather than a defect: AnimalTFDB assessed *C. elegans* cofactors and nobody has
published its transcription factors.

:func:`cofactor_table` answers ``None`` for a species no table ships for — the raw
absence, and the one place ``None`` is how it is said, because this is the layer
below the one a caller touches. Everything above turns it into an error naming the
species that do have a table.

This module is pure. It reads shipped package resources and nothing else — never the
**Data dir**, never the network — and the files themselves are built by
``scripts/build_tf_cofactor.py``, which lives outside the wheel. A shipped file that
cannot be trusted never answers: it is validated as it is read and raises
:class:`CofactorTableError` naming the file and the repair, since a table that ships
broken is a defect in this package rather than anything a caller did.

Examples
--------
>>> from genome.tf.cofactor import cofactor_species, cofactor_table
>>> "mus_musculus" in cofactor_species()
True
>>> mouse = cofactor_table("Mus musculus")
>>> len(mouse), mouse.columns[:4]
(970, ('gene_id_stem', 'symbol', 'is_cofactor', 'source'))
>>> mouse.provenance.sources[0].publisher
'AnimalTFDB'
>>> len(cofactor_table("Caenorhabditis elegans"))
317
>>> human = cofactor_table("Homo sapiens")
>>> len(human), len(human.provenance.sources)
(1466, 3)
>>> cofactor_table("Danio rerio") is None
True
"""

from __future__ import annotations

import gzip
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any

import pandas as pd

# Imported rather than written a third time. ``species_slug`` is the file-naming
# convention every shipped-data directory here uses and belongs to none of them in
# particular; promoting it to a shared home is a one-line refactor whenever a third
# caller wants it. The flag spellings come along for the same reason: one spelling of
# *yes* across every table this package ships.
from genome.tf.gene.census import FALSE_CELL, TRUE_CELL, species_slug

#: Directory inside the package holding one **Cofactor table** per species, plus the
#: two provenance tables beside them.
COFACTOR_SUBDIR = "data/tf_cofactor"

#: What one of those files is called: the species slug, then this. The slug is what a
#: table is looked up by, so the two halves are never spelled apart.
COFACTOR_SUFFIX = ".cofactor_table.tsv.gz"

#: The provenance table keyed by species: which file, whose taxid, and the checksum of
#: the unpacked bytes. **Plain**, where the tables themselves are gzipped — bulk gzipped,
#: small metadata plain, the convention every shipped-data directory here follows.
COFACTOR_METADATA_RESOURCE = f"{COFACTOR_SUBDIR}/cofactor_metadata.tsv"

#: The provenance table keyed by species **and** source: who published that part of the
#: table, which release, and what to cite for it. Two tables and not one because a table
#: built from three publishers cannot be described by one row, and joining publishers,
#: versions and PubMed ids positionally inside a cell is the shape that breaks quietly.
#: Deliberately ragged: one row per species today, three for a species built from three.
COFACTOR_SOURCE_METADATA_RESOURCE = f"{COFACTOR_SUBDIR}/cofactor_source_metadata.tsv"

#: The four columns every **Cofactor table** carries under the same name and in the same
#: place: the **Gene id stem**, the symbol, the cofactor flag and the source. Uniform in
#: position; everything after them is one publisher's own vocabulary under a namespaced
#: name and is never compared with another's (ADR-0014).
UNIFORM_COLUMNS: tuple[str, ...] = ("gene_id_stem", "symbol", "is_cofactor", "source")

#: The three legal values of the uniform ``source`` column, and the whole vocabulary.
#: ``both`` says two publishers listed the gene — agreement on **membership only**, never
#: on how either of them classified it.
ANIMALTFDB, EPIFACTORS, BOTH = "animaltfdb", "epifactors", "both"
SOURCES: tuple[str, ...] = (ANIMALTFDB, EPIFACTORS, BOTH)

#: A source that supplies identifiers and no membership: it makes a stem readable and
#: lists nobody, so it never spells a row's ``source`` and is only ever cited.
HGNC = "hgnc"

#: What the *provenance* table's own ``source`` column may say — every source a table
#: was built from and owes a citation to. A different list from :data:`SOURCES` on
#: purpose, in both directions: ``both`` is a fact about a row of a table and describes
#: no publisher, and ``hgnc`` describes a publisher that listed no gene. Human's 442
#: EpiFactors-only stems exist only because HGNC said so, which is why identifiers earn
#: a citation of their own rather than passing as an implementation detail.
CITED_SOURCES: tuple[str, ...] = (ANIMALTFDB, EPIFACTORS, HGNC)

#: Where the cofactor flag and the source sit in every table, which the uniform four fix.
_IS_COFACTOR = UNIFORM_COLUMNS.index("is_cofactor")
_SOURCE = UNIFORM_COLUMNS.index("source")

#: The species-keyed provenance table's columns, in table order.
_METADATA_COLUMNS: tuple[str, ...] = ("species", "ncbi_taxid", "file", "sha256")

#: The source-keyed provenance table's columns, in table order.
_SOURCE_METADATA_COLUMNS: tuple[str, ...] = (
    "species",
    "source",
    "publisher",
    "version",
    "pubmed_id",
    "source_url",
)

#: Provenance columns holding a number rather than text.
_NUMERIC_METADATA_COLUMNS = frozenset({"ncbi_taxid", "pubmed_id"})

#: What to do about any of it, named in every message this module raises.
_REBUILD = "scripts/build_tf_cofactor.py"


class CofactorTableError(ValueError):
    r"""A shipped **Cofactor table** cannot be read, so it is not allowed to answer.

    A **packaging defect** and not a caller error: these files ship inside the package
    and are written by a generator, so a header that does not lead with the uniform
    four, a repeated **Gene id stem**, a source spelled a way the vocabulary does not,
    or a table with no row in either provenance table are all faults in what was
    committed here. A :class:`ValueError`, because a file that says something the
    format does not is a bad value rather than a broken program.

    The message names the file and the repair, since regenerating that file is the only
    thing anyone can do about it.

    Examples
    --------
    >>> provenance = cofactor_metadata()[0]
    >>> try:
    ...     parse_cofactor_table("gene_id\n", provenance=provenance)
    ... except CofactorTableError as error:
    ...     print(provenance.file in str(error))
    True
    """


@dataclass(frozen=True)
class CofactorSource:
    """One publisher's contribution to one species' table (one row of the ragged table).

    The single declaration of what a source's provenance consists of. Every column is
    required — a source nobody can cite is one this package may not redistribute — and
    a species has as many of these as it has publishers, which is why they live in a
    table of their own rather than joined into one cell of another.

    Attributes
    ----------
    species : str
        The species, as the assembly metadata table spells it.
    source : str
        Which source this row is about, spelled as the uniform ``source`` column spells
        one where the two vocabularies meet. One of :data:`CITED_SOURCES`: never
        ``both``, which is a fact about a *row* of the table rather than a publisher,
        and possibly ``hgnc``, which lists no gene and supplies identifiers only.
    publisher : str
        Who published it, and who is to be cited for it.
    version : str
        The publisher's own release identifier, e.g. ``"4.0"``.
    pubmed_id : int
        PubMed id of the paper to cite.
    source_url : str
        Where the publisher's own file was downloaded from.

    Examples
    --------
    >>> cofactor_table("Mus musculus").provenance.sources[0].pubmed_id
    36268869
    """

    species: str
    source: str
    publisher: str
    version: str
    pubmed_id: int
    source_url: str

    @classmethod
    def from_row(cls, row: Mapping[str, str], *, origin: str) -> CofactorSource:
        """Build a record from one row of the source-keyed provenance table.

        Parameters
        ----------
        row : mapping of str to str
            Column name to cell, as the shipped TSV spells one.
        origin : str
            Where the row came from; named in every message, since fixing that file is
            the only repair.

        Returns
        -------
        CofactorSource
            The record the row spells.

        Raises
        ------
        CofactorTableError
            If a column is missing or blank, a numeric column holds something
            :class:`int` cannot read, or the source is outside :data:`CITED_SOURCES`.
            The message names the column.

        Examples
        --------
        >>> row = dict(
        ...     species="Tiny beast",
        ...     source="animaltfdb",
        ...     publisher="Someone et al. 1999",
        ...     version="v1",
        ...     pubmed_id="2",
        ...     source_url="https://example.org/beast",
        ... )
        >>> CofactorSource.from_row(row, origin="cofactor_source_metadata.tsv").version
        'v1'
        """
        record = cls(
            **{name: _parse_cell(name, row, origin=origin) for name in _SOURCE_METADATA_COLUMNS}
        )
        if record.source not in CITED_SOURCES:
            raise CofactorTableError(
                f"{origin} names the source {record.source!r} for {record.species!r}, and the "
                f"vocabulary a provenance row is spelled from is {list(CITED_SOURCES)}. A row "
                f"here describes one source and what to cite for it, so {BOTH!r} — which is a "
                f"fact about a row of the table and about no publisher — is not one of them. "
                f"Re-run {_REBUILD} for that species, which writes the table and its provenance "
                f"rows together."
            )
        return record

    def attribution(self) -> str:
        """Return the one line to print for this publisher's part of a table.

        Returns
        -------
        str
            Publisher, version, PubMed id and source URL.

        Examples
        --------
        >>> print(cofactor_table("Mus musculus").provenance.sources[0].attribution())
        AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Mus_musculus_Cof
        """
        return (
            f"{self.publisher} {self.version} (PMID {self.pubmed_id}) \N{EM DASH} {self.source_url}"
        )


@dataclass(frozen=True)
class CofactorProvenance:
    """Where one species' **Cofactor table** came from, and who to cite for it.

    Two provenance tables read as one record: the species-keyed row says which file and
    what its bytes hash to, and the ragged source-keyed rows say who published each part
    of it. A species has at least one source and may have several; they arrive in the
    provenance table's own row order.

    Attributes
    ----------
    species : str
        The species, as the assembly metadata table spells it — ``"Mus musculus"``. Its
        slug names the table's file.
    ncbi_taxid : int
        NCBI taxonomy id for that species.
    file : str
        The table's file name within :data:`COFACTOR_SUBDIR`.
    sha256 : str
        Digest of the **unpacked** table — the TSV inside the gzip, not the gzip bytes,
        so a copy recompressed elsewhere still matches (ADR-0006).
    sources : tuple of CofactorSource
        One record per source this species' table was built from and owes a citation
        to — every publisher that listed genes, and any that only made a stem
        readable. Human has three where mouse and worm have one.

    Examples
    --------
    >>> provenance = cofactor_table("Caenorhabditis elegans").provenance
    >>> provenance.ncbi_taxid, len(provenance.sources)
    (6239, 1)
    >>> [source.source for source in cofactor_table("Homo sapiens").provenance.sources]
    ['animaltfdb', 'epifactors', 'hgnc']
    """

    species: str
    ncbi_taxid: int
    file: str
    sha256: str
    sources: tuple[CofactorSource, ...]

    @classmethod
    def from_row(
        cls, row: Mapping[str, str], *, sources: tuple[CofactorSource, ...], origin: str
    ) -> CofactorProvenance:
        """Build a record from one row of the species-keyed table and its source rows.

        Parameters
        ----------
        row : mapping of str to str
            Column name to cell, as the shipped TSV spells one.
        sources : tuple of CofactorSource
            The source-keyed rows for the same species, in that table's row order.
        origin : str
            Where the row came from; named in every message, since fixing that file is
            the only repair.

        Returns
        -------
        CofactorProvenance
            The record the row spells.

        Raises
        ------
        CofactorTableError
            If a column is missing or blank, a numeric column holds something
            :class:`int` cannot read, or no source is given. The message names the
            column.

        Examples
        --------
        >>> row = dict(
        ...     species="Tiny beast",
        ...     ncbi_taxid="1",
        ...     file="tiny_beast.cofactor_table.tsv.gz",
        ...     sha256="0" * 64,
        ... )
        >>> source = CofactorSource.from_row(
        ...     dict(
        ...         species="Tiny beast",
        ...         source="animaltfdb",
        ...         publisher="Someone et al. 1999",
        ...         version="v1",
        ...         pubmed_id="2",
        ...         source_url="https://example.org/beast",
        ...     ),
        ...     origin="cofactor_source_metadata.tsv",
        ... )
        >>> CofactorProvenance.from_row(
        ...     row, sources=(source,), origin="cofactor_metadata.tsv"
        ... ).ncbi_taxid
        1
        """
        if not sources:
            raise CofactorTableError(
                f"{origin} names {row.get('species')!r} and {COFACTOR_SOURCE_METADATA_RESOURCE} "
                f"gives it no source, so nothing says who published that table or what to cite "
                f"for it — and citing the publisher is the condition on redistributing one here. "
                f"Re-run {_REBUILD} for that species, which writes both provenance rows."
            )
        fields = {name: _parse_cell(name, row, origin=origin) for name in _METADATA_COLUMNS}
        return cls(**fields, sources=sources)

    def attribution(self) -> str:
        """Return the one line to print beside anything this species' table answered.

        What a caller owes the publishers, rendered once here so the CLI, a notebook and
        an error message all say it the same way. One line per species however many
        publishers contributed: each source's own line, joined.

        Returns
        -------
        str
            Every source's publisher, version, PubMed id and source URL.

        Examples
        --------
        >>> print(cofactor_table("Mus musculus").provenance.attribution())
        AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Mus_musculus_Cof
        """
        return " \N{MIDDLE DOT} ".join(source.attribution() for source in self.sources)


@dataclass(frozen=True)
class CofactorTable:
    """One species' cofactors as shipped: every gene a publisher listed, and what it said.

    What one shipped file says, read back and frozen. It carries the whole table rather
    than the listed-positive part of it, because ``is_cofactor`` is a column and not the
    file's mere existence — :attr:`cofactor_stems` is the filter.

    Attributes
    ----------
    species : str
        The species, as its provenance row spells it.
    provenance : CofactorProvenance
        Where the table came from and who to cite for it.
    columns : tuple of str
        The table's columns, in file order. The first four are always
        :data:`UNIFORM_COLUMNS`; the rest are one publisher's own, namespaced.
    rows : tuple of tuple of (str or None)
        One tuple per gene, in the publisher's own row order, parallel to
        :attr:`columns`. A blank cell — the publisher recorded nothing there — is
        ``None``.
    gene_id_stems : tuple of str
        Every **Gene id stem** the table lists, in row order. Unique within a table;
        never blank.
    cofactor_stems : tuple of str
        The stems the table says are cofactors, in row order.

    Examples
    --------
    >>> mouse = cofactor_table("Mus musculus")
    >>> mouse.columns
    ('gene_id_stem', 'symbol', 'is_cofactor', 'source', 'animaltfdb_family', 'animaltfdb_category')
    >>> len(mouse), len(mouse.cofactor_stems)
    (970, 970)
    >>> mouse.frame().loc[0, "animaltfdb_category"]
    'Other Cofactors'
    """

    species: str
    provenance: CofactorProvenance
    columns: tuple[str, ...]
    rows: tuple[tuple[str | None, ...], ...]

    def __len__(self) -> int:
        """Return how many genes the table lists.

        Examples
        --------
        >>> len(cofactor_table("Caenorhabditis elegans"))
        317
        """
        return len(self.rows)

    @property
    def gene_id_stems(self) -> tuple[str, ...]:
        """Return every **Gene id stem** the table lists, in row order.

        Examples
        --------
        >>> cofactor_table("Caenorhabditis elegans").gene_id_stems[0]
        'WBGene00000064'
        """
        # A stem is never blank — that is checked as the file is read — so ``or ""``
        # only narrows the type of a cell that is always text.
        return tuple(row[0] or "" for row in self.rows)

    @property
    def cofactor_stems(self) -> tuple[str, ...]:
        """Return the stems this table says are **Transcription cofactor**s, in row order.

        The default a **TF cofactor list** is built from, and the peer of a census's
        assessed-positive genes. It is every stem in the table today, because no
        publisher here releases a rejected set — which is that publisher's shape and not
        a promise of this format: a source that did record one would ship ``no`` rows,
        and they would be excluded here rather than needing a new column.

        Examples
        --------
        >>> len(cofactor_table("Mus musculus").cofactor_stems)
        970
        """
        return tuple(row[0] or "" for row in self.rows if row[_IS_COFACTOR] == TRUE_CELL)

    def frame(self) -> pd.DataFrame:
        """Return the table as a fresh :class:`~pandas.DataFrame`, one row per gene.

        Built for the caller each time it is asked for, so mutating it cannot reach the
        cached table behind it. Every column is text with ``None`` for a blank cell,
        except the uniform cofactor flag, which reads back as a boolean.

        Returns
        -------
        pandas.DataFrame
            The table's columns in file order, indexed from zero in row order.

        Examples
        --------
        >>> frame = cofactor_table("Mus musculus").frame()
        >>> int(frame["is_cofactor"].sum()), frame["source"].unique().tolist()
        (970, ['animaltfdb'])
        >>> frame["animaltfdb_category"].nunique()
        6
        """
        frame = pd.DataFrame(list(self.rows), columns=pd.Index(self.columns), dtype=object)
        frame["is_cofactor"] = [row[_IS_COFACTOR] == TRUE_CELL for row in self.rows]
        return frame


@cache
def cofactor_metadata() -> tuple[CofactorProvenance, ...]:
    """Return the provenance of every **Cofactor table** the shipped tables record.

    Both provenance tables read as one answer, in the species-keyed table's own row
    order: where each table came from, who published each part of it and what to cite.
    Read once and cached; the records are frozen, so the tuple is safe to hold on to.

    Returns
    -------
    tuple of CofactorProvenance
        One record per species, each carrying its publishers' records.

    Raises
    ------
    CofactorTableError
        If either table cannot be read, a species has no source row, or a source row
        names a species the other table does not. The message names the column or the
        species.

    Examples
    --------
    >>> {record.species for record in cofactor_metadata()} >= {"Mus musculus"}
    True
    >>> print(cofactor_metadata()[0].attribution())
    AnimalTFDB 4.0 (PMID 36268869) — https://guolab.wchscu.cn/AnimalTFDB4_static/download/Cof_list_final/Caenorhabditis_elegans_Cof
    """
    sources = files("genome").joinpath(COFACTOR_SOURCE_METADATA_RESOURCE)
    species = files("genome").joinpath(COFACTOR_METADATA_RESOURCE)
    return _read_metadata(
        species.read_text(encoding="utf-8"),
        sources=_read_source_metadata(sources.read_text(encoding="utf-8"), origin=str(sources)),
        origin=str(species),
    )


@cache
def cofactor_species() -> tuple[str, ...]:
    """Return every species a **Cofactor table** ships for, as its file name spells it.

    What can be asked about at all, and the answer an error names when a species cannot
    be. The directory is enumerated rather than any list of species being kept in code,
    so adding a species is dropping in a file. The names are slugs —
    ``mus_musculus`` — and :func:`cofactor_table` takes either those or the species as
    the assembly metadata table spells it.

    Returns
    -------
    tuple of str
        The species slugs, sorted. Empty only if the package ships no cofactor table.

    Examples
    --------
    >>> cofactor_species()
    ('caenorhabditis_elegans', 'homo_sapiens', 'mus_musculus')
    """
    directory = files("genome").joinpath(COFACTOR_SUBDIR)
    if not directory.is_dir():
        return ()
    return tuple(
        sorted(
            entry.name[: -len(COFACTOR_SUFFIX)]
            for entry in directory.iterdir()
            if entry.name.endswith(COFACTOR_SUFFIX)
        )
    )


@cache
def cofactor_table(species: str) -> CofactorTable | None:
    """Return the **Cofactor table** shipped for ``species``, or ``None``.

    The raw absence, and the only place ``None`` is an acceptable way to say it: this is
    the layer below the one a caller touches, and everything above it turns the ``None``
    into an error naming the species that do have a table, so that *nobody has published
    one for this species* can never be read as *this species has no cofactors*.

    The species is slugged and then looked up among what :func:`cofactor_species` found,
    rather than joined onto the resource directory, so a name shaped like a path finds
    nothing instead of walking out of it. The file is validated as it is read. Read once
    per species and cached; everything it returns is frozen.

    Parameters
    ----------
    species : str
        The species to look up, either as the assembly metadata table spells it
        (``"Mus musculus"``) or as its slug (``"mus_musculus"``).

    Returns
    -------
    CofactorTable or None
        The table, or ``None`` when none ships for that species. ``None`` is legal and
        ordinary — most species have no published cofactor list.

    Raises
    ------
    CofactorTableError
        If a table ships for that species and cannot be read, or ships with no row in
        the provenance table; the message names the file.

    Examples
    --------
    >>> cofactor_table("Mus musculus").species
    'Mus musculus'
    >>> cofactor_table("mus_musculus") == cofactor_table("Mus musculus")
    True
    >>> cofactor_table("Danio rerio") is None
    True
    """
    slug = species_slug(species)
    if slug not in cofactor_species():
        return None
    resource = files("genome").joinpath(COFACTOR_SUBDIR, f"{slug}{COFACTOR_SUFFIX}")
    origin = str(resource)
    provenance = next(
        (record for record in cofactor_metadata() if species_slug(record.species) == slug), None
    )
    if provenance is None:
        raise CofactorTableError(
            f"{origin} ships with no row in {COFACTOR_METADATA_RESOURCE}, so nothing says whose "
            f"table it is or what to cite for it — and citing the publisher is the condition on "
            f"redistributing one here. Re-run {_REBUILD} for {slug!r}, which writes the table "
            f"and its provenance rows together."
        )
    # These are tens of kilobytes of shipped rows, so the seam is between the bytes and
    # the format: unpacking happens here and the parse below is a pure function of text.
    return parse_cofactor_table(
        gzip.decompress(resource.read_bytes()).decode("utf-8"),
        provenance=provenance,
        origin=origin,
    )


def parse_cofactor_table(
    text: str, *, provenance: CofactorProvenance, origin: str | None = None
) -> CofactorTable:
    r"""Read one **Cofactor table**'s text, holding it to what a shipped table promises.

    A pure function from text to a table: it opens nothing, downloads nothing and
    decompresses nothing. Public and separate from the resource it came out of, as
    :func:`~genome.tf.link.parse_motif_link_table` is, so every way a shipped file can
    be malformed is reachable without writing a broken one into the package.

    Parameters
    ----------
    text : str
        The whole table, unpacked. These are tens of kilobytes, so they are read whole.
    provenance : CofactorProvenance
        Where the text came from and who to cite for it. Passed in rather than read off
        the rows, because no row carries it: a table says which publisher listed each
        gene and the provenance says who that publisher is.
    origin : str or None, optional
        Where the text came from, named in every message since regenerating that file is
        the only repair. Defaults to ``provenance.file``.

    Returns
    -------
    CofactorTable
        The table the text spells, in file order.

    Raises
    ------
    CofactorTableError
        If the header does not lead with :data:`UNIFORM_COLUMNS`, a column is named
        twice, a row holds the wrong number of cells, the cofactor flag or the source is
        spelled a way no table spells one, a **Gene id stem** is blank or repeated, or
        the file declares columns and no genes.

    Examples
    --------
    >>> header = "\t".join(UNIFORM_COLUMNS)
    >>> row = "ENSTEST0001\tAbc1\tyes\tanimaltfdb"
    >>> table = parse_cofactor_table(
    ...     f"{header}\n{row}\n", provenance=cofactor_metadata()[0], origin="tiny.tsv"
    ... )
    >>> table.gene_id_stems, table.cofactor_stems
    (('ENSTEST0001',), ('ENSTEST0001',))
    """
    where = provenance.file if origin is None else origin
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise CofactorTableError(
            f"{where} is empty. A cofactor table carries a header and at least one gene; absence "
            f"is spelled by shipping no file at all, and an empty one is the second spelling of "
            f"it that would read as *this species has no cofactors*. Re-run {_REBUILD} for that "
            f"species, or remove the file."
        )
    columns = tuple(lines[0].split("\t"))
    _check_columns(columns, origin=where)
    rows = tuple(
        _read_row(line, columns, number, origin=where)
        for number, line in enumerate(lines[1:], start=2)
    )
    if not rows:
        raise CofactorTableError(
            f"{where} declares columns and no genes. A table listing nothing says no more than "
            f"an absent file does — re-run {_REBUILD} for that species, or remove the file."
        )
    _check_stems(rows, origin=where)
    return CofactorTable(
        species=provenance.species, provenance=provenance, columns=columns, rows=rows
    )


def _read_metadata(
    text: str, *, sources: tuple[CofactorSource, ...], origin: str
) -> tuple[CofactorProvenance, ...]:
    """Read the species-keyed provenance table, joining each species' source rows onto it."""
    rows = _read_metadata_rows(text, _METADATA_COLUMNS, origin=origin)
    named = {row["species"] for row in rows}
    orphan = sorted({source.species for source in sources} - named)
    if orphan:
        raise CofactorTableError(
            f"{COFACTOR_SOURCE_METADATA_RESOURCE} names the species {orphan}, and {origin} does "
            f"not, so nothing says which file those sources describe. The two provenance tables "
            f"are keyed on the same spelling of a species — re-run {_REBUILD}, which writes both."
        )
    return tuple(
        CofactorProvenance.from_row(
            row,
            sources=tuple(source for source in sources if source.species == row["species"]),
            origin=origin,
        )
        for row in rows
    )


def _read_source_metadata(text: str, *, origin: str) -> tuple[CofactorSource, ...]:
    """Read the source-keyed provenance table, refusing a species-and-source named twice."""
    rows = _read_metadata_rows(text, _SOURCE_METADATA_COLUMNS, origin=origin)
    records = tuple(CofactorSource.from_row(row, origin=origin) for row in rows)
    keys = [(record.species, record.source) for record in records]
    if len(set(keys)) != len(keys):
        repeated = sorted({key for key in keys if keys.count(key) > 1})
        raise CofactorTableError(
            f"{origin} names these species and sources more than once: {repeated}. One row says "
            f"what to cite for one publisher's part of one species' table, so two would let a "
            f"caller cite either. Re-run {_REBUILD} for that species."
        )
    return records


def _read_metadata_rows(
    text: str, columns: tuple[str, ...], *, origin: str
) -> list[dict[str, str]]:
    """Read one provenance table's rows as text, validating its header as it goes."""
    lines = text.splitlines()
    if not lines:
        raise CofactorTableError(
            f"{origin} is empty, and a cofactor table with no provenance is one nobody can cite. "
            f"Re-run {_REBUILD}, which writes it."
        )
    header = tuple(lines[0].split("\t"))
    if header != columns:
        raise CofactorTableError(
            f"{origin} carries the columns {list(header)} where that provenance table's are "
            f"{list(columns)}. Re-run {_REBUILD}, which writes both tables with the columns "
            f"this reader expects."
        )
    rows = []
    for number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        cells = line.split("\t")
        if len(cells) != len(header):
            raise CofactorTableError(
                f"{origin} line {number} holds {len(cells)} cells where the header declares "
                f"{len(header)}. Re-run {_REBUILD}."
            )
        rows.append(dict(zip(header, cells, strict=True)))
    return rows


def _parse_cell(name: str, row: Mapping[str, str], *, origin: str) -> Any:
    """Return one provenance cell, parsed by its column and never blank."""
    text = row.get(name, "").strip()
    if not text:
        raise CofactorTableError(
            f"{origin} leaves the {name!r} column blank for {row.get('species')!r}. Every "
            f"provenance column is required: a table nobody can cite is one this package may "
            f"not redistribute. Fill that cell in, or re-run {_REBUILD}."
        )
    if name not in _NUMERIC_METADATA_COLUMNS:
        return text
    try:
        return int(text)
    except ValueError as error:
        raise CofactorTableError(
            f"{origin} holds {text!r} in the {name!r} column, which is not a number. Fix that cell."
        ) from error


def _check_columns(columns: tuple[str, ...], *, origin: str) -> None:
    """Hold a table's header to the uniform four in front and distinct names after."""
    if columns[: len(UNIFORM_COLUMNS)] != UNIFORM_COLUMNS:
        raise CofactorTableError(
            f"{origin} leads with the columns {list(columns[: len(UNIFORM_COLUMNS)])} where "
            f"every cofactor table leads with {list(UNIFORM_COLUMNS)}. Those four are the only "
            f"columns one table shares with another, so a file without them cannot be read as "
            f"one. Re-run {_REBUILD} for that species."
        )
    if len(set(columns)) != len(columns):
        raise CofactorTableError(
            f"{origin} names a column twice: {list(columns)}. Each publisher's columns are "
            f"namespaced so that two of them never collide — re-run {_REBUILD}."
        )


def _read_row(
    line: str, columns: tuple[str, ...], number: int, *, origin: str
) -> tuple[str | None, ...]:
    """Read one gene's row, blank cells becoming ``None`` and the uniform flags being checked."""
    cells = line.split("\t")
    if len(cells) != len(columns):
        raise CofactorTableError(
            f"{origin} line {number} holds {len(cells)} cells where its header declares "
            f"{len(columns)}. A cofactor table is a plain TSV with no quoting, so a cell "
            f"carrying a tab is a defect in the generator rather than something to parse "
            f"around. Re-run {_REBUILD} for that species."
        )
    if cells[_IS_COFACTOR] not in (TRUE_CELL, FALSE_CELL):
        raise CofactorTableError(
            f"{origin} line {number} spells its cofactor flag {cells[_IS_COFACTOR]!r}, and a "
            f"table spells it {TRUE_CELL!r} or {FALSE_CELL!r}. The flag is one of the four "
            f"uniform columns, so its spelling is this package's and not the publisher's — "
            f"re-run {_REBUILD} for that species."
        )
    if cells[_SOURCE] not in SOURCES:
        raise CofactorTableError(
            f"{origin} line {number} names the source {cells[_SOURCE]!r}, and the vocabulary is "
            f"{list(SOURCES)}. It is closed on purpose: it says which publisher listed the gene, "
            f"and a value outside it names nobody the provenance table can be asked about. "
            f"Re-run {_REBUILD} for that species."
        )
    return tuple(cell if cell else None for cell in cells)


def _check_stems(rows: tuple[tuple[str | None, ...], ...], *, origin: str) -> None:
    """Hold every **Gene id stem** to being present and naming its gene once."""
    seen: set[str] = set()
    for number, row in enumerate(rows, start=2):
        stem = row[0]
        if stem is None:
            raise CofactorTableError(
                f"{origin} line {number} carries no gene id stem. A cofactor table is keyed by "
                f"stem, so a row without one cannot be looked up or resolved against an "
                f"annotation. Re-run {_REBUILD} for that species."
            )
        if stem in seen:
            raise CofactorTableError(
                f"{origin} names the gene id stem {stem!r} more than once. One row per gene, so "
                f"two rows for one stem would let a caller read either — reconcile the "
                f"publishers' files and re-run {_REBUILD}."
            )
        seen.add(stem)
