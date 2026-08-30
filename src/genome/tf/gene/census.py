"""The shipped censuses — which genes one published census judges transcription factors.

**Attribution.** The human census is Lambert *et al.*, "The Human Transcription
Factors", *Cell* 172(4):650-665, 2018 (PMID 29425488), database extract v_1.01,
redistributed here from https://humantfs.ccbr.utoronto.ca/. The mouse census is
AnimalTFDB 4.0, Shen *et al.*, *Nucleic Acids Research* 51(D1):D39-D45, 2023
(PMID 36268869), redistributed from https://guolab.wchscu.cn/. Every verdict in
either is its publisher's and none of it is this package's. Cite the publisher whose
census you used; :meth:`CensusProvenance.attribution` renders the line to print, and
the provenance table beside the data carries the same facts for every census that
ships.

One **TF gene table** per species ships inside the package under
``data/tf_gene/<species>.tf_gene_table.tsv.gz``, and they are found by enumerating
that directory rather than by any list of species in code — so adding a species is
dropping in a file. The name is the species **slug**: the assembly metadata table's
own spelling, lower-cased with each run of non-alphanumerics turned into one
underscore, so ``Homo sapiens`` is ``homo_sapiens``. That is the convention because
species is the key an assembly's metadata already carries (ADR-0003), and
:func:`tf_gene_table` accepts either spelling. The taxid would have been the other
choice and was not taken: a directory listing of ``9606.tf_gene_table.tsv.gz`` says
nothing to a reader, and the taxid is recorded in the provenance table anyway.

Four columns are uniform across every census and lead every file — the **Gene id
stem**, the symbol, the TF flag and the **DBD family** (:data:`UNIFORM_COLUMNS`).
Everything after them is its publisher's own column under a snake_case spelling of
its published name, and is never compared with another publisher's: the family
vocabularies in particular are deliberately not crosswalked (ADR-0014), so
Lambert's ``ARID/BRIGHT`` and AnimalTFDB's ``ARID`` are two publishers' spellings
and not one family asserted twice. How many columns follow the four is the
publisher's business too — Lambert publishes judgements to spare and AnimalTFDB
publishes none beyond the family, so mouse ships the uniform four and nothing else.
A blank cell is the publisher recording nothing and reads back as ``None``, the
reading every other table in this package gives a blank cell.

**Absence is not emptiness**, and this module is where the distinction starts. A
gene the census never assessed and a gene it assessed and rejected are different
answers, so the rejected genes ship too: Lambert assessed 2,765 genes and judged
1,639 of them TFs, and the other 1,126 are a verdict rather than a silence. Whether
a census *has* rejections is its publisher's shape and not this package's — every
gene AnimalTFDB lists is one it accepts, so mouse's 1,611 rows are all positive and
the genes it left out get no verdict at all rather than a fabricated ``no``.
:func:`tf_gene_table` answers ``None`` for a species no census ships for — the raw
absence, and the one place ``None`` is how it is said, exactly as
:func:`genome.gene_list.curated_gene_list` says it, because this is the layer below
the one a caller touches. Everything above turns it into an error naming the
species that do have a census.

This module is pure. It reads shipped package resources and nothing else — never
the **Data dir**, never the network — and the files themselves are built by
``scripts/build_tf_census.py``, which lives outside the wheel. A shipped file that
cannot be trusted never answers: it is validated as it is read and raises
:class:`TFGeneTableError` naming the file and what is wrong with it, since a census
that ships broken is a defect in this package rather than anything a caller did.

Examples
--------
>>> from genome.tf.gene import census_species, tf_gene_table
>>> "homo_sapiens" in census_species()
True
>>> census = tf_gene_table("Homo sapiens")
>>> len(census), len(census.assessed_positive)
(2765, 1639)
>>> census.provenance.publisher
'Lambert et al. 2018'
>>> mouse = tf_gene_table("Mus musculus")
>>> len(mouse), len(mouse.assessed_positive)
(1611, 1611)
>>> tf_gene_table("Caenorhabditis elegans") is None
True
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from typing import Any, get_type_hints

import pandas as pd

# Re-exported rather than defined here. Naming a file after a species and spelling a flag
# `yes` are conventions every **Shipped table** in this package shares and none of them
# owns, so they live with the reader that holds every such table to them. Kept importable
# from here because this is where they used to be.
from genome.shipped import (
    FALSE_CELL,
    TRUE_CELL,
    ShippedTable,
    ShippedTableError,
    species_slug,
)

__all__ = [
    "CENSUS_FORMAT",
    "CENSUS_METADATA_FORMAT",
    "CENSUS_METADATA_RESOURCE",
    "CENSUS_SUBDIR",
    "CENSUS_SUFFIX",
    "FALSE_CELL",
    "TRUE_CELL",
    "UNIFORM_COLUMNS",
    "CensusProvenance",
    "TFGeneTable",
    "TFGeneTableError",
    "census_metadata",
    "census_species",
    "species_slug",
    "tf_gene_table",
]

#: Directory inside the package holding one **TF gene table** per species, plus the
#: provenance table beside them.
CENSUS_SUBDIR = "data/tf_gene"

#: What one of those files is called: the species slug, then this. The slug is what a
#: census is looked up by, so the two halves are never spelled apart.
CENSUS_SUFFIX = ".tf_gene_table.tsv.gz"

#: The provenance table beside the censuses, in the shape ``assembly_metadata.tsv``
#: and ``annotation_metadata.tsv`` already use — and **plain**, where the censuses are
#: gzipped. That is the convention for shipped data here: bulk gzipped, small metadata
#: plain, because a file short enough to be read and edited by hand is worth more as a
#: reviewable diff than as the bytes gzip would save.
CENSUS_METADATA_RESOURCE = f"{CENSUS_SUBDIR}/census_metadata.tsv"

#: The four columns every census carries under the same name and in the same place:
#: the **Gene id stem**, the symbol, the TF flag and the **DBD family**. Uniform in
#: position and deliberately not in content — two publishers' family vocabularies are
#: not crosswalked (ADR-0014).
UNIFORM_COLUMNS: tuple[str, ...] = ("gene_id_stem", "symbol", "is_tf", "dbd_family")

#: Where the TF flag sits in every census, which the uniform four fix.
_IS_TF = UNIFORM_COLUMNS.index("is_tf")

#: The provenance table's columns, in table order.
_METADATA_COLUMNS: tuple[str, ...] = (
    "species",
    "ncbi_taxid",
    "file",
    "publisher",
    "version",
    "pubmed_id",
    "family_column",
    "source_url",
    "sha256",
)

#: What to do about any of it, named in every message this module raises.
_REBUILD = "re-run scripts/build_tf_census.py"


class TFGeneTableError(ShippedTableError):
    r"""A shipped **TF gene table** cannot be read, so it is not allowed to answer.

    A **packaging defect** and not a caller error: these files ship inside the
    package and are written by a generator, so a header that does not lead with the
    uniform four, a repeated **Gene id stem**, a TF flag spelled a way no census
    spells one, or a census with no row in the provenance table are all faults in
    what was committed here. A :class:`ValueError`, because a file that says
    something the format does not is a bad value rather than a broken program.

    The message names the file and what is wrong with it, since regenerating that
    file is the only thing anyone can do about it.

    Examples
    --------
    >>> try:
    ...     _read_metadata("species\n", origin="census_metadata.tsv")
    ... except TFGeneTableError as error:
    ...     print("census_metadata.tsv" in str(error))
    True
    """


@dataclass(frozen=True)
class CensusProvenance:
    """Where one census came from (one row of the provenance table).

    The single declaration of what a census's provenance consists of: the table is
    parsed through these fields, column by column, exactly as
    :class:`~genome.metadata.AssemblyMetadata` parses its own. Every column is
    required — a census with nothing said about its origin is one nobody can cite,
    and citing it is the condition on redistributing it here.

    Attributes
    ----------
    species : str
        The species, as the assembly metadata table spells it — ``"Homo sapiens"``.
        Its slug names the census's file.
    ncbi_taxid : int
        NCBI taxonomy id for that species.
    file : str
        The census's file name within :data:`CENSUS_SUBDIR`.
    publisher : str
        Who published the census, and who is to be cited for it.
    version : str
        The publisher's own release identifier, e.g. ``"v_1.01"``.
    pubmed_id : int
        PubMed id of the paper to cite.
    family_column : str
        The publisher's own name for the column the **DBD family** was taken from —
        ``DBD`` for Lambert. Recorded because the vocabularies are deliberately not
        crosswalked (ADR-0014), so a reader grouping by family has to know whose
        family names they are.
    source_url : str
        Where the publisher's own file was downloaded from.
    sha256 : str
        Digest of the **unpacked** census — the TSV inside the gzip, not the gzip
        bytes, so a copy recompressed elsewhere still matches (ADR-0006).

    Examples
    --------
    >>> tf_gene_table("Homo sapiens").provenance.pubmed_id
    29425488
    """

    species: str
    ncbi_taxid: int
    file: str
    publisher: str
    version: str
    pubmed_id: int
    family_column: str
    source_url: str
    sha256: str

    @classmethod
    def from_row(cls, row: Mapping[str, str], *, origin: str) -> CensusProvenance:
        """Build a record from one row of the provenance table.

        Parameters
        ----------
        row : mapping of str to str
            Column name to cell, as the shipped TSV spells one.
        origin : str
            Where the row came from; named in every message, since fixing that file
            is the only repair.

        Returns
        -------
        CensusProvenance
            The record the row spells.

        Raises
        ------
        TFGeneTableError
            If a column is missing or blank, or a numeric column holds something
            :class:`int` cannot read. The message names the column.

        Examples
        --------
        >>> row = dict(
        ...     species="Homo sapiens",
        ...     ncbi_taxid="9606",
        ...     file="homo_sapiens.tf_gene_table.tsv.gz",
        ...     publisher="Lambert et al. 2018",
        ...     version="v_1.01",
        ...     pubmed_id="29425488",
        ...     family_column="DBD",
        ...     source_url="https://example.org/x.csv",
        ...     sha256="0" * 64,
        ... )
        >>> CensusProvenance.from_row(row, origin="census_metadata.tsv").ncbi_taxid
        9606
        """
        return cls(**CENSUS_METADATA_FORMAT.record(row, _PROVENANCE_TYPES, origin=origin))

    def attribution(self) -> str:
        """Return the one line to print beside anything this census answered.

        What a caller owes the publisher, rendered once here so the CLI, a notebook
        and an error message all say it the same way.

        Returns
        -------
        str
            Publisher, version, PubMed id and source URL.

        Examples
        --------
        >>> print(tf_gene_table("Homo sapiens").provenance.attribution())
        Lambert et al. 2018 v_1.01 (PMID 29425488) — https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
        """
        return (
            f"{self.publisher} {self.version} (PMID {self.pubmed_id}) \N{EM DASH} {self.source_url}"
        )


@dataclass(frozen=True)
class TFGeneTable:
    """One species' census as shipped: every gene it assessed, and what it said.

    What one shipped file says, read back and frozen. It carries the whole census
    and not the assessed-positive part of it, because a gene the census rejected is
    a verdict a caller may want — :attr:`assessed_positive` is the filter, not the
    file.

    Attributes
    ----------
    species : str
        The species, as its provenance row spells it.
    provenance : CensusProvenance
        Where the census came from and who to cite for it.
    columns : tuple of str
        The census's columns, in file order. The first four are always
        :data:`UNIFORM_COLUMNS`; the rest are the publisher's own.
    rows : tuple of tuple of (str or None)
        One tuple per gene, in the publisher's own row order, parallel to
        :attr:`columns`. A blank cell — the publisher recorded nothing there — is
        ``None``.
    gene_id_stems : tuple of str
        Every **Gene id stem** the census assessed, in row order. Unique within a
        census; never blank.
    assessed_positive : tuple of str
        The stems the census judged transcription factors, in row order. A strict
        subset of :attr:`gene_id_stems` for a census that records rejections, and
        all of it for one that lists only the genes it accepts.

    Examples
    --------
    >>> census = tf_gene_table("Homo sapiens")
    >>> census.columns[:4]
    ('gene_id_stem', 'symbol', 'is_tf', 'dbd_family')
    >>> len(census), len(census.assessed_positive)
    (2765, 1639)
    >>> census.frame().loc[0, "symbol"]
    'TFAP2A'
    """

    species: str
    provenance: CensusProvenance
    columns: tuple[str, ...]
    rows: tuple[tuple[str | None, ...], ...]

    def __len__(self) -> int:
        """Return how many genes the census assessed.

        Examples
        --------
        >>> len(tf_gene_table("Homo sapiens"))
        2765
        """
        return len(self.rows)

    @property
    def gene_id_stems(self) -> tuple[str, ...]:
        """Return every **Gene id stem** the census assessed, in row order.

        Examples
        --------
        >>> tf_gene_table("Homo sapiens").gene_id_stems[0]
        'ENSG00000137203'
        """
        # A stem is never blank — that is checked as the file is read — so ``or ""``
        # only narrows the type of a cell that is always text.
        return tuple(row[0] or "" for row in self.rows)

    @property
    def assessed_positive(self) -> tuple[str, ...]:
        """Return the stems this census judged transcription factors, in row order.

        The default a **TF gene list** is built from. Reading it as *the genes that
        are TFs* is the census speaking and never this package: a stem absent from
        :attr:`gene_id_stems` was not assessed at all, which is a different answer
        from one that is here and not in this tuple.

        Examples
        --------
        >>> len(tf_gene_table("Homo sapiens").assessed_positive)
        1639
        """
        return tuple(row[0] or "" for row in self.rows if row[_IS_TF] == TRUE_CELL)

    def frame(self) -> pd.DataFrame:
        """Return the census as a fresh :class:`~pandas.DataFrame`, one row per gene.

        Built for the caller each time it is asked for, so mutating it cannot reach
        the cached census behind it. Every column is text with ``None`` for a blank
        cell, except the uniform TF flag, which reads back as a boolean.

        Returns
        -------
        pandas.DataFrame
            The census's columns in file order, indexed from zero in row order.

        Examples
        --------
        >>> frame = tf_gene_table("Homo sapiens").frame()
        >>> list(frame.columns[:4])
        ['gene_id_stem', 'symbol', 'is_tf', 'dbd_family']
        >>> int(frame["is_tf"].sum())
        1639
        """
        frame = pd.DataFrame(list(self.rows), columns=pd.Index(self.columns), dtype=object)
        frame["is_tf"] = [row[_IS_TF] == TRUE_CELL for row in self.rows]
        return frame


#: Each provenance field's declared type, which parses that field's column of the table.
_PROVENANCE_TYPES: dict[str, Any] = get_type_hints(CensusProvenance)

#: The provenance table as a **Shipped table**: where it lives, what its header is, what it
#: is called and what repairs it. Every check the file is held to lives in
#: :mod:`genome.shipped`. One row per species, which is also the row the generator replaces
#: when it rebuilds one census and leaves the others alone (:mod:`genome.shipped_writer`).
CENSUS_METADATA_FORMAT = ShippedTable(
    resource=CENSUS_METADATA_RESOURCE,
    columns=_METADATA_COLUMNS,
    noun="census provenance table",
    repair=_REBUILD,
    error=TFGeneTableError,
    key=("species",),
    because=(
        "Every provenance column is required: a census nobody can cite is one this package "
        "may not redistribute"
    ),
    identify=("file", "species"),
)

#: One census as a **Shipped table**, one per species. The uniform four are the header's
#: leading columns and the publisher's own follow them; the **Gene id stem** is the key, so
#: a census naming one twice is refused as the two verdicts it would let a caller read.
#: ``scripts/build_tf_census.py`` writes a census through this same declaration, so the file
#: it produces is held to what this module will hold it to before it reaches disk.
CENSUS_FORMAT = ShippedTable(
    resource=f"{CENSUS_SUBDIR}/{{slug}}{CENSUS_SUFFIX}",
    columns=UNIFORM_COLUMNS,
    noun="census",
    repair=_REBUILD,
    error=TFGeneTableError,
    unit="gene",
    absence="this species has no transcription factors",
    leading=True,
    key=("gene_id_stem",),
    required=("gene_id_stem",),
    flags=("is_tf",),
    because=(
        "A census is keyed by gene id stem, so a row without one cannot be looked up or "
        "resolved against an annotation"
    ),
)


@cache
def census_metadata() -> tuple[CensusProvenance, ...]:
    """Return the provenance of every census the shipped table records, in table order.

    Where each **TF gene table** came from, who published it and what to cite. Read
    once and cached; the records are frozen, so the tuple is safe to hold on to.

    Returns
    -------
    tuple of CensusProvenance
        One record per row of the provenance table beside the censuses.

    Raises
    ------
    TFGeneTableError
        If a row cannot be read; the message names the column.

    Examples
    --------
    >>> {record.species for record in census_metadata()} >= {"Homo sapiens"}
    True
    """
    return _read_metadata(CENSUS_METADATA_FORMAT.text(), origin=CENSUS_METADATA_FORMAT.origin())


@cache
def census_species() -> tuple[str, ...]:
    """Return every species a **TF gene table** ships for, as its file name spells it.

    What can be asked about at all, and the answer an error names when a species
    cannot be. The directory is enumerated rather than any list of species being kept
    in code, so adding a species is dropping in a file. The names are slugs —
    ``homo_sapiens`` — and :func:`tf_gene_table` takes either those or the species as
    the assembly metadata table spells it.

    Returns
    -------
    tuple of str
        The species slugs, sorted. Empty only if the package ships no census.

    Examples
    --------
    >>> "homo_sapiens" in census_species()
    True
    """
    directory = files("genome").joinpath(CENSUS_SUBDIR)
    if not directory.is_dir():
        return ()
    return tuple(
        sorted(
            entry.name[: -len(CENSUS_SUFFIX)]
            for entry in directory.iterdir()
            if entry.name.endswith(CENSUS_SUFFIX)
        )
    )


@cache
def tf_gene_table(species: str) -> TFGeneTable | None:
    """Return the **TF gene table** shipped for ``species``, or ``None``.

    The raw absence, and the only place ``None`` is an acceptable way to say it: this
    is the layer below the one a caller touches, and everything above it turns the
    ``None`` into an error naming the species that do have a census, so that
    *nobody has published one for this species* can never be read as *this species
    has no transcription factors*.

    The species is slugged and then looked up among what :func:`census_species`
    found, rather than joined onto the resource directory, so a name shaped like a
    path finds nothing instead of walking out of it. The file is validated as it is
    read. Read once per species and cached; everything it returns is frozen.

    Parameters
    ----------
    species : str
        The species to look up, either as the assembly metadata table spells it
        (``"Homo sapiens"``) or as its slug (``"homo_sapiens"``).

    Returns
    -------
    TFGeneTable or None
        The census, or ``None`` when none ships for that species. ``None`` is legal
        and ordinary — most species have no published census.

    Raises
    ------
    TFGeneTableError
        If a census ships for that species and cannot be read, or ships with no row
        in the provenance table; the message names the file.

    Examples
    --------
    >>> tf_gene_table("Homo sapiens").provenance.family_column
    'DBD'
    >>> tf_gene_table("homo_sapiens") == tf_gene_table("Homo sapiens")
    True
    >>> tf_gene_table("Danio rerio") is None
    True
    """
    slug = species_slug(species)
    if slug not in census_species():
        return None
    origin = CENSUS_FORMAT.origin(slug=slug)
    provenance = next(
        (record for record in census_metadata() if species_slug(record.species) == slug), None
    )
    if provenance is None:
        raise TFGeneTableError(
            f"{origin} ships with no row in {CENSUS_METADATA_RESOURCE}, so nothing says who "
            f"published it or what to cite for it — and citing the publisher is the condition "
            f"on redistributing a census here. Re-run scripts/build_tf_census.py for "
            f"{slug!r}, which writes the census and its provenance row together."
        )
    # These are hundreds of kilobytes of shipped rows, so the seam is between the bytes and
    # the format: unpacking happens at the resource boundary and the parse below is a pure
    # function of text.
    return _read_census(CENSUS_FORMAT.text(slug=slug), provenance=provenance, origin=origin)


def _read_metadata(text: str, *, origin: str) -> tuple[CensusProvenance, ...]:
    """Read the provenance table from ``text``, holding it to what the shared reader checks."""
    return tuple(
        CensusProvenance.from_row(row, origin=origin)
        for row in CENSUS_METADATA_FORMAT.parse(text, origin=origin).mappings()
    )


def _read_census(text: str, *, provenance: CensusProvenance, origin: str) -> TFGeneTable:
    """Read one census from ``text``, holding it to what a shipped census promises.

    Separate from the resource it came out of, so every way a file can be wrong is
    reachable without writing a broken one into the package. ``origin`` is where the
    text came from and is named in every message, since regenerating that file is the
    only repair.
    """
    read = CENSUS_FORMAT.parse(text, origin=origin)
    return TFGeneTable(
        species=provenance.species, provenance=provenance, columns=read.columns, rows=read.rows
    )
