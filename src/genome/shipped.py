r"""Reading a **Shipped table** — the one loader every curated table in the wheel goes through.

A **Shipped table** is a small tab-separated file that ships inside the package: the curated
**Assembly metadata** and **Annotation metadata** tables, the **Xref source** table, the
provenance tables beside each shipped-data directory, the censuses, the **Cofactor table**s and
the **Motif link** tables. Six modules read one, and each of them declares what its own table is
— its resource path, its columns, the row type they become and the command that repairs a broken
one — while everything else happens here: resource lookup, gzip, header validation, cell parsing,
the blank-cell rules, duplicate-key refusal, and the shape of the error a broken file raises.

**Six failures recur, and they are checked here so that they reach every table**: an empty file, a
header that is not the declared columns, a row with the wrong cell count, a blank cell in a
required column, a flag spelled a way no table spells one, and a repeated key. Only the last is
declared per table rather than assumed — a **Motif link** table carries many rows per **Gene id
stem** by design, so it names no key and the refusal does not run.

**The messages are the tables' own.** A shipped file that cannot be read is a packaging defect and
not anything a caller did, so the only useful message names the file and the repair. Each table
supplies both halves — the noun it is called by (:attr:`ShippedTable.noun`) and the command that
rebuilds it (:attr:`ShippedTable.repair`) — and the templates here compose them, so ``re-run
scripts/build_tf_cofactor.py for that species`` reads the same whichever of the six failures a
file hit.

**The rule against reading a whole genomic file into memory does not reach here.** The largest
shipped table is the human census at 495 KB unpacked and 52 KB on disk, and a table is read whole
because validating it as it is read is the point: a file that cannot be trusted must never answer,
and half a file cannot be held to a duplicate-key rule. The rule is about the FASTA, the GTF and
the aligner inputs, which are three to six orders of magnitude larger and are streamed everywhere
they are touched.

**The shipped gene lists are not Shipped tables**, which is why :mod:`genome.gene_list` does not
read through this module. Those are JSON keyed by **Registered name**, nested rather than tabular:
no header, no columns, no rows, and a duplicate rule about a gene id appearing in two categories
rather than about a repeated key. The only thing they share with a table here is being found by
enumerating a directory on a file-name suffix, and that is not the reader.

Examples
--------
>>> from genome.shipped import ShippedTable
>>> tiny = ShippedTable(resource="", columns=("id", "name"), noun="tiny table")
>>> tiny.parse("id\tname\nx\tA\n", origin="tiny.tsv").rows
(('x', 'A'),)
>>> try:
...     tiny.parse("id\n", origin="tiny.tsv")
... except ShippedTableError as error:
...     print("tiny.tsv" in str(error))
True
"""

from __future__ import annotations

import gzip
import math
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.resources import files
from typing import Any, get_args

#: How every **Shipped table** spells a flag: one pair of spellings across the censuses, the
#: **Cofactor table**s and the **Motif link** tables, so a third spelling is a defect in what
#: wrote the file rather than a publisher's own convention. The curated metadata tables are
#: looser — see :func:`parse_cell`, where a blank flag cell is the real answer *no*.
TRUE_CELL, FALSE_CELL = "yes", "no"

#: How one cell spells more than one value — a **Motif link**'s partners and tax ids, a
#: **Cofactor table**'s functions and complexes, a census's Interpro ids. One separator across
#: every table this package ships, so a caller never has to remember which column uses which,
#: and a semicolon rather than a comma or a tab: the files carry no quoting and are read by
#: splitting on the tab. A published value that already contains it is refused by whatever
#: writes the cell, since writing it unchanged would make one value read back as two.
VALUE_SEPARATOR = ";"

#: Cell spellings a curated table's boolean column accepts, lower-cased. Anything else in one is
#: a typo in a hand-maintained table and says so rather than reading as ``False``.
_TRUE_CELLS = frozenset({"yes", "true", "1"})
_FALSE_CELLS = frozenset({"no", "false", "0"})

#: Why a column may not be left blank, where a table gives no reason of its own: the one
#: that is true of every table here, since blank is how each of them spells *unknown*.
_UNKNOWN = (
    "A blank cell reads back as unknown, and a column whose type declares no unknown has no "
    "way to say so"
)


class ShippedTableError(ValueError):
    """A **Shipped table** cannot be read, so it is not allowed to answer.

    The base of every table's own error class, so a caller may catch one kind of defect across
    all of them and still tell a census from a **Cofactor table** by the class that was raised.
    A :class:`ValueError`, because a file that says something the format does not is a bad value
    rather than a broken program — and never a :class:`LookupError`, which is what absence is
    spelled with, so catching an absent table cannot swallow a broken one.

    Examples
    --------
    >>> from genome.tf.gene import TFGeneTableError
    >>> issubclass(TFGeneTableError, ShippedTableError)
    True
    """


class MetadataRowError(ShippedTableError):
    """A row cannot be read as a record, and the message names the column that refused.

    Raised by :func:`parse_cell` for the curated tables that declare their columns as a frozen
    dataclass — the **Assembly metadata** and **Annotation metadata** tables and the **Xref
    source** table — for a cell no column's type can read and for a blank cell in a column that
    has no unknown. Re-exported from :mod:`genome.metadata`, which is where those tables live.

    Examples
    --------
    >>> from genome.metadata import AssemblyMetadata
    >>> try:
    ...     AssemblyMetadata.from_row({"assembly_name": "tiny", "ncbi_taxid": "many"})
    ... except MetadataRowError as error:
    ...     print("ncbi_taxid" in str(error))
    True
    """


def species_slug(species: str) -> str:
    """Return the file-name spelling of ``species``.

    Lower case, with each run of anything that is not a letter or a digit collapsed to one
    underscore. It is the one place a curated table's spelling of a species and a file name are
    reconciled, so neither has to be written the other's way — and it lives with the reader
    because naming a file after a species is what every shipped-data directory in this package
    does and none of them owns.

    Parameters
    ----------
    species : str
        A species name, in any spelling.

    Returns
    -------
    str
        Its slug.

    Examples
    --------
    >>> species_slug("Homo sapiens")
    'homo_sapiens'
    >>> species_slug("Escherichia coli HT115")
    'escherichia_coli_ht115'
    >>> species_slug("homo_sapiens")
    'homo_sapiens'
    """
    kept = "".join(character if character.isalnum() else " " for character in species.lower())
    return "_".join(kept.split())


@dataclass(frozen=True)
class ShippedRows:
    r"""One **Shipped table** read and validated: its header, its rows and where it came from.

    What :meth:`ShippedTable.parse` hands back, and the whole of what the shared reader knows
    about a table — turning these into a record is the declaring module's own business, since
    only it knows what a row means.

    Attributes
    ----------
    columns : tuple of str
        The header as the file spells it. For a table whose columns are merely the *leading*
        ones this carries the publisher's own columns after them, in file order.
    rows : tuple of tuple of (str or None)
        One tuple per row, parallel to :attr:`columns`. A blank cell is ``None``, which is the
        reading every table in this package gives one.
    origin : str
        Where the text came from, as every message about it names it.

    Examples
    --------
    >>> rows = ShippedTable(resource="", columns=("id", "name"), noun="tiny table").parse(
    ...     "id\tname\nx\t\n", origin="tiny.tsv"
    ... )
    >>> rows.rows
    (('x', None),)
    >>> rows.mappings()
    ({'id': 'x', 'name': ''},)
    """

    columns: tuple[str, ...]
    rows: tuple[tuple[str | None, ...], ...]
    origin: str

    def mappings(self) -> tuple[dict[str, str], ...]:
        """Return each row keyed by column name, with a blank cell as the empty string.

        The shape a record's ``from_row`` takes, and the one :func:`parse_cell` reads.

        Returns
        -------
        tuple of dict of str to str
            One mapping per row, in file order.

        Examples
        --------
        >>> ShippedRows(("id",), (("x",),), "tiny.tsv").mappings()
        ({'id': 'x'},)
        """
        return tuple(
            {column: cell or "" for column, cell in zip(self.columns, row, strict=True)}
            for row in self.rows
        )


@dataclass(frozen=True)
class ShippedTable:
    r"""What one **Shipped table** declares about itself, which is all its module has to say.

    A declaration and never a reader: the module that owns a table states where it lives, what
    its columns are, what it is called, what repairs it and which of the optional rules apply,
    and every check runs here. Frozen, so a declaration is a module-level constant.

    Attributes
    ----------
    resource : str
        The table's path inside the package, as a template: ``"data/tf_gene/{slug}.tf_gene_table
        .tsv.gz"`` for a table there is one of per species, and a plain path for a table there is
        one of. Empty for a table only ever read from text, which is how the malformed cases are
        reached without writing a broken file into the package. A name ending ``.gz`` is unpacked
        on the way in.
    columns : tuple of str
        The header the table carries, in order — the whole of it, or, under ``leading``, the
        columns it must begin with.
    noun : str
        What the table is called, as its own messages call it: ``"census"``, ``"cofactor
        table"``, ``"link table"``.
    repair : str
        The imperative that fixes a broken one, lower-cased and without its full stop —
        ``"re-run scripts/build_tf_cofactor.py for that species"``. Empty for a hand-maintained
        table with no generator, which drops the clause rather than inventing a command.
    error : type of ShippedTableError
        The class every refusal about this table is raised as.
    unit : str
        What one row is: ``"gene"``, ``"link"``, ``"row"``.
    absence : str
        The reading an empty file would wrongly invite — ``"this species has no transcription
        factors"``. Setting it also declares that the table carries at least one row, since that
        reading is the whole reason a file with none is refused.
    leading : bool
        Whether ``columns`` are the leading columns rather than the whole header. A census and a
        **Cofactor table** carry the publisher's own columns after the uniform ones, which are
        held only to being named once each.
    key : tuple of str
        The columns one row's identity is, which no two rows may repeat. Empty opts out, which
        is what a **Motif link** table does: it carries many rows per **Gene id stem** by design.
    required : tuple of str
        Columns no row may leave blank. Every column of a table whose record declares its fields
        is required by that record's own types instead — see :func:`parse_cell`.
    flags : tuple of str
        Columns spelled :data:`TRUE_CELL` or :data:`FALSE_CELL`, and nothing else.
    because : str
        Why a required column may not be blank, as a sentence without its full stop. The one
        piece of prose a shared refusal cannot compose, because it is the table's own reason;
        the default is the reason true of every table here.
    identify : tuple of str
        The columns that name a row in a message, best first. Only the non-blank ones are
        printed, so a row that is malformed *because* one of them is blank still says which row.

    Examples
    --------
    >>> table = ShippedTable(
    ...     resource="",
    ...     columns=("gene_id_stem", "symbol"),
    ...     noun="tiny table",
    ...     repair="re-run scripts/build_tiny.py",
    ...     key=("gene_id_stem",),
    ... )
    >>> table.parse("gene_id_stem\tsymbol\ng1\tA\n", origin="tiny.tsv").rows
    (('g1', 'A'),)
    >>> try:
    ...     table.parse("gene_id_stem\tsymbol\ng1\tA\ng1\tB\n", origin="tiny.tsv")
    ... except ShippedTableError as error:
    ...     print("re-run scripts/build_tiny.py" in str(error))
    True
    """

    resource: str
    columns: tuple[str, ...]
    noun: str
    repair: str = ""
    error: type[ShippedTableError] = ShippedTableError
    unit: str = "row"
    absence: str = ""
    leading: bool = False
    key: tuple[str, ...] = ()
    required: tuple[str, ...] = ()
    flags: tuple[str, ...] = ()
    because: str = _UNKNOWN
    identify: tuple[str, ...] = ()

    # -- where the bytes are ------------------------------------------------------------

    def origin(self, **keys: str) -> str:
        """Return the shipped file's path, which every message about it names.

        Parameters
        ----------
        **keys : str
            What :attr:`resource`'s template names — ``slug=``, ``release=``.

        Returns
        -------
        str
            The resource's path as the installed package spells it.

        Examples
        --------
        >>> table = ShippedTable(
        ...     resource="data/tf_gene/{slug}.tf_gene_table.tsv.gz",
        ...     columns=("gene_id_stem",),
        ...     noun="census",
        ... )
        >>> table.origin(slug="homo_sapiens").endswith("homo_sapiens.tf_gene_table.tsv.gz")
        True
        """
        return str(files("genome").joinpath(self.resource.format(**keys)))

    def text(self, **keys: str) -> str:
        r"""Return the shipped file's text, unpacking it when its name says gzip.

        The one place a shipped resource is opened. Bulk ships gzipped and small metadata ships
        plain — the convention every shipped-data directory here follows — so the suffix is what
        decides, and no caller writes the decompression again.

        Parameters
        ----------
        **keys : str
            What :attr:`resource`'s template names.

        Returns
        -------
        str
            The whole file, decoded as UTF-8.

        Examples
        --------
        >>> table = ShippedTable(
        ...     resource="data/tf_gene/{slug}.tf_gene_table.tsv.gz",
        ...     columns=("gene_id_stem",),
        ...     noun="census",
        ... )
        >>> table.text(slug="homo_sapiens").startswith("gene_id_stem\t")
        True
        """
        path = self.resource.format(**keys)
        resource = files("genome").joinpath(path)
        if path.endswith(".gz"):
            return gzip.decompress(resource.read_bytes()).decode("utf-8")
        return resource.read_text(encoding="utf-8")

    def read(self, **keys: str) -> ShippedRows:
        """Return the shipped file, read and held to everything this table declares.

        Parameters
        ----------
        **keys : str
            What :attr:`resource`'s template names.

        Returns
        -------
        ShippedRows
            The header, the rows and the file they came from.

        Raises
        ------
        ShippedTableError
            Of :attr:`error`'s class, for any of the failures :meth:`parse` checks.

        Examples
        --------
        >>> table = ShippedTable(
        ...     resource="data/tf_cofactor/cofactor_metadata.tsv",
        ...     columns=("species", "ncbi_taxid", "file", "sha256"),
        ...     noun="provenance table",
        ... )
        >>> table.read().columns[0]
        'species'
        """
        return self.parse(self.text(**keys), origin=self.origin(**keys))

    # -- what a table is held to --------------------------------------------------------

    def parse(self, text: str, *, origin: str) -> ShippedRows:
        r"""Read ``text`` as this table, refusing it the moment it is not one.

        Separate from the resource it came out of, so every way a shipped file can be wrong is
        reachable without writing a broken one into the package.

        Parameters
        ----------
        text : str
            The whole table, tab separated, header first.
        origin : str
            Where the text came from, named in every message.

        Returns
        -------
        ShippedRows
            The header and the rows, blank cells as ``None``.

        Raises
        ------
        ShippedTableError
            Of :attr:`error`'s class: for an empty file, a header that is not
            :attr:`columns`, a column named twice, a row with the wrong cell count, a blank cell
            in a :attr:`required` column, a :attr:`flags` cell spelled a way no table spells one,
            a repeated :attr:`key`, or — for a table declaring an :attr:`absence` — a file that
            declares columns and no rows.

        Examples
        --------
        >>> table = ShippedTable(resource="", columns=("id",), noun="tiny table")
        >>> table.parse("id\nx\n", origin="tiny.tsv").rows
        (('x',),)
        """
        lines = text.split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        if not lines:
            raise self.error(
                (
                    f"{origin} is empty. {self._an.capitalize()} {self.noun} carries a header and at "
                    f"least one {self.unit}; absence is spelled by shipping no file at all, and an empty "
                    f"one is the second spelling of it that would read as *{self.absence}*."
                    f"{self._then(', or remove the file')}"
                )
                if self.absence
                else (
                    f"{origin} is empty, and {self._an} {self.noun} that says nothing is one nothing "
                    f"can be read from.{self._then('')}"
                )
            )
        columns = tuple(lines[0].split("\t"))
        self._check_header(columns, origin=origin)
        rows = tuple(
            self._read_row(line, columns, number, origin=origin)
            for number, line in enumerate(lines[1:], start=2)
        )
        if self.absence and not rows:
            raise self.error(
                f"{origin} declares columns and no {self.unit}s. {self._an.capitalize()} {self.noun} "
                f"carrying nothing says no more than an absent file does."
                f"{self._then(', or remove the file')}"
            )
        self._check_key(columns, rows, origin=origin)
        return ShippedRows(columns=columns, rows=rows, origin=origin)

    def record(
        self, row: Mapping[str, object], types: Mapping[str, Any], *, origin: str
    ) -> dict[str, Any]:
        """Parse one row into a record's fields, each cell by its own field's declared type.

        The reader half of a table whose record declares its columns as a frozen dataclass: the
        fields are the columns, and their types are what parses them.

        Parameters
        ----------
        row : mapping of str to object
            Column name to cell, as :meth:`ShippedRows.mappings` spells one.
        types : mapping of str to object
            Each field's declared type, normally :func:`typing.get_type_hints` of the record.
        origin : str
            Where the row came from, named in every message.

        Returns
        -------
        dict of str to object
            Field name to parsed value, ready to build the record with.

        Raises
        ------
        ShippedTableError
            Of :attr:`error`'s class, as :func:`parse_cell` raises it.

        Examples
        --------
        >>> table = ShippedTable(resource="", columns=("pubmed_id",), noun="provenance table")
        >>> table.record({"pubmed_id": "26896847"}, {"pubmed_id": int}, origin="test")
        {'pubmed_id': 26896847}
        """
        return {name: parse_cell(name, row, types, table=self, origin=origin) for name in types}

    # -- the six failures ---------------------------------------------------------------

    def _check_header(self, columns: tuple[str, ...], *, origin: str) -> None:
        """Hold a header to the columns declared, whole or leading, and to distinct names."""
        if not self.leading:
            if columns != self.columns:
                raise self.error(
                    f"{origin} carries the columns {list(columns)} where every {self.noun} "
                    f"carries {list(self.columns)}. The header is the format and not a hint: "
                    f"fix it, keeping the columns in that order{self.repair_clause}"
                )
            return
        wide = len(self.columns)
        if columns[:wide] != self.columns:
            raise self.error(
                f"{origin} leads with the columns {list(columns[:wide])} where every "
                f"{self.noun} leads with {list(self.columns)}. Those {wide} are the only "
                f"columns one {self.noun} shares with another, so a file without them cannot be "
                f"read as one{self.repair_clause}"
            )
        if len(set(columns)) != len(columns):
            raise self.error(
                f"{origin} names a column twice: {list(columns)}. Each column is named once, so "
                f"two of one name would let a reader take either{self.repair_clause}"
            )

    def _read_row(
        self, line: str, columns: tuple[str, ...], number: int, *, origin: str
    ) -> tuple[str | None, ...]:
        """Read one row, blank cells becoming ``None`` and every declared rule holding."""
        cells = line.split("\t")
        if len(cells) != len(columns):
            raise self.error(
                f"{origin} line {number} holds {len(cells)} cells where the header declares "
                f"{len(columns)}. {self._an.capitalize()} {self.noun} is a plain TSV with no "
                f"quoting, so a cell carrying a tab is a defect in what wrote it rather than "
                f"something to parse around{self.repair_clause}"
            )
        row = dict(zip(columns, cells, strict=True))
        for column in self.flags:
            if row[column] not in (TRUE_CELL, FALSE_CELL):
                raise self.error(
                    f"{origin} line {number} spells its {column!r} cell {row[column]!r}, and "
                    f"{self._an} {self.noun} spells a flag {TRUE_CELL!r} or {FALSE_CELL!r}. One "
                    f"spelling of *yes* across every table this package ships, so a third is a "
                    f"defect in what wrote the file and not a new kind of answer{self.repair_clause}"
                )
        for column in self.required:
            if not row[column]:
                raise self.error(
                    f"{origin} line {number} leaves the {column!r} column blank."
                    f"{self._reason}{self.repair_clause}"
                )
        return tuple(cell if cell else None for cell in cells)

    def _check_key(
        self, columns: tuple[str, ...], rows: tuple[tuple[str | None, ...], ...], *, origin: str
    ) -> None:
        """Hold every row's key to naming its subject once — unless the table declares none."""
        if not self.key:
            return
        places = [columns.index(column) for column in self.key]
        words = " and ".join(column.replace("_", " ") for column in self.key)
        seen: set[tuple[str | None, ...]] = set()
        for row in rows:
            identity = tuple(row[place] for place in places)
            if identity in seen:
                spelled = identity[0] if len(identity) == 1 else identity
                raise self.error(
                    f"{origin} names the {words} {spelled!r} more than once. "
                    f"{self._an.capitalize()} {self.noun} carries one row per {words}, so two "
                    f"rows for one would let a caller read either{self.repair_clause}"
                )
            seen.add(identity)

    # -- the two halves of a message, composed ------------------------------------------

    @property
    def _reason(self) -> str:
        """The table's own reason a required column may not be blank, as its own sentence."""
        return f" {self.because}." if self.because else ""

    @property
    def _an(self) -> str:
        """``a`` or ``an``, whichever this table's noun takes."""
        return "an" if self.noun[:1].lower() in "aeiou" else "a"

    @property
    def repair_clause(self) -> str:
        """Return the repair as the trailing clause a refusal about this table closes with.

        The half of a message a table declares and a shared refusal composes, so that
        whatever writes a table can close its own refusals in the same words the reader
        closes its own with. A table with no generator has no command to name and closes
        the sentence with a full stop instead.

        Returns
        -------
        str
            ``" — re-run scripts/build_tf_census.py."``, or ``"."``.

        Examples
        --------
        >>> ShippedTable(resource="", columns=("id",), noun="tiny table").repair_clause
        '.'
        >>> ShippedTable(
        ...     resource="", columns=("id",), noun="tiny table", repair="re-run the build"
        ... ).repair_clause
        ' — re-run the build.'
        """
        return f" — {self.repair}." if self.repair else "."

    def _then(self, tail: str) -> str:
        """Return the repair as a sentence of its own, ``tail`` before its full stop."""
        if not self.repair:
            return ""
        return f" {self.repair[0].upper()}{self.repair[1:]}{tail}."

    def _names(self, row: Mapping[str, object]) -> str:
        """Return the clause saying which row a message is about, or nothing when none says."""
        found = [text for column in self.identify if (text := str(row.get(column) or "").strip())]
        return f" for {'/'.join(repr(text) for text in found)}" if found else ""


def parse_cell(
    name: str,
    row: Mapping[str, object],
    types: Mapping[str, Any],
    *,
    table: ShippedTable | None = None,
    origin: str | None = None,
) -> Any:
    """Parse the ``name`` column of ``row`` with that field's own declared type.

    How every **Shipped table** whose record declares its columns as a frozen dataclass turns a
    cell into a field: the fields are the columns, and their types are the parser. The curated
    **Assembly metadata**, **Annotation metadata** and **Xref source** tables read this way, and
    so do the provenance tables beside the shipped data.

    A field declared optional (``T | None``) is parsed by ``T`` when its cell carries text, and
    is ``None`` when the cell is blank or its column is absent — a union is not callable, so the
    type inside it does the parsing. A field declared ``bool`` is a flag column, where an empty
    cell is a real answer *no* rather than an unknown. Any other field is parsed by its declared
    type and has no unknown, so a blank cell there is a malformed row.

    Parameters
    ----------
    name : str
        The column to read, which is also the field it fills.
    row : mapping of str to object
        Column name to cell, as the shipped TSV spells one. An absent column is blank, and so is
        the NaN a caller reading the table with pandas gets for one.
    types : mapping of str to object
        Each field's declared type, normally :func:`typing.get_type_hints` of the record.
    table : ShippedTable, optional
        The table the row came out of, when there is one. It supplies the error class, the noun
        and the repair, so a defect in a shipped file names the command that fixes it. Omitted
        for a row a caller built by hand, where none of those is true and the message says only
        which column refused.
    origin : str, optional
        Where the row came from, named in every message when there is one to name.

    Returns
    -------
    object
        The cell parsed by its column's own type, or ``None`` for a blank cell in a column that
        has an unknown.

    Raises
    ------
    ShippedTableError
        Of ``table``'s class, or :class:`MetadataRowError` with no table: if the cell cannot be
        read as its column's type, or the column has no unknown and the cell is blank. The
        message names the column.

    Examples
    --------
    >>> parse_cell("ncbi_taxid", {"ncbi_taxid": "9606"}, {"ncbi_taxid": int})
    9606
    >>> parse_cell("species", {}, {"species": str | None}) is None
    True
    """
    declared = types[name]
    if declared is bool:
        return _parse_flag(name, row, table=table, origin=origin)
    inside = [arg for arg in get_args(declared) if arg is not type(None)]
    text = _cell_text(name, row)
    if inside:
        return _parse_text(name, text, inside[0], table=table, origin=origin) if text else None
    if not text:
        raise _blank(name, row, table=table, origin=origin)
    return _parse_text(name, text, declared, table=table, origin=origin)


def _cell_text(name: str, row: Mapping[str, object]) -> str:
    """Return the ``name`` cell of ``row`` as stripped text, empty when it is blank.

    Blank is a missing column, ``None``, the NaN pandas reads a blank cell as, or whitespace.
    Anything else is that cell's own text, so a row that spells a value in its column's own
    type — an ``int`` taxid, a ``bool`` flag — reads back the same as the table's text for it.
    """
    cell = row.get(name)
    if cell is None or (isinstance(cell, float) and math.isnan(cell)):
        return ""
    return str(cell).strip()


def _parse_flag(
    name: str,
    row: Mapping[str, object],
    *,
    table: ShippedTable | None,
    origin: str | None,
) -> bool:
    """Parse a curated table's boolean column — blank is ``False``, a spelling nobody uses raises."""
    text = _cell_text(name, row).lower()
    if not text or text in _FALSE_CELLS:
        return False
    if text in _TRUE_CELLS:
        return True
    accepted = ", ".join(sorted(_TRUE_CELLS | _FALSE_CELLS))
    raise _raised(table)(
        f"{_says(name, text, origin)}, which is not a flag. Fix that cell to one of: {accepted} "
        f"— or leave it blank, which reads as false."
    )


def _parse_text(
    name: str,
    text: str,
    declared: Any,
    *,
    table: ShippedTable | None,
    origin: str | None,
) -> Any:
    """Parse one cell's text with its column's type, or say which cell refused."""
    try:
        return declared(text)
    except (TypeError, ValueError) as error:
        raise _raised(table)(
            f"{_says(name, text, origin)}, which {declared.__name__} cannot read. Fix that cell "
            f"to a value {declared.__name__} accepts{'.' if table is None else table.repair_clause}"
        ) from error


def _blank(
    name: str,
    row: Mapping[str, object],
    *,
    table: ShippedTable | None,
    origin: str | None,
) -> ShippedTableError:
    """Return the refusal a blank cell in a column with no unknown raises."""
    if table is None:
        return MetadataRowError(
            f"the {name!r} column is blank, and it is one no row may leave blank. Fill that "
            f"cell in: a blank cell reads back as unknown, and {name!r} has none."
        )
    where = f"{origin} leaves" if origin else "the table leaves"
    return table.error(
        f"{where} the {name!r} column blank{table._names(row)}.{table._reason} Fill that cell "
        f"in{table.repair_clause}"
    )


def _says(name: str, text: str, origin: str | None) -> str:
    """Return how a cell-level message opens — the file where there is one, the cell always."""
    if origin:
        return f"{origin} holds {text!r} in the {name!r} column"
    return f"the {name!r} column holds {text!r}"


def _raised(table: ShippedTable | None) -> type[ShippedTableError]:
    """Return the class a refusal about ``table`` is raised as."""
    return MetadataRowError if table is None else table.error
