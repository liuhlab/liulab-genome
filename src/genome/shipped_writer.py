r"""Writing a **Shipped table** — the one writer the generators in ``scripts/`` go through.

A **Shipped table** is written by exactly one kind of thing: a build script run by hand
against a publisher's download, never at install time and never in CI. Three of them write
one today — ``build_tf_census.py``, ``build_tf_cofactor.py`` and ``build_tf_links.py`` — and
each used to carry the same renderer, the same gzip call and the same provenance merge under
the same names. They are one implementation here, and what each script keeps is its
publisher's recipe: which file, which columns, which cleaning, which checks are that
publisher's.

**The writer and the reader are handed the same declaration.** A file is held to its own
:class:`~genome.shipped.ShippedTable` *before it reaches disk* — the header it must carry, the
columns no row may leave blank, the flag spellings, the key no two rows may repeat — so the
two halves of a format agree because they read one object, not because two people kept two
copies of the rules in step. A generator supplies its own error class and its own repair, by
:func:`dataclasses.replace` on that declaration, so a refusal names the file that would have
been wrong and the recipe to fix rather than telling whoever ran the build to re-run it.

**Two runs write the same bytes.** gzip is given ``compresslevel=9`` and ``mtime=0``, so
nothing about *when* a build ran reaches the file: rebuilding from an unchanged input produces
no diff at all, which is what makes a shipped ``.tsv.gz`` reviewable and its pinned digest
worth pinning. The digest is over the **unpacked** bytes (ADR-0006), so a copy recompressed
elsewhere still matches. Row order is not this module's business: a generator hands rows in
whatever deterministic order its publisher's file or its own sort gives them, and this writes
them down in that order.

**Bulk ships gzipped and small metadata plain**, and the file name is what decides — the same
rule :meth:`~genome.shipped.ShippedTable.text` unpacks by, read off the same declaration, so
the two halves cannot pick differently. The name itself comes from the declaration's resource
template too (:func:`shipped_name`), which is why no script spells a file suffix of its own.

**It sits beside the reader rather than under** ``io/``. A shipped table is a package resource
and not anything under the **Data dir**, so the two halves of one format belong in one place;
and the gate is the deciding reason — ``scripts/`` is linted but neither typechecked nor
tested, so a writer left there would have been the one module in this package that nothing
holds to anything.

Examples
--------
>>> from pathlib import Path
>>> from tempfile import TemporaryDirectory
>>> from genome.shipped import ShippedTable
>>> tiny = ShippedTable(
...     resource="data/tiny/{slug}.tiny_table.tsv.gz",
...     columns=("gene_id_stem", "symbol"),
...     noun="tiny table",
...     key=("gene_id_stem",),
... )
>>> with TemporaryDirectory() as directory:
...     written = write_table(tiny, Path(directory), tiny.columns, [("g1", "A")], slug="worm")
...     (written.path.name, written.unpacked)
('worm.tiny_table.tsv.gz', b'gene_id_stem\tsymbol\ng1\tA\n')
"""

from __future__ import annotations

import gzip
import hashlib
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from genome.shipped import ShippedTable

__all__ = ["WrittenTable", "merge_rows", "render", "shipped_name", "write_table"]

#: What gzip is asked for, and the whole of why: level 9 because these files are committed and
#: read far more often than they are built, and ``mtime=0`` because the alternative is a header
#: that changes on every rebuild. Together they are the promise two runs of a generator agree
#: byte for byte on — stated here and nowhere else, so a fourth generator inherits it rather
#: than having to notice it.
_COMPRESSION = 9

#: The file-name suffix that says a table's bytes are gzipped, as
#: :meth:`~genome.shipped.ShippedTable.text` unpacks by.
_PACKED_SUFFIX = ".gz"

#: What no cell of a plain TSV may contain, and what a refusal calls each. The file carries no
#: quoting and is read by splitting on the tab and the newline, so a cell holding one is a
#: defect in what built it rather than something to escape around.
_UNWRITABLE = (("\t", "a tab"), ("\n", "a newline"), ("\r", "a carriage return"))


@dataclass(frozen=True)
class WrittenTable:
    r"""One **Shipped table** as it was just written: the file, and the bytes inside it.

    What :func:`write_table` and :func:`merge_rows` hand back. The unpacked bytes travel with
    the path because they are what a provenance row's digest is over (ADR-0006) and what a
    generator reports its size from, and reading them back off a gzipped file to say so would
    be doing the work twice.

    Attributes
    ----------
    path : pathlib.Path
        The file written.
    unpacked : bytes
        The TSV as it reads — before compression, whether or not the file is compressed.

    Examples
    --------
    >>> WrittenTable(Path("tiny.tsv"), b"id\nx\n").sha256[:8]
    '346af203'
    """

    path: Path
    unpacked: bytes

    @property
    def sha256(self) -> str:
        r"""Return the digest of the unpacked bytes, which is the one a provenance row pins.

        Over the TSV and never the gzip around it (ADR-0006), so a copy recompressed
        elsewhere still matches what the table says about itself.

        Returns
        -------
        str
            The SHA-256 digest, hex.

        Examples
        --------
        >>> len(WrittenTable(Path("tiny.tsv"), b"id\nx\n").sha256)
        64
        """
        return hashlib.sha256(self.unpacked).hexdigest()

    @property
    def packed(self) -> int:
        """Return the size of the file on disk, in bytes.

        Returns
        -------
        int
            What the file occupies — the compressed size for a gzipped table.

        Examples
        --------
        >>> from tempfile import TemporaryDirectory
        >>> from genome.shipped import ShippedTable
        >>> tiny = ShippedTable(resource="tiny.tsv", columns=("id",), noun="tiny table")
        >>> with TemporaryDirectory() as directory:
        ...     write_table(tiny, Path(directory), tiny.columns, [("x",)]).packed
        5
        """
        return self.path.stat().st_size


def shipped_name(table: ShippedTable, **keys: str) -> str:
    """Return what one of ``table``'s files is called, from the declaration's own template.

    The file name is declared once — in the resource path the reader looks the table up by —
    so a generator writing a table never spells a suffix, a release prefix or a species slug
    into a name of its own.

    Parameters
    ----------
    table : genome.shipped.ShippedTable
        The declaration of the format being written.
    **keys : str
        What the declaration's resource template names — ``slug=``, ``release=``.

    Returns
    -------
    str
        The bare file name, with no directory in front of it.

    Examples
    --------
    >>> from genome.tf.gene import CENSUS_FORMAT
    >>> shipped_name(CENSUS_FORMAT, slug="homo_sapiens")
    'homo_sapiens.tf_gene_table.tsv.gz'
    """
    return Path(table.resource.format(**keys)).name


def render(columns: Sequence[str], rows: Iterable[Sequence[str]]) -> bytes:
    r"""Return one table as the bytes it ships as: the header, then one line per row.

    Tab separated, UTF-8, ``\n`` throughout and **no quoting** — every cell as it stands, so
    what a reader splits on is what was written. A pure function of what it is handed: the
    same columns and rows render to the same byte string every time, which is the first half
    of two runs of a generator agreeing.

    Parameters
    ----------
    columns : sequence of str
        The header, in table order.
    rows : iterable of sequence of str
        One sequence of cells per row, parallel to ``columns``, in the order they ship in.

    Returns
    -------
    bytes
        The whole table, ending in a newline.

    Examples
    --------
    >>> render(("id", "name"), [("g1", "A"), ("g2", "B")])
    b'id\tname\ng1\tA\ng2\tB\n'
    """
    lines = ["\t".join(columns)]
    lines.extend("\t".join(row) for row in rows)
    return ("\n".join(lines) + "\n").encode("utf-8")


def write_table(
    table: ShippedTable,
    directory: Path,
    columns: Sequence[str],
    rows: Iterable[Sequence[str]],
    **keys: str,
) -> WrittenTable:
    r"""Write one **Shipped table** into ``directory``, refusing it unless it reads back.

    The whole of writing one: the cells are held to what a plain TSV can carry, the table is
    rendered, the rendering is parsed back through ``table``'s own declaration — so the file
    is held to everything :mod:`genome.shipped` will hold it to, before it exists — and only
    then is it compressed and written. Gzipped when the declared name says so and plain
    otherwise.

    Parameters
    ----------
    table : genome.shipped.ShippedTable
        The declaration of the format being written, carrying the error class and the repair
        its refusals are raised with. A generator supplies its own of both with
        :func:`dataclasses.replace`, since what repairs a file being built is the recipe that
        built it and not the command that runs it.
    directory : pathlib.Path
        Where the file goes. Created if it is not there.
    columns : sequence of str
        The header, in table order. The declared columns for a table whose header is fixed,
        and those followed by the publisher's own for a table whose leading columns are what
        is uniform.
    rows : iterable of sequence of str
        One sequence of cells per row, in the order they ship in.
    **keys : str
        What the declaration's resource template names — ``slug=``, ``release=``.

    Returns
    -------
    WrittenTable
        The file written and the unpacked bytes, which is what the digest is over.

    Raises
    ------
    ShippedTableError
        Of ``table``'s class: if a cell carries a tab, a newline or a carriage return or is
        not text at all, or if what was rendered is not a file the reader would accept — a
        header that is not the declared columns, a row with the wrong cell count, a blank cell
        in a required column, a flag spelled a way no table spells one, or a repeated key.

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> from genome.shipped import ShippedTable
    >>> tiny = ShippedTable(
    ...     resource="{slug}.tiny_table.tsv", columns=("id", "flag"), noun="tiny table"
    ... )
    >>> with TemporaryDirectory() as directory:
    ...     write_table(tiny, Path(directory), tiny.columns, [("g1", "yes")], slug="worm").unpacked
    b'id\tflag\ng1\tyes\n'
    """
    path = directory / shipped_name(table, **keys)
    origin = str(path)
    unpacked = render(columns, _checked(table, columns, rows, origin=origin))
    table.parse(unpacked.decode("utf-8"), origin=origin)
    directory.mkdir(parents=True, exist_ok=True)
    path.write_bytes(
        gzip.compress(unpacked, compresslevel=_COMPRESSION, mtime=0)
        if path.name.endswith(_PACKED_SUFFIX)
        else unpacked
    )
    return WrittenTable(path=path, unpacked=unpacked)


def merge_rows(
    table: ShippedTable,
    directory: Path,
    rows: Sequence[Mapping[str, str]],
    **keys: str,
) -> WrittenTable:
    r"""Rewrite a **Shipped table** with ``rows`` in place of every row they key over.

    How a generator that rebuilds one species leaves the others alone: the file is read, the
    rows sharing a key with one being written are dropped, the new ones take their place, and
    the whole file is written back sorted by that key — so a rebuild is a diff of the rows
    that changed and nothing else. The key is the one the declaration already names as the
    identity no two rows may repeat, so what makes a row the same row is stated once.

    Parameters
    ----------
    table : genome.shipped.ShippedTable
        The declaration of the table being merged into. It must name a key.
    directory : pathlib.Path
        Where the file is. Created, and the file written from nothing, if neither is there.
    rows : sequence of mapping of str to str
        The rows to write, each keyed by column name and carrying every declared column.
    **keys : str
        What the declaration's resource template names.

    Returns
    -------
    WrittenTable
        The file written and the bytes in it.

    Raises
    ------
    ValueError
        If the declaration names no key or declares only its leading columns, so there is
        nothing to replace a row by or no whole header to write back, or if a row leaves out a
        column the table declares.
    ShippedTableError
        Of ``table``'s class: if the file already there is not one this table can read, or the
        merged file would not be.

    Examples
    --------
    >>> from tempfile import TemporaryDirectory
    >>> from genome.shipped import ShippedTable
    >>> provenance = ShippedTable(
    ...     resource="tiny_metadata.tsv",
    ...     columns=("species", "file"),
    ...     noun="tiny provenance table",
    ...     key=("species",),
    ... )
    >>> with TemporaryDirectory() as directory:
    ...     mouse = merge_rows(provenance, Path(directory), [{"species": "Mus musculus", "file": "m"}])
    ...     human = merge_rows(provenance, Path(directory), [{"species": "Homo sapiens", "file": "h"}])
    ...     human.unpacked
    b'species\tfile\nHomo sapiens\th\nMus musculus\tm\n'
    """
    path = directory / shipped_name(table, **keys)
    if not table.key or table.leading:
        raise ValueError(
            f"{path.name} is not a table a merge can write: a merge replaces a row by the key "
            f"its declaration names, and writes back exactly the columns it declares. A table "
            f"naming no key, or declaring only the columns its header leads with, is written "
            f"whole with write_table instead."
        )
    for row in rows:
        missing = [column for column in table.columns if column not in row]
        if missing:
            raise ValueError(
                f"a row for {path.name} leaves out the columns {missing}, which that table "
                f"declares. Every column is written as it stands, so fill them in — a blank "
                f"cell is a value and never a column nobody supplied."
            )
    written = {tuple(row[column] for column in table.key) for row in rows}
    kept = [
        dict(existing)
        for existing in _existing(table, path)
        if tuple(existing[column] for column in table.key) not in written
    ]
    kept.extend(dict(row) for row in rows)
    kept.sort(key=lambda row: tuple(row[column] for column in table.key))
    return write_table(
        table,
        directory,
        table.columns,
        [tuple(row[column] for column in table.columns) for row in kept],
        **keys,
    )


def _existing(table: ShippedTable, path: Path) -> tuple[dict[str, str], ...]:
    """Return the rows already in ``path``, or none at all when no file is there yet."""
    if not path.is_file():
        return ()
    return table.parse(path.read_text(encoding="utf-8"), origin=str(path)).mappings()


def _checked(
    table: ShippedTable,
    columns: Sequence[str],
    rows: Iterable[Sequence[str]],
    *,
    origin: str,
) -> list[Sequence[str]]:
    """Return ``rows``, held to what one cell of a plain unquoted TSV can carry."""
    checked = list(rows)
    for number, row in enumerate(checked, start=1):
        for column, cell in zip(columns, row, strict=False):
            if not isinstance(cell, str):
                raise table.error(
                    f"{origin} would carry {cell!r} in the {column!r} column of {table.unit} "
                    f"{number}, which is not text. Every cell is written as it stands, so a "
                    f"value has to be rendered to text before it is written"
                    f"{table.repair_clause}"
                )
            for character, name in _UNWRITABLE:
                if character in cell:
                    raise table.error(
                        f"{origin} would carry {name} inside the {column!r} column of "
                        f"{table.unit} {number}: {cell!r}. Every {table.noun} is a plain TSV "
                        f"with no quoting, so a cell carrying one is a defect in what built it "
                        f"rather than something to escape around{table.repair_clause}"
                    )
    return checked
