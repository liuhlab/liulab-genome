#!/usr/bin/env python
r"""Rebuild one shipped **Motif link** table from one census and one JASPAR **Release**.

Run by hand, never at install time and never in CI — the release is a download CI
cannot make, which is why the shipped tables are guarded by ``tests/test_tf_link.py``
instead. This script lives outside the package for the same reason its sibling
``build_tf_census.py`` does: the wheel carries the tables and not the tooling that
made them.

**No test proves a shipped table still matches JASPAR, and none can.** Regenerating
one needs a download and CI has no network, so the guard tests pin counts instead:
drift shows up as a loud failure rather than as a table that quietly stopped
agreeing with the release it names. That is the most that is available, and it is
the cost ADR-0015 accepts.

Attribution
-----------
The profiles are JASPAR's — <https://jaspar.elixir.no/> — read from the SQLite dump
of one release, ``https://jaspar.elixir.no/download/database/JASPAR<year>.sqlite``.
Cite the JASPAR release you used. The genes are the census's, and the census carries
its own attribution: see ``src/genome/data/tf_gene/ATTRIBUTION.md``.

What it does
------------
Reads one **Release**'s SQLite dump and one shipped census, joins them by name, and
writes a gzipped TSV, as ``build_tf_census.py`` writes a census (ADR-0015).

*The join.* Upper-case the **Motif name**, split it on ``::``, and match each part
against the census's own symbol column. Only assessed-positive genes receive links.
The parts the census does not spell that way are handled by the shipped alias table
beside the data, which is keyed on **Gene id stem** and not on a symbol: JASPAR's
``SCAND3`` is Lambert's ``ZBED9``, and guessing at symbol history got that one wrong
twice. A profile naming no gene at all — the oncogenic fusion ``EWSR1-FLI1`` — stays
unlinked by design rather than asserting one.

*The species.* Per-profile species comes from the dump's ``MATRIX_SPECIES`` table,
which is the whole reason this reads SQLite: the transfac file the package's own
loader reads carries no species, so **the loader needs no change**. A profile
measured on another vertebrate is kept and marked rather than dropped (ADR-0013).

*The number.* The total **Information content** on each row is
:attr:`genome.tf.motif.motif.Motif.information_content` summed over the matrix's
columns — the package's own quantity, computed here from the dump's ``MATRIX_DATA``
rather than reimplemented. Those counts and the transfac file's agree to within
1e-6 bits on every profile of both releases, so the number is the same one whichever
serialization it is read from. **No quality score is computed or shipped**: JASPAR
publishes none and this package invents none.

Two properties are the point, as they are in ``build_tf_census.py``:

*It fails loudly on an unrecognised column.* Every table and column this reads out of
the dump is listed in :data:`JASPAR_TABLES`, and the census columns it reads are the
uniform four. A re-spelling is one missing name and one unexpected name, and both are
named in the error. The profile count is checked against the count the package's own
loader pins for that release, so a dump that is not the release it claims to be is an
error rather than a short table.

*It writes rows in a total order.* They are sorted by **Gene id stem** and then by
rank, and the total information content is rendered to a fixed number of decimals, so
what this hands over is already the same on every machine. Everything after that, the
rendering and the compression, is :mod:`genome.shipped_writer`'s.

Usage
-----
``python scripts/build_tf_links.py <species> <release> <the release's SQLite dump>``,
where ``<species>`` is a species slug a census ships for and ``<release>`` is one of
the releases the package prepares.

Add ``--data-dir`` to write somewhere other than ``src/genome/data/tf_link``.
"""

from __future__ import annotations

import argparse
import sqlite3
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType

import numpy as np

# The format is the package's and this script writes through it: the columns, the file
# name, the two **Role**s, the flag spellings and the separator are declared once, beside
# the reader, and the writer holds what is about to be written to exactly them.
from genome.shipped_writer import shipped_name, write_table
from genome.tf import COMPLEX, LINK_COLUMNS, LINK_FORMAT, LINK_SUBDIR, MONOMER, VALUE_SEPARATOR
from genome.tf.gene import (
    FALSE_CELL,
    TRUE_CELL,
    UNIFORM_COLUMNS,
    TFGeneTable,
    census_species,
    tf_gene_table,
)
from genome.tf.motif.jaspar import DEFAULT_TAX_GROUP, JASPAR_RELEASES, MOTIF_COUNTS
from genome.tf.motif.motif import BASES, Motif

#: Where the shipped link tables live, relative to the repository root — the package's
#: own resource directory, under the source tree it is packaged from.
DATA_SUBDIR = Path("src") / "genome" / LINK_SUBDIR

#: The alias table beside them: the **Motif name** parts no census spells that way. Three
#: rows, and **plain**, as every small metadata table in this package's data is — the
#: gzip is for the bulk tables, and a three-row file curated by hand is worth more as a
#: readable diff than as 200 saved bytes. It is read here and by nothing in the wheel:
#: it is an input to the join and never an answer the package gives.
ALIAS_FILE = "motif_name_alias.tsv"

#: The alias table's columns, in table order. ``gene_id_stem`` is the key: a symbol is
#: what moved, so keying on one is keying on the thing that changed.
ALIAS_COLUMNS = ("species", "motif_name_part", "gene_id_stem", "census_symbol", "note")

#: How the **Motif name** spells the join between the genes a complex names.
NAME_SEPARATOR = "::"

#: Decimals the total **Information content** is written to. The rank is computed from
#: the written value rather than from the full float, so re-sorting the shipped table on
#: the shipped columns reproduces the shipped rank exactly.
IC_DECIMALS = 4

#: Every table and column this script reads out of a release's SQLite dump. **This is the
#: fail-loud surface**: a dump that re-spells one of these is refused by name rather than
#: having the column silently drop out of a regenerated table.
JASPAR_TABLES: Mapping[str, tuple[str, ...]] = MappingProxyType(
    {
        "MATRIX": ("ID", "COLLECTION", "BASE_ID", "VERSION", "NAME"),
        "MATRIX_ANNOTATION": ("LOCAL_ID", "ID", "TAG", "VAL"),
        "MATRIX_DATA": ("ID", "row", "col", "val"),
        "MATRIX_SPECIES": ("ID", "TAX_ID"),
    }
)

#: JASPAR's collection this reads, and the only one: CORE is the curated, non-redundant
#: set the package's own loader ships.
COLLECTION = "CORE"

#: The annotation tag carrying a profile's **Tax group**, and the value kept.
TAX_GROUP_TAG = "tax_group"


class MotifLinkSourceError(ValueError):
    """A link table cannot be built from what it was handed.

    Either the release's dump is not the file this script was written against, the
    census is not the shape a census is, the alias table says something no census
    confirms, or what the join would write is not a file a **Motif link** reader would
    accept back. Every message names what was wrong and what to fix.
    """


#: What repairs a refused build: never *re-run this*, which is what the reader's own
#: refusals say, because re-running is what has just failed.
_REPAIR = (
    "fix the join in scripts/build_tf_links.py, or reconcile the census and the release's "
    "dump before building from them"
)

#: A **Motif link** table as this script writes one: the package's own declaration of the
#: format, refusing in this script's words. Everything about what a link table *is* comes
#: from :mod:`genome.tf.link`; only the error class and the repair are the build's.
WRITTEN_LINKS = replace(LINK_FORMAT, error=MotifLinkSourceError, repair=_REPAIR)


@dataclass(frozen=True)
class Profile:
    """One JASPAR profile of one **Release**, as much of it as a **Motif link** needs.

    Attributes
    ----------
    motif_id : str
        The **Motif id** — ``BASE_ID.VERSION``, ``"MA0139.2"``.
    motif_name : str
        The **Motif name** JASPAR publishes, in JASPAR's own spelling.
    tax_ids : tuple of str
        The NCBI taxonomy ids the dump files this profile under, ascending. Empty for
        a profile the dump records no species for.
    total_information_content : float
        The matrix's **Information content** summed over its columns, in bits,
        rounded to :data:`IC_DECIMALS`.
    """

    motif_id: str
    motif_name: str
    tax_ids: tuple[str, ...]
    total_information_content: float

    @property
    def name_parts(self) -> tuple[str, ...]:
        """Return the upper-cased **Motif name** split on ``::`` — the genes it names."""
        return tuple(self.motif_name.upper().split(NAME_SEPARATOR))

    @property
    def role(self) -> str:
        """Return :data:`MONOMER` where this profile names one gene, :data:`COMPLEX` otherwise."""
        return MONOMER if len(self.name_parts) == 1 else COMPLEX


@dataclass(frozen=True)
class Alias:
    """One row of the alias table: a **Motif name** part the census spells another way.

    Attributes
    ----------
    species : str
        The species whose census this row speaks about, as that census spells it.
    motif_name_part : str
        The upper-cased **Motif name** part JASPAR publishes — ``"SCAND3"``.
    gene_id_stem : str
        The **Gene id stem** it names. The key: a symbol is what moved.
    census_symbol : str
        What the census calls that gene — ``"ZBED9"``. Checked against the census, so
        a stale row is an error rather than a comment nobody reads.
    note : str
        Why the row is here, for whoever reviews the diff.
    """

    species: str
    motif_name_part: str
    gene_id_stem: str
    census_symbol: str
    note: str


@dataclass(frozen=True)
class Link:
    """One **Motif link**: one profile answering for one **TF gene**.

    Attributes
    ----------
    gene_id_stem : str
        The gene the census keys this verdict by.
    symbol : str
        The census's own symbol for that gene, which is not always the **Motif name**.
    profile : Profile
        The profile answering for it.
    partners : tuple of str
        The other genes the **Motif name** names, upper-cased as the name spells them.
        Empty exactly when the **Role** is :data:`MONOMER`.
    is_cross_species : bool
        Whether the profile was measured on a vertebrate other than this gene's own
        species — a **Cross-species link**, kept and marked (ADR-0013).
    """

    gene_id_stem: str
    symbol: str
    profile: Profile
    partners: tuple[str, ...]
    is_cross_species: bool

    def specificity_key(self) -> tuple[int, int, float, str]:
        """Return this link's place under **Attribution specificity**, lowest first.

        Four keys, so the order is total and stable across machines and **Release**s:
        **Role** ``monomer`` before ``complex``, species-matched before a
        **Cross-species link**, then higher total **Information content**, then
        **Motif id**. It states what a matrix is attributable to and explicitly not
        which motif is better — the canonical AP-1 matrix is a complex and describes
        JUN's binding better than any JUN monomer does.
        """
        return (
            0 if self.profile.role == MONOMER else 1,
            1 if self.is_cross_species else 0,
            -self.profile.total_information_content,
            self.profile.motif_id,
        )


def read_profiles(path: Path, release: str) -> tuple[Profile, ...]:
    """Return one **Release**'s CORE vertebrate profiles, read from its SQLite dump.

    The non-redundant release is reconstructed the way JASPAR publishes it: within the
    CORE collection, the highest version of each base id, filed under the
    ``vertebrates`` **Tax group**. The result is counted against what the package's own
    loader pins for that release, so a dump that is not that release is refused.

    Parameters
    ----------
    path : pathlib.Path
        The release's SQLite dump.
    release : str
        Which **Release** the dump is meant to be.

    Returns
    -------
    tuple of Profile
        Every profile of that release, ordered by **Motif id**.

    Raises
    ------
    MotifLinkSourceError
        If the dump re-spells a table or a column, holds a matrix this cannot read, or
        yields a profile count the loader does not pin for that release.
    """
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        check_schema(connection, origin=str(path))
        rows = connection.execute(
            "SELECT ID, BASE_ID, VERSION, NAME FROM MATRIX WHERE COLLECTION = ?", (COLLECTION,)
        ).fetchall()
        newest: dict[str, tuple[int, int, str]] = {}
        for identifier, base, version, name in rows:
            if base not in newest or int(version) > newest[base][1]:
                newest[base] = (int(identifier), int(version), str(name))
        groups = dict(
            connection.execute(
                "SELECT ID, VAL FROM MATRIX_ANNOTATION WHERE TAG = ?", (TAX_GROUP_TAG,)
            ).fetchall()
        )
        kept = {
            base: record
            for base, record in newest.items()
            if groups.get(record[0]) == DEFAULT_TAX_GROUP
        }
        check_release_size(len(kept), release, origin=str(path))
        matrices = read_matrices(
            connection, {record[0] for record in kept.values()}, origin=str(path)
        )
        species = read_species(connection)
    finally:
        connection.close()
    profiles = [
        Profile(
            motif_id=f"{base}.{version}",
            motif_name=name,
            tax_ids=species.get(identifier, ()),
            total_information_content=total_information_content(
                f"{base}.{version}", name, matrices[identifier], origin=str(path)
            ),
        )
        for base, (identifier, version, name) in kept.items()
    ]
    return tuple(sorted(profiles, key=lambda profile: profile.motif_id))


def total_information_content(
    motif_id: str, motif_name: str, counts: np.ndarray, *, origin: str
) -> float:
    """Return one matrix's **Information content** summed over its columns, in bits.

    The package's own quantity and not a second implementation of it: the counts are
    handed to a :class:`~genome.tf.motif.motif.Motif` and its per-position information
    content is summed, so the number shipped is the number a caller would compute from
    the same release.

    Raises
    ------
    MotifLinkSourceError
        If the dump's counts are not a matrix a motif can be made from — a negative
        cell, or a column of zeros. Nothing can be said about such a profile, and
        shipping a row with a made-up number for it would be worse than refusing.
    """
    try:
        motif = Motif(motif_id, motif_name, counts)
    except ValueError as error:
        raise MotifLinkSourceError(
            f"{origin} holds a count matrix for {motif_id} that no motif can be made from: "
            f"{error}. Re-download the dump; if it reads the same way, that profile is broken "
            f"upstream and JASPAR is who to tell."
        ) from error
    return round(float(motif.information_content.sum()), IC_DECIMALS)


def check_schema(connection: sqlite3.Connection, *, origin: str) -> None:
    """Hold the dump to the tables and columns this script reads, naming every difference.

    Raises
    ------
    MotifLinkSourceError
        If a table is absent, or one carries columns other than the ones
        :data:`JASPAR_TABLES` names. A release that re-spells a column would otherwise
        have it silently drop out of a regenerated table, which is what this prevents.
    """
    faults: list[str] = []
    for table, expected in JASPAR_TABLES.items():
        seen = tuple(
            str(row[1]) for row in connection.execute(f"PRAGMA table_info('{table}')").fetchall()
        )
        if not seen:
            faults.append(f"{table} is not in the dump at all")
            continue
        missing = [name for name in expected if name not in seen]
        unexpected = [name for name in seen if name not in expected]
        if missing:
            faults.append(f"{table} does not carry {missing}")
        if unexpected:
            faults.append(f"{table} carries {unexpected}, which this script does not name")
    if faults:
        raise MotifLinkSourceError(
            f"{origin} is not the dump this script was written against — "
            + "; ".join(faults)
            + ". Update JASPAR_TABLES in this script to match the release you are building "
            "from, and say in the commit which columns moved."
        )


def check_release_size(count: int, release: str, *, origin: str) -> None:
    """Hold the selected profiles to the count the package's own loader pins for that release.

    Raises
    ------
    MotifLinkSourceError
        If the two disagree. The loader counts the release's published transfac file and
        this counts its SQLite dump, so a disagreement means one of them is not that
        release — a truncated download, or a selection rule that has stopped
        reconstructing the non-redundant set.
    """
    pinned = MOTIF_COUNTS.get((release, DEFAULT_TAX_GROUP))
    if pinned is None:
        raise MotifLinkSourceError(
            f"nothing pins how many {DEFAULT_TAX_GROUP} profiles the {release} release holds, so "
            f"there is nothing to check {origin} against. Add the count to MOTIF_COUNTS first."
        )
    if count != pinned:
        raise MotifLinkSourceError(
            f"{origin} yields {count} CORE {DEFAULT_TAX_GROUP} profiles where the {release} "
            f"release ships {pinned}. Either the dump is not that release — re-download it from "
            f"https://jaspar.elixir.no/download/database/JASPAR{release}.sqlite — or the rule "
            f"that reconstructs the non-redundant set from the dump no longer does."
        )


def read_matrices(
    connection: sqlite3.Connection, wanted: set[int], *, origin: str
) -> dict[int, np.ndarray]:
    """Return one 4 x L **Count matrix** per wanted profile, in :data:`BASES` order.

    Raises
    ------
    MotifLinkSourceError
        If a profile has no matrix, names a base outside :data:`BASES`, or is missing a
        cell — a matrix that cannot be read is a dump that cannot be trusted.
    """
    cells: dict[int, dict[tuple[str, int], float]] = {}
    for identifier, base, column, value in connection.execute(
        "SELECT ID, row, col, val FROM MATRIX_DATA"
    ):
        if int(identifier) not in wanted:
            continue
        if str(base) not in BASES:
            raise MotifLinkSourceError(
                f"{origin} gives matrix {identifier} a row named {base!r}, and a **Count matrix** "
                f"has one row per base in {''.join(BASES)}. This dump is not the shape this "
                f"script reads."
            )
        cells.setdefault(int(identifier), {})[(str(base), int(column))] = float(value)
    matrices: dict[int, np.ndarray] = {}
    for identifier in sorted(wanted):
        table = cells.get(identifier)
        if not table:
            raise MotifLinkSourceError(
                f"{origin} carries no matrix data for profile {identifier}, so nothing can be "
                f"said about its information content. Re-download the dump."
            )
        length = max(column for _, column in table)
        try:
            matrices[identifier] = np.array(
                [[table[(base, column)] for column in range(1, length + 1)] for base in BASES],
                dtype=np.float64,
            )
        except KeyError as error:
            raise MotifLinkSourceError(
                f"{origin} leaves matrix {identifier} a cell short at {error.args[0]}. A count "
                f"matrix with a hole in it cannot be read — re-download the dump."
            ) from error
    return matrices


def read_species(connection: sqlite3.Connection) -> dict[int, tuple[str, ...]]:
    """Return each profile's NCBI taxonomy ids, ascending.

    This table is the whole reason a link table is built from the dump rather than from
    the transfac file the package's own loader reads: that file carries no per-profile
    species, so the loader needs no change to make this join possible.
    """
    collected: dict[int, set[str]] = {}
    for identifier, tax_id in connection.execute("SELECT ID, TAX_ID FROM MATRIX_SPECIES"):
        text = str(tax_id).strip() if tax_id is not None else ""
        if text:
            collected.setdefault(int(identifier), set()).add(text)
    return {identifier: tuple(sorted(values, key=int)) for identifier, values in collected.items()}


def read_aliases(path: Path) -> tuple[Alias, ...]:
    """Return the shipped alias table, or an empty tuple when none ships yet.

    Raises
    ------
    MotifLinkSourceError
        If the header is not :data:`ALIAS_COLUMNS`, a row has the wrong number of cells,
        or a row leaves a cell blank. Every column is required: a row nobody can check
        is a claim nobody made.
    """
    if not path.is_file():
        return ()
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    header = tuple(lines[0].split("\t")) if lines else ()
    if header != ALIAS_COLUMNS:
        raise MotifLinkSourceError(
            f"{path} carries the columns {list(header)} where the alias table's are "
            f"{list(ALIAS_COLUMNS)}. Fix the file's header."
        )
    aliases = []
    for number, line in enumerate(lines[1:], start=2):
        cells = line.split("\t")
        if len(cells) != len(ALIAS_COLUMNS):
            raise MotifLinkSourceError(
                f"{path} line {number} holds {len(cells)} cells where its header declares "
                f"{len(ALIAS_COLUMNS)}."
            )
        if not all(cell.strip() for cell in cells):
            raise MotifLinkSourceError(
                f"{path} line {number} leaves a cell blank, and every alias column is required: "
                f"a row nobody can check against a census is a claim nobody made."
            )
        aliases.append(Alias(*(cell.strip() for cell in cells)))
    return tuple(aliases)


def symbol_lookup(census: TFGeneTable, aliases: Sequence[Alias], *, origin: str) -> dict[str, str]:
    """Return upper-cased symbol to **Gene id stem**, for the genes the census assessed positive.

    The census's own symbols first, then the aliases on top of them. Only
    assessed-positive genes are in it, so a gene the census assessed and turned down
    receives no link however a profile is named.

    Raises
    ------
    MotifLinkSourceError
        If the census does not lead with the uniform four, if two assessed-positive
        genes share a symbol, or if an alias names a gene the census does not assess
        positive, contradicts the census's own symbol, is not upper-cased, or shadows a
        symbol the census already spells that way.
    """
    if census.columns[: len(UNIFORM_COLUMNS)] != UNIFORM_COLUMNS:
        raise MotifLinkSourceError(
            f"{origin} leads with the columns {list(census.columns[: len(UNIFORM_COLUMNS)])} "
            f"where every census leads with {list(UNIFORM_COLUMNS)}. This script reads the "
            f"symbol and the TF flag off those four, so a census without them cannot be joined."
        )
    genes = assessed_genes(census)
    lookup: dict[str, str] = {}
    for stem, (symbol, positive) in genes.items():
        if not positive or not symbol:
            continue
        if symbol.upper() in lookup:
            raise MotifLinkSourceError(
                f"{origin} spells the symbol {symbol.upper()!r} for both "
                f"{lookup[symbol.upper()]} and {stem}, so a profile named that way would link "
                f"to whichever came first. Reconcile the census before building from it."
            )
        lookup[symbol.upper()] = stem
    for alias in aliases:
        check_alias(alias, lookup, genes, origin=origin)
        lookup[alias.motif_name_part] = alias.gene_id_stem
    return lookup


def assessed_genes(census: TFGeneTable) -> dict[str, tuple[str, bool]]:
    """Return each **Gene id stem** the census assessed, with its symbol and its verdict."""
    symbol_at = UNIFORM_COLUMNS.index("symbol")
    is_tf = UNIFORM_COLUMNS.index("is_tf")
    return {row[0] or "": (row[symbol_at] or "", row[is_tf] == TRUE_CELL) for row in census.rows}


def check_alias(
    alias: Alias,
    lookup: Mapping[str, str],
    genes: Mapping[str, tuple[str, bool]],
    *,
    origin: str,
) -> None:
    """Hold one alias row to the census it speaks about.

    Raises
    ------
    MotifLinkSourceError
        If the row is not upper-cased, shadows a symbol the census already spells that
        way, names a gene the census does not carry or does not assess positive, or
        disagrees with the census about that gene's symbol.
    """
    if alias.motif_name_part != alias.motif_name_part.upper():
        raise MotifLinkSourceError(
            f"the alias {alias.motif_name_part!r} is not upper-cased, and the join upper-cases "
            f"every **Motif name** before splitting it — so this row could never match."
        )
    if alias.motif_name_part in lookup:
        raise MotifLinkSourceError(
            f"the alias {alias.motif_name_part!r} names a symbol {origin} already spells that "
            f"way, for {lookup[alias.motif_name_part]}. An alias exists for the parts no census "
            f"spells that way; this one would silently shadow a gene the census does carry."
        )
    if alias.gene_id_stem not in genes:
        raise MotifLinkSourceError(
            f"the alias {alias.motif_name_part!r} names the gene {alias.gene_id_stem}, which "
            f"{origin} does not assess at all. An alias resolves to a gene the census carries, "
            f"or it resolves to nothing anyone can check."
        )
    published, positive = genes[alias.gene_id_stem]
    if not positive:
        raise MotifLinkSourceError(
            f"the alias {alias.motif_name_part!r} names the gene {alias.gene_id_stem}, which "
            f"{origin} assessed and did not judge a transcription factor. Only assessed-positive "
            f"genes receive links, so this row would assert a verdict the census did not reach."
        )
    if published.upper() != alias.census_symbol.upper():
        raise MotifLinkSourceError(
            f"the alias for {alias.gene_id_stem} records the census symbol "
            f"{alias.census_symbol!r} where {origin} spells it {published!r}. The alias table is "
            f"keyed on the gene id precisely because symbols move — so fix the recorded symbol, "
            f"and check the gene id is still the gene you meant."
        )


def build_links(
    profiles: Sequence[Profile], lookup: Mapping[str, str], symbols: Mapping[str, str], taxid: str
) -> list[Link]:
    """Return every **Motif link** between these profiles and this census's genes.

    One link per (profile, gene) pair, so a ``FOS::JUN`` profile contributes two — and a
    gene whose only motifs are complexes is linked rather than reported motif-less.
    """
    links: list[Link] = []
    for profile in profiles:
        parts = profile.name_parts
        for index, part in enumerate(parts):
            stem = lookup.get(part)
            if stem is None:
                continue
            links.append(
                Link(
                    gene_id_stem=stem,
                    symbol=symbols.get(stem, ""),
                    profile=profile,
                    partners=tuple(
                        other for position, other in enumerate(parts) if position != index
                    ),
                    is_cross_species=taxid not in profile.tax_ids,
                )
            )
    return links


def rank_links(links: Sequence[Link]) -> list[tuple[Link, int]]:
    """Return every link with its rank, dense from one within a gene, in table order.

    Ranked under **Attribution specificity** and then laid out by **Gene id stem** and
    rank, which is a total order over the whole table and the reason two runs agree byte
    for byte.
    """
    by_gene: dict[str, list[Link]] = {}
    for link in links:
        by_gene.setdefault(link.gene_id_stem, []).append(link)
    ranked: list[tuple[Link, int]] = []
    for stem in sorted(by_gene):
        ordered = sorted(by_gene[stem], key=lambda link: link.specificity_key())
        ranked.extend((link, rank) for rank, link in enumerate(ordered, start=1))
    return ranked


def build_table(
    ranked: Sequence[tuple[Link, int]], *, release: str, species: str
) -> list[tuple[str, ...]]:
    """Return the shipped rows as text, one tuple per **Motif link**, in table order."""
    return [
        (
            release,
            species,
            link.gene_id_stem,
            link.symbol,
            link.profile.motif_id,
            link.profile.motif_name,
            link.profile.role,
            VALUE_SEPARATOR.join(link.partners),
            VALUE_SEPARATOR.join(link.profile.tax_ids),
            TRUE_CELL if link.is_cross_species else FALSE_CELL,
            f"{link.profile.total_information_content:.{IC_DECIMALS}f}",
            str(rank),
        )
        for link, rank in ranked
    ]


def check_table(rows: Sequence[Sequence[str]], *, origin: str) -> None:
    """Hold the built table to what a shipped link table promises, or refuse to write it.

    What is checked here is what a link *means*, since the shape a **Shipped table**
    has — the header, the width of a row, what a cell may carry, the columns no row may
    leave blank — is checked as the file is written, against the package's own
    declaration of the format.

    Raises
    ------
    MotifLinkSourceError
        If the table is empty, a **Role** is neither of the two declared, or a complex
        names no partner or a monomer names one — a role that disagrees with its
        partners is the confusion **Role** exists to end.
    """
    if not rows:
        raise MotifLinkSourceError(
            f"{origin} would carry no links at all. A census with no profile answering for any "
            f"of its genes is a join that did not happen — check the release and the census "
            f"name each other's genes before shipping an empty table."
        )
    role_at = LINK_COLUMNS.index("role")
    partners_at = LINK_COLUMNS.index("partners")
    for number, row in enumerate(rows, start=2):
        if row[role_at] not in (MONOMER, COMPLEX):
            raise MotifLinkSourceError(
                f"{origin} line {number} would spell its role {row[role_at]!r}, and a link is "
                f"{MONOMER!r} or {COMPLEX!r} and nothing else."
            )
        if (row[role_at] == COMPLEX) != bool(row[partners_at]):
            raise MotifLinkSourceError(
                f"{origin} line {number} spells its role {row[role_at]!r} with partners "
                f"{row[partners_at]!r}. A complex names at least one partner and a monomer names "
                f"none — that is what keeps a heterodimer matrix from being read as a monomer's."
            )


def report_residuals(
    profiles: Sequence[Profile], lookup: Mapping[str, str], taxid: str, aliases: Sequence[Alias]
) -> None:
    """Print what this species' own profiles resolved to, and what did not resolve.

    The denominator is the profiles the dump files under this species' taxid, which is
    the only set where a profile failing to name one of this census's genes is worth
    printing: a profile measured on another vertebrate naming a gene this census never
    assessed is ordinary rather than a residual.

    **An unresolved profile is not a pending alias.** It is as often a profile naming a
    gene the census never assessed at all, which no alias can reach and which is a
    correct answer — mouse ``MA0611.3 Dux`` is one, and ``ATTRIBUTION.md`` beside the
    data records why. This list is a report, not a queue.
    """
    tagged = [profile for profile in profiles if taxid in profile.tax_ids]
    unresolved = [
        profile for profile in tagged if not any(part in lookup for part in profile.name_parts)
    ]
    print(f"  species-tagged profiles ({taxid}): {len(tagged)}")
    print(f"    resolved to at least one assessed-positive gene: {len(tagged) - len(unresolved)}")
    for profile in unresolved:
        print(f"    unresolved: {profile.motif_id} {profile.motif_name}")
    if unresolved:
        print("      (a profile naming a gene the census never assessed is a correct answer)")
    named = {part for profile in profiles for part in profile.name_parts}
    unused = sorted(
        alias.motif_name_part for alias in aliases if alias.motif_name_part not in named
    )
    print(
        f"  aliases: {len(aliases) - len(unused)} of {len(aliases)} name a profile in this release"
    )
    for part in unused:
        print(f"    unused here: {part}")


def build(species: str, release: str, dump: Path, data_dir: Path) -> None:
    """Build one link table from one census and one release, and report what was written."""
    census = tf_gene_table(species)
    if census is None:
        raise MotifLinkSourceError(
            f"no census ships for {species!r}, so there are no genes to link profiles to. "
            f"The censuses that ship are {list(census_species())} — build one with "
            f"scripts/build_tf_census.py first."
        )
    origin = census.provenance.file
    aliases = [
        alias for alias in read_aliases(data_dir / ALIAS_FILE) if alias.species == census.species
    ]
    lookup = symbol_lookup(census, aliases, origin=origin)
    symbols = {stem: symbol for stem, (symbol, _) in assessed_genes(census).items()}
    taxid = str(census.provenance.ncbi_taxid)
    profiles = read_profiles(dump, release)
    ranked = rank_links(build_links(profiles, lookup, symbols, taxid))
    file_name = shipped_name(LINK_FORMAT, slug=species, release=release)
    rows = build_table(ranked, release=release, species=census.species)
    check_table(rows, origin=file_name)
    written = write_table(
        WRITTEN_LINKS, data_dir, LINK_COLUMNS, rows, slug=species, release=release
    )
    genes = {link.gene_id_stem for link, _ in ranked}
    cross = sum(1 for link, _ in ranked if link.is_cross_species)
    monomer = sum(1 for link, _ in ranked if link.profile.role == MONOMER)
    print(f"{file_name}: {len(genes)} genes, {len(rows)} links")
    print(f"  {monomer} monomer, {len(rows) - monomer} complex; {cross} cross-species")
    print(f"  {len(profiles)} CORE {DEFAULT_TAX_GROUP} profiles in the {release} release")
    report_residuals(profiles, lookup, taxid, aliases)
    print(f"  unpacked {len(written.unpacked)} bytes")
    print(f"  packed   {written.packed} bytes ({len(written.unpacked) / written.packed:.1f}x)")


def main(argv: Sequence[str] | None = None) -> None:
    """Parse the command line and build the link table it names."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "species", choices=sorted(census_species()), help="which census's genes to link"
    )
    parser.add_argument(
        "release", choices=JASPAR_RELEASES, help="which JASPAR release to link it to"
    )
    parser.add_argument("dump", type=Path, help="that release's SQLite dump")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / DATA_SUBDIR,
        help="where the link table is written",
    )
    arguments = parser.parse_args(argv)
    try:
        build(arguments.species, arguments.release, arguments.dump, arguments.data_dir)
    except MotifLinkSourceError as error:
        # The message already says what was wrong and what to do about it; a traceback on
        # top of it would only bury that.
        raise SystemExit(str(error)) from None


if __name__ == "__main__":
    main()
