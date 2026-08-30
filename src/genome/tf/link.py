"""The **Motif link** join — which JASPAR motifs answer for one **TF gene**.

The join neither half owns. It imports the gene half, :mod:`genome.tf.gene`, and the
motif half, :mod:`genome.tf.motif`, and neither of them imports it: the two are keyed
differently — one by gene, one by motif — and the mapping between them is many-to-many,
so it belongs beside both rather than inside either.

**The mapping is data and not a rule** (ADR-0015). One gzipped TSV ships per species per
**Release**, CORE ``vertebrates`` only, under ``data/tf_link/<species
slug>.jaspar<release>.motif_link_table.tsv.gz``; this module reads them and nothing more.
It is a leaf: it belongs to no **Assembly**, needs no **Data dir** and reaches no network,
and the tables it reads are readable in R or a shell by collaborators who never import
this package — which was the reason for shipping them, and gzip is what every such tool
already opens. ``ATTRIBUTION.md`` beside the files records how the join was made, what it
deliberately leaves unlinked, and the counts that are pinned.

**Bulk ships gzipped; small metadata ships plain.** That is the convention across every
data directory in this package, and this one is where both halves of it are visible: the
four link tables are gzipped and the three-row ``motif_name_alias.tsv`` beside them is
not, exactly as ``census_metadata.tsv`` and the assembly and annotation metadata tables
are not. A file small enough to be curated by hand is worth more as a readable diff than
as the bytes gzip would save; a few hundred kilobytes of generated rows is not.

**Order is Attribution specificity, and it is not a quality score.** A gene's links come
back **Role** ``monomer`` before ``complex``, species-matched before a **Cross-species
link**, then higher total **Information content**, then **Motif id** — four keys, so the
order is total and stable and "the motif for this factor" means the same thing on two
machines and in two releases. It states what a matrix is *attributable to* and explicitly
**not** which motif is better: JASPAR's canonical AP-1 matrix is the complex
``MA0099.4 FOS::JUN``, it describes JUN's binding better than any JUN monomer does, and
it still ranks below both JUN monomers here because it is a motif *of* a complex. No
quality score is computed or shipped anywhere — JASPAR publishes none, and matrix depth
is normalised per assay, so ranking on depth ranks the assay. A caller who disagrees
re-sorts on the attributes every link already carries. The shipped ``rank`` column is
where this order lives; this module honours it rather than re-deriving a second ordering
that could disagree with the file.

**Absence is not emptiness**, in the two layers the curated gene lists and the censuses
already set. :func:`motif_link_table` answers ``None`` for a species, **Release** or
**Tax group** no table ships for — the raw absence, and the one place ``None`` is how it
is said, because this is the layer below the one a caller touches.
:func:`motif_links` turns that into :class:`NoMotifLinkTableError`, and a gene no census
assessed into :class:`GeneNotAssessedError` — narrowed to
:class:`TranscriptionCofactorError` where a publisher lists that gene as a **Transcription
cofactor**, which is the same absence with a reason attached to it. None of them is ever an
empty collection. A gene that *is* assessed and has no motif is a real answer and comes
back with no links — :attr:`MotifLinks.is_tf` is what tells a gene its census turned down
from one JASPAR has no profile for.

Examples
--------
>>> from genome.tf import motif_links
>>> jun = motif_links("JUN", "Homo sapiens")
>>> jun.release, jun.tax_group
('2026', 'vertebrates')
>>> [(link.motif_id, link.role) for link in jun][:3]
[('MA0488.2', 'monomer'), ('MA0489.3', 'monomer'), ('MA1131.2', 'complex')]
>>> jun[2].partners
('FOSL2',)
>>> mouse = motif_links("Ctcf", "Mus musculus", cross_species=False)
>>> len(mouse), mouse.symbol, mouse.is_tf
(0, 'Ctcf', True)
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator, Mapping
from dataclasses import dataclass
from functools import cache
from importlib.resources import files
from types import MappingProxyType

import pandas as pd

from genome.metadata import species_slug
from genome.tf.cofactor import BOTH, SOURCES, UNIFORM_COLUMNS, CofactorTable, cofactor_table
from genome.tf.gene import FALSE_CELL, TRUE_CELL, TFGeneTable, tf_gene_table
from genome.tf.motif.jaspar import DEFAULT_RELEASE, DEFAULT_TAX_GROUP

#: Directory inside the package holding one **Motif link** table per species per
#: **Release**, plus the alias table the generator used and ``ATTRIBUTION.md``.
LINK_SUBDIR = "data/tf_link"

#: What names one of those files: the species slug, this prefix and the **Release**, then
#: :data:`LINK_SUFFIX` — ``homo_sapiens.jaspar2026.motif_link_table.tsv.gz``. Two keys
#: name one table, so both are in the name and neither is enumerated in code: a new
#: species or a new release is a file dropped in, exactly as it is for a census.
RELEASE_PREFIX = "jaspar"

#: The end of every link table's name, and what enumerating the directory matches on. The
#: tables are gzipped, as the censuses beside them are, and the decompression happens at
#: the resource boundary so :func:`parse_motif_link_table` stays a pure function of text.
LINK_SUFFIX = ".motif_link_table.tsv.gz"

#: The one **Tax group** the shipped tables cover, and the reason it is not in a file
#: name: JASPAR's CORE ``vertebrates``, which is where every human and mouse profile is.
#: Asking for another is answered rather than returning nothing — the tables for it do
#: not exist, which is a different fact from a gene having no motifs. Declared here in
#: full rather than aliased to the motif half's default, because what ships is a property
#: of these files and not of what a scan defaults to.
LINK_TAX_GROUP = "vertebrates"

#: Every link table's columns, in table order. Identical across tables, which is what
#: lets two of them concatenate into one frame that still says what each row came from.
LINK_COLUMNS: tuple[str, ...] = (
    "release",
    "species",
    "gene_id_stem",
    "symbol",
    "motif_id",
    "motif_name",
    "role",
    "partners",
    "motif_tax_ids",
    "is_cross_species",
    "total_information_content",
    "rank",
)

#: The two **Role**s, and the only two. ``monomer`` where the profile names one gene,
#: ``complex`` otherwise, so a heterodimer matrix is never read as a monomer's.
MONOMER, COMPLEX = "monomer", "complex"

#: What separates one value from the next inside a cell — a complex's partners, and a
#: profile's tax ids. A semicolon and never a tab: the file carries no quoting.
VALUE_SEPARATOR = ";"

#: The shape of a versioned gene id, tried only after a **Gene id stem** and a symbol
#: have both been looked for. Lambert's census carries clone-style symbols — ``AC023509.3``
#: is one — so a name of this shape is a symbol first and a versioned id second.
_VERSIONED_GENE_ID = re.compile(r"(?P<stem>[^.]+)\.\d+\w*")

#: Where a **Cofactor table** keeps the flag and the publisher, among the uniform columns
#: every such table leads with. Read off that tuple rather than written down again, so the
#: two modules cannot drift apart: the flag says whether the publisher listed the gene as a
#: cofactor or recorded a rejection, and the source says who listed it.
_COFACTOR_FLAG = UNIFORM_COLUMNS.index("is_cofactor")
_COFACTOR_SOURCE = UNIFORM_COLUMNS.index("source")


class MotifLinkTableError(ValueError):
    r"""A shipped **Motif link** table cannot be read, so it is not allowed to answer.

    A **packaging defect** and not a caller error: these files ship inside the package and
    are written by ``scripts/build_tf_links.py``, so a header that is not the twelve
    columns, a row with the wrong number of cells, a **Role** nothing spells that way or
    two releases inside one file are faults in what was committed here. A
    :class:`ValueError`, because a file that says something the format does not is a bad
    value rather than a broken program.

    The message names the file and what is wrong with it, since regenerating that file is
    the only thing anyone can do about it.

    Examples
    --------
    >>> try:
    ...     parse_motif_link_table("release\tspecies\n", source="broken.tsv")
    ... except MotifLinkTableError as error:
    ...     print("broken.tsv" in str(error))
    True
    """


class NoMotifLinkTableError(LookupError):
    """No **Motif link** table ships for that species, **Release** or **Tax group**.

    The first of the two absences, and the one a caller must never read as *no motif
    answers for this gene*: no table was ever built for what was asked, so the question
    was not answered rather than answered in the negative. A :class:`LookupError`, as the
    curated gene lists' own pair of absences are, so a caller may catch that and still
    tell this from :class:`GeneNotAssessedError`.

    The message names the species, releases and **Tax group** that do ship.

    Examples
    --------
    >>> try:
    ...     motif_links("CTCF", "Homo sapiens", tax_group="plants")
    ... except NoMotifLinkTableError as error:
    ...     print("vertebrates" in str(error))
    True
    """


class GeneNotAssessedError(LookupError):
    """The census for that species never assessed that gene, so nothing can answer for it.

    The second absence. A gene the census assessed and turned down is *not* this: it has a
    verdict, and comes back with no links and :attr:`MotifLinks.is_tf` ``False``. This is
    the gene the census never looked at — a symbol it does not spell, a **Gene id stem**
    of another species, or a typo — and answering it emptily would read as *this gene has
    no motifs*. A gene a publisher lists as a **Transcription cofactor** raises the narrower
    :class:`TranscriptionCofactorError` instead, which is this absence with a reason.

    The message names the species and the census that speaks for it.

    Examples
    --------
    >>> try:
    ...     motif_links("ENSMUSG00000005698", "Homo sapiens")
    ... except GeneNotAssessedError as error:
    ...     print("Lambert et al. 2018" in str(error))
    True
    """


class TranscriptionCofactorError(GeneNotAssessedError):
    """No census assessed that gene, and a publisher lists it as a **Transcription cofactor**.

    The same absence with something known in its place: a cofactor acts on transcription
    without recognising a sequence of its own, so there is **no motif to look for** rather
    than one nobody has found yet, and the census's plain silence would read as *nothing here
    knows this gene*.

    It subclasses :class:`GeneNotAssessedError` because that is literally true — no TF census
    assessed this gene — so an ``except`` clause written before this error existed keeps
    covering every gene it covered. The is-a is about the censuses and not about biology: a
    cofactor is not a kind of transcription factor, and being one never suppresses the motifs
    a census already reached, which is why the census is asked first (ADR-0016).

    The message names the census that did not assess the gene, the publisher that lists it as
    a cofactor, and that there is no motif here to look for.

    Examples
    --------
    >>> try:
    ...     motif_links("WDR5", "Homo sapiens")
    ... except TranscriptionCofactorError as error:
    ...     print("transcription cofactor" in str(error))
    True
    """


class VersionedGeneIdError(ValueError):
    """A versioned gene id was passed where a **Gene id stem** is the key.

    Not an absence — the gene is assessed and its links are here — so it is a
    :class:`ValueError` and a caller catching :class:`LookupError` for a missing gene does
    not swallow it. A stem may name more than one gene id in one **Annotation**:
    ``ENSG00000182378.14`` and ``ENSG00000182378.14_PAR_Y`` are two genes of one stem in
    ``gencode_v50lift37``, and the census reached one verdict for the stem. Answering a
    versioned id would therefore answer for the stem — which names a gene the caller did
    not — so it is refused, in the same spirit as
    :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`, which answers a stem with
    *every* gene id it names and never picks one.

    The message names the stem to pass instead.

    Examples
    --------
    >>> try:
    ...     motif_links("ENSG00000177606.6", "Homo sapiens")
    ... except VersionedGeneIdError as error:
    ...     print("ENSG00000177606" in str(error))
    True
    """


@dataclass(frozen=True)
class MotifLink:
    """One **Motif link**: one JASPAR profile that answers for one **TF gene**.

    One row of a shipped table, read back and frozen. It says what the matrix is a motif
    *of* — this gene alone, or a complex and which partners — and carries everything a
    caller needs to re-sort the answer on: the profile's tax ids, its **Cross-species
    link** flag and its total **Information content**.

    Attributes
    ----------
    release : str
        The JASPAR **Release** this link was built from. On the row rather than left to
        the file name, so two tables concatenate into one frame that still says which
        release each row came from.
    species : str
        The gene's species, as the assembly metadata table and its census spell it.
    gene_id_stem : str
        The **Gene id stem** the census is keyed by, and what this table is keyed by too.
    symbol : str
        The **census's** own symbol for the gene, never JASPAR's. The two differ exactly
        where the shipped alias table says they do — Lambert's ``T`` is JASPAR's ``TBXT``
        — and a row that mixed them would be unreadable.
    motif_id : str
        The **Motif id**, versioned: ``MA0099.4``.
    motif_name : str
        The **Motif name** JASPAR publishes, in JASPAR's own spelling and case, with
        ``::`` between the genes a complex names.
    role : str
        :data:`MONOMER` where the profile names one gene, :data:`COMPLEX` otherwise.
    partners : tuple of str
        The other genes the **Motif name** names, upper-cased as the name spells them.
        Empty for a monomer, and never empty for a complex.
    motif_tax_ids : tuple of str
        The NCBI taxonomy ids JASPAR files this profile under, ascending. Empty for a
        profile the release records no species for — ``MA0108``, TBP, is the one.
    is_cross_species : bool
        Whether the matrix was measured on a vertebrate other than the gene's own species
        (ADR-0013). A profile with no recorded species is marked ``True``: the row cannot
        claim a species match it has no evidence for.
    total_information_content : float
        The matrix's **Information content** summed over its columns, in bits. The third
        key of **Attribution specificity**, and not a quality score.
    rank : int
        This link's place under **Attribution specificity** among *all* of this gene's
        links, dense from 1 in the shipped table. It is the row's own number and not the
        answer's position, so filtering out **Cross-species link**s leaves gaps in it —
        deliberately, since a rank that renumbered itself would no longer say how
        specific the attribution was.

    Examples
    --------
    >>> link = motif_links("JUN", "Homo sapiens")[6]
    >>> link.motif_id, link.motif_name, link.role, link.partners
    ('MA0099.4', 'FOS::JUN', 'complex', ('FOS',))
    >>> link.is_cross_species, link.rank
    (False, 7)
    """

    release: str
    species: str
    gene_id_stem: str
    symbol: str
    motif_id: str
    motif_name: str
    role: str
    partners: tuple[str, ...]
    motif_tax_ids: tuple[str, ...]
    is_cross_species: bool
    total_information_content: float
    rank: int

    @property
    def is_complex(self) -> bool:
        """Return whether this profile is a motif of a complex rather than of this gene alone.

        Examples
        --------
        >>> [link.is_complex for link in motif_links("CTCF", "Homo sapiens")]
        [False, False, False]
        """
        return self.role == COMPLEX


@dataclass(frozen=True)
class MotifLinks:
    """One gene's **Motif link**s, in **Attribution specificity** order.

    :func:`motif_links`' answer. It is a sequence of links — iterate it, index it, take
    its length — that also says where it was cut from: the species, the **Release**, the
    **Tax group** and the shipped file. **That provenance is captured before any
    filtering**, for the reason
    :class:`~genome.tf.motif.jaspar.JasparDatabase` hands back a plain
    :class:`~genome.tf.motif.motif.MotifSet` when it is filtered: what comes out of a
    filter is no longer the release it came from, so unless it was written down first
    there is nothing left to say which release it was. Drop every **Cross-species link**
    from mouse ``Ctcf`` and no link survives to carry the release on its row — and this
    still says which release found nothing.

    **An empty answer is a real answer**, and it is never how absence is spelled: a
    species or release with no table raises :class:`NoMotifLinkTableError` and a gene no
    census assessed raises :class:`GeneNotAssessedError`. What is left is a gene with no
    links, and :attr:`is_tf` says which of the two kinds it is — a gene the census turned
    down, which receives no links by design, or one it judged a transcription factor that
    this **Release** has no profile for. 763 of human's 1,639 assessed-positive genes are
    the second kind on the 2026 release.

    Attributes
    ----------
    species : str
        The species, as its census spells it.
    release : str
        The JASPAR **Release** the table was built from.
    tax_group : str
        The **Tax group** the table covers — :data:`LINK_TAX_GROUP`. Recorded here
        because it is the one key of the three that no row carries.
    source : str
        The shipped file these links were read from.
    gene_id_stem : str
        The **Gene id stem** the gene was resolved to, whichever spelling was asked for.
    symbol : str
        The census's own symbol for that gene.
    is_tf : bool
        The census's own verdict. ``True`` for every gene that receives links, since only
        assessed-positive genes do; ``False`` says an empty answer is empty by design.
    links : tuple of MotifLink
        The links, most specifically attributable first, after any filtering asked for.

    Examples
    --------
    >>> jun = motif_links("JUN", "Homo sapiens", release="2024")
    >>> jun.species, jun.release, jun.gene_id_stem
    ('Homo sapiens', '2024', 'ENSG00000177606')
    >>> jun.motif_ids[:2]
    ('MA0488.2', 'MA0489.3')
    >>> matched = motif_links("JUN", "Homo sapiens", release="2024", cross_species=False)
    >>> matched.motif_ids[:2]
    ('MA0488.2', 'MA1131.2')
    >>> smad2 = motif_links("SMAD2", "Homo sapiens")
    >>> len(smad2), smad2.is_tf
    (0, False)
    """

    species: str
    release: str
    tax_group: str
    source: str
    gene_id_stem: str
    symbol: str
    is_tf: bool
    links: tuple[MotifLink, ...]

    def __len__(self) -> int:
        """Return how many links this answer holds.

        Examples
        --------
        >>> len(motif_links("TP53", "Homo sapiens"))
        1
        """
        return len(self.links)

    def __iter__(self) -> Iterator[MotifLink]:
        """Iterate the links, most specifically attributable first.

        Examples
        --------
        >>> [link.role for link in motif_links("AHR", "Homo sapiens")]
        ['complex']
        """
        return iter(self.links)

    def __getitem__(self, index: int) -> MotifLink:
        """Return one link by its position in this answer.

        Examples
        --------
        >>> motif_links("CTCF", "Homo sapiens")[0].motif_id
        'MA1930.2'
        """
        return self.links[index]

    @property
    def motif_ids(self) -> tuple[str, ...]:
        """Return every **Motif id** here, in **Attribution specificity** order.

        What is handed to a scan: the matrices to pull out of a
        :class:`~genome.tf.motif.jaspar.JasparDatabase` of the same **Release**.

        Examples
        --------
        >>> motif_links("AHR", "Homo sapiens").motif_ids
        ('MA0006.2',)
        """
        return tuple(link.motif_id for link in self.links)


@dataclass(frozen=True)
class MotifLinkTable:
    """One shipped **Motif link** table: one species, one **Release**, every gene in it.

    What one file says, read back and frozen. It is the whole join for that species and
    release — the layer :meth:`links_for` cuts one gene's answer out of — and the thing
    that knows which release and **Tax group** it is, which is why a
    :class:`MotifLinks` copies that down before it filters.

    Attributes
    ----------
    species : str
        The species, as every row of the file names it.
    release : str
        The JASPAR **Release**, as every row of the file names it.
    tax_group : str
        The **Tax group** the file covers — :data:`LINK_TAX_GROUP`, which no row carries
        because one value is all that ships.
    source : str
        Where the bytes came from.
    links : tuple of MotifLink
        Every link in the file, in file order: by **Gene id stem**, then by rank.

    Examples
    --------
    >>> table = motif_link_table("Homo sapiens", "2026")
    >>> table.release, table.tax_group, len(table)
    ('2026', 'vertebrates', 1085)
    >>> len(table.gene_id_stems)
    876
    >>> table.links_for("TP53").motif_ids
    ('MA0106.3',)
    """

    species: str
    release: str
    tax_group: str
    source: str
    links: tuple[MotifLink, ...]

    def __len__(self) -> int:
        """Return how many links the table holds.

        Examples
        --------
        >>> len(motif_link_table("Mus musculus", "2026"))
        896
        """
        return len(self.links)

    @property
    def gene_id_stems(self) -> tuple[str, ...]:
        """Return every **Gene id stem** the table links, once each, in file order.

        The genes this release has a motif for, which is a strict subset of the genes the
        census judged transcription factors — the rest are assessed positive and have no
        JASPAR profile.

        Examples
        --------
        >>> motif_link_table("Homo sapiens", "2026").gene_id_stems[0]
        'ENSG00000001167'
        """
        return tuple(dict.fromkeys(link.gene_id_stem for link in self.links))

    def frame(self) -> pd.DataFrame:
        """Return the table as a fresh :class:`~pandas.DataFrame`, one row per link.

        Built for the caller each time it is asked for, so mutating it cannot reach the
        cached table behind it. The columns are the file's own twelve in file order, with
        the multi-value cells as tuples and the flag as a boolean.

        Returns
        -------
        pandas.DataFrame
            The links, indexed from zero in file order.

        Examples
        --------
        >>> frame = motif_link_table("Homo sapiens", "2026").frame()
        >>> list(frame.columns) == list(LINK_COLUMNS)
        True
        >>> int(frame["is_cross_species"].sum())
        162
        """
        rows = [[getattr(link, column) for column in LINK_COLUMNS] for link in self.links]
        return pd.DataFrame(rows, columns=pd.Index(LINK_COLUMNS))

    def links_for(self, gene: str, *, cross_species: bool = True) -> MotifLinks:
        """Return one gene's **Motif link**s, most specifically attributable first.

        The gene is named by its **Gene id stem** or by the symbol its own census
        publishes, and a versioned gene id is refused rather than stemmed — see
        :func:`motif_links`, which is this method with the table looked up for you and is
        what a caller normally holds.

        Parameters
        ----------
        gene : str
            A **Gene id stem** — ``"ENSG00000177606"`` — or the census's own symbol for
            it, in any case: ``"JUN"``, ``"Jun"``.
        cross_species : bool, default True
            Whether to keep links whose profile was measured on another vertebrate. Pass
            ``False`` for a question that demands species-matched profiles; it can empty
            the answer, and for mouse it usually thins it (ADR-0013).

        Returns
        -------
        MotifLinks
            The links, in **Attribution specificity** order, carrying this table's
            provenance — copied down *before* the filter, since what a filter returns is
            no longer the release it came from.

        Raises
        ------
        GeneNotAssessedError
            If this species' census never assessed that gene.
        TranscriptionCofactorError
            If it never assessed that gene and a publisher lists it as a **Transcription
            cofactor**; a narrowing of the error above.
        VersionedGeneIdError
            If ``gene`` is a versioned gene id whose stem the census does assess.
        MotifLinkTableError
            If no census ships for this table's species, which no shipped table can be in
            — a link table is built from a census.

        Examples
        --------
        >>> table = motif_link_table("Mus musculus", "2026")
        >>> table.links_for("Ctcf").motif_ids
        ('MA1930.2', 'MA1929.2', 'MA0139.2')
        >>> table.links_for("Ctcf", cross_species=False).motif_ids
        ()
        """
        census = tf_gene_table(self.species)
        if census is None:
            raise MotifLinkTableError(
                f"{self.source} links genes of {self.species!r} and no census ships for that "
                f"species, so nothing says which of its genes are transcription factors. A link "
                f"table is built from a census — re-run scripts/build_tf_links.py, or drop the "
                f"table."
            )
        stem, symbol, is_tf = _resolve_gene(census, gene)
        # Provenance first, and the filter second. A filtered answer records no release or
        # tax group of its own, so what it was cut from is written down before anything is
        # dropped — the move `JasparDatabase` makes by handing back a plain `MotifSet`.
        found = sorted(
            (link for link in self.links if link.gene_id_stem == stem), key=lambda link: link.rank
        )
        if not cross_species:
            found = [link for link in found if not link.is_cross_species]
        return MotifLinks(
            species=self.species,
            release=self.release,
            tax_group=self.tax_group,
            source=self.source,
            gene_id_stem=stem,
            symbol=symbol,
            is_tf=is_tf,
            links=tuple(found),
        )


@cache
def shipped_link_tables() -> tuple[tuple[str, str], ...]:
    """Return every **Motif link** table that ships, as ``(species slug, release)``.

    What can be asked about at all, and what an error names when something cannot be. The
    directory is enumerated and the two keys read out of each file name, so neither the
    species nor the releases are listed in code and adding either is dropping a file in.

    Returns
    -------
    tuple of (str, str)
        One pair per shipped table, sorted. Empty only if the package ships no table.

    Examples
    --------
    >>> ("homo_sapiens", "2026") in shipped_link_tables()
    True
    """
    directory = files("genome").joinpath(LINK_SUBDIR)
    if not directory.is_dir():
        return ()
    found = []
    for entry in directory.iterdir():
        if not entry.name.endswith(LINK_SUFFIX):
            continue
        slug, _, release = entry.name[: -len(LINK_SUFFIX)].rpartition(".")
        if slug and release.startswith(RELEASE_PREFIX):
            found.append((slug, release.removeprefix(RELEASE_PREFIX)))
    return tuple(sorted(found))


@cache
def motif_link_table(
    species: str, release: str = DEFAULT_RELEASE, tax_group: str = DEFAULT_TAX_GROUP
) -> MotifLinkTable | None:
    """Return the **Motif link** table shipped for one species and **Release**, or ``None``.

    The raw absence, and the only place ``None`` is an acceptable way to say it: this is
    the layer below the one a caller touches, and :func:`motif_links` above it turns the
    ``None`` into an error naming what does ship, so that *no table was built for this*
    can never be read as *this gene has no motifs*.

    The species is slugged and then looked up among what :func:`shipped_link_tables`
    found, rather than joined onto the resource directory, so a name shaped like a path
    finds nothing instead of walking out of it. Read once per table and cached; everything
    it returns is frozen.

    Parameters
    ----------
    species : str
        The species, either as the assembly metadata table spells it (``"Homo sapiens"``)
        or as its slug (``"homo_sapiens"``).
    release : str, default ``"2026"``
        The JASPAR **Release** to read the links of. The motif half's default, so a new
        analysis links against the same release a fresh scan loads.
    tax_group : str, default ``"vertebrates"``
        The **Tax group**. Only :data:`LINK_TAX_GROUP` ships; any other answers ``None``,
        which is absence and not emptiness.

    Returns
    -------
    MotifLinkTable or None
        The table, or ``None`` when none ships for that species, release and tax group.

    Raises
    ------
    MotifLinkTableError
        If a table ships and cannot be read, or names a species or release its file name
        does not; the message names the file.

    Examples
    --------
    >>> motif_link_table("Homo sapiens", "2024").release
    '2024'
    >>> motif_link_table("homo_sapiens") == motif_link_table("Homo sapiens")
    True
    >>> motif_link_table("Danio rerio") is None
    True
    >>> motif_link_table("Homo sapiens", tax_group="plants") is None
    True
    """
    slug = species_slug(species)
    if tax_group != LINK_TAX_GROUP or (slug, release) not in shipped_link_tables():
        return None
    resource = files("genome").joinpath(
        LINK_SUBDIR, f"{slug}.{RELEASE_PREFIX}{release}{LINK_SUFFIX}"
    )
    origin = str(resource)
    table = parse_motif_link_table(
        gzip.decompress(resource.read_bytes()).decode("utf-8"), source=origin
    )
    if species_slug(table.species) != slug or table.release != release:
        raise MotifLinkTableError(
            f"{origin} is named for {slug!r} and release {release!r} and its rows say "
            f"{species_slug(table.species)!r} and {table.release!r}. The file name is what a "
            f"table is found by, so a file whose rows disagree with it would answer a question "
            f"nobody asked — re-run scripts/build_tf_links.py for that species and release."
        )
    return table


def parse_motif_link_table(text: str, *, source: str) -> MotifLinkTable:
    r"""Read one **Motif link** table's text into links, holding it to what a table promises.

    A pure function from text to links: it opens nothing, downloads nothing and
    decompresses nothing. Separate from the resource it came out of, as
    :func:`~genome.tf.motif.jaspar.parse_transfac` is, so every way a file can be wrong is
    reachable without writing a broken one into the package. The shipped tables are
    gzipped and :func:`motif_link_table` unpacks them at the resource boundary, which is
    where :func:`~genome.tf.gene.tf_gene_table` unpacks a census too — the seam is between
    the bytes and the format, and it does not move because the bytes are compressed.

    The **Release** and the species are read off the rows rather than passed in, and every
    row must agree about both — they are on each row so that two tables concatenate into
    one frame that still says where each row came from, which is a promise only if they
    are uniform within a file. The **Tax group** is the one key no row carries, since
    :data:`LINK_TAX_GROUP` is all that ships.

    Parameters
    ----------
    text : str
        The whole table, unpacked. These are a few hundred kilobytes, so they are read
        whole.
    source : str
        Where the text came from; named in every message, since regenerating that file is
        the only repair.

    Returns
    -------
    MotifLinkTable
        The table the text spells, in file order.

    Raises
    ------
    MotifLinkTableError
        If the header is not :data:`LINK_COLUMNS`, a row holds the wrong number of cells,
        a **Role** or a flag is spelled a way no table spells one, a number is not one, a
        key cell is blank, the file declares no links, or two rows name different releases
        or species.

    Examples
    --------
    >>> header = "\t".join(LINK_COLUMNS)
    >>> row = "2026\tHomo sapiens\tENSG00000141510\tTP53\tMA0106.3\tTP53\tmonomer\t\t9606\tno\t20.6607\t1"
    >>> table = parse_motif_link_table(f"{header}\n{row}\n", source="one.tsv")
    >>> table.species, table.release, table.links[0].motif_name
    ('Homo sapiens', '2026', 'TP53')
    """
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise MotifLinkTableError(
            f"{source} is empty. A link table carries a header and at least one link; absence "
            f"is spelled by shipping no file at all, and an empty one is the second spelling of "
            f"it that would read as *no motif answers for any of these genes*. Re-run "
            f"scripts/build_tf_links.py for that species and release, or remove the file."
        )
    header = tuple(lines[0].split("\t"))
    if header != LINK_COLUMNS:
        raise MotifLinkTableError(
            f"{source} carries the columns {list(header)} where every link table carries "
            f"{list(LINK_COLUMNS)}. Identical headers are what let two tables concatenate into "
            f"one frame — re-run scripts/build_tf_links.py, which writes them."
        )
    links = tuple(
        _read_link(line, number, source=source)
        for number, line in enumerate(lines[1:], start=2)
        if line
    )
    if not links:
        raise MotifLinkTableError(
            f"{source} declares columns and no links. A table linking nothing says no more than "
            f"an absent file does — re-run scripts/build_tf_links.py for that species and "
            f"release, or remove the file."
        )
    return MotifLinkTable(
        species=_one_value(links, "species", source=source),
        release=_one_value(links, "release", source=source),
        tax_group=LINK_TAX_GROUP,
        source=source,
        links=links,
    )


def motif_links(
    gene: str,
    species: str,
    *,
    release: str = DEFAULT_RELEASE,
    tax_group: str = DEFAULT_TAX_GROUP,
    cross_species: bool = True,
) -> MotifLinks:
    """Return the JASPAR motifs that answer for one **TF gene**, most specific first.

    The entry point. It reads the shipped table for that species and **Release** and cuts
    one gene's links out of it — nothing is downloaded, nothing on disk is touched, and no
    **Assembly** is involved.

    **What names a gene here.** A **Gene id stem**, which is what the tables and the
    censuses are keyed by, or the symbol the gene's *own census* publishes, in any case —
    the two censuses spell one factor ``JUN`` and ``Jun``, and each spelling is unique
    within its own census. A versioned gene id is **refused** rather than stemmed: a stem
    may name more than one gene id in one **Annotation**, so answering ``ENSG00000182378.14``
    would answer for a stem that also names ``ENSG00000182378.14_PAR_Y``, and this package
    never picks a gene the caller did not name (see
    :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`, which crosses that gap in
    the other direction). Pass the stem, and the error says which one.

    **The species is passed and never inferred.** A table is named by a species and a
    release, and a **Gene id stem**'s prefix is not a claim about which species it belongs
    to — deriving one from the other is the guess ADR-0003 exists to forbid. A caller
    holding an assembly has its species already, in the assembly's own metadata row.

    **The census is asked first, and the order is what keeps this correct.** A gene the
    census assessed is answered whatever else is known about it — the 151 human genes that
    are both a **TF gene** and a **Transcription cofactor**, TBP and KMT2A and DNMT1 among
    them, come back exactly as they always did, because a second table must never suppress
    an answer the census already reached. Only then is that species' **Cofactor table**
    asked, and a gene it lists raises :class:`TranscriptionCofactorError`: a cofactor
    recognises no sequence of its own, so *no motif to look for* is a truer answer than the
    census's silence. A gene neither knows raises :class:`GeneNotAssessedError`, as it
    always has, and so does every gene of a species that ships no cofactor table.

    Parameters
    ----------
    gene : str
        A **Gene id stem** or the census's own symbol for the gene.
    species : str
        The species, as the assembly metadata table spells it or as its slug.
    release : str, default ``"2026"``
        The JASPAR **Release** to link against. Both releases the package prepares have
        tables; asking for another raises and names the ones that ship.
    tax_group : str, default ``"vertebrates"``
        The **Tax group**. Only :data:`LINK_TAX_GROUP` ships, and asking for another
        raises rather than answering emptily.
    cross_species : bool, default True
        Whether to keep links whose profile was measured on another vertebrate (ADR-0013).
        ``False`` is how a question that demands species-matched profiles asks for them.

    Returns
    -------
    MotifLinks
        The gene's links in **Attribution specificity** order, carrying the release, tax
        group and file they were cut from.

    Raises
    ------
    NoMotifLinkTableError
        If no table ships for that species, release or tax group; the message names what
        does.
    GeneNotAssessedError
        If that species' census never assessed that gene.
    TranscriptionCofactorError
        If it never assessed that gene and a publisher lists it as a **Transcription
        cofactor**; a narrowing of the error above, so an ``except`` clause written for
        that one catches this too.
    VersionedGeneIdError
        If ``gene`` is a versioned gene id whose stem the census does assess.

    Examples
    --------
    >>> ctcf = motif_links("CTCF", "Homo sapiens")
    >>> ctcf.motif_ids
    ('MA1930.2', 'MA1929.2', 'MA0139.2')
    >>> jun = motif_links("Jun", "Mus musculus")
    >>> [(link.motif_id, link.is_cross_species) for link in jun][:2]
    [('MA0489.3', False), ('MA0488.2', True)]
    >>> motif_links("Jun", "Mus musculus", cross_species=False).motif_ids
    ('MA0489.3',)
    >>> motif_links("T", "Homo sapiens").links[0].motif_name
    'TBXT'
    >>> len(motif_links("FOXM1", "Homo sapiens", release="2024"))
    0
    """
    table = motif_link_table(species, release, tax_group)
    if table is None:
        raise NoMotifLinkTableError(_nothing_ships(species, release, tax_group))
    return table.links_for(gene, cross_species=cross_species)


def _nothing_ships(species: str, release: str, tax_group: str) -> str:
    """Return the message for a species, release or tax group no table was built for."""
    shipped = shipped_link_tables()
    species_names = ", ".join(sorted({slug for slug, _ in shipped})) or "none"
    releases = ", ".join(sorted({found for _, found in shipped})) or "none"
    return (
        f"no motif link table ships for {species_slug(species)!r}, release {release!r}, tax "
        f"group {tax_group!r}. The tables ship for {species_names}, for the releases "
        f"{releases}, and for the tax group {LINK_TAX_GROUP!r} alone — JASPAR's CORE "
        f"{LINK_TAX_GROUP} is where every censused species' profiles are. Ask for one of "
        f"those, and read *no table* as nobody having built one rather than as this gene "
        f"having no motifs."
    )


def _read_link(line: str, number: int, *, source: str) -> MotifLink:
    """Read one link's row, holding every cell to what a shipped table promises."""
    cells = line.split("\t")
    if len(cells) != len(LINK_COLUMNS):
        raise MotifLinkTableError(
            f"{source} line {number} holds {len(cells)} cells where the header declares "
            f"{len(LINK_COLUMNS)}. A link table is a plain TSV with no quoting, so a cell "
            f"carrying a tab is a defect in the generator rather than something to parse "
            f"around — re-run scripts/build_tf_links.py."
        )
    row = dict(zip(LINK_COLUMNS, cells, strict=True))
    for column in ("release", "species", "gene_id_stem", "symbol", "motif_id", "motif_name"):
        if not row[column]:
            raise MotifLinkTableError(
                f"{source} line {number} leaves the {column!r} column blank. Every link names "
                f"the gene it answers for, the profile that answers, and the release and "
                f"species it was built from — re-run scripts/build_tf_links.py."
            )
    if row["role"] not in (MONOMER, COMPLEX):
        raise MotifLinkTableError(
            f"{source} line {number} spells its role {row['role']!r}, and a link is {MONOMER!r} "
            f"or {COMPLEX!r}. The role is what keeps a heterodimer matrix from being read as a "
            f"monomer's, so a third spelling is a defect and not a new kind of link — re-run "
            f"scripts/build_tf_links.py."
        )
    partners = _split(row["partners"])
    if (row["role"] == COMPLEX) != bool(partners):
        raise MotifLinkTableError(
            f"{source} line {number} is a {row['role']} naming the partners {list(partners)}. A "
            f"complex names at least one partner and a monomer names none — a complex with none "
            f"beside it is the monomer reading with a label on it. Re-run "
            f"scripts/build_tf_links.py."
        )
    if row["is_cross_species"] not in (TRUE_CELL, FALSE_CELL):
        raise MotifLinkTableError(
            f"{source} line {number} spells its cross-species flag "
            f"{row['is_cross_species']!r}, and a flag is spelled {TRUE_CELL!r} or "
            f"{FALSE_CELL!r}, as a census spells its own. That flag is the only thing a caller "
            f"needing species-matched profiles can filter on — re-run scripts/build_tf_links.py."
        )
    return MotifLink(
        release=row["release"],
        species=row["species"],
        gene_id_stem=row["gene_id_stem"],
        symbol=row["symbol"],
        motif_id=row["motif_id"],
        motif_name=row["motif_name"],
        role=row["role"],
        partners=partners,
        motif_tax_ids=_split(row["motif_tax_ids"]),
        is_cross_species=row["is_cross_species"] == TRUE_CELL,
        total_information_content=_bits(row["total_information_content"], number, source=source),
        rank=_rank(row["rank"], number, source=source),
    )


def _split(cell: str) -> tuple[str, ...]:
    """Split one multi-value cell on its separator, dropping what says nothing."""
    return tuple(part for part in cell.split(VALUE_SEPARATOR) if part)


def _bits(cell: str, line: int, *, source: str) -> float:
    """Read one total **Information content** cell, or say which line is not a number."""
    try:
        return float(cell)
    except ValueError as error:
        raise _not_a_number(cell, "total_information_content", line, source) from error


def _rank(cell: str, line: int, *, source: str) -> int:
    """Read one rank cell, or say which line is not a whole number."""
    try:
        return int(cell)
    except ValueError as error:
        raise _not_a_number(cell, "rank", line, source) from error


def _not_a_number(cell: str, column: str, line: int, source: str) -> MotifLinkTableError:
    """Return the error a numeric cell that is not a number raises."""
    return MotifLinkTableError(
        f"{source} line {line} holds {cell!r} in the {column!r} column, which is not a number. "
        f"Re-run scripts/build_tf_links.py, which writes it."
    )


def _one_value(links: tuple[MotifLink, ...], field: str, *, source: str) -> str:
    """Return the one value every link carries for ``field``, or say which two disagree."""
    values = {str(getattr(link, field)) for link in links}
    if len(values) != 1:
        raise MotifLinkTableError(
            f"{source} names {len(values)} values in its {field!r} column: "
            f"{sorted(values)}. One file is one species and one release — the column is on "
            f"every row so that two tables concatenate and still say where each row came from, "
            f"which promises nothing unless it is uniform within a file. Re-run "
            f"scripts/build_tf_links.py for the species and release the file is named for."
        )
    return values.pop()


#: One census indexed by **Gene id stem** and by upper-cased symbol, per species. Built on
#: first use and kept, because a caller asking a few hundred genes for their links would
#: otherwise walk the whole census once per gene. Both indexes are read-only.
_CENSUS_INDEXES: dict[str, tuple[Mapping[str, tuple[str, bool]], Mapping[str, str]]] = {}


def _census_index(census: TFGeneTable) -> tuple[Mapping[str, tuple[str, bool]], Mapping[str, str]]:
    """Return one census indexed by **Gene id stem** and by upper-cased symbol."""
    slug = species_slug(census.species)
    index = _CENSUS_INDEXES.get(slug)
    if index is None:
        assessed = [(row[0] or "", row[1] or "", row[2] == TRUE_CELL) for row in census.rows]
        index = (
            MappingProxyType({stem: (symbol, flag) for stem, symbol, flag in assessed}),
            MappingProxyType({symbol.upper(): stem for stem, symbol, _ in assessed if symbol}),
        )
        _CENSUS_INDEXES[slug] = index
    return index


def _resolve_gene(census: TFGeneTable, gene: str) -> tuple[str, str, bool]:
    """Return the stem, the census's symbol and the census's verdict for one named gene."""
    by_stem, by_symbol = _census_index(census)
    stem = gene if gene in by_stem else by_symbol.get(gene.upper())
    if stem is not None:
        symbol, is_tf = by_stem[stem]
        return stem, symbol, is_tf
    versioned = _VERSIONED_GENE_ID.fullmatch(gene)
    if versioned is not None and versioned["stem"] in by_stem:
        raise VersionedGeneIdError(
            f"{gene!r} is a versioned gene id and this is keyed by gene id stem. Pass "
            f"{versioned['stem']!r}: a stem may name more than one gene id in one annotation — "
            f"eight of gencode_v50lift37's nine are a pseudoautosomal gene and its _PAR_Y copy — "
            f"so answering the versioned id would answer for the stem, which names a gene you "
            f"did not."
        )
    known = versioned["stem"] if versioned is not None else gene
    # The census is asked first and this second, so a gene it assessed is answered whatever
    # any other table says about it — the 151 human genes that are both a TF gene and a
    # cofactor never reach here. A species with no cofactor table finds nobody and falls
    # through to the error below, exactly as every species did before this table existed.
    publishers = _cofactor_publishers(census.species, known)
    if publishers is not None:
        raise TranscriptionCofactorError(
            f"{census.provenance.publisher} never assessed {known!r} as a transcription factor, "
            f"and it is listed as a transcription cofactor by {publishers}. A cofactor acts on "
            f"transcription without recognising a sequence of its own, so it has no motif and "
            f"none is missing here — which is a different answer from a gene nobody has "
            f"assessed. Stop looking for a motif, and see "
            f"genome.tf.cofactor.cofactor_table({census.species!r}) for the class it is listed "
            f"under."
        )
    raise GeneNotAssessedError(
        f"{census.provenance.publisher} never assessed {known!r}, so nothing here says whether "
        f"it is a transcription factor and no motif can answer for it. That census speaks for "
        f"{census.species!r} and is keyed by gene id stem, with its own symbol beside each — "
        f"check the species you named, and see genome.tf.gene.tf_gene_table({census.species!r}) "
        f"for what it does assess. Absent from a census is not the same as judged not a TF."
    )


def _cofactor_publishers(species: str, gene: str) -> str | None:
    """Return who lists that gene as a **Transcription cofactor**, or ``None`` if nobody does."""
    table = cofactor_table(species)
    if table is None:
        return None
    wanted = gene.upper()
    for row in table.rows:
        named = row[0] == gene or (row[1] or "").upper() == wanted
        # A publisher recording a rejection is not a publisher listing the gene, so the flag
        # is read rather than the row's presence. It says `yes` on every row that ships today.
        if named and row[_COFACTOR_FLAG] == TRUE_CELL:
            return _publishers_of(table, row[_COFACTOR_SOURCE] or "")
    return None


def _publishers_of(table: CofactorTable, source: str) -> str:
    """Return the publishers one row's source cell names, as a phrase to print."""
    # A row spelled `both` names every publisher of that table that lists genes, which is
    # every source it cites but one supplying identifiers alone — HGNC lists nobody.
    listed = [
        entry.publisher
        for entry in table.provenance.sources
        if entry.source == source or (source == BOTH and entry.source in SOURCES)
    ]
    return " and ".join(listed) or source
