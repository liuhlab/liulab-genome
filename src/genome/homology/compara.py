"""One Ensembl Compara **Release** on disk, sliced to a species pair you can query.

I/O boundary module: name two species and the per-species dump that holds that pair is
fetched once into the ``homology/`` subtree of the **Data dir**, verified against the
publisher's own md5, sliced to the pair and recorded with a **Completion marker**, then
read back as a :class:`HomologySet`. A second construction re-reads what is there and
fetches nothing — the shape :class:`~genome.tf.motif.jaspar.JasparDatabase` established.
The types it hands back are :mod:`genome.io.results`'s, which read no file.

The slice is stored as a plain gzipped TSV carrying the publisher's own header and the
publisher's own rows, unedited, so a collaborator can read it in R or in a shell without
this package. Every value in a **Homology link** is a cell of that file and none of them
is computed here: this package publishes no quality score, no ranking and no "best
ortholog" of its own, and it never uses homology to build a table it ships (ADR-0019).

**The partition guard is why this module exists.** Compara's per-species files are a
de-duplicated partition *at the pair level* — the README says each one holds "an arbitrary
subset of orthologies involving the given genome" — so a pair lives in exactly one of its
two files and which one is not promised stable across releases. Measured on release 116:
the human file holds **0** human↔mouse rows and 23,982 human↔worm rows; the mouse file
holds 23,764 human rows and 25,006 worm rows; the worm file holds neither pair. Which file
holds a pair is written down in the shipped provenance table as a *measurement*, and it is
measured again on every prepare: a slice that comes back empty raises
:class:`ComparaPartitionError` naming the other file rather than answering nothing. A pair
is never *partially* present, which is what makes zero a trustworthy signal.

**On paralogy, and a measurement that matters.** A **Paralogy link** is kept and marked
rather than excluded, and ``paralogs=True`` is what returns one — but **release 116
publishes none for these pairs, so the switch changes nothing on it.** Counted over the
whole human dump (4.0 M rows, ~200 partner species): zero ``between_species_paralog``, and
every ``other_paralog`` (128,020), ``within_species_paralog`` (13,144) and ``gene_split``
(9) row relates two genes of *one* species. A **Homology link** relates two species, so
those same-species rows are not this pair's and are not in its set — the fixtures carry
real ones precisely to hold that boundary. The switch is the place a cross-species
duplication label would land the release Compara publishes one, and it is here rather than
added later so that *not an ortholog* has somewhere to be distinguishable from *absent*
(ADR-0013).

**On the quality scores.** Both are null on 100% of the rows of *either* worm pairing, not
only human↔worm. Which columns a set holds nothing in is therefore measured over the
prepared slice rather than listed against a pair, and it rides on every answer.

**On a fetch that lies.** A resumed download of one of these gzips has been observed to
pass ``gzip -t`` — "decompression OK, trailing garbage ignored" — with the wrong md5. The
publisher's checksum is checked against the bytes as they are fetched, and the slice's own
sha256 is recorded and re-checked on every read, because opening cleanly is not evidence.

Examples
--------
>>> from genome.homology.compara import HomologySet, compara_url
>>> compara_url("Homo sapiens", "116")                       # doctest: +ELLIPSIS
'https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/homo_sapiens/...'
>>> homologs = HomologySet("Homo sapiens", "Mus musculus")   # doctest: +SKIP
>>> homologs.homologs(["ENSG00000141510"]).links[0]          # doctest: +SKIP
HomologyLink(gene_id_stem='ENSG00000141510', ...)
"""

from __future__ import annotations

import gzip
import hashlib
import re
import shlex
from collections.abc import Iterable, Mapping
from pathlib import Path
from types import MappingProxyType

from genome.homology import metadata as metadata_mod
from genome.homology.metadata import HomologyMetadata
from genome.io import fetch
from genome.io.completion import (
    CompletionRecord,
    build_record,
    check_registration,
    clear_work_dir,
    work_dir,
    write_record,
)
from genome.io.registration import liulab_data_dir
from genome.io.results import HomologyAnswer, HomologyLink

# Imported rather than written a second time: the file-naming convention every
# shipped-data directory here uses belongs to none of them in particular, and for the
# three species this prepares it is also Ensembl's own genome name. Promoting it to a
# shared home is a one-line refactor whenever that is worth doing.
from genome.tf.gene.census import species_slug

#: Where Ensembl publishes its releases. The FTP tree is used and no REST or BioMart API
#: is: everything here is a bulk file fetched once and read locally, so a run is
#: reproducible and nothing fails intermittently against a remote service.
ENSEMBL_BASE_URL = "https://ftp.ensembl.org/pub"

#: Which of Compara's gene-tree dumps is read: the protein-coding trees of the ``default``
#: collection. The ncRNA dumps and the whole-clade collection dumps are deliberately not
#: read — they are different objects with different membership, not more rows of this one.
COMPARA_DUMP = "protein_default"

#: Subdirectory of the **Data dir** holding homology data, a *sibling* of the assembly
#: tree beside ``motif/``: a **Homology set** belongs to no **Assembly**, so there is no
#: assembly directory it could go under.
HOMOLOGY_SUBDIR = "homology"

#: Subdirectory of the homology tree holding Ensembl Compara's sets. One publisher, one
#: directory: a second homology source would get its own beside this.
COMPARA_SUBDIR = "ensembl_compara"

#: The **Release** used when none is named — the newest this package pins, so a new
#: analysis starts on current data and reproducing an old one is the call that says which.
DEFAULT_RELEASE = "116"

#: The columns of Compara's homology TSV, in the order the publisher writes them. The
#: slice keeps every one of them verbatim, so the file on disk is the publisher's rows
#: and nothing here has re-derived a cell.
COMPARA_COLUMNS: tuple[str, ...] = (
    "gene_stable_id",
    "protein_stable_id",
    "species",
    "identity",
    "homology_type",
    "homology_gene_stable_id",
    "homology_protein_stable_id",
    "homology_species",
    "homology_identity",
    "dn",
    "ds",
    "goc_score",
    "wga_coverage",
    "is_high_confidence",
    "homology_id",
)

#: The two confidence columns that carry a **score**, as opposed to Compara's boolean
#: high-confidence flag. Which of them a set holds no value in *anywhere* is measured when
#: it is prepared and reported on every answer, because both are null on **100% of** the
#: rows of *either* worm pairing — 23,982 human↔worm and 25,006 mouse↔worm on release 116 —
#: and a filter written against one would empty itself in silence. Measured rather than
#: listed by pair, so a release that starts scoring worm needs no code change.
QUALITY_SCORE_COLUMNS: tuple[str, ...] = ("goc_score", "wga_coverage")

#: How Compara spells a cell it recorded nothing in.
NULL_CELL = "NULL"

#: How Compara spells its high-confidence flag.
_TRUE_CELL, _FALSE_CELL = "1", "0"

#: What a **Completion marker** written here calls what it recorded.
RECORD_KIND = "homology"

#: A gene id carrying a version, and the stem inside it. Compara keys its dumps by stem,
#: so a versioned id would match nothing — the quietest failure in this landscape, and the
#: one refused by name rather than answered with an empty result.
_VERSIONED_GENE_ID = re.compile(r"(?P<stem>[^.]+)\.\d+\w*")

#: Where the query gene and the homologous gene sit in a Compara row.
_GENE, _SPECIES = COMPARA_COLUMNS.index("gene_stable_id"), COMPARA_COLUMNS.index("species")
_HOMOLOG = COMPARA_COLUMNS.index("homology_gene_stable_id")
_HOMOLOG_SPECIES = COMPARA_COLUMNS.index("homology_species")
_TYPE = COMPARA_COLUMNS.index("homology_type")
_HIGH_CONFIDENCE = COMPARA_COLUMNS.index("is_high_confidence")

#: Where each quality column sits in a row, looked up once rather than per cell.
_QUALITY_AT: Mapping[str, int] = MappingProxyType(
    {name: COMPARA_COLUMNS.index(name) for name in QUALITY_SCORE_COLUMNS}
)


class UnknownHomologySpeciesError(LookupError):
    """A species no **Homology set** is prepared for, so the question cannot be answered.

    A :class:`LookupError` and never an empty answer: *nobody pinned this species* must
    never read as *this species has no homologs*. The message names the species that do.

    Examples
    --------
    >>> try:
    ...     check_species("Danio rerio")
    ... except UnknownHomologySpeciesError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """


class NoHomologyPairError(LookupError):
    """No shipped row pins this species pair in this **Release**.

    Distinct from :class:`UnknownHomologySpeciesError`: both species are prepared, but not
    together in the release asked for. The message names the pairs that are.

    Examples
    --------
    >>> issubclass(NoHomologyPairError, LookupError)
    True
    """


class ComparaPartitionError(RuntimeError):
    """The file recorded as holding a species pair holds none of its rows.

    **The trap this module exists for.** Compara's per-species dumps are a de-duplicated
    partition at the pair level, so a pair lives in exactly one of its two files and the
    assignment is arbitrary and unstable across releases. A slice that comes back empty
    therefore means the partition moved, not that the two species share no homologs — and
    a pair is never *partially* present, which is what makes zero trustworthy. The message
    names the other file, which is where the pair now is.

    A :class:`RuntimeError` rather than a :class:`LookupError`: nothing the caller asked
    for is missing, the shipped provenance row is out of date.
    """


class HomologySetNotDownloadedError(RuntimeError):
    """The set is not prepared here, and the publisher's dump could not be fetched.

    A :class:`RuntimeError` rather than a :class:`LookupError`: nothing the caller asked
    for is missing, this machine could not reach the publisher. Fetching is the one step
    in this package that needs the network and the lab's CPU cluster compute nodes have
    none, so the message names the call to make on a login node instead of reporting a
    transport error and stopping there.
    :class:`genome.xref.xref.XrefSetNotDownloadedError` says the same thing for an **Xref
    set**; they are two classes because they are two contexts.

    Examples
    --------
    >>> issubclass(HomologySetNotDownloadedError, RuntimeError)
    True
    """


class ComparaFileError(ValueError):
    """A file read as Compara's is not one, or a prepared slice is not what was recorded.

    One class for both because they mean the same thing to a caller — what is on disk is
    not what it should be — and the repair is the same: delete the set's directory and
    construct it again. The message names the file and the repair.
    """


class VersionedGeneIdError(ValueError):
    """A versioned gene id was passed where a **Gene id stem** is the key.

    Compara keys its dumps by stem, so ``ENSG00000141510.18`` matches nothing and would
    ride back in ``unresolved`` looking exactly like a gene Compara never placed in a tree.
    Joining a versioned id to a bare one returning zero matches *in silence* is the most
    error-prone detail in this landscape, so it is refused by name instead. The message
    names the stem to pass.

    :class:`genome.tf.link.VersionedGeneIdError` refuses the same thing for the same
    reason in the TF context; they are two classes because they are two contexts, and
    nothing is expected to catch both.

    Examples
    --------
    >>> try:
    ...     check_stem("ENSG00000141510.18")
    ... except VersionedGeneIdError as error:
    ...     print("ENSG00000141510" in str(error))
    True
    """


def homology_data_dir() -> Path:
    """Return the directory holding homology data, which belongs to no assembly.

    ``<liulab_data>/homology/``, a **sibling** of the assembly tree beside ``motif/``: a
    **Homology set** is anchored to a species pair and a **Release**, so it names no
    **Assembly** and there is no per-assembly directory it could be filed under. Shared by
    every project on the machine. Nothing is created by asking.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/homology``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> homology_data_dir()
    PosixPath('/scratch/liulab/homology')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / HOMOLOGY_SUBDIR


def compara_url(species: str, release: str) -> str:
    """Return the URL of one **Release**'s homology dump for one species.

    Parameters
    ----------
    species : str
        The species whose per-species dump is wanted, in any spelling.
    release : str
        The Compara **Release**.

    Returns
    -------
    str
        The published file's URL.

    Notes
    -----
    **This builds a message, never a pin.** The shipped provenance table's own
    ``source_url`` is the authority for a set that is actually fetched; this is what names
    the *other* file of a pair when the partition has moved, where nothing is pinned to
    name. It is not a template a release can be added through: release 113 ships these
    dumps **uncompressed**, so its file has no ``.gz``, and only releases 90 and 116
    publish an ``MD5SUM`` at all — 91 to 112 publish no checksum of any kind. A new release
    is a measured row, not a formatted string.

    Ensembl's genome name is the species slug for each of the three species prepared here.
    It is not for every genome Ensembl carries — several are named for a subspecies or a
    collection — so a fourth species would bring its own name rather than this deriving
    one.

    Examples
    --------
    >>> compara_url("Mus musculus", "116")
    'https://ftp.ensembl.org/pub/release-116/tsv/ensembl-compara/homologies/mus_musculus/Compara.116.protein_default.homologies.tsv.gz'
    """
    genome = species_slug(species)
    return (
        f"{ENSEMBL_BASE_URL}/release-{release}/tsv/ensembl-compara/homologies/{genome}/"
        f"Compara.{release}.{COMPARA_DUMP}.homologies.tsv.gz"
    )


def homology_prepare_command(species: str, other_species: str, release: str) -> str:
    r"""Return the call that prepares one **Homology set**, for an error message to quote.

    One spelling of it, so a renamed entry point is renamed once. Quoted by the error a
    caller repairs by fetching the set on a machine with internet.

    Parameters
    ----------
    species : str
        One species of the pair, as the shipped table spells it.
    other_species : str
        The other.
    release : str
        The pinned **Release**.

    Returns
    -------
    str
        A shell command, unquoted and unfenced — the caller decides how to set it.

    Examples
    --------
    >>> homology_prepare_command("Homo sapiens", "Mus musculus", "116")
    'python -c "from genome.homology import HomologySet; HomologySet(\'Homo sapiens\', \'Mus musculus\', \'116\')"'
    """
    call = (
        f"from genome.homology import HomologySet; "
        f"HomologySet({species!r}, {other_species!r}, {release!r})"
    )
    return f'python -c "{call}"'


def pair_name(species: str, other_species: str) -> str:
    """Return the directory name one species pair is filed under, order-independent.

    The two slugs, sorted and joined with a hyphen, so the pair asked for either way round
    reaches one directory and one download.

    Parameters
    ----------
    species : str
        One species of the pair.
    other_species : str
        The other.

    Returns
    -------
    str
        The pair's directory name.

    Examples
    --------
    >>> pair_name("Mus musculus", "Homo sapiens")
    'homo_sapiens-mus_musculus'
    """
    return "-".join(sorted((species_slug(species), species_slug(other_species))))


def slice_filename(row: HomologyMetadata) -> str:
    """Return the name one pair's stored slice is cached under.

    Compara's own naming with the pair in place of the member type, so a directory listing
    says which release and which two species a file is without opening it.

    Parameters
    ----------
    row : genome.homology.metadata.HomologyMetadata
        The provenance row for the pair.

    Returns
    -------
    str
        The local file name.
    """
    return f"Compara.{row.release}.{pair_name(*row.pair)}.homologies.tsv.gz"


def source_filename(row: HomologyMetadata) -> str:
    """Return the name the publisher's own dump is downloaded under, inside the work area.

    Parameters
    ----------
    row : genome.homology.metadata.HomologyMetadata
        The provenance row for the pair.

    Returns
    -------
    str
        The download's file name, carrying the release and the species whose dump it is.
    """
    return (
        f"Compara.{row.release}.{species_slug(row.holding_species)}."
        f"{COMPARA_DUMP}.homologies.tsv.gz"
    )


def set_dir(root: Path, row: HomologyMetadata) -> Path:
    """Return the directory one pair's set is prepared in, under ``root``.

    ``<root>/ensembl_compara/<release>/<pair>/``: one directory per set, because each
    carries a **Completion marker** of its own, and releases sit side by side so holding
    two is not a re-download.

    Parameters
    ----------
    root : pathlib.Path
        The homology root — :func:`homology_data_dir`, or a ``cache_dir`` override.
    row : genome.homology.metadata.HomologyMetadata
        The provenance row for the pair.

    Returns
    -------
    pathlib.Path
        The set's own directory.
    """
    return root / COMPARA_SUBDIR / row.release / pair_name(*row.pair)


def check_species(species: str) -> str:
    """Return ``species`` as the shipped table spells it, else say which it prepares.

    Either spelling is accepted — the assembly metadata table's ``"Homo sapiens"`` or the
    slug ``"homo_sapiens"`` — and the table's own spelling comes back, so an answer names
    a species one way however it was asked for.

    Parameters
    ----------
    species : str
        A species name, in either spelling.

    Returns
    -------
    str
        The species as the shipped provenance table spells it.

    Raises
    ------
    UnknownHomologySpeciesError
        If no shipped row names that species.

    Examples
    --------
    >>> check_species("homo_sapiens")
    'Homo sapiens'
    """
    known = metadata_mod.homology_species()
    slug = species_slug(species)
    for name in known:
        if species_slug(name) == slug:
            return name
    raise UnknownHomologySpeciesError(
        f"no Homology set for {species!r}: this package prepares {', '.join(known)}. Nothing "
        f"here answers a cross-species question about another species rather than answering "
        f"it with a translated guess (ADR-0019) — ask about one of those, or add a row to "
        f"data/homology/homology_metadata.tsv measured against the release you want."
    )


def check_release(release: str) -> str:
    """Return ``release`` if the shipped table pins it, else say which it does.

    Examples
    --------
    >>> check_release("116")
    '116'
    """
    pinned = metadata_mod.homology_releases()
    if release not in pinned:
        raise ValueError(
            f"no Ensembl Compara release {release!r} is pinned here: this package prepares "
            f"{', '.join(pinned)}. Pass one of those as a string, or add its rows to "
            f"data/homology/homology_metadata.tsv with the publisher's own md5 beside each."
        )
    return release


def check_pair(species: str, other_species: str, release: str) -> HomologyMetadata:
    """Return the shipped provenance row for one species pair and **Release**.

    Parameters
    ----------
    species : str
        One species of the pair, in either spelling.
    other_species : str
        The other.
    release : str
        The Compara **Release**.

    Returns
    -------
    genome.homology.metadata.HomologyMetadata
        The row, with the species as the shipped table spells them.

    Raises
    ------
    UnknownHomologySpeciesError
        If either species is not one this package prepares.
    ValueError
        If the release is not pinned, or if the two species are the same one — a
        **Homology set** relates two species, and a gene's paralogs within one species are
        a different question this does not answer.
    NoHomologyPairError
        If both species are prepared but not together in that release.

    Examples
    --------
    >>> check_pair("homo_sapiens", "Mus musculus", "116").holding_species
    'Mus musculus'
    """
    first, second = check_species(species), check_species(other_species)
    pinned = check_release(release)
    if first == second:
        raise ValueError(
            f"a Homology set relates two different species and both of these are {first!r}. "
            f"Which genes of *another* species a gene is homologous to is the question this "
            f"answers; a gene's paralogs within one species is a different one, and Compara "
            f"publishes it in a file this does not read."
        )
    row = metadata_mod.homology_metadata(first, second, pinned)
    if row is None:
        pairs = ", ".join(
            f"{a}/{b}" for a, b in sorted({r.pair for r in metadata_mod.homology_table()})
        )
        raise NoHomologyPairError(
            f"no Homology set is pinned for {first!r} and {second!r} in Compara release "
            f"{pinned}: the pairs that ship are {pairs}. Add a row to "
            f"data/homology/homology_metadata.tsv naming the file that holds this pair — "
            f"counted, not guessed, since which of the two files holds it is arbitrary."
        )
    return row


def check_stem(stem: str) -> str:
    """Return ``stem`` if it is a **Gene id stem**, else refuse the versioned id it is.

    Examples
    --------
    >>> check_stem("WBGene00020462")
    'WBGene00020462'
    """
    versioned = _VERSIONED_GENE_ID.fullmatch(stem)
    if versioned is not None:
        raise VersionedGeneIdError(
            f"{stem!r} is a versioned gene id and a Homology set is keyed by gene id stem. "
            f"Pass {versioned['stem']!r}: Compara writes its gene ids bare, so the versioned "
            f"spelling would match nothing and come back in `unresolved` looking exactly like "
            f"a gene the publisher never placed in a tree."
        )
    return stem


class HomologySet:
    """One Ensembl Compara **Release**, sliced to one species pair and read into an index.

    **Constructing one prepares it**, as opening a :class:`~genome.genome.Genome` does: the
    per-species dump that holds the pair is fetched into :func:`homology_data_dir` on the
    first construction — verified against the publisher's own md5 as it arrives — sliced to
    the pair, and recorded with a **Completion marker**; every construction after re-reads
    what is there and fetches nothing. The lab's CPU cluster compute nodes have no
    internet, so the first construction of a set must happen on a login node, exactly as
    :class:`~genome.tf.motif.jaspar.JasparDatabase` already documents.

    The pair is unordered as far as the download goes and ordered as far as the *question*
    goes: ``HomologySet("Homo sapiens", "Mus musculus")`` and its reverse read one prepared
    file, and each answers about the species it was named with first.

    **Nothing here is computed.** Every field of every **Homology link** is a cell of the
    publisher's own file, the **Homology type** most of all: it is Compara's tree-derived
    label and is never recomputed, not after a filter and not after resolution into an
    **Annotation** (ADR-0020). No quality score, ranking or "best ortholog" of this
    package's own exists, and no table this package publishes is derived through homology
    (ADR-0019).

    Parameters
    ----------
    species : str
        The species whose genes are asked about, in either spelling.
    other_species : str
        The species whose homologous genes come back.
    release : str, default ``"116"``
        The Compara **Release**, one of those the shipped table pins.
    cache_dir : str or pathlib.Path, optional
        The homology root to prepare under, overriding :func:`homology_data_dir`. The
        directory itself, not a root to file ``homology/`` under; the
        ``ensembl_compara/<release>/<pair>/`` layout is still applied beneath it, since a
        set carries a **Completion marker** and needs a directory of its own.
    progressbar : bool, default True
        Show the download's progress bar. Nothing is drawn when the set is already there.

    Attributes
    ----------
    species : str
        The species asked about, as the shipped table spells it.
    other_species : str
        The species answered with, likewise.
    release : str
        The **Release** this is.
    path : pathlib.Path
        The stored slice these links were read from — a plain gzipped TSV carrying the
        publisher's own header and rows.
    source_url : str
        Where those rows came from: the per-species dump that holds this pair.
    provenance : genome.homology.metadata.HomologyMetadata
        The shipped row behind all of the above.
        :meth:`~genome.homology.metadata.HomologyMetadata.attribution` renders the line to
        print beside anything this answered.
    null_quality_scores : tuple of str
        Which of :data:`QUALITY_SCORE_COLUMNS` this set holds no value in *anywhere*,
        measured over the prepared slice. Both, for either worm pairing; empty for a pair
        Compara scored — so a caller filtering on one is told rather than left to discover
        it when the filter empties.

    Raises
    ------
    UnknownHomologySpeciesError
        If either species is not one this package prepares.
    NoHomologyPairError
        If the pair is not pinned in that release.
    ValueError
        If the release is not pinned, or both species are the same one.
    ComparaPartitionError
        If the recorded file holds none of the pair's rows — the partition moved.
    ComparaFileError
        If the fetched dump is not Compara's, or a prepared slice disagrees with its
        record.

    Examples
    --------
    >>> from genome.homology import HomologySet
    >>> worms = HomologySet("Homo sapiens", "Caenorhabditis elegans")   # doctest: +SKIP
    >>> len(worms)                                                      # doctest: +SKIP
    23982
    >>> worms.null_quality_scores                                       # doctest: +SKIP
    ('goc_score', 'wga_coverage')
    >>> worms.homologs(["ENSG00000152670"]).homolog_gene_id_stems       # doctest: +SKIP
    ['WBGene00001598', 'WBGene00001599', 'WBGene00001600']
    """

    def __init__(
        self,
        species: str,
        other_species: str,
        release: str = DEFAULT_RELEASE,
        *,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> None:
        row = check_pair(species, other_species, release)
        self.provenance = row
        self.release = row.release
        self.species = check_species(species)
        self.other_species = check_species(other_species)
        self.source_url = row.source_url
        root = Path(cache_dir).expanduser() if cache_dir is not None else homology_data_dir()
        directory = set_dir(root, row)
        self.path, record = _prepare(row, directory=directory, progressbar=progressbar)
        links, nulls = _read_slice(self.path, row=row, species=self.species, record=record)
        self._links = links
        self.null_quality_scores = nulls

    def __len__(self) -> int:
        """Return how many **Homology link**s this set holds, of every type."""
        return sum(len(links) for links in self._links.values())

    def __repr__(self) -> str:
        """Return which pair and release this is and how many links it holds."""
        return (
            f"HomologySet(species={self.species!r}, other_species={self.other_species!r}, "
            f"release={self.release!r}, links={len(self)})"
        )

    def homologs(self, stems: Iterable[str], *, paralogs: bool = False) -> HomologyAnswer:
        """Return the other species' genes homologous to each **Gene id stem** asked about.

        The one question a **Homology set** answers. Every stem that named at least one
        link maps to all of them, in the order the stems were asked about, and no value is
        ever empty — a stem this set names no homolog for is in
        :attr:`~genome.io.results.HomologyAnswer.unresolved` instead, so what your list
        holds and this release does not is visible rather than dropped.

        **Orthologs are the default and paralogs come back only on request**, so the common
        question stays the easy one. A **Paralogy link** is kept in the set and marked by
        its own **Homology type** rather than excluded, which is what keeps *not an
        ortholog* distinguishable from *absent* — the stance ADR-0013 takes for a
        **Cross-species link**. Release 116 publishes none for these pairs; see this
        module's own documentation for the count.

        Whatever a filter removed is counted in
        :attr:`~genome.io.results.HomologyAnswer.dropped_partners` rather than silently
        gone, and the **Homology type** on a link that survived is untouched by it
        (ADR-0020): a view can look one-to-one and still be labelled ``ortholog_one2many``.

        Parameters
        ----------
        stems : iterable of str
            The **Gene id stem**s to ask about, in the order they should come back.
            Repeats are asked once. Compara writes its gene ids bare, so a versioned id is
            refused rather than answered emptily.
        paralogs : bool, default False
            Return every link the publisher wrote for these genes, rather than only the
            ones its label calls a speciation event.

        Returns
        -------
        genome.io.results.HomologyAnswer
            The stems that named homologs, mapped to every link each names, and the stems
            that named none.

        Raises
        ------
        VersionedGeneIdError
            If a stem carries a version; the message names the stem to pass.

        Examples
        --------
        >>> homologs = HomologySet("Mus musculus", "Homo sapiens")      # doctest: +SKIP
        >>> answer = homologs.homologs(["ENSMUSG00000059552"])          # doctest: +SKIP
        >>> answer.resolved["ENSMUSG00000059552"][0].homology_type      # doctest: +SKIP
        'ortholog_one2one'
        """
        asked = tuple(dict.fromkeys(check_stem(stem) for stem in stems))
        resolved: dict[str, tuple[HomologyLink, ...]] = {}
        for stem in asked:
            found = self._links.get(stem, ())
            kept = found if paralogs else tuple(link for link in found if link.is_ortholog)
            if kept:
                resolved[stem] = kept
        named = {link.homolog_gene_id_stem for stem in asked for link in self._links.get(stem, ())}
        kept_partners = {link.homolog_gene_id_stem for links in resolved.values() for link in links}
        return HomologyAnswer(
            species=self.species,
            other_species=self.other_species,
            release=self.release,
            resolved=resolved,
            unresolved=tuple(stem for stem in asked if stem not in resolved),
            dropped_partners=tuple(sorted(named - kept_partners)),
            null_quality_scores=self.null_quality_scores,
        )


def _prepare(
    row: HomologyMetadata, *, directory: Path, progressbar: bool
) -> tuple[Path, CompletionRecord]:
    """Return the stored slice and its record, preparing the set once if it is not there.

    Four steps, in the order every build in this package takes them: fetch into the working
    area, slice, place under the final name, write the record. The record is written last,
    after the file it claims exists, so an interrupted run reads as unfinished rather than
    present — and the working area survives an interruption, so repairing one costs no
    second download of a 110 MB dump.

    The partition guard fires *before* anything is placed: a slice holding no rows is never
    written, so a set that would answer empty does not exist on disk to be re-read.
    """
    path = directory / slice_filename(row)
    repair = f"rm -rf {shlex.quote(str(directory))}"
    record = check_registration(directory, repair=repair)
    if record is not None:
        return path, record

    try:
        source = fetch.fetch_url(
            row.source_url,
            work_dir(directory),
            known_hash=f"md5:{row.md5}",
            fname=source_filename(row),
            progressbar=progressbar,
        )
    except OSError as error:
        # `OSError` alone, and not the `ValueError` a checksum mismatch raises: a pin that
        # does not match names bytes that arrived, which is a different fact needing a
        # different next action, and re-running the fetch would not repair it.
        raise HomologySetNotDownloadedError(
            f"the Compara release {row.release} homology set for {row.species!r} and "
            f"{row.other_species!r} is not prepared here and {row.source_url} could not be "
            f"fetched: {error}. Nothing else in this package needs the network, so this is "
            f"the one step that does. Prepare it on a machine with internet — a login node, "
            f"since the lab's compute nodes have none — with "
            f"`{homology_prepare_command(row.species, row.other_species, row.release)}`, after "
            f"which it is read from the Data dir and shared by every project on the machine."
        ) from error
    staged = work_dir(directory) / f"{path.name}.part"
    rows, digest, nulls = _slice(source, staged, row=row)
    if rows == 0:
        staged.unlink(missing_ok=True)
        raise ComparaPartitionError(
            f"{row.source_url} holds no {row.species}/{row.other_species} rows at all, so "
            f"Compara release {row.release} has filed this pair in the other file of the two: "
            f"{compara_url(row.other, row.release)}. Compara says so itself — each per-species "
            f"dump carries 'an arbitrary subset of orthologies involving the given genome', and "
            f"its README tells you to take the files of both genomes — so zero means the pair "
            f"moved rather than that these species share no homologs, and a pair is never "
            f"partially present. Update this pair's row in "
            f"data/homology/homology_metadata.tsv to name {row.other!r} as its "
            f"holding_species, with that file's own URL and md5 read from the release's own "
            f"listing rather than built from this one's."
        )
    directory.mkdir(parents=True, exist_ok=True)
    staged.replace(path)
    written = build_record(
        directory,
        kind=RECORD_KIND,
        name=pair_name(*row.pair),
        files=[path],
        source_url=row.source_url,
        sha256=digest,
        details={
            "publisher": row.publisher,
            "release": row.release,
            "species": list(row.pair),
            "holding_species": row.holding_species,
            "source_md5": row.md5,
            "links": rows,
            "null_quality_scores": list(nulls),
        },
    )
    write_record(directory, written)
    clear_work_dir(directory)
    return path, written


def _slice(
    source: Path, target: Path, *, row: HomologyMetadata
) -> tuple[int, str, tuple[str, ...]]:
    """Write the pair's rows out of ``source`` into ``target``, and say what they were.

    Streams a line at a time: the published dumps run to millions of rows and are never
    held in memory. Rows are copied byte for byte — the **Homology type** and both quality
    scores reach disk exactly as the publisher wrote them.

    Returns the row count, the sha256 of the **unpacked** slice (ADR-0006) and which
    quality columns held nothing anywhere in it.
    """
    wanted = {species_slug(name) for name in row.pair}
    digest = hashlib.sha256()
    written = 0
    scored: set[str] = set()
    with gzip.open(source, "rt", encoding="utf-8") as reading:
        header = reading.readline().rstrip("\n")
        found = tuple(header.split("\t"))
        if found != COMPARA_COLUMNS:
            raise ComparaFileError(
                f"{source} leads with the columns {list(found)} where a Compara homology dump "
                f"leads with {list(COMPARA_COLUMNS)}. This is not the file {row.source_url} "
                f"publishes — delete it and construct the set again."
            )
        target.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(target, "wt", encoding="utf-8") as writing:
            first = header + "\n"
            writing.write(first)
            digest.update(first.encode("utf-8"))
            for line in reading:
                cells = line.rstrip("\n").split("\t")
                if len(cells) != len(COMPARA_COLUMNS):
                    continue
                here, there = cells[_SPECIES], cells[_HOMOLOG_SPECIES]
                if here == there or {here, there} != wanted:
                    continue
                kept = line if line.endswith("\n") else line + "\n"
                writing.write(kept)
                digest.update(kept.encode("utf-8"))
                written += 1
                scored.update(_scored_columns(cells))
    return written, digest.hexdigest(), tuple(n for n in QUALITY_SCORE_COLUMNS if n not in scored)


def _scored_columns(cells: list[str]) -> set[str]:
    """Return which of the quality columns this row records an actual value in."""
    return {name for name in QUALITY_SCORE_COLUMNS if cells[_QUALITY_AT[name]] != NULL_CELL}


def _read_slice(
    path: Path, *, row: HomologyMetadata, species: str, record: CompletionRecord
) -> tuple[Mapping[str, tuple[HomologyLink, ...]], tuple[str, ...]]:
    """Read a stored slice into an index keyed by ``species``'s own gene id stems.

    The slice keeps the publisher's orientation — its ``species`` column is the holding
    species' — so a set asked about the *other* species reads each row the other way round.
    Nothing about the row changes: the **Homology type** is a property of the pair of genes
    and not of which of them was asked about.

    The unpacked bytes are hashed as they are read and held to the digest the **Completion
    marker** recorded, so a slice edited or truncated after it was prepared raises instead
    of answering short. A pair's slice is tens of thousands of rows rather than a genomic
    file, so it is indexed in memory; the multi-million-row dump it was cut from never is.
    """
    flip = species_slug(species) != species_slug(row.holding_species)
    digest = hashlib.sha256()
    index: dict[str, list[HomologyLink]] = {}
    scored: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8") as reading:
        header = reading.readline()
        digest.update(header.encode("utf-8"))
        if tuple(header.rstrip("\n").split("\t")) != COMPARA_COLUMNS:
            raise ComparaFileError(
                f"{path} does not lead with Compara's own columns, so it is not a slice this "
                f"package wrote. Delete {path.parent} and construct the set again."
            )
        for line in reading:
            digest.update(line.encode("utf-8"))
            cells = line.rstrip("\n").split("\t")
            if len(cells) != len(COMPARA_COLUMNS):
                continue
            stem = cells[_HOMOLOG if flip else _GENE]
            index.setdefault(stem, []).append(
                HomologyLink(
                    gene_id_stem=stem,
                    homolog_gene_id_stem=cells[_GENE if flip else _HOMOLOG],
                    homology_type=cells[_TYPE],
                    is_high_confidence=_flag(cells[_HIGH_CONFIDENCE]),
                    goc_score=_whole(cells[_QUALITY_AT["goc_score"]]),
                    wga_coverage=_fraction(cells[_QUALITY_AT["wga_coverage"]]),
                )
            )
            scored.update(_scored_columns(cells))
    if record.sha256 is not None and digest.hexdigest() != record.sha256:
        raise ComparaFileError(
            f"{path} does not hash to what its completion record claims, so it changed after "
            f"it was prepared. Nothing here can be trusted as the release it is filed under — "
            f"delete {path.parent} and construct the set again."
        )
    links = {
        stem: tuple(sorted(found, key=lambda link: link.homolog_gene_id_stem))
        for stem, found in index.items()
    }
    return MappingProxyType(links), tuple(n for n in QUALITY_SCORE_COLUMNS if n not in scored)


def _flag(cell: str) -> bool | None:
    """Read Compara's high-confidence flag — ``None`` where it recorded nothing."""
    if cell == _TRUE_CELL:
        return True
    if cell == _FALSE_CELL:
        return False
    return None


def _whole(cell: str) -> int | None:
    """Read a whole-number quality score — ``None`` where Compara recorded nothing."""
    if cell == NULL_CELL or not cell:
        return None
    try:
        return int(cell)
    except ValueError:
        return None


def _fraction(cell: str) -> float | None:
    """Read a fractional quality score — ``None`` where Compara recorded nothing."""
    if cell == NULL_CELL or not cell:
        return None
    try:
        return float(cell)
    except ValueError:
        return None
