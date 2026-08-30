"""One **Xref set** on disk, and the two verbs it answers.

I/O boundary module: name a species and the **Xref source**'s file is fetched once, sliced
to that species, written into the ``xref/`` subtree of the **Data dir** as a plain gzipped
TSV, and read back as an :class:`XrefSet`. A second construction re-reads what is there
and fetches nothing. The set belongs to no **Assembly** and opens no **Genome** — an
identifier is a name and not a place — so it is filed beside the assembly tree rather than
inside it, the way a **Motif set** already is.

**Two directions and only two**, to the hub and from it (ADR-0017).
:meth:`XrefSet.to_stems` turns foreign ids into **Gene id stem**s,
:meth:`XrefSet.from_stems` turns stems back into foreign ids, and there is no verb that
turns an Entrez id into an HGNC id in one call: a caller wanting that makes both calls and
owns the join, which keeps the hop visible in their code rather than invisible in ours. A
query reads exactly one set, so the **Xref source** is a property of the whole answer
rather than a column on any row, and two publishers are two answers rather than one merged
one.

**A symbol travels one of those directions under its own verb**, because the two are not
mirror images. Away from the hub a stem yields the authority's single current approved
symbol through :meth:`XrefSet.from_stems`, which is labelling a plot. Toward it,
:meth:`XrefSet.match_symbols` matches approved, previous and alias spellings, answers with
every stem any of them names and says which kind each match was —
:meth:`XrefSet.to_stems` refuses the symbol **Namespace** rather than answering it on
approved spellings alone, which would drop the 31 EpiFactors rows this exists for.

**The stored form is a plain gzipped TSV** — three columns, ``namespace``, ``xref_id`` and
``gene_id_stem``, sorted and unique — so a collaborator who does not use Python reads it in
R or in a shell with no library at all. It is a *derived slice* and not the publisher's
bytes, which is why the **Completion marker** beside it carries two checksums: the
publisher's own over its unpacked file, as provenance, and this slice's own sha256, as the
integrity check. A marker that disagrees with what is on disk means the set is unfinished
rather than present, and says so instead of answering a query half-way.

The lab's CPU cluster compute nodes have no internet, so the first construction of a set
must happen on a login node — exactly as :class:`~genome.tf.motif.jaspar.JasparDatabase`
already documents. A construction that cannot fetch raises
:class:`XrefSetNotDownloadedError` naming the call to make there.

Examples
--------
>>> from genome.xref import XrefSet, xref_data_dir
>>> import os
>>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
>>> xref_data_dir()
PosixPath('/scratch/liulab/xref')
>>> del os.environ["LIULAB_DATA"]
>>> worms = XrefSet("Caenorhabditis elegans")                     # doctest: +SKIP
>>> worms.namespaces                                              # doctest: +SKIP
('ensembl', 'entrez', 'uniprot', 'wormbase')
>>> worms.to_stems(["G5EDP9"], "uniprot").resolved                # doctest: +SKIP
{'G5EDP9': ('WBGene00000001',)}
"""

from __future__ import annotations

import gzip
import hashlib
import shlex
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Protocol
from urllib.parse import urlparse

from genome.io import fetch
from genome.io.completion import (
    CompletionRecord,
    RegistrationMismatchError,
    build_record,
    check_registration,
    clear_work_dir,
    work_dir,
    write_record,
)
from genome.io.registration import xref_data_dir
from genome.io.results import ResolvedStems, ResolvedSymbols, ResolvedXrefIds, SymbolMatch
from genome.metadata import species_slug
from genome.xref.alliance import ALLIANCE, read_alliance
from genome.xref.bgi import ALLIANCE_BGI, BGI_SYMBOL_LIMIT, read_bgi
from genome.xref.ensembl import ENSEMBL_TSV, read_ensembl
from genome.xref.evidence import normalise_evidence
from genome.xref.hgnc import HGNC_ARCHIVE, HGNC_SYMBOL_LIMIT, read_hgnc
from genome.xref.ids import ENSEMBL, NAMESPACES, SYMBOL, normalise_id
from genome.xref.metadata import NoXrefSetError, XrefMetadata, lookup_xref
from genome.xref.symbols import (
    KIND_NAMESPACES,
    SYMBOL_KINDS,
    SYMBOL_NAMESPACES,
    SymbolDirectionError,
    fold_symbol,
    normalise_symbol,
)

if TYPE_CHECKING:  # pragma: no cover - typeshed's name for a live hash object
    from hashlib import _Hash as Hash

#: What a **Completion marker** written here calls the thing it recorded, beside the
#: ``genome``, ``annotation`` and ``index`` kinds the assembly tree writes.
XREF_KIND = "xref"

#: The stored slice's columns, in file order. Three and no more: what the answer is keyed
#: by, what it answers with, and which **Namespace** the key belongs to. The **Xref
#: source** and the **Release** are *not* columns — they are properties of the whole set,
#: which is what makes merging two publishers into one file inexpressible (ADR-0017).
SLICE_COLUMNS: tuple[str, ...] = ("namespace", "xref_id", "gene_id_stem")

#: Every value the slice's ``namespace`` column may hold: the **Namespace**s a caller may
#: name, plus the two stored spellings that carry a **Symbol match**'s kind
#: (:data:`~genome.xref.symbols.KIND_NAMESPACES`). A kind is written into this column
#: rather than into a fourth one, because a column saying *what sort of row this is* is the
#: level discriminator this design refuses everywhere else.
STORED_NAMESPACES: tuple[str, ...] = NAMESPACES + tuple(
    namespace for namespace in SYMBOL_NAMESPACES if namespace not in NAMESPACES
)

#: What one stored slice is called, after the species slug. Gzipped, because it is bulk;
#: plain TSV inside, because a collaborator reads it in R.
SLICE_SUFFIX = ".xref_table.tsv.gz"


class XrefReader(Protocol):
    """Reads one **Xref source**'s published file into ``(namespace, id, stem)`` triples.

    What a source *is*, beside its row in the curated table: a pure function from the
    publisher's lines to this package's triples, which is what makes adding a source data
    plus a reader rather than a refactor. It opens nothing and downloads nothing —
    :func:`_prepare` has already fetched, unpacked and verified the bytes by the time one
    of these is called.
    """

    def __call__(
        self,
        lines: Iterable[str],
        *,
        ncbi_taxid: int,
        origin: str,
        evidence: tuple[str, ...] = (),
    ) -> tuple[tuple[str, str, str], ...]:
        """Return this species' triples, sorted and unique, from ``lines``.

        ``evidence`` is the grading filter, empty for none. Whether a source records a
        grading at all is the reader's own fact — a file with no such column meets a filter
        with :class:`~genome.xref.evidence.EvidenceNotRecordedError` rather than ignoring
        it, since a filter silently dropped is a quality claim nobody made.
        """
        ...


#: Every **Xref source** with a reader, keyed by the name its curated rows use. Adding one
#: is a module and an entry here; nothing else in this file changes.
_READERS: Mapping[str, XrefReader] = MappingProxyType(
    {
        ALLIANCE: read_alliance,
        ENSEMBL_TSV: read_ensembl,
        HGNC_ARCHIVE: read_hgnc,
        ALLIANCE_BGI: read_bgi,
    }
)

#: What each **Xref source** that carries symbols cannot match, and why — ``None`` for one
#: that publishes all three kinds. Read onto every answer :meth:`XrefSet.match_symbols`
#: gives, so the explanation is part of the behaviour rather than a comment: *this gene is
#: not in the release* and *this source does not publish the spelling you used* are
#: different answers and must not both be silence. A source absent from here carries no
#: symbols at all and says so through :class:`NamespaceNotCarriedError`.
_SYMBOL_LIMITS: Mapping[str, str | None] = MappingProxyType(
    {HGNC_ARCHIVE: HGNC_SYMBOL_LIMIT, ALLIANCE_BGI: BGI_SYMBOL_LIMIT}
)


class NamespaceNotCarriedError(LookupError):
    """The set carries no such **Namespace**, and the message names the ones it does.

    A :class:`LookupError` and not a :class:`ValueError`: the namespace is a name that this
    release resolves nothing under, which is the same kind of miss as an unknown species.
    The three species carry three different authorities, so a mouse set asked for ``hgnc``
    lands here rather than answering nothing — the failure that would otherwise look like
    a gene list with no matches.

    Examples
    --------
    >>> mouse = XrefSet("Mus musculus")                            # doctest: +SKIP
    >>> mouse.to_stems(["HGNC:11998"], "hgnc")                     # doctest: +SKIP
    Traceback (most recent call last):
    NamespaceNotCarriedError: ...
    """


class XrefSetNotDownloadedError(RuntimeError):
    """The set is not on disk and could not be fetched, so nothing can answer.

    A :class:`RuntimeError`, because nothing about the call was wrong: the bytes are simply
    not here and this machine could not go and get them. The lab's compute nodes have no
    internet, so the message names the call to make on a login node instead.

    Examples
    --------
    >>> XrefSet("Homo sapiens")                                    # doctest: +SKIP
    Traceback (most recent call last):
    XrefSetNotDownloadedError: ...
    """


class XrefTableError(ValueError):
    r"""A table read here is not the shape it must be, so it is not allowed to answer.

    Covers both files this module reads: the publisher's, when what arrived does not match
    the checksum the curated row pins or carries no row for the species asked for; and the
    stored slice, when its header, its columns or a **Namespace** in it is not what this
    package writes. A :class:`ValueError`, because a file that says something the format
    does not is a bad value rather than a broken program, and the message names the file
    and the repair.

    Examples
    --------
    >>> try:
    ...     parse_slice("wrong\theader\n", origin="example.tsv")
    ... except XrefTableError as error:
    ...     print("namespace" in str(error))
    True
    """


def xref_set_dir(species: str, source: str, release: str, *, evidence: Sequence[str] = ()) -> Path:
    """Return the directory one **Xref set** is prepared in, whether or not it exists.

    ``<liulab_data>/xref/<source>/<release>/<species slug>/``. Source and release above
    species, so two releases of one publisher sit side by side and neither is *the* xref
    directory — holding two releases at once is the whole point of pinning one.

    A set built under an evidence filter is a **different set** and gets a directory of its
    own beside the unfiltered one, since a filter that changes which rows are read changes
    what the stored slice holds. Two callers who named the same types in different orders
    land on one directory, the filter being sorted on the way in.

    Parameters
    ----------
    species : str
        The species, in either the curated table's spelling or its slug.
    source : str
        The **Xref source**.
    release : str
        The pinned **Release**.
    evidence : sequence of str, optional
        The evidence filter the set was built under, as
        :func:`~genome.xref.evidence.normalise_evidence` spells it. Empty for none.

    Returns
    -------
    pathlib.Path
        The set's own directory. Nothing is created by asking.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> xref_set_dir("Homo sapiens", "alliance", "9.0.0")
    PosixPath('/scratch/liulab/xref/alliance/9.0.0/homo_sapiens')
    >>> xref_set_dir("Homo sapiens", "ensembl", "116", evidence=("DEPENDENT",))
    PosixPath('/scratch/liulab/xref/ensembl/116/homo_sapiens.evidence-dependent')
    >>> del os.environ["LIULAB_DATA"]
    """
    return xref_data_dir() / source / release / f"{species_slug(species)}{_suffix(evidence)}"


def _suffix(evidence: Sequence[str]) -> str:
    """Return what an evidence filter adds to a set's directory name, or nothing."""
    return f".evidence-{'-'.join(kind.lower() for kind in evidence)}" if evidence else ""


def xref_slice_name(species: str) -> str:
    """Return the file name one species' stored slice is written under.

    Parameters
    ----------
    species : str
        The species, in either the curated table's spelling or its slug.

    Returns
    -------
    str
        ``<species slug>.xref_table.tsv.gz``.

    Examples
    --------
    >>> xref_slice_name("Caenorhabditis elegans")
    'caenorhabditis_elegans.xref_table.tsv.gz'
    """
    return f"{species_slug(species)}{SLICE_SUFFIX}"


def xref_prepare_command(
    species: str, source: str, release: str, *, evidence: Sequence[str] = ()
) -> str:
    r"""Return the call that prepares one **Xref set**, for an error message to quote.

    One spelling of it, so a renamed entry point is renamed once. Quoted by every error
    here that a caller repairs by fetching the set on a machine with internet — which is
    why an evidence filter travels in it: repairing a filtered set by preparing the
    unfiltered one would leave the caller exactly where they started.

    Parameters
    ----------
    species : str
        The species, as the curated table spells it.
    source : str
        The **Xref source**.
    release : str
        The pinned **Release**.
    evidence : sequence of str, optional
        The evidence filter the set was built under. Omitted from the command when empty.

    Returns
    -------
    str
        A shell command, unquoted and unfenced — the caller decides how to set it.

    Examples
    --------
    >>> xref_prepare_command("Homo sapiens", "alliance", "9.0.0")
    'python -c "from genome.xref import XrefSet; XrefSet(\'Homo sapiens\', \'alliance\', \'9.0.0\')"'
    >>> "evidence=('DEPENDENT',)" in xref_prepare_command(
    ...     "Homo sapiens", "ensembl", "116", evidence=("DEPENDENT",)
    ... )
    True
    """
    filtered = f", evidence={tuple(evidence)!r}" if evidence else ""
    call = (
        f"from genome.xref import XrefSet; XrefSet({species!r}, {source!r}, {release!r}{filtered})"
    )
    return f'python -c "{call}"'


def parse_slice(text: str, *, origin: str) -> tuple[tuple[str, str, str], ...]:
    r"""Read a stored slice's text into ``(namespace, id, stem)`` triples, in file order.

    A pure function from text to triples: it opens nothing and holds no opinion about which
    **Release** the text came from. Public because it is what says the stored form *is* a
    plain TSV — anything that can read three tab-separated columns can read one of these.

    Parameters
    ----------
    text : str
        The whole slice, header line included.
    origin : str
        The file the text came from; named in every message, since deleting that file and
        constructing the set again is the repair.

    Returns
    -------
    tuple of tuple of str
        One triple per data row, in file order.

    Raises
    ------
    XrefTableError
        If the header is not :data:`SLICE_COLUMNS`, a row is not three fields, a cell is
        blank, or a **Namespace** is one this package does not know.

    Examples
    --------
    >>> parse_slice("namespace\txref_id\tgene_id_stem\nentrez\t7157\tENSG00000141510\n",
    ...             origin="example.tsv")
    (('entrez', '7157', 'ENSG00000141510'),)
    """
    lines = text.splitlines()
    if not lines or tuple(lines[0].split("\t")) != SLICE_COLUMNS:
        found = lines[0] if lines else "<empty file>"
        raise XrefTableError(
            f"{origin} opens with {found!r} where a stored xref slice opens with "
            f"{'/'.join(SLICE_COLUMNS)}. It is not one of this package's slices, or it was "
            f"written by a version that spelled them differently. Delete it and construct "
            f"the set again."
        )
    triples: list[tuple[str, str, str]] = []
    for number, line in enumerate(lines[1:], start=2):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != len(SLICE_COLUMNS) or not all(fields):
            raise XrefTableError(
                f"{origin} line {number} is {line!r}, and every row of a stored xref slice "
                f"is {len(SLICE_COLUMNS)} non-empty tab-separated fields. Delete the file "
                f"and construct the set again."
            )
        if fields[0] not in STORED_NAMESPACES:
            raise XrefTableError(
                f"{origin} line {number} names the namespace {fields[0]!r}, and the ones "
                f"this package knows are {', '.join(STORED_NAMESPACES)}. Delete the file "
                f"and construct the set again."
            )
        triples.append((fields[0], fields[1], fields[2]))
    return tuple(triples)


class XrefSet:
    """One species, one **Xref source**, one pinned **Release**, prepared on disk.

    **Constructing one prepares it**, as opening a :class:`~genome.genome.Genome` does: the
    publisher's file is fetched on the first construction, sliced to this species and
    written under :func:`xref_set_dir`, and every construction after re-reads what is there
    and fetches nothing. It answers with no genome open and belongs to no assembly.

    It answers two questions and only two — :meth:`to_stems` and :meth:`from_stems`, with
    :meth:`match_symbols` the first of them asked about a symbol, where the answer carries
    the kind of spelling that matched — and never converts one foreign **Namespace**
    directly into another. Gene level only: a gene, a transcript and a protein have
    different keys and different sources and are three objects rather than one table with a
    level column.

    Parameters
    ----------
    species : str
        The species, in either the curated table's spelling (``"Homo sapiens"``) or its
        slug (``"homo_sapiens"``).
    source : str, optional
        The **Xref source**. Omitted, the species' **Default xref source** answers — which
        is a default and not a recommendation: naming one is how the scientific choice gets
        made deliberately, and two publishers disagreeing are two answers rather than one
        merged one.

        **This constructor fills in the identifier default, and a set that carries no
        symbol keeps refusing to match one** — it holds one publisher's bytes and answering
        from another's is exactly what one query reading one set forbids (ADR-0017).
        :meth:`for_symbols` is the constructor with the symbol question named, and is what
        fills in the source that carries them; :meth:`for_namespace` is the same fill-in for
        a caller holding a **Namespace** rather than a verb.

        **The sources are not equals, and the choice is nearly half the answer.** Measured
        on human release 116 against NCBI's own file, Ensembl and NCBI agree on only
        **57.6%** of the gene-level (GeneID, ENSG) pairs they assert between them. The
        cause is method rather than release skew: NCBI's mapping is a sequence match at a
        published overlap threshold and is all but one-to-one, while **Ensembl's fans out
        to 72 stems for one GeneID** (``79166``) **and 208 GeneIDs for one stem**
        (``ENSG00000278233``). Ask ``ensembl`` for an id and expect a wider answer than
        ``alliance`` gives for the same one; the width is the publisher's assertion, and
        nothing here narrows it or reconciles the two.
    release : str, optional
        The pinned **Release**. Omitted, the newest the curated table lists. Each source
        pins its own numbering and they do not correspond — ``9.0.0`` is Alliance's and
        ``116`` is Ensembl's.
    evidence : str or iterable of str, optional
        Keep only the rows the publisher graded with one of these ``info_type``s, or
        ``None`` for every row. A capability of the source rather than of every set: a
        publisher whose file grades nothing raises
        :class:`~genome.xref.evidence.EvidenceNotRecordedError` rather than ignoring the
        filter. **A filter that keeps nothing raises too**, because every human
        ``EntrezGene`` row Ensembl release 116 publishes is ``DEPENDENT`` and not one is
        ``DIRECT``, so the intuitive quality filter empties the set rather than narrowing
        it. A filtered set is prepared beside the unfiltered one and never over it.
    cache_dir : str or pathlib.Path, optional
        The directory to prepare in, overriding :func:`xref_set_dir`. The directory itself,
        not a root to file under.
    progressbar : bool, default True
        Show the download's progress bar. Nothing is drawn when the set is already there.

    Attributes
    ----------
    species : str
        The species, as the curated table spells it.
    source : str
        The **Xref source** whose assertions this carries.
    release : str
        The pinned **Release**.
    evidence : tuple of str
        The evidence filter this set was built under, empty for none.
    path : pathlib.Path
        The stored slice these mappings were read from — a plain gzipped TSV.
    source_url : str
        Where the publisher's own file was fetched from.
    provenance : genome.xref.metadata.XrefMetadata
        The curated row this set actually resolved to, defaults filled in — who published
        it, which release, and the paper to cite. Read off the set rather than looked up
        again, so what is cited is what answered.
    namespaces : tuple of str
        The **Namespace**s this set actually carries, read off the slice rather than
        declared, in :data:`~genome.xref.ids.NAMESPACES` order. One that is not here raises
        rather than answering nothing.
    symbol_kinds : tuple of str
        Which kinds of **Symbol match** this set can make — ``approved``, ``previous``,
        ``alias`` — read off the slice the same way, and empty for a source that carries no
        symbols at all.
    symbol_limits : str or None
        Why the kinds not in :attr:`symbol_kinds` are missing, or ``None`` when all three
        are there or none is. It rides back on every :meth:`match_symbols` answer, because
        *this gene is not in the release* and *this source does not publish the spelling
        you used* are different answers and must not both be silence.

    Raises
    ------
    genome.xref.metadata.NoXrefSetError
        If no set exists for that species, source or release. The message names what does.
    XrefSetNotDownloadedError
        If the set is not on disk and could not be fetched.
    XrefTableError
        If the publisher's file or the stored slice is not the shape it must be.
    genome.xref.evidence.EvidenceNotRecordedError
        If an evidence filter is named and this source's file grades nothing.
    genome.xref.evidence.EmptyEvidenceFilterError
        If an evidence filter is named and it keeps none of the release's rows.
    genome.io.completion.RegistrationMismatchError
        If the **Completion marker** disagrees with what is on disk, either about the
        slice's size or about its checksum — both mean unfinished rather than present.

    Examples
    --------
    >>> from genome.xref import XrefSet
    >>> human = XrefSet("Homo sapiens")                            # doctest: +SKIP
    >>> human.namespaces                                           # doctest: +SKIP
    ('ensembl', 'entrez', 'uniprot', 'hgnc')
    >>> human.to_stems(["7157"], "entrez").resolved                # doctest: +SKIP
    {'7157': ('ENSG00000141510',)}
    >>> human.from_stems(["ENSG00000141510.18"], "hgnc").resolved  # doctest: +SKIP
    {'ENSG00000141510.18': ('HGNC:11998',)}
    >>> print(human.provenance.attribution())                      # doctest: +SKIP
    Alliance of Genome Resources 9.0.0 (PMID 38552170) — https://download.alliancegenome.org/...
    """

    def __init__(
        self,
        species: str,
        source: str | None = None,
        release: str | None = None,
        *,
        evidence: str | Iterable[str] | None = None,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> None:
        row = lookup_xref(species, source, release)
        self.species = row.species
        self.source = row.source
        self.release = row.release
        self.source_url = row.url
        self.provenance = row
        self.evidence = normalise_evidence(evidence)
        directory = (
            Path(cache_dir).expanduser()
            if cache_dir is not None
            else xref_set_dir(row.species, row.source, row.release, evidence=self.evidence)
        )
        self.path = directory / xref_slice_name(row.species)
        record = _prepare(
            row,
            directory=directory,
            path=self.path,
            evidence=self.evidence,
            progressbar=progressbar,
        )
        triples, digest = _read_slice(self.path)
        if record.sha256 != digest:
            raise RegistrationMismatchError(
                f"{self.path} hashes to {digest} where its record claims {record.sha256}, so "
                f"the slice on disk is not the one that was prepared. Something rewrote it "
                f"after the record was written; nothing here can be trusted as complete. "
                f"Re-prepare it with `rm -rf {shlex.quote(str(directory))} && "
                f"{xref_prepare_command(*_names(row), evidence=self.evidence)}`."
            )
        self._to_stems, self._from_stems = _index(triples)
        self.namespaces: tuple[str, ...] = tuple(
            namespace for namespace in NAMESPACES if namespace in self._to_stems
        )
        self.symbol_kinds: tuple[str, ...] = tuple(
            kind for kind in SYMBOL_KINDS if KIND_NAMESPACES[kind] in self._to_stems
        )
        self.symbol_limits: str | None = (
            _SYMBOL_LIMITS.get(self.source) if self.symbol_kinds else None
        )
        self._exact, self._folded = _symbol_index(self._to_stems)

    @classmethod
    def for_namespace(
        cls,
        species: str,
        namespace: str,
        source: str | None = None,
        release: str | None = None,
        *,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> XrefSet:
        """Return the set that answers a question about ``namespace`` when none is named.

        **A Default xref source is per species and per question** (ADR-0021), and this is
        where a caller holding a **Namespace** rather than a verb names the question: the
        symbol one fills in the species' symbol-carrying default, every other one fills in
        its identifier default. It is :meth:`for_symbols` generalised to the shape a caller
        who read a namespace off a flag already has, so that surface has nothing left to
        decide — the choosing happens here and in no second place.

        Which namespace a set actually carries is still the set's to say: this fills a
        source in and does not check, so a namespace the resolved set does not carry raises
        on the verb, naming the ones it does.

        Parameters
        ----------
        species : str
            The species, in either the curated table's spelling or its slug.
        namespace : str
            The **Namespace** the question is about, which is what picks the default. Case
            is not significant, as it is not on the verbs.
        source : str, optional
            The **Xref source**. Omitted, the default this namespace's question implies
            answers. Named, it is honoured whatever the namespace is — naming one is the
            deliberate scientific choice and is never swapped for a flagged row.
        release : str, optional
            The pinned **Release**. Omitted, the newest that source has.
        cache_dir : str or pathlib.Path, optional
            The directory to prepare in, as on the ordinary constructor.
        progressbar : bool, default True
            Show the download's progress bar.

        Returns
        -------
        XrefSet
            The prepared set, which is an ordinary one in every respect — every verb answers
            on it and nothing about it remembers which question filled its source in.

        Raises
        ------
        genome.xref.metadata.NoXrefSetError
            If no set exists for that species or release, or if the namespace is the symbol
            one and no row for the species is flagged to answer symbols.

        Examples
        --------
        >>> from genome.xref import XrefSet
        >>> XrefSet.for_namespace("Homo sapiens", "symbol").source   # doctest: +SKIP
        'hgnc'
        >>> XrefSet.for_namespace("Homo sapiens", "entrez").source   # doctest: +SKIP
        'alliance'
        """
        row = lookup_xref(species, source, release, for_symbols=namespace.strip().lower() == SYMBOL)
        return cls(
            species,
            row.source,
            row.release,
            cache_dir=cache_dir,
            progressbar=progressbar,
        )

    @classmethod
    def for_symbols(
        cls,
        species: str,
        source: str | None = None,
        release: str | None = None,
        *,
        cache_dir: str | Path | None = None,
        progressbar: bool = True,
    ) -> XrefSet:
        """Return the set that answers **symbols** for a species when none is named.

        The ordinary constructor with the question named, and nothing else: a **Default
        xref source** is per species *and* per question (ADR-0021), because the publisher
        carrying a species' identifiers is usually not the one carrying its symbols — human
        ids default to the Alliance, whose cross-reference file publishes no human symbol at
        all, and HGNC's quarterly archive is what does. Mouse and worm reach a third source
        again, ``alliance_bgi``.

        **A source named here is honoured exactly as it is anywhere else.** Naming one is
        how the scientific choice gets made deliberately, so this fills in a default and
        never overrides a choice: ``XrefSet.for_symbols(species, "alliance")`` is
        ``XrefSet(species, "alliance")``, symbols and all — which for human means a set that
        matches none and says so.

        This is :meth:`for_namespace` with the namespace fixed, so the two cannot answer
        differently: one fill-in, named by a verb here and by a namespace there.

        Parameters
        ----------
        species : str
            The species, in either the curated table's spelling or its slug.
        source : str, optional
            The **Xref source**. Omitted, the species' symbol-carrying default answers.
        release : str, optional
            The pinned **Release**. Omitted, the newest that source has. Honoured against
            whichever source was filled in, exactly as on the ordinary constructor.
        cache_dir : str or pathlib.Path, optional
            The directory to prepare in, as on the ordinary constructor.
        progressbar : bool, default True
            Show the download's progress bar.

        Returns
        -------
        XrefSet
            The prepared set, which is an ordinary one in every respect — every verb answers
            on it and nothing about it remembers which question filled its source in.

        Raises
        ------
        genome.xref.metadata.NoXrefSetError
            If no set exists for that species or release, or if no row for the species is
            flagged to answer symbols. The message names what does exist.

        Examples
        --------
        >>> from genome.xref import XrefSet
        >>> human = XrefSet.for_symbols("Homo sapiens")             # doctest: +SKIP
        >>> human.source                                           # doctest: +SKIP
        'hgnc'
        >>> human.match_symbols(["ARNTL"]).gene_id_stems           # doctest: +SKIP
        ['ENSG00000133794']
        """
        return cls.for_namespace(
            species, SYMBOL, source, release, cache_dir=cache_dir, progressbar=progressbar
        )

    def __len__(self) -> int:
        """Return how many **Gene id stem**s this set carries.

        Examples
        --------
        >>> len(XrefSet("Caenorhabditis elegans"))                 # doctest: +SKIP
        46926
        """
        return len(self._from_stems[ENSEMBL])

    def __repr__(self) -> str:
        """Return which set this is and how many stems it holds.

        The evidence filter appears only when there is one, so a filtered set never reads
        as the unfiltered set it is not.
        """
        filtered = f", evidence={self.evidence!r}" if self.evidence else ""
        return (
            f"XrefSet(species={self.species!r}, source={self.source!r}, "
            f"release={self.release!r}{filtered}, stems={len(self)})"
        )

    def to_stems(self, ids: Iterable[str], namespace: str) -> ResolvedStems:
        """Return the **Gene id stem**s this release says each foreign id names.

        The hop *toward* the hub. Every id is reduced to one spelling on the way in —
        version dropped, the namespace's own CURIE prefix accepted whether or not it is
        written — so ``ENSG00000141510.18`` and ``ENSG00000141510``, or ``HGNC:11998`` and
        ``11998``, are one identifier and resolve identically. Joining a versioned id to a
        bare one otherwise returns zero matches and says nothing, which is the most
        error-prone detail in this landscape.

        **Every stem, and never a chosen one.** A foreign id naming two stems answers with
        both: 2,535 of 40,665 human genes carry more than one Ensembl cross-reference in
        Alliance 9.0.0, so 6.2% of HGNC ids are ambiguous and nothing here picks a side.

        **Nothing is dropped.** Ids this release names no stem for come back in
        :attr:`~genome.io.results.ResolvedStems.unresolved`, in ask order, so what a list
        holds and this release does not is visible rather than silently shorter.

        Parameters
        ----------
        ids : iterable of str
            The foreign ids, in the order they should come back. Repeats are asked once,
            on the caller's own spelling, so a versioned and an unversioned spelling of one
            id are two entries with identical values and the answer still zips against the
            caller's table row for row.
        namespace : str
            The **Namespace** those ids belong to, one of :attr:`namespaces`. Named rather
            than sniffed, because the string does not say. Case is not significant.

        Returns
        -------
        genome.io.results.ResolvedStems
            The ids that named stems, mapped to every stem each names, and the ids that
            named none — with the species, source, release and namespace that answered.

        Raises
        ------
        NamespaceNotCarriedError
            If this set carries no such namespace. The message names the ones it does.
        genome.xref.symbols.SymbolDirectionError
            If the namespace is the symbol one. A symbol is not answered like an id — it
            matches previous and alias spellings too, and each match carries which kind it
            was — so the message names :meth:`match_symbols` rather than answering here on
            approved spellings alone.

        Examples
        --------
        >>> human = XrefSet("Homo sapiens")                        # doctest: +SKIP
        >>> human.to_stems(["7157", "999999999"], "entrez")        # doctest: +SKIP
        ResolvedStems(species='Homo sapiens', ...)
        """
        if namespace.strip().lower() == SYMBOL:
            raise SymbolDirectionError(
                f"a {SYMBOL!r} namespace is not asked toward the hub with to_stems: it "
                f"would match this release's approved spellings and nothing else, so a "
                f"table spelling a gene the way the authority used to would come back "
                f"unresolved rather than matched — which is what happens to 31 of "
                f"EpiFactors' 801 rows. Call match_symbols(symbols) instead, which matches "
                f"approved, previous and alias spellings and says on each match which kind "
                f"it was. The other direction, from_stems(stems, {SYMBOL!r}), is the "
                f"labelling one and is answered here."
            )
        checked = self._checked(namespace)
        index = self._to_stems[checked]
        asked = tuple(dict.fromkeys(ids))
        found = {key: index[normalise_id(key, checked)] for key in asked}
        return ResolvedStems(
            species=self.species,
            source=self.source,
            release=self.release,
            namespace=checked,
            resolved={key: stems for key, stems in found.items() if stems},
            unresolved=tuple(key for key in asked if not found[key]),
        )

    def from_stems(self, stems: Iterable[str], namespace: str) -> ResolvedXrefIds:
        """Return the foreign ids this release says each **Gene id stem** names.

        The hop *away* from the hub, and :meth:`to_stems`'s mirror in every respect: the
        same normalisation on the way in, every id and never a chosen one, and the stems
        that named nothing riding back in ask order. A stem this release never carried and
        a stem it carries with no id in *this* namespace are one bucket, since no id
        history is held that could tell a retirement from an absence (ADR-0017).

        Parameters
        ----------
        stems : iterable of str
            The stems, in the order they should come back. A versioned gene id is accepted
            and reduced to its stem, so an annotation's own ids may be passed straight in.
        namespace : str
            The **Namespace** to answer in, one of :attr:`namespaces`. Case is not
            significant.

        Returns
        -------
        genome.io.results.ResolvedXrefIds
            The stems that named ids, mapped to every id each names, and the stems that
            named none — with the species, source, release and namespace that answered.

        Raises
        ------
        NamespaceNotCarriedError
            If this set carries no such namespace. The message names the ones it does, and
            for the symbol namespace it names this species' symbol source and
            :meth:`for_symbols` too — labelling is a symbol question like matching one, and
            misses on a set carrying none for the same reason.

        Examples
        --------
        >>> worms = XrefSet("Caenorhabditis elegans")              # doctest: +SKIP
        >>> worms.from_stems(["WBGene00000001"], "uniprot").xref_ids   # doctest: +SKIP
        ['G5EDP9']
        """
        checked = self._checked(namespace)
        index = self._from_stems[checked]
        asked = tuple(dict.fromkeys(stems))
        found = {key: index[normalise_id(key, ENSEMBL)] for key in asked}
        return ResolvedXrefIds(
            species=self.species,
            source=self.source,
            release=self.release,
            namespace=checked,
            resolved={key: ids for key, ids in found.items() if ids},
            unresolved=tuple(key for key in asked if not found[key]),
        )

    def match_symbols(
        self, symbols: Iterable[str], *, case_insensitive: bool = False
    ) -> ResolvedSymbols:
        """Return every gene this release says each symbol names, and how each matched.

        The hop toward the hub from the one **Namespace** answered unlike the rest, and
        deliberately **not** :meth:`from_stems`'s mirror. A symbol is matched against
        approved, previous **and** alias spellings, answers with every **Gene id stem** any
        of them names, and each hit says which kind of spelling it was — so ambiguity is
        the return type and not an edge case. ``ADCY3`` is HGNC's approved symbol for one
        gene and a symbol it retired from another, and both come back.

        **Exact by default.** The species is fixed by the set, so ``Brca1`` asked of a human
        set is a mouse spelling asked of the wrong authority and matches nothing rather than
        half-working. ``case_insensitive=True`` folds both sides and **still answers with
        every gene matched** rather than picking one.

        **What this source could not have matched rides back on the answer**, in
        :attr:`~genome.io.results.ResolvedSymbols.kinds` and
        :attr:`~genome.io.results.ResolvedSymbols.limits`: mouse and worm match approved
        spellings only, their authorities' typed previous and alias spellings belonging to
        publishers that cannot be pinned or cannot be fetched (ADR-0018), and an answer that
        did not say so would look exactly like a gene that is not in the release.

        Parameters
        ----------
        symbols : iterable of str
            The symbols, in the order they should come back. Surrounding whitespace goes;
            case does not, unless ``case_insensitive`` is set. Repeats are asked once, on
            the caller's own spelling, so the answer still zips against their table.
        case_insensitive : bool, default False
            Fold case on both sides — the caller's spelling and the authority's.

        Returns
        -------
        genome.io.results.ResolvedSymbols
            The symbols that matched, mapped to every match each made, and the symbols that
            matched nothing — with the species, source and release that answered, which
            kinds it could match and why the others are missing.

        Raises
        ------
        NamespaceNotCarriedError
            If this set carries no symbols at all — which the species' identifier default
            does not, for human. It raises rather than reaching for another publisher's
            bytes, and the message names the one source that answers this species' symbols
            and the constructor that fills it in, :meth:`for_symbols`.

        Examples
        --------
        >>> human = XrefSet("Homo sapiens", "hgnc")                # doctest: +SKIP
        >>> human.match_symbols(["ARNTL"]).resolved["ARNTL"]       # doctest: +SKIP
        (SymbolMatch(symbol='ARNTL', gene_id_stem='ENSG00000133794', kind='previous'),)
        >>> XrefSet.for_symbols("Homo sapiens").source             # doctest: +SKIP
        'hgnc'
        """
        if not self.symbol_kinds:
            raise self._no_symbols(SYMBOL)
        index = self._folded if case_insensitive else self._exact
        key = fold_symbol if case_insensitive else normalise_symbol
        asked = tuple(dict.fromkeys(symbols))
        found = {symbol: index.get(key(symbol), ()) for symbol in asked}
        return ResolvedSymbols(
            species=self.species,
            source=self.source,
            release=self.release,
            case_insensitive=case_insensitive,
            kinds=self.symbol_kinds,
            limits=self.symbol_limits,
            resolved={symbol: matches for symbol, matches in found.items() if matches},
            unresolved=tuple(symbol for symbol in asked if not found[symbol]),
        )

    def _checked(self, namespace: str) -> str:
        """Return ``namespace`` if this set carries it, else say which ones it does.

        The symbol namespace misses differently from every other one: which source carries
        a species' symbols is a question the curated table answers, so the message routes
        there rather than stopping at the list of what this set holds.
        """
        wanted = namespace.strip().lower()
        if wanted in self.namespaces:
            return wanted
        if wanted == SYMBOL:
            raise self._no_symbols(namespace)
        raise NamespaceNotCarriedError(
            f"the {self.source} {self.release} set for {self.species!r} carries no "
            f"{namespace!r} namespace: it carries {', '.join(self.namespaces)}. Ask in "
            f"one of those — the three species have three different authorities, so a "
            f"namespace that answers for one is not one that answers for another."
        )

    def _no_symbols(self, namespace: str) -> NamespaceNotCarriedError:
        """Return the error both symbol questions raise on a set that carries none.

        One message for both directions: matching a symbol and labelling a stem with one
        miss on this set for the same reason and have the same next action, so they say the
        same thing rather than one of them naming nowhere to go.
        """
        return NamespaceNotCarriedError(
            f"the {self.source} {self.release} set for {self.species!r} carries no "
            f"{namespace!r} namespace: it carries {', '.join(self.namespaces)}. "
            f"{_symbol_source_hint(self.species)} This set is not asked on another "
            f"publisher's behalf — two publishers are two answers and one query reads "
            f"exactly one set (ADR-0017)."
        )


def _names(row: XrefMetadata) -> tuple[str, str, str]:
    """Return the three strings that name one set, for a command to be built from."""
    return row.species, row.source, row.release


def _symbol_source_hint(species: str) -> str:
    """Return the sentence naming the source that answers this species' symbols.

    Read off the curated table rather than listed from the readers, so the hint is *this
    species'* one source and not every source that publishes symbols for anybody — a mouse
    set told to try ``hgnc`` sends the caller to a publisher with no mouse row in it.
    """
    try:
        row = lookup_xref(species, for_symbols=True)
    except NoXrefSetError:
        carried = ", ".join(sorted(_SYMBOL_LIMITS))
        return f"The sources that publish symbols are {carried} — construct the set naming one."
    return (
        f"Symbols for {species!r} come from {row.source} — construct the set as "
        f"XrefSet({species!r}, {row.source!r}), or as XrefSet.for_symbols({species!r}), "
        f"which fills that in for you."
    )


def _prepare(
    row: XrefMetadata,
    *,
    directory: Path,
    path: Path,
    evidence: tuple[str, ...],
    progressbar: bool,
) -> CompletionRecord:
    """Return the finished set's record, fetching and slicing it once if it is not there.

    The publisher's file lands in the directory's working area, is verified against the
    checksum the curated row pins, is sliced to this species and written under a name of
    its own, and is renamed into place only once the whole slice exists. The **Completion
    marker** is written last, after the file it claims, and the working area is emptied
    once it is — so the set is finished or it is not, and a run killed anywhere in between
    leaves a directory that reads as unfinished rather than as present.

    An evidence filter is applied by the reader, as the publisher's rows go past, and a
    filter that keeps nothing raises there — before anything is written, so a set that
    could only answer nothing is never left on disk to be re-read as present.
    """
    existing = check_registration(
        directory,
        repair=(
            f"rm -rf {shlex.quote(str(directory))} && "
            f"{xref_prepare_command(*_names(row), evidence=evidence)}"
        ),
    )
    if existing is not None:
        return existing
    reader = _reader(row)
    clear_work_dir(directory)
    work = work_dir(directory)
    fetched = _fetch(row, work=work, evidence=evidence, progressbar=progressbar)
    digest = hashlib.new(_algorithm(row))
    triples = reader(
        _unpacked_lines(fetched, digest),
        ncbi_taxid=row.ncbi_taxid,
        origin=str(fetched),
        evidence=evidence,
    )
    _check_source_checksum(row, digest.hexdigest(), path=fetched)
    if not triples:
        raise XrefTableError(
            f"{fetched} carries no row for {row.species!r} (NCBITaxon:{row.ncbi_taxid}), so "
            f"the slice would be empty and every query would answer nothing. The file at "
            f"{row.url} is not the one this row pins, or the publisher dropped the species. "
            f"Fix the row in the curated xref metadata table."
        )
    staged = work / path.name
    sha256 = _write_slice(staged, triples)
    directory.mkdir(parents=True, exist_ok=True)
    staged.replace(path)
    record = build_record(
        directory,
        kind=XREF_KIND,
        name=f"{row.source}/{row.release}/{species_slug(row.species)}{_suffix(evidence)}",
        files=[path],
        source_url=row.url,
        sha256=sha256,
        details={
            "species": row.species,
            "ncbi_taxid": row.ncbi_taxid,
            "source": row.source,
            "release": row.release,
            "publisher": row.publisher,
            "version": row.version,
            # The publisher's own checksum of its own file, kept as *provenance*: what is
            # stored here is a derived slice, so `sha256` above — this slice's own digest —
            # is what the set is held to, and this says which bytes it was cut from.
            "source_checksum": row.source_checksum,
            "evidence": list(evidence),
            "namespaces": sorted({namespace for namespace, _id, _stem in triples}),
            "symbol_kinds": [
                kind
                for kind in SYMBOL_KINDS
                if KIND_NAMESPACES[kind] in {namespace for namespace, _id, _stem in triples}
            ],
            "rows": len(triples),
        },
    )
    write_record(directory, record)
    clear_work_dir(directory)
    return record


def _reader(row: XrefMetadata) -> XrefReader:
    """Return the reader for this row's **Xref source**, or say which sources have one."""
    reader = _READERS.get(row.source)
    if reader is None:
        raise NoXrefSetError(
            f"the curated table lists the source {row.source!r} and this package ships no "
            f"reader for it: the sources it can read are {', '.join(sorted(_READERS))}. "
            f"Adding one is a module and an entry in the reader table."
        )
    return reader


def _fetch(row: XrefMetadata, *, work: Path, evidence: tuple[str, ...], progressbar: bool) -> Path:
    """Download the publisher's file into the working area, or say what to do instead.

    Nothing is verified here, and ``known_hash`` is deliberately not passed. **The checksum
    a row pins is over the publisher's unpacked bytes** (ADR-0006) while what lands is the
    compressed file, so pooch would compare the pin against the wrong bytes and reject
    every download. It is checked instead as the file is read — see :func:`_unpacked_lines`
    and :func:`_check_source_checksum` — which costs one pass rather than two.
    """
    name = Path(urlparse(row.url).path).name or f"{row.source}-{row.release}"
    try:
        return fetch.fetch_url(row.url, work, fname=name, progressbar=progressbar)
    except (OSError, ValueError) as error:
        raise XrefSetNotDownloadedError(
            f"the {row.source} {row.release} xref set for {row.species!r} is not prepared "
            f"here and {row.url} could not be fetched: {error}. Nothing else in this "
            f"package needs the network, so this is the one step that does. Prepare it on a "
            f"machine with internet — a login node, since the lab's compute nodes have "
            f"none — with `{xref_prepare_command(*_names(row), evidence=evidence)}`."
        ) from error


def _algorithm(row: XrefMetadata) -> str:
    """Return the hash algorithm the row's pinned checksum names, e.g. ``md5``."""
    algorithm, separator, _digest = row.source_checksum.partition(":")
    if not separator:
        raise XrefTableError(
            f"the curated xref row for {row.species!r} pins the checksum "
            f"{row.source_checksum!r}, and a checksum here is spelled "
            f"'<algorithm>:<hexdigest>' so that the algorithm travels with it. Fix that "
            f"cell in the xref metadata table."
        )
    return algorithm


def _unpacked_lines(path: Path, digest: Hash) -> Iterator[str]:
    """Yield the file's unpacked lines, hashing the bytes exactly as the publisher did.

    Streamed a line at a time and never held whole: the publisher's file unpacks to over
    half a gigabyte. Gzip is undone here rather than in a reader, so a source that ships
    plain text needs no branch of its own.
    """
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as packed:
            yield from _hashed(packed, digest)
    else:
        with path.open("rb") as plain:
            yield from _hashed(plain, digest)


def _hashed(lines: Iterable[bytes], digest: Hash) -> Iterator[str]:
    """Decode each line, feeding the bytes to ``digest`` exactly as they were read."""
    for raw in lines:
        digest.update(raw)
        yield raw.decode("utf-8")


def _check_source_checksum(row: XrefMetadata, found: str, *, path: Path) -> None:
    """Hold the publisher's file to the digest the curated row pins for it."""
    _algorithm_name, _, expected = row.source_checksum.partition(":")
    if found != expected:
        raise XrefTableError(
            f"{path} hashes to {found} where the curated xref row for {row.species!r} pins "
            f"{expected}. A truncated download is not a smaller release, and slicing it "
            f"would answer a query with silently fewer genes — delete {path} and construct "
            f"the set again to fetch it afresh. If {row.url} has genuinely been "
            f"re-published under the same name, that source cannot be pinned and does not "
            f"belong in the table (ADR-0018)."
        )


def _write_slice(path: Path, triples: tuple[tuple[str, str, str], ...]) -> str:
    """Write the slice as a plain gzipped TSV and return its **unpacked** sha256.

    The gzip is stamped with no modification time, so one release sliced on two machines
    produces byte-identical files. The digest is of what is *inside* the gzip (ADR-0006),
    so a slice recompressed elsewhere still matches its record.
    """
    payload = "".join(
        ["\t".join(SLICE_COLUMNS) + "\n", *(f"{a}\t{b}\t{c}\n" for a, b, c in triples)]
    ).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
        out.write(payload)
    return hashlib.sha256(payload).hexdigest()


def _read_slice(path: Path) -> tuple[tuple[tuple[str, str, str], ...], str]:
    """Return the stored slice's triples and the sha256 of its **unpacked** bytes.

    One pass: the digest that holds the file to its record costs nothing extra, since the
    bytes are being read anyway.
    """
    with gzip.open(path, "rb") as handle:
        payload = handle.read()
    return parse_slice(payload.decode("utf-8"), origin=str(path)), hashlib.sha256(
        payload
    ).hexdigest()


def _index(
    triples: tuple[tuple[str, str, str], ...],
) -> tuple[dict[str, _Index], dict[str, _Index]]:
    """Return the two directions of the hop, one index per **Namespace**, ids ascending.

    Both are built together off one pass, and both default to an empty tuple for a key they
    do not hold — so a verb reads the index once and never asks whether a key is in it.
    """
    to_stems: dict[str, dict[str, list[str]]] = {}
    from_stems: dict[str, dict[str, list[str]]] = {}
    for namespace, identifier, stem in triples:
        to_stems.setdefault(namespace, {}).setdefault(identifier, []).append(stem)
        from_stems.setdefault(namespace, {}).setdefault(stem, []).append(identifier)
    return (
        {namespace: _Index(index) for namespace, index in to_stems.items()},
        {namespace: _Index(index) for namespace, index in from_stems.items()},
    )


def _symbol_index(
    to_stems: Mapping[str, _Index],
) -> tuple[dict[str, tuple[SymbolMatch, ...]], dict[str, tuple[SymbolMatch, ...]]]:
    """Return the exact and the case-folded symbol lookups, matches in kind order.

    Both are built from the three stored symbol **Namespace**s in one pass, and both hold
    the same matches: folding is a property of the *lookup* and never of what was stored,
    so the insensitive path answers with every gene matched rather than a folded set having
    quietly merged two spellings into one entry on the way in.

    Empty for a set whose source carries no symbols, which costs such a set nothing.
    """
    exact: dict[str, list[SymbolMatch]] = {}
    folded: dict[str, list[SymbolMatch]] = {}
    for kind in SYMBOL_KINDS:
        index = to_stems.get(KIND_NAMESPACES[kind])
        if index is None:
            continue
        for spelling in index:
            for stem in index[spelling]:
                match = SymbolMatch(symbol=spelling, gene_id_stem=stem, kind=kind)
                exact.setdefault(spelling, []).append(match)
                folded.setdefault(fold_symbol(spelling), []).append(match)
    return (
        {key: tuple(matches) for key, matches in exact.items()},
        {key: tuple(matches) for key, matches in folded.items()},
    )


class _Index(Mapping[str, tuple[str, ...]]):
    """One namespace's lookups, answering an unknown key with an empty tuple.

    A mapping rather than a plain ``dict`` so that *no answer* has one spelling: a verb
    reads ``index[key]`` and decides between :attr:`resolved` and :attr:`unresolved` on
    whether what came back is empty, instead of branching on membership first and
    remembering to keep the two branches in step.
    """

    def __init__(self, values: dict[str, list[str]]) -> None:
        self._values = {key: tuple(sorted(set(found))) for key, found in values.items()}

    def __getitem__(self, key: str) -> tuple[str, ...]:
        """Return what ``key`` names, ascending — empty when it names nothing."""
        return self._values.get(key, ())

    def __iter__(self) -> Iterator[str]:
        """Iterate the keys this index holds."""
        return iter(self._values)

    def __len__(self) -> int:
        """Return how many keys this index holds."""
        return len(self._values)
