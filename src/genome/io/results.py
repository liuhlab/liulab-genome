"""What the API answers with: the values it returns and ``--json`` writes.

Registering something in this package — an **Assembly**, an annotation — or asking it a
question about one ends by returning a value that says what it found, and every surface
reads that same value: a script keeps it, the CLI prints it, ``--json`` serializes it.
Those values live here, together, because they change for one reason: *what an answer must
be able to say*. The modules that do the work import them; nothing here does any of it,
reads a directory, or reaches the network.

**A shipped table's answer is not one of them.** A **TF gene list** and a **TF cofactor
list** say what one published file holds, so they live with the code that reads that file,
under :mod:`genome.tf`, and change when a publisher's columns do rather than when a
registration's report does. Keeping them out is what leaves this module importing nothing
of that half: registering an annotation must not cost a caller the census reader, the
cofactor reader and the whole motif tree behind them.

Two of them carry a **Completion marker** whole rather than copying it out field by
field — :class:`RegisteredAssembly` and :class:`RegisteredAnnotation` — so every later
question is answered off the record in hand instead of by opening the directory again. The
rest are reports assembled from more than one source: :class:`VerifiedAssembly` sets a
digest against what pinned it, :class:`AnnotationStatus` sets what the annotation table
offers against what this machine holds, one :class:`AnnotationStatusRow` per name,
:class:`GeneList` sets one **Gene category**'s genes against the **Curated gene list**s
that contributed them, one :class:`GeneListSource` apiece, :class:`ResolvedGeneIds`
sets the **Gene id stem**s a caller asked about against the gene ids one annotation
actually carries — including the stems it carries none for. :class:`ResolvedStems` and
:class:`ResolvedXrefIds` are the two directions of an **Xref set**'s hop, and they are
:class:`ResolvedGeneIds`'s shape again with one difference: what produced the answer is a
publisher and a **Release** rather than an assembly and an annotation.
:class:`ResolvedSymbols` is the third direction of that hop and the one that is not a
mirror of the other two: a symbol matches approved, previous and alias spellings, so each
hit is a :class:`SymbolMatch` carrying the kind it matched rather than a bare string.

**An answer knows how it reads.** :attr:`AnnotationStatusRow.state` and
:attr:`AnnotationStatus.default_summary` are here rather than in the surface printing them,
so the precedence between *broken* and *registered* is stated once — where the invariant
lives — and the command a report names is built by :func:`annotation_register_command`
rather than concatenated wherever it is wanted. ``as_json`` is the same rule for the
machine-readable half: the keys are the ones written on disk, in the order the record
declares them, and no surface respells them.

Examples
--------
>>> from genome.io.results import chromosome_check_summary
>>> chromosome_check_summary({"chromosomes_checked": True})
'chromosomes checked — every name the GTF uses is one the assembly carries'
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from genome.io.completion import CompletionRecord
from genome.io.components import ChimeraDetails

#: What answered *which digest should this FASTA have?* — the assembly's curated
#: metadata row, or the completion record its own registration wrote. Reported by
#: :func:`~genome.io.download.verify_assembly` so that being held to a pin and being held
#: only to what this machine last produced are never read as the same result. Public
#: because the CLI keys its two sentences on them: a surface that spelled the strings
#: again would print the raw status the day one of these was renamed, rather than failing.
EXPECTED_FROM_TABLE = "table"
EXPECTED_FROM_RECORD = "record"

#: What every one of Compara's speciation labels begins with, and nothing else does:
#: ``ortholog_one2one``, ``ortholog_one2many``, ``ortholog_many2many``. It is how a
#: **Homology link** tells an ortholog from a **Paralogy link** without this package
#: keeping a list of the publisher's labels that would go stale the release it added one
#: — the prefix is Compara's own naming and a duplication label never carries it.
ORTHOLOG_TYPE_PREFIX = "ortholog_"

#: What ``details["chromosomes_unchecked_because"]`` says when the caller stood the check
#: down — ``check_chromosomes=False``, or ``--no-check-chromosomes`` from a shell. There is
#: no advice to give about it: the assembly may be registered and the names deliberately
#: accepted, so all a surface can say is what the record therefore does not vouch for.
UNCHECKED_CALLER_OVERRIDE = "caller-override"

#: …and when the check was asked for but had nothing to run against, the assembly having
#: no ``chrom.sizes`` yet. Registering the assembly is what makes the check possible, so
#: this is the one of the two states where saying so is useful advice rather than noise.
UNCHECKED_NO_CHROM_SIZES = "no-chrom-sizes"

#: What each state of the chromosome check reads as, one sentence apiece — including the
#: one where it ran and passed, since a surface that says nothing about it reads as a pass.
#: Keyed by ``details["chromosomes_unchecked_because"]``; ``None`` is the check that ran.
_CHECK_SUMMARIES = {
    None: "chromosomes checked — every name the GTF uses is one the assembly carries",
    UNCHECKED_NO_CHROM_SIZES: (
        "chromosomes not checked — nothing to check against; register the assembly first "
        "to verify them"
    ),
    UNCHECKED_CALLER_OVERRIDE: (
        "chromosomes not checked — the check was stood down, so the record does not vouch "
        "for the names"
    ),
}

#: What a record written before the reason was recorded reads as. Its bare ``False`` was
#: written for either reason and nothing on disk says which, so it is reported as neither.
_UNKNOWN_REASON_SUMMARY = (
    "chromosomes not checked — this record does not say why, so whether the names match "
    "the assembly is unknown"
)


def annotation_register_command(assembly: str, name: str) -> str:
    """Return the command that registers the annotation ``name`` for ``assembly``.

    The one spelling of it. Errors quote it, the repair adds ``--force`` to it, and
    :attr:`AnnotationStatus.default_summary` names it for a **Default annotation** nobody
    has fetched yet — so a renamed command is renamed once.

    Parameters
    ----------
    assembly : str
        The **Assembly** the annotation belongs to, e.g. ``"hg38"``.
    name : str
        The **Registered name** to address it by.

    Returns
    -------
    str
        A shell command, unquoted and unfenced — the caller decides how to set it.

    Examples
    --------
    >>> annotation_register_command("hg38", "gencode_v50")
    'genome register-annotation hg38 gencode_v50'
    """
    return f"genome register-annotation {assembly} {name}"


def chromosome_check_summary(details: Mapping[str, Any]) -> str:
    """Return the one line a surface prints about an annotation's chromosome-name check.

    Four states, four sentences, and one of them is always returned: the check ran and
    passed; it had nothing to run against, and registering the assembly is what fixes
    that; the caller stood it down, which is not something to advise about; or the record
    does not say which, and none of the three may be claimed. Silence is not a fifth
    state — a surface that prints nothing about the check reads as one that passed.

    ``details`` is a registration record's ``details``; a caller holding what a
    registration answered with asks :attr:`RegisteredAnnotation.chromosome_check` instead
    and never spells the two fields. Those are ``chromosomes_checked`` — the check ran and
    the GTF's names were all the assembly's — and ``chromosomes_unchecked_because``, which
    says which of the two reasons it did not, and is ``None`` when it did.

    A record written before the second field existed carries a bare
    ``chromosomes_checked: false`` that was written for either reason, and nothing on disk
    says which. It reads as *unknown* rather than as either one, and rather than raising:
    the reason is a fact that was never gathered, which is what an absent entry in
    ``tool_versions`` means too.

    Parameters
    ----------
    details : mapping of str to object
        A registration record's ``details``. Anything else it holds is ignored, and a
        mapping holding neither field reads as unknown.

    Returns
    -------
    str
        One sentence, with no trailing punctuation and no leading indent — the caller
        decides how to set it.

    Examples
    --------
    >>> chromosome_check_summary({"chromosomes_checked": True})
    'chromosomes checked — every name the GTF uses is one the assembly carries'
    >>> print(chromosome_check_summary({"chromosomes_unchecked_because": "caller-override"}))
    chromosomes not checked — the check was stood down, so the record does not vouch for the names
    """
    if details.get("chromosomes_checked") is True:
        return _CHECK_SUMMARIES[None]
    # Anything else — the field absent, or a reason a later version writes and this one
    # has never heard of — is a reason that cannot be reported, which is the unknown.
    because = details.get("chromosomes_unchecked_because")
    if isinstance(because, str) and because in _CHECK_SUMMARIES:
        return _CHECK_SUMMARIES[because]
    return _UNKNOWN_REASON_SUMMARY


@dataclass(frozen=True)
class RegisteredAssembly:
    """What preparing an assembly on disk produced: its record, and where that landed.

    :func:`~genome.io.download.register_assembly`'s answer — what ``genome register``
    prints, and what its ``--json`` serializes. The **Completion marker** the run wrote
    *is* the answer, so it is carried whole rather than copied out field by field, and the
    two questions a surface then asks — which files are claimed, and is this a **Chimera**
    — are answered from that one record instead of by reading the directory again.

    Attributes
    ----------
    assembly : str
        The **Assembly** that was registered, under the name the caller asked for.
    directory : pathlib.Path
        Its **Assembly dir** — where those files and that record are.
    record : genome.io.completion.CompletionRecord
        The record the registration wrote, read back.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io.completion import CompletionRecord
    >>> registered = RegisteredAssembly(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38"),
    ...     record=CompletionRecord(
    ...         kind="genome",
    ...         name="hg38",
    ...         files={"hg38.fa.fai": 21, "hg38.fa": 12},
    ...         source_url="https://example.org/hg38.fa.gz",
    ...         sha256="1a2b3c",
    ...         tool_versions={},
    ...         package_version="2026.8.0",
    ...         completed_at="2026-08-12T09:00:00+00:00",
    ...         details={},
    ...     ),
    ... )
    >>> registered.file_names
    ['hg38.fa', 'hg38.fa.fai']
    >>> registered.chimera is None
    True
    >>> registered.as_json()["directory"]
    '/data/genome/hg38'
    """

    assembly: str
    directory: Path
    record: CompletionRecord

    @property
    def source_url(self) -> str | None:
        """Where the bytes were fetched from, or ``None`` when nothing was — a chimera's."""
        return self.record.source_url

    @property
    def sha256(self) -> str | None:
        """Digest of the unpacked FASTA, or ``None`` when none was computed."""
        return self.record.sha256

    @property
    def file_names(self) -> list[str]:
        """Every file the record claims, sorted — a fresh list each call."""
        return sorted(self.record.files)

    @property
    def chimera(self) -> ChimeraDetails | None:
        """What the build recorded about its components, or ``None`` for anything else.

        The record is what says an assembly is a **Chimera**, here as everywhere else —
        and the record is already in hand, so a surface reporting the registration that
        just happened never reads the same file a second time to find out.
        """
        return ChimeraDetails.from_record(self.record)

    def as_json(self) -> dict[str, Any]:
        """Return this registration as ``--json`` serializes it.

        The record's own fields under the record's own names, then the ``assembly`` asked
        for and the ``directory`` it landed in — the two facts a record does not hold
        about itself. The names are the ones written on disk and are never respelled here.

        Returns
        -------
        dict
            The record's fields, followed by ``assembly`` and ``directory``.
        """
        return {**asdict(self.record), "assembly": self.assembly, "directory": str(self.directory)}


@dataclass(frozen=True)
class VerifiedAssembly:
    """What re-reading a FASTA proved: its digest, what that was held to, and the components.

    :func:`~genome.io.download.verify_assembly`'s answer, and three results a caller must
    be able to tell apart, so each is a field of its own: the digest computed, *what
    supplied* the digest it was held to — being held to the lab's pin and being held to
    what this machine last produced are different answers — and, for a **Chimera**, what
    comparing its components settled. A digest that disagreed raises rather than arriving
    here, so this is what nothing refused.

    Attributes
    ----------
    assembly : str
        The **Assembly** whose row supplied the digest to check against.
    fasta : pathlib.Path
        The file that was read.
    sha256 : str
        The digest computed over it.
    expected : str or None
        The digest it was held to, or ``None`` when nothing pinned one.
    expected_from : str or None
        What answered with ``expected`` — :data:`EXPECTED_FROM_TABLE`,
        :data:`EXPECTED_FROM_RECORD`, or ``None`` when nothing did.
    components : str or None
        :data:`~genome.io.components.COMPONENTS_UNCHANGED` or
        :data:`~genome.io.components.COMPONENTS_UNKNOWN` for a chimera, and ``None`` for
        anything else — including every ``fasta`` checked on its own.

    Examples
    --------
    >>> from pathlib import Path
    >>> checked = VerifiedAssembly(
    ...     assembly="sacCer3",
    ...     fasta=Path("/data/genome/sacCer3/sacCer3.fa"),
    ...     sha256="6ff72f07",
    ...     expected="6ff72f07",
    ...     expected_from=EXPECTED_FROM_TABLE,
    ...     components=None,
    ... )
    >>> checked.verified
    True
    >>> checked.as_json()["expected_from"]
    'table'
    """

    assembly: str
    fasta: Path
    sha256: str
    expected: str | None
    expected_from: str | None
    components: str | None

    @property
    def verified(self) -> bool:
        """Whether there was a digest to check against at all, rather than merely one computed."""
        return self.expected is not None

    def as_json(self) -> dict[str, Any]:
        """Return this verification as ``--json`` serializes it.

        Returns
        -------
        dict
            Every attribute above, with ``fasta`` rendered as text and :attr:`verified`
            written out beside the fields it is read from.
        """
        return {
            "assembly": self.assembly,
            "fasta": str(self.fasta),
            "sha256": self.sha256,
            "expected": self.expected,
            "expected_from": self.expected_from,
            "verified": self.verified,
            "components": self.components,
        }


@dataclass(frozen=True)
class RegisteredAnnotation:
    """What registering one annotation produced: its record, and where that landed.

    :func:`~genome.io.gtf.register_annotation`'s answer and
    :func:`~genome.io.gtf.register_gtf`'s — what ``genome register-annotation`` and
    ``genome register-gtf`` print, and what their ``--json`` serializes. A :class:`~genome.io.gtf.GtfAnnotation` says where an annotation's two
    files are; this says what the run that wrote them did, which is the **Completion
    marker** itself, carried whole. Every question a surface then asks — the digest, the
    source, the files claimed, whether the chromosome names were actually checked — is
    answered from that one record rather than by reading the directory again.

    Attributes
    ----------
    assembly : str
        The **Assembly** the annotation belongs to. It is not in the record, which names
        the annotation rather than what it annotates.
    directory : pathlib.Path
        The annotation's own directory, ``<assembly dir>/gtf/<name>/``.
    record : genome.io.completion.CompletionRecord
        The record the registration wrote, read back.

    Examples
    --------
    >>> from pathlib import Path
    >>> from genome.io.completion import CompletionRecord
    >>> registered = RegisteredAnnotation(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38/gtf/gencode_v50"),
    ...     record=CompletionRecord(
    ...         kind="annotation",
    ...         name="gencode_v50",
    ...         files={"gencode_v50.gtf": 12, "gencode_v50.db": 34},
    ...         source_url="https://example.org/gencode_v50.gtf.gz",
    ...         sha256="1a2b3c",
    ...         tool_versions={},
    ...         package_version="2026.8.0",
    ...         completed_at="2026-08-12T09:00:00+00:00",
    ...         details={"chromosomes_checked": True},
    ...     ),
    ... )
    >>> registered.name, registered.file_names
    ('gencode_v50', ['gencode_v50.db', 'gencode_v50.gtf'])
    >>> print(registered.chromosome_check)
    chromosomes checked — every name the GTF uses is one the assembly carries
    """

    assembly: str
    directory: Path
    record: CompletionRecord

    @property
    def name(self) -> str:
        """The **Registered name** it is addressed by — the record's own name."""
        return self.record.name

    @property
    def source_url(self) -> str | None:
        """The URL fetched, or the path a GTF was handed over at; ``None`` for a merge."""
        return self.record.source_url

    @property
    def sha256(self) -> str | None:
        """Digest of the placed GTF, or ``None`` when none was computed."""
        return self.record.sha256

    @property
    def file_names(self) -> list[str]:
        """Every file the record claims, sorted — a fresh list each call."""
        return sorted(self.record.files)

    @property
    def chromosome_check(self) -> str:
        """The one line saying what the chromosome-name check settled for this annotation.

        :func:`chromosome_check_summary` over the record this registration wrote, so the
        surface that prints it never reads the record's own keys. Always a sentence:
        silence would read as a pass.
        """
        return chromosome_check_summary(self.record.details)

    def as_json(self) -> dict[str, Any]:
        """Return this registration as ``--json`` serializes it.

        The record's own fields under the record's own names, then the ``assembly`` it
        belongs to and the ``directory`` it landed in — the two facts a record does not
        hold about itself. The names are the ones written on disk and are never respelled
        here.

        Returns
        -------
        dict
            The record's fields, followed by ``assembly`` and ``directory``.
        """
        return {**asdict(self.record), "assembly": self.assembly, "directory": str(self.directory)}


@dataclass(frozen=True)
class GeneListSource:
    """One **Curated gene list** that contributed to an answer, and what it contributed.

    What makes a **Merged annotation**'s genes attributable: one of these per contributing
    annotation, so a caller counting worm ribosomal RNA can drop the *E. coli* entry
    rather than being handed one number it cannot take apart. An annotation that is not a
    merge has exactly one, whose ``component`` is ``None``.

    ``description`` and ``source`` travel with the ids rather than being looked up
    separately, because they are what says whether these ids mean what the caller's metric
    needs: two annotations spelling a category the same way need not have curated it the
    same way.

    Attributes
    ----------
    component : str or None
        The **Component** assembly whose genes these are, for a contributor to a **Merged
        annotation**; ``None`` for anything else.
    annotation : str
        The **Registered name** of the contributing annotation.
    description : str
        What membership in this category means for that annotation.
    source : str
        Where that membership came from, and the caveats on using it.
    gene_ids : tuple of str
        The gene ids it contributed, in the order its curated list lists them.

    Examples
    --------
    >>> source = GeneListSource(
    ...     component="ce11",
    ...     annotation="wormbase_ws298",
    ...     description="the mature ribosomal RNA genes",
    ...     source="WormBase WS298 gene_biotype",
    ...     gene_ids=("WBGene00004512", "WBGene00004513"),
    ... )
    >>> source.as_json()["component"]
    'ce11'
    >>> len(source.gene_ids)
    2
    """

    component: str | None
    annotation: str
    description: str
    source: str
    gene_ids: tuple[str, ...]

    def as_json(self) -> dict[str, Any]:
        """Return this contribution as ``--json`` serializes it: every attribute, in order.

        Returns
        -------
        dict
            The fields above, under their own names, with ``gene_ids`` as a list.
        """
        return asdict(self)


@dataclass(frozen=True)
class GeneList:
    """The genes one annotation puts in one **Gene category**, attributed to their sources.

    :meth:`~genome.io.gtf.AnnotationRegistry.gene_list`'s answer — what ``genome
    gene-list`` prints and what its ``--json`` serializes. There is no empty one: an
    annotation that declares no categories, and one that declares categories but not this
    one, each raise a :class:`LookupError` of their own rather than answering with
    nothing, so holding one of these means the category was really declared and really has
    genes in it.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** asked about — the merged name for a **Merged annotation**,
        whose contributors are named in :attr:`sources`.
    category : str
        The **Gene category**, as the curated lists spell it.
    sources : tuple of GeneListSource
        One entry per contributing **Curated gene list**, in contributor order. Never
        empty. A contributor that does not declare this category is simply absent — a
        bacterium has no mitochondria, and that is not a failure to report.

    Examples
    --------
    >>> genes = GeneList(
    ...     assembly="ce11",
    ...     annotation="wormbase_ws298",
    ...     category="rRNA",
    ...     sources=(
    ...         GeneListSource(None, "wormbase_ws298", "rRNA genes", "WormBase", ("a", "b")),
    ...     ),
    ... )
    >>> genes.gene_ids
    ['a', 'b']
    >>> genes.as_json()["category"]
    'rRNA'
    """

    assembly: str
    annotation: str
    category: str
    sources: tuple[GeneListSource, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every source's gene ids, concatenated in source order — a fresh list each call.

        **Concatenated and not de-duplicated.** A merge rewrites only the seqname and
        never the ``gene_id``, so two components carrying the same id would be a real
        ambiguity in the merged annotation, and collapsing it here would hide exactly that
        — a caller summing over these would silently under-count one of the two. Where
        that matters, :attr:`sources` says which contributor each id came from.
        """
        return [gene_id for source in self.sources for gene_id in source.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``category``, the concatenated ``gene_ids``, and
            ``sources`` as a list of :meth:`GeneListSource.as_json` entries.
            :attr:`gene_ids` is written out beside the sources it is read from rather than
            left to the reader: assembling it is where a reader would reach for a set and
            de-duplicate, which is the one thing this answer must not do.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "category": self.category,
            "gene_ids": self.gene_ids,
            "sources": [source.as_json() for source in self.sources],
        }


@dataclass(frozen=True)
class ResolvedGeneIds:
    """The gene ids one **Annotation** carries for the **Gene id stem**s it was asked about.

    :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`'s answer. A stem is a gene
    id with its version dropped, and inside one annotation it may name more than one gene
    id — nine do in ``gencode_v50lift37``, eight of them pseudoautosomal-Y — so every stem
    answers with **all** of them and nothing here picks one. Two stems never name the same
    gene id, since an id has exactly one stem.

    **What was asked about and is not there rides back on the answer.** A caller holding a
    few thousand stems gets the ones this annotation carries no gene for in
    :attr:`unresolved` rather than a shorter list than it passed, so what the thing it was
    holding contains and this annotation does not is visible instead of dropped.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    resolved : mapping of str to tuple of str
        Every stem that named at least one gene id, in the order the stems were asked
        about, to the ids it names, in ascending order. No value is ever an empty tuple —
        a stem that named nothing is in :attr:`unresolved` instead.
    unresolved : tuple of str
        The stems no gene id in the annotation is of, in the order they were asked about.

    Examples
    --------
    >>> answer = ResolvedGeneIds(
    ...     assembly="hg19",
    ...     annotation="gencode_v50lift37",
    ...     resolved={
    ...         "ENSG00000182378": ("ENSG00000182378.14", "ENSG00000182378.14_PAR_Y"),
    ...         "ENSG00000141510": ("ENSG00000141510.18",),
    ...     },
    ...     unresolved=("ENSG00000288541",),
    ... )
    >>> answer.gene_ids
    ['ENSG00000182378.14', 'ENSG00000182378.14_PAR_Y', 'ENSG00000141510.18']
    >>> answer.as_json()["unresolved"]
    ['ENSG00000288541']
    """

    assembly: str
    annotation: str
    resolved: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id resolved, stem order and then id order — a fresh list each call.

        **Every** id, not one per stem. Flattening is exactly where a reader would take
        the first id of each stem and lose the second, which is the pseudoautosomal gene
        this answer's shape exists to keep; :attr:`resolved` is what says which stem an id
        came from.
        """
        return [gene_id for ids in self.resolved.values() for gene_id in ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``resolved`` as a mapping of stem to a list of
            gene ids, ``unresolved`` as a list, and the flattened ``gene_ids``. The last
            is written out beside the mapping it is read from for the reason
            :attr:`gene_ids` gives: a reader assembling it is a reader who might take one
            id per stem.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "resolved": {stem: list(ids) for stem, ids in self.resolved.items()},
            "unresolved": list(self.unresolved),
            "gene_ids": self.gene_ids,
        }


@dataclass(frozen=True)
class ResolvedStems:
    """The **Gene id stem**s one **Xref set** says a foreign id names.

    :meth:`~genome.xref.xref.XrefSet.to_stems`'s answer — the hop *toward* the hub. Every
    field before :attr:`resolved` says what produced it, because one publisher's
    assertions are not another's: NCBI and Ensembl agree on 57.6% of human gene-level
    (GeneID, ENSG) pairs, so an answer that did not name its **Xref source** and
    **Release** would be unreproducible a year later. A query reads exactly one set
    (ADR-0017), which is why the source is one field here rather than a column on every
    row.

    **A foreign id naming two stems answers with both**, and nothing picks one — the same
    guarantee :class:`ResolvedGeneIds` gives for a stem naming two gene ids. **What named
    nothing rides back** in :attr:`unresolved` rather than shortening the answer.

    The keys are the caller's **own spelling** of the ids it asked about, so a versioned
    and an unversioned spelling of one id are two keys with identical values and the
    answer still zips against the caller's table row for row.

    Attributes
    ----------
    species : str
        The species this set is for, as the curated metadata table spells it.
    source : str
        The **Xref source** whose assertions these are.
    release : str
        The pinned **Release** of that source.
    namespace : str
        The **Namespace** the ids asked about belong to.
    resolved : mapping of str to tuple of str
        Every id that named at least one stem, in the order they were asked about, to the
        stems it names, in ascending order. No value is ever an empty tuple — an id that
        named nothing is in :attr:`unresolved` instead.
    unresolved : tuple of str
        The ids this release names no stem for, in the order they were asked about.

    Examples
    --------
    >>> answer = ResolvedStems(
    ...     species="Homo sapiens",
    ...     source="alliance",
    ...     release="8.4.0",
    ...     namespace="entrez",
    ...     resolved={"7157": ("ENSG00000141510",)},
    ...     unresolved=("999999999",),
    ... )
    >>> answer.gene_id_stems
    ['ENSG00000141510']
    >>> answer.as_json()["source"]
    'alliance'
    """

    species: str
    source: str
    release: str
    namespace: str
    resolved: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]

    @property
    def gene_id_stems(self) -> list[str]:
        """Every stem resolved, ask order and then stem order — a fresh list each call.

        **Every** stem, not one per id. Flattening loses which id named which stem, and
        with it the fact that an id named more than one: a reader taking the first stem of
        each id would silently pick one of two genes a **Namespace** is ambiguous
        between. :attr:`resolved` is what says which id a stem came from. It also loses
        the ask order *of the ids*, since one id contributing two stems contributes two
        entries here.
        """
        return [stem for stems in self.resolved.values() for stem in stems]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``species``, ``source``, ``release`` and ``namespace``, ``resolved`` as a
            mapping of id to a list of stems, ``unresolved`` as a list, and the flattened
            ``gene_id_stems``. The last is written out beside the mapping it is read from
            for the reason :attr:`gene_id_stems` gives.
        """
        return {
            "species": self.species,
            "source": self.source,
            "release": self.release,
            "namespace": self.namespace,
            "resolved": {asked: list(stems) for asked, stems in self.resolved.items()},
            "unresolved": list(self.unresolved),
            "gene_id_stems": self.gene_id_stems,
        }


@dataclass(frozen=True)
class ResolvedXrefIds:
    """The foreign ids one **Xref set** says a **Gene id stem** names.

    :meth:`~genome.xref.xref.XrefSet.from_stems`'s answer — the hop *away* from the hub,
    and :class:`ResolvedStems`'s mirror in every respect: the same four provenance fields,
    the same ask order, the same never-empty resolved value, and the same tuple of what
    named nothing. Two verbs and only two, so a caller wanting one **Namespace** from
    another makes both calls and owns the join (ADR-0017).

    Attributes
    ----------
    species : str
        The species this set is for, as the curated metadata table spells it.
    source : str
        The **Xref source** whose assertions these are.
    release : str
        The pinned **Release** of that source.
    namespace : str
        The **Namespace** the answering ids belong to.
    resolved : mapping of str to tuple of str
        Every stem that named at least one id in that namespace, in the order the stems
        were asked about, to the ids it names, in ascending order. No value is ever an
        empty tuple.
    unresolved : tuple of str
        The stems this release gives no id in that namespace, in the order they were asked
        about. One bucket and not two: a stem this release never carried and a stem it
        carries with no id in *this* namespace are both *this set answers nothing*, and no
        id history is held that could tell a retired stem from an unknown one (ADR-0017).

    Examples
    --------
    >>> answer = ResolvedXrefIds(
    ...     species="Homo sapiens",
    ...     source="alliance",
    ...     release="8.4.0",
    ...     namespace="hgnc",
    ...     resolved={"ENSG00000141510": ("HGNC:11998",)},
    ...     unresolved=("ENSG00000288541",),
    ... )
    >>> answer.xref_ids
    ['HGNC:11998']
    >>> answer.as_json()["namespace"]
    'hgnc'
    """

    species: str
    source: str
    release: str
    namespace: str
    resolved: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]

    @property
    def xref_ids(self) -> list[str]:
        """Every foreign id resolved, ask order and then id order — a fresh list each call.

        **Every** id, not one per stem. Flattening loses which stem named which id, so a
        reader taking the first id of each stem would hand a collaborator one of two
        accessions a gene genuinely has without saying it had chosen; and a stem naming two
        ids contributes two entries, so the flattened list no longer runs parallel to the
        stems asked about. :attr:`resolved` is what says which stem an id came from.
        """
        return [xref_id for ids in self.resolved.values() for xref_id in ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``species``, ``source``, ``release`` and ``namespace``, ``resolved`` as a
            mapping of stem to a list of ids, ``unresolved`` as a list, and the flattened
            ``xref_ids``, written out for the reason :attr:`xref_ids` gives.
        """
        return {
            "species": self.species,
            "source": self.source,
            "release": self.release,
            "namespace": self.namespace,
            "resolved": {stem: list(ids) for stem, ids in self.resolved.items()},
            "unresolved": list(self.unresolved),
            "xref_ids": self.xref_ids,
        }


@dataclass(frozen=True)
class SymbolMatch:
    """One hit of a gene symbol against an **Xref set**, and which spelling matched it.

    The kind rides on the match rather than being filtered away on the way out, because a
    table that spells a gene the way it was spelled five years ago is otherwise dropped
    without a word — the failure that would have hit 31 of EpiFactors' 801 rows.

    Attributes
    ----------
    symbol : str
        The **authority's own** spelling that matched, which is not always the one asked
        about: on the case-insensitive path ``brca1`` matches and this says ``BRCA1``.
    gene_id_stem : str
        The **Gene id stem** that spelling names.
    kind : str
        ``approved``, ``previous`` or ``alias`` — see :mod:`genome.xref.symbols`.

    Examples
    --------
    >>> match = SymbolMatch(symbol="ARNTL", gene_id_stem="ENSG00000133794", kind="previous")
    >>> match.as_json()
    {'symbol': 'ARNTL', 'gene_id_stem': 'ENSG00000133794', 'kind': 'previous'}
    """

    symbol: str
    gene_id_stem: str
    kind: str

    def as_json(self) -> dict[str, Any]:
        """Return this match as ``--json`` serializes it, in field order.

        Returns
        -------
        dict
            ``symbol``, ``gene_id_stem`` and ``kind``.
        """
        return asdict(self)


@dataclass(frozen=True)
class ResolvedSymbols:
    """The genes one **Xref set** says each gene symbol names, and how each one matched.

    :meth:`~genome.xref.xref.XrefSet.match_symbols`'s answer — the hop *toward* the hub
    from the one **Namespace** that is not answered like an identifier.
    :class:`ResolvedStems`'s shape in every respect a caller relies on — ask order, no
    empty resolved value, what named nothing riding back — with one difference: a value is
    a tuple of :class:`SymbolMatch` rather than of stems, because **ambiguity is the return
    type here and not an edge case**. A symbol naming several genes answers with all of
    them and nothing picks one, and each says whether it matched an approved, a previous or
    an alias spelling so the caller can judge the ambiguity themselves.

    **What the set could not have matched is on the answer too.** :attr:`kinds` says which
    kinds of spelling this **Xref source** publishes and :attr:`limits` says why the others
    are missing, so *this gene is not in the release* and *this source cannot match the way
    you spelled it* are distinguishable rather than both being silence.

    Attributes
    ----------
    species : str
        The species this set is for, as the curated metadata table spells it.
    source : str
        The **Xref source** whose assertions these are.
    release : str
        The pinned **Release** of that source.
    case_insensitive : bool
        Whether case was ignored. ``False`` is the default: the species is fixed by the
        set, so a mouse-cased spelling asked of a human set is the wrong authority's and
        matches nothing rather than half-working.
    kinds : tuple of str
        The kinds of **Symbol match** this set could make, in
        :data:`~genome.xref.symbols.SYMBOL_KINDS` order.
    limits : str or None
        Why the kinds not in :attr:`kinds` are missing, or ``None`` when all three are
        there.
    resolved : mapping of str to tuple of SymbolMatch
        Every symbol that matched at least one gene, in the order they were asked about,
        to every match it made — approved first, then previous, then alias, and by stem
        within a kind. No value is ever an empty tuple.
    unresolved : tuple of str
        The symbols this release matched nothing for, in the order they were asked about.

    Examples
    --------
    >>> answer = ResolvedSymbols(
    ...     species="Homo sapiens",
    ...     source="hgnc",
    ...     release="2026-07-07",
    ...     case_insensitive=False,
    ...     kinds=("approved", "previous", "alias"),
    ...     limits=None,
    ...     resolved={
    ...         "ADCY3": (
    ...             SymbolMatch("ADCY3", "ENSG00000138031", "approved"),
    ...             SymbolMatch("ADCY3", "ENSG00000155897", "previous"),
    ...         )
    ...     },
    ...     unresolved=("Brca1",),
    ... )
    >>> answer.gene_id_stems
    ['ENSG00000138031', 'ENSG00000155897']
    >>> answer.as_json()["resolved"]["ADCY3"][1]["kind"]
    'previous'
    """

    species: str
    source: str
    release: str
    case_insensitive: bool
    kinds: tuple[str, ...]
    limits: str | None
    resolved: Mapping[str, tuple[SymbolMatch, ...]]
    unresolved: tuple[str, ...]

    @property
    def gene_id_stems(self) -> list[str]:
        """Every stem matched, ask order and then match order — a fresh list each call.

        **Every** stem, not one per symbol, and it may repeat: one gene answering a symbol
        on both an approved and an alias spelling contributes two matches and so two
        entries. Flattening loses the two things this answer exists to carry — which
        symbol named which gene, and which kind of spelling each match was on — so a
        reader who takes this list has thrown away the means of judging the ambiguity.
        :attr:`resolved` is what keeps both.
        """
        return [match.gene_id_stem for matches in self.resolved.values() for match in matches]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``species``, ``source``, ``release``, ``case_insensitive``, ``kinds`` and
            ``limits``, ``resolved`` as a mapping of symbol to a list of match objects,
            ``unresolved`` as a list, and the flattened ``gene_id_stems`` — written out
            beside the mapping it is read from for the reason :attr:`gene_id_stems` gives.
        """
        return {
            "species": self.species,
            "source": self.source,
            "release": self.release,
            "case_insensitive": self.case_insensitive,
            "kinds": list(self.kinds),
            "limits": self.limits,
            "resolved": {
                asked: [match.as_json() for match in matches]
                for asked, matches in self.resolved.items()
            },
            "unresolved": list(self.unresolved),
            "gene_id_stems": self.gene_id_stems,
        }


@dataclass(frozen=True)
class HomologyLink:
    """One row of a **Homology set**: two genes in two species, and what relates them.

    The publisher's assertion, carried unchanged. :attr:`homology_type` is Compara's own
    tree-derived label and is **never recomputed** — not after a filter, not after
    resolution into an **Annotation**, not after a caller slices the answer (ADR-0020) —
    so a link can read one-to-one in a view and still be labelled ``ortholog_one2many``.
    What a view removed is counted separately, as a **Dropped partner**, rather than the
    label being quietly corrected.

    The three confidence fields ride through exactly as the publisher wrote them, so a
    caller can filter on them. **Two of them are null for a whole species pair rather
    than for a row**: Compara records no ``goc_score`` and no ``wga_coverage`` on any link
    of either worm pairing, and a filter written against one would silently empty. Which
    fields a set holds nothing in is measured when it is prepared and said out loud on the
    answer — see :attr:`HomologyAnswer.null_quality_scores`.

    Attributes
    ----------
    gene_id_stem : str
        The **Gene id stem** asked about, in the set's first species.
    homolog_gene_id_stem : str
        The **Gene id stem** of the homologous gene, in the other species.
    homology_type : str
        Compara's own ``homology_type``, verbatim: ``ortholog_one2one``,
        ``ortholog_one2many``, ``ortholog_many2many`` for speciation, and its duplication
        labels for a **Paralogy link**.
    is_high_confidence : bool or None
        Compara's high-confidence flag, or ``None`` where it recorded none.
    goc_score : int or None
        Gene order conservation score, or ``None`` where Compara recorded none.
    wga_coverage : float or None
        Whole-genome-alignment coverage, or ``None`` where Compara recorded none.

    Examples
    --------
    >>> link = HomologyLink(
    ...     gene_id_stem="ENSG00000141510",
    ...     homolog_gene_id_stem="ENSMUSG00000059552",
    ...     homology_type="ortholog_one2one",
    ...     is_high_confidence=True,
    ...     goc_score=100,
    ...     wga_coverage=96.79,
    ... )
    >>> link.is_ortholog
    True
    >>> link.as_json()["homology_type"]
    'ortholog_one2one'
    """

    gene_id_stem: str
    homolog_gene_id_stem: str
    homology_type: str
    is_high_confidence: bool | None
    goc_score: int | None
    wga_coverage: float | None

    @property
    def is_ortholog(self) -> bool:
        """Whether the publisher's label is a speciation one rather than a duplication one.

        Read off :attr:`homology_type` and never off a count of rows: this is what the
        publisher said about the gene tree, which is why it survives every filter
        (ADR-0020). A link this is ``False`` for is a **Paralogy link**, kept and marked
        rather than dropped, so *not an ortholog* stays distinguishable from *absent*.
        """
        return self.homology_type.startswith(ORTHOLOG_TYPE_PREFIX)

    def as_json(self) -> dict[str, Any]:
        """Return this link as ``--json`` serializes it.

        Returns
        -------
        dict
            The six fields above under their own names, plus ``is_ortholog``. A null
            confidence field stays ``None`` rather than being filled in.
        """
        return {
            "gene_id_stem": self.gene_id_stem,
            "homolog_gene_id_stem": self.homolog_gene_id_stem,
            "homology_type": self.homology_type,
            "is_ortholog": self.is_ortholog,
            "is_high_confidence": self.is_high_confidence,
            "goc_score": self.goc_score,
            "wga_coverage": self.wga_coverage,
        }


@dataclass(frozen=True)
class HomologyAnswer:
    """The homologous genes one **Homology set** names for the stems it was asked about.

    :meth:`~genome.homology.compara.HomologySet.homologs`'s answer, in the shape
    :class:`ResolvedGeneIds` establishes: every stem that named at least one **Homology
    link** maps to all of them, in ask order, and no value is ever empty — a stem that
    named none is in :attr:`unresolved` instead. Nothing here picks a "best" homolog, and
    this package computes no ranking or quality score of its own.

    **What a filter removed is counted rather than dropped.** Asking for orthologs — the
    default — removes every **Paralogy link**, and the partner genes the answer therefore
    no longer names at all are counted in :attr:`dropped_partners`. Ask again with
    ``paralogs=True`` to see them.

    **A quality score that is null for the whole set says so here.** Compara records
    neither ``goc_score`` nor ``wga_coverage`` on any link of either worm pairing, so a
    caller filtering on one would get an empty result and no reason for it.
    :attr:`null_quality_scores` names the fields the set holds no value in anywhere,
    measured over the prepared slice rather than listed against a pair.

    Attributes
    ----------
    species : str
        The species the stems asked about belong to.
    other_species : str
        The species the homologous genes belong to.
    release : str
        The Ensembl Compara **Release** that asserted these links.
    resolved : mapping of str to tuple of HomologyLink
        Every stem that named at least one link, in the order the stems were asked about,
        to its links, ordered by the partner's stem. No value is ever an empty tuple.
    unresolved : tuple of str
        The stems this set names no homolog for, in the order they were asked about.
    dropped_partners : tuple of str
        The **Dropped partner**s: every partner **Gene id stem** this answer would have
        named had nothing been filtered out, and now names nowhere, in ascending order.
    null_quality_scores : tuple of str
        The names of the confidence fields the whole set holds no value in —
        ``("goc_score", "wga_coverage")`` for either worm pairing, empty for a pair
        Compara scored.

    Examples
    --------
    >>> link = HomologyLink(
    ...     "ENSG00000141510", "ENSMUSG00000059552", "ortholog_one2one", True, 100, 96.79
    ... )
    >>> answer = HomologyAnswer(
    ...     species="Homo sapiens",
    ...     other_species="Mus musculus",
    ...     release="116",
    ...     resolved={"ENSG00000141510": (link,)},
    ...     unresolved=("ENSG00000288541",),
    ...     dropped_partners=(),
    ...     null_quality_scores=(),
    ... )
    >>> answer.homolog_gene_id_stems
    ['ENSMUSG00000059552']
    >>> answer.as_json()["unresolved"]
    ['ENSG00000288541']
    """

    species: str
    other_species: str
    release: str
    resolved: Mapping[str, tuple[HomologyLink, ...]]
    unresolved: tuple[str, ...]
    dropped_partners: tuple[str, ...]
    null_quality_scores: tuple[str, ...]

    @property
    def links(self) -> list[HomologyLink]:
        """Every link, stem order and then partner order — a fresh list each call."""
        return [link for links in self.resolved.values() for link in links]

    @property
    def homolog_gene_id_stems(self) -> list[str]:
        """Every homologous **Gene id stem** named, stem order then partner order.

        **Flattening loses the two things the answer exists to carry.** It drops the
        **Homology type** and the confidence fields, so a one-to-one and a many-to-many
        partner become the same string; and it drops which asked stem each partner came
        from, so a partner two asked stems both name appears twice rather than once. Read
        :attr:`resolved` for either. This is here because a caller assembling the list
        themselves is a caller who might take one partner per stem and lose the rest,
        which is the ``ortholog_one2many`` case this shape exists to keep.
        """
        return [link.homolog_gene_id_stem for links in self.resolved.values() for link in links]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``species``, ``other_species``, ``release``, ``resolved`` as a mapping of stem
            to a list of :meth:`HomologyLink.as_json` links, ``unresolved``,
            ``dropped_partners`` and ``null_quality_scores`` as lists, and the flattened
            ``homolog_gene_id_stems``. The last is written out beside the mapping it is
            read from for the reason :attr:`homolog_gene_id_stems` gives.
        """
        return {
            "species": self.species,
            "other_species": self.other_species,
            "release": self.release,
            "resolved": {
                stem: [link.as_json() for link in links] for stem, links in self.resolved.items()
            },
            "unresolved": list(self.unresolved),
            "dropped_partners": list(self.dropped_partners),
            "null_quality_scores": list(self.null_quality_scores),
            "homolog_gene_id_stems": self.homolog_gene_id_stems,
        }


@dataclass(frozen=True)
class ResolvedHomologs:
    """A :class:`HomologyAnswer` put into one registered **Annotation**'s own gene ids.

    :func:`~genome.homology.annotation.resolve_homologs`'s answer, and the crossing a
    caller makes once they have homologs and want to join them to their own counts
    matrix. The hop itself is
    :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`, used unchanged.

    **The Homology type is the publisher's and stands** (ADR-0020). An annotation that
    spells one gene of an ``ortholog_one2many`` link and not the other leaves a view that
    looks one-to-one, and the label still reads ``ortholog_one2many``; what the crossing
    removed is in :attr:`dropped_partners` rather than folded into the label.

    **Both qualifications the answer carried ride through.** A **Dropped partner** counts
    partners lost to either step — a **Homology type** filter before the crossing, or an
    annotation missing the gene during it — so the count a caller reads is what the whole
    path cost rather than what its last step did. And
    :attr:`~HomologyAnswer.null_quality_scores` is a fact about the set the answer came
    from, not about the crossing, so it is repeated here for a caller who filters on
    ``goc_score`` after resolving.

    Attributes
    ----------
    species : str
        The species the stems asked about belong to.
    other_species : str
        The species the homologous genes belong to, and the one the annotation annotates.
    release : str
        The Ensembl Compara **Release** that asserted these links.
    assembly : str
        The **Assembly** whose annotation these gene ids belong to.
    annotation : str
        The **Registered name** whose own gene ids these are.
    resolved : mapping of str to tuple of HomologyLink
        Every asked stem that still names at least one link, in ask order, to the links
        whose partner this annotation carries a gene for. No value is ever empty.
    gene_ids : mapping of str to tuple of str
        Every partner **Gene id stem** that survived, to the gene ids this annotation
        spells it with, ascending. Keyed by partner and not by asked stem, because two
        asked stems may name one partner and its ids are the same ids.
    unresolved : tuple of str
        The asked stems left naming nothing here: first those the crossing emptied — every
        partner missing from this annotation — in ask order, then those the set already
        named no homolog for, in ask order. Two groups rather than one interleaved list,
        because *this annotation is missing every partner* and *this release knows no
        homolog* are different facts about a gene.
    dropped_partners : tuple of str
        The **Dropped partner**s: every partner **Gene id stem** this answer no longer
        names, ascending — those a **Homology type** filter removed before the crossing
        and those this annotation carries no gene for, in one count, since the definition
        covers both and a caller wants what the whole path cost.
    null_quality_scores : tuple of str
        The names of the confidence fields the **Homology set** behind this holds no value
        in, carried through from :attr:`HomologyAnswer.null_quality_scores` unchanged. The
        crossing neither adds a score nor removes one.

    Examples
    --------
    >>> link = HomologyLink(
    ...     "ENSG00000141510", "ENSMUSG00000059552", "ortholog_one2many", True, 100, 96.79
    ... )
    >>> crossed = ResolvedHomologs(
    ...     species="Homo sapiens",
    ...     other_species="Mus musculus",
    ...     release="116",
    ...     assembly="mm39",
    ...     annotation="gencode_vM39",
    ...     resolved={"ENSG00000141510": (link,)},
    ...     gene_ids={"ENSMUSG00000059552": ("ENSMUSG00000059552.5",)},
    ...     unresolved=(),
    ...     dropped_partners=("ENSMUSG00000000001",),
    ...     null_quality_scores=(),
    ... )
    >>> crossed.homolog_gene_ids
    ['ENSMUSG00000059552.5']
    >>> crossed.resolved["ENSG00000141510"][0].homology_type
    'ortholog_one2many'
    """

    species: str
    other_species: str
    release: str
    assembly: str
    annotation: str
    resolved: Mapping[str, tuple[HomologyLink, ...]]
    gene_ids: Mapping[str, tuple[str, ...]]
    unresolved: tuple[str, ...]
    dropped_partners: tuple[str, ...]
    null_quality_scores: tuple[str, ...]

    @property
    def homolog_gene_ids(self) -> list[str]:
        """Every homologous gene id named, partner order then id order — a fresh list.

        **Flattening loses what the mapping carries**: which asked stem reached the gene,
        and under what **Homology type**. It keeps every id rather than one per partner,
        since one stem may be spelled by two gene ids — the pseudoautosomal case
        :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids` answers with both of.
        """
        return [gene_id for ids in self.gene_ids.values() for gene_id in ids]

    def as_json(self) -> dict[str, Any]:
        """Return this crossing as ``--json`` serializes it.

        Returns
        -------
        dict
            The species pair, the ``release``, the ``assembly`` and ``annotation``,
            ``resolved`` as a mapping of stem to a list of :meth:`HomologyLink.as_json`
            links, ``gene_ids`` as a plain mapping, ``unresolved``, ``dropped_partners``
            and ``null_quality_scores`` as lists, and the flattened ``homolog_gene_ids``.
        """
        return {
            "species": self.species,
            "other_species": self.other_species,
            "release": self.release,
            "assembly": self.assembly,
            "annotation": self.annotation,
            "resolved": {
                stem: [link.as_json() for link in links] for stem, links in self.resolved.items()
            },
            "gene_ids": {stem: list(ids) for stem, ids in self.gene_ids.items()},
            "unresolved": list(self.unresolved),
            "dropped_partners": list(self.dropped_partners),
            "null_quality_scores": list(self.null_quality_scores),
            "homolog_gene_ids": self.homolog_gene_ids,
        }


@dataclass(frozen=True)
class AnnotationStatusRow:
    """One annotation, in whichever of its states it is: offered, registered, broken.

    One shape for all of them, so a reader never has to ask which fields a row has — a
    name the table does not list carries the table's columns as ``None``, and one nothing
    is wrong with carries the broken columns as ``None``. :attr:`registered` and
    :attr:`broken` are never both true: a registration nothing vouches for is not one, and
    :attr:`state` is that invariant said in one word.

    Attributes
    ----------
    name : str
        The **Registered name** this row is about.
    offered : bool
        Whether the annotation table lists it for this assembly.
    registered : bool
        Whether a record here vouches for it.
    broken : bool
        Whether its directory is here and cannot be trusted.
    default : bool
        The table's own default flag, ``False`` for a name no row lists.
    provider : str or None
        Who publishes it, from the table's row; ``None`` for an unlisted one.
    version : str or None
        The provider's release identifier; ``None`` for an unlisted one.
    url : str or None
        Where the table says its GTF is fetched from; ``None`` for an unlisted one.
    sha256 : str or None
        The digest the table pins; ``None`` when it pins none, and for an unlisted one.
    path : str or None
        The registered GTF's path, or ``None`` when it is not registered here.
    problem : str or None
        What is wrong, when :attr:`broken`; ``None`` otherwise.
    repair : str or None
        The command that registers it again from scratch, when :attr:`broken`.

    Examples
    --------
    >>> row = AnnotationStatusRow(
    ...     name="gencode_v50",
    ...     offered=True,
    ...     registered=False,
    ...     broken=False,
    ...     default=True,
    ...     provider="GENCODE",
    ...     version="v50",
    ...     url="https://example.org/gencode_v50.gtf.gz",
    ...     sha256=None,
    ...     path=None,
    ...     problem=None,
    ...     repair=None,
    ... )
    >>> row.as_json()["offered"]
    True
    >>> row.state
    'offered, not registered'
    """

    name: str
    offered: bool
    registered: bool
    broken: bool
    default: bool
    provider: str | None
    version: str | None
    url: str | None
    sha256: str | None
    path: str | None
    problem: str | None
    repair: str | None

    @property
    def state(self) -> str:
        """Which of its states this row is in, in the words a surface prints.

        ``broken`` first, because the three fields it is read from are not independent:
        a broken annotation is not registered — no record vouches for it — so answering
        with the absence of one would be true and useless, and it is the state that needs
        acting on. A row nothing offers is one this disk holds and the table does not, so
        the only thing left to say about it is that.

        Returns
        -------
        str
            One of ``broken``, ``registered, not offered``, ``registered``, or
            ``offered, not registered``.
        """
        if self.broken:
            return "broken"
        if not self.offered:
            return "registered, not offered"
        return "registered" if self.registered else "offered, not registered"

    def as_json(self) -> dict[str, Any]:
        """Return this row as ``--json`` serializes it: every attribute above, in order.

        :attr:`state` is not among them: it is read from :attr:`broken`, :attr:`offered`
        and :attr:`registered`, which are all here, so writing it out too would be a
        second spelling of the same rule for a reader to disagree with.

        Returns
        -------
        dict
            The row's fields, under their own names.
        """
        return asdict(self)


@dataclass(frozen=True)
class AnnotationStatus:
    """What one assembly's table offers, set against what is registered on this machine.

    :meth:`~genome.io.gtf.AnnotationRegistry.status`'s answer, and what ``genome
    annotations`` prints. Two questions joined for one reader, with a third riding along
    because this is where anyone would look for it: a directory that cannot be trusted is
    ``broken`` rather than registered, and reporting one is the point — nothing here
    raises.

    Attributes
    ----------
    assembly : str
        The **Assembly** reported on.
    directory : pathlib.Path
        Its **Assembly dir**, whether or not anything is there.
    default_annotation : str or None
        The **Default annotation**'s name, or ``None`` when nothing decides one. It may
        name one nobody has registered here, which is a fresh machine's ordinary state.
    annotations : tuple of AnnotationStatusRow
        One row per name: the offered ones in table order, then anything on this disk
        that no row lists.

    Examples
    --------
    >>> from pathlib import Path
    >>> status = AnnotationStatus(
    ...     assembly="hg38",
    ...     directory=Path("/data/genome/hg38"),
    ...     default_annotation=None,
    ...     annotations=(),
    ... )
    >>> status.default_row is None
    True
    >>> status.default_summary
    'default: (none)'
    >>> status.as_json()["directory"]
    '/data/genome/hg38'
    """

    assembly: str
    directory: Path
    default_annotation: str | None
    annotations: tuple[AnnotationStatusRow, ...]

    @property
    def default_row(self) -> AnnotationStatusRow | None:
        """The **Default annotation**'s own row, or ``None`` when no row is about it.

        ``None`` covers both of the ways that happens, and a caller wanting to tell them
        apart reads :attr:`default_annotation` beside this: nothing decided a default, or
        one is decided and the table lists it under a name this disk knows nothing about.
        """
        return next((row for row in self.annotations if row.name == self.default_annotation), None)

    @property
    def default_summary(self) -> str:
        """The closing sentence naming the **Default annotation**, and how to get it.

        Four answers, and three of them tell the reader what to do next: nothing decided a
        default; one is decided and registered, which needs no advice; one is decided and
        broken, which is repaired by the command its own row already carries; or one is
        decided and absent, which is the ordinary state of a fresh machine and is
        registered by :func:`annotation_register_command`. Both commands come off an
        interface rather than being assembled here, so the two halves of this sentence
        cannot drift apart.

        Returns
        -------
        str
            One line, beginning ``default: `` — the whole of what ``genome annotations``
            prints last.
        """
        default = self.default_annotation
        if default is None:
            return "default: (none)"
        row = self.default_row
        if row is not None and row.broken:
            return f"default: {default} — broken here; repair it with `{row.repair}`"
        if row is not None and row.registered:
            return f"default: {default}"
        return (
            f"default: {default} — not registered here; register it with "
            f"`{annotation_register_command(self.assembly, default)}`"
        )

    def as_json(self) -> dict[str, Any]:
        """Return this report as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, the ``directory`` as text, the ``default_annotation`` name, and
            ``annotations`` as a list of :meth:`AnnotationStatusRow.as_json` rows.
        """
        return {
            "assembly": self.assembly,
            "directory": str(self.directory),
            "default_annotation": self.default_annotation,
            "annotations": [row.as_json() for row in self.annotations],
        }
