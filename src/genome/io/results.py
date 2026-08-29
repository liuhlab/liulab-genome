"""What the API answers with: the values it returns and ``--json`` writes.

Registering something in this package — an **Assembly**, an annotation — or asking it a
question about one ends by returning a value that says what it found, and every surface
reads that same value: a script keeps it, the CLI prints it, ``--json`` serializes it.
Those values live here, together, because they change for one reason: *what an answer must
be able to say*. The modules that do the work import them; nothing here does any of it,
reads a directory, or reaches the network.

Two of them carry a **Completion marker** whole rather than copying it out field by
field — :class:`RegisteredAssembly` and :class:`RegisteredAnnotation` — so every later
question is answered off the record in hand instead of by opening the directory again. The
rest are reports assembled from more than one source: :class:`VerifiedAssembly` sets a
digest against what pinned it, :class:`AnnotationStatus` sets what the annotation table
offers against what this machine holds, one :class:`AnnotationStatusRow` per name,
:class:`GeneList` sets one **Gene category**'s genes against the **Curated gene list**s
that contributed them, one :class:`GeneListSource` apiece, :class:`ResolvedGeneIds`
sets the **Gene id stem**s a caller asked about against the gene ids one annotation
actually carries — including the stems it carries none for — and :class:`TFGeneList` sets
one published census against one annotation, one :class:`TFGene` per gene it named here,
carrying the census's own provenance so the verdict never travels without it.
:class:`TFCofactorList` and :class:`TFCofactor` are that pair's counterpart for the genes
a publisher lists as a **Transcription cofactor**, in the same shape and with the same
rules, because a caller who has read one answer has read both.

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
from genome.tf.cofactor import CofactorProvenance
from genome.tf.gene import CensusProvenance

#: What answered *which digest should this FASTA have?* — the assembly's curated
#: metadata row, or the completion record its own registration wrote. Reported by
#: :func:`~genome.io.download.verify_assembly` so that being held to a pin and being held
#: only to what this machine last produced are never read as the same result. Public
#: because the CLI keys its two sentences on them: a surface that spelled the strings
#: again would print the raw status the day one of these was renamed, rather than failing.
EXPECTED_FROM_TABLE = "table"
EXPECTED_FROM_RECORD = "record"

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
class TFGene:
    """One gene a census assessed, named in one **Annotation**'s own gene ids.

    An entry of a :class:`TFGeneList`. The census's four uniform columns are fields of
    their own — the **Gene id stem** it is keyed by, the symbol, the TF flag and the **DBD
    family** — and everything the publisher recorded beyond them stays under the
    publisher's own name in :attr:`judgements`, because beyond those four no two censuses
    carry the same columns and nothing here compares one publisher's with another's
    (ADR-0014).

    ``gene_ids`` is a tuple because one stem may name more than one gene id in one
    annotation and this never picks one — see
    :meth:`~genome.io.gtf.AnnotationRegistry.resolve_gene_ids`. It is never empty: a stem
    the annotation carries no gene for is in :attr:`TFGeneList.unresolved` instead of here.

    Attributes
    ----------
    gene_id_stem : str
        The **Gene id stem** the census is keyed by.
    gene_ids : tuple of str
        Every gene id this annotation spells that stem with, ascending.
    symbol : str or None
        The gene symbol the census records, or ``None`` where it records none.
    is_tf : bool
        The census's own TF flag. ``True`` for every gene of a **TF gene list** left at
        its default, and the field that tells a rejected gene from an accepted one when a
        caller widened to carry both.
    dbd_family : str or None
        The **DBD family** this census classifies the gene under, in the publisher's own
        vocabulary — group by it within a species and never across two.
    judgements : mapping of str to (str or None)
        Every other column the census records for this gene, under the census's own
        snake_case spelling of its published name: the **TF assessment** a caller tightens
        or loosens on, the binding mode, the motif status, the KRAB flag and the
        third-party votes, for a census that records them. ``None`` is a cell the
        publisher left blank.

    Examples
    --------
    >>> gene = TFGene(
    ...     gene_id_stem="ENSG00000214717",
    ...     gene_ids=("ENSG00000214717.13", "ENSG00000214717.13_PAR_Y"),
    ...     symbol="ZBED1",
    ...     is_tf=True,
    ...     dbd_family="BED ZF",
    ...     judgements={"tf_assessment": "Known motif"},
    ... )
    >>> gene.judgements["tf_assessment"]
    'Known motif'
    >>> gene.as_json()["gene_ids"]
    ['ENSG00000214717.13', 'ENSG00000214717.13_PAR_Y']
    """

    gene_id_stem: str
    gene_ids: tuple[str, ...]
    symbol: str | None
    is_tf: bool
    dbd_family: str | None
    judgements: Mapping[str, str | None]

    def as_json(self) -> dict[str, Any]:
        """Return this gene as ``--json`` serializes it.

        Returns
        -------
        dict
            The fields above under their own names, with ``gene_ids`` as a list and
            ``judgements`` as a plain mapping under the census's own column names.
        """
        return {
            "gene_id_stem": self.gene_id_stem,
            "gene_ids": list(self.gene_ids),
            "symbol": self.symbol,
            "is_tf": self.is_tf,
            "dbd_family": self.dbd_family,
            "judgements": dict(self.judgements),
        }


@dataclass(frozen=True)
class TFGeneList:
    """One **Assembly**'s **TF gene**s, in its registered annotation's own gene ids.

    :meth:`~genome.io.gtf.AnnotationRegistry.tf_gene_list`'s answer, and what a
    ``--json`` surface over it serializes. The census's **Gene id stem**s
    resolved against one annotation, so the ids join to a counts matrix with nothing left
    to normalise, and assessed-positive by default: the common case is not 2,765 rows to
    filter down to 1,639.

    **Nothing here decides what a transcription factor is.** Every verdict is the census's
    and travels with :attr:`provenance`, which names the publisher to cite. Two censuses
    that classify one factor differently are two answers rather than a contradiction, and
    this says which one is speaking.

    **What the census holds and this annotation does not is visible.** A stem no gene id
    here is of comes back in :attr:`unresolved` rather than being dropped, so a caller can
    count what the crossing cost instead of wondering.

    There is no empty one for the reasons an absent census would give: an assembly whose
    species has no census, and one nothing names a species for, each raise a
    :class:`LookupError` of their own.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    species : str
        The species the assembly's own metadata row names, which is what selected the
        census. Never passed in by a caller, so asking for one species' transcription
        factors while holding another species' assembly is not expressible (ADR-0003).
    provenance : genome.tf.gene.census.CensusProvenance
        Where the census came from: publisher, version, PubMed id, source URL and digest.
        :meth:`~genome.tf.gene.census.CensusProvenance.attribution` renders the line to
        print beside anything it answered.
    genes : tuple of TFGene
        One entry per **Gene id stem** that named at least one gene id here, in the
        census's own row order.
    unresolved : tuple of str
        The stems this annotation carries no gene for, in census row order.

    Examples
    --------
    >>> from genome.tf.gene import tf_gene_table
    >>> answer = TFGeneList(
    ...     assembly="hg38",
    ...     annotation="gencode_v50",
    ...     species="Homo sapiens",
    ...     provenance=tf_gene_table("Homo sapiens").provenance,
    ...     genes=(
    ...         TFGene("ENSG00000137203", ("ENSG00000137203.12",), "TFAP2A", True, "AP-2", {}),
    ...     ),
    ...     unresolved=("ENSG00000214717",),
    ... )
    >>> answer.gene_ids
    ['ENSG00000137203.12']
    >>> answer.provenance.publisher
    'Lambert et al. 2018'
    >>> answer.as_json()["unresolved"]
    ['ENSG00000214717']
    """

    assembly: str
    annotation: str
    species: str
    provenance: CensusProvenance
    genes: tuple[TFGene, ...]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id, gene order then id order — a fresh list each call.

        **Every** id, not one per gene, for the reason
        :attr:`ResolvedGeneIds.gene_ids` gives: flattening is where a reader would take
        the first id of a stem that names two and lose the other. :attr:`genes` is what
        says which gene an id came from, and what the census said about it.
        """
        return [gene_id for gene in self.genes for gene_id in gene.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``species``, the census's ``provenance`` under
            its own field names, ``genes`` as a list of :meth:`TFGene.as_json` entries,
            the flattened ``gene_ids``, and ``unresolved`` as a list. The ids are written
            out beside the genes they are read from for the reason :attr:`gene_ids` gives.
        """
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "species": self.species,
            "provenance": asdict(self.provenance),
            "genes": [gene.as_json() for gene in self.genes],
            "gene_ids": self.gene_ids,
            "unresolved": list(self.unresolved),
        }


@dataclass(frozen=True)
class TFCofactor:
    """One gene a publisher lists as a cofactor, named in one **Annotation**'s gene ids.

    An entry of a :class:`TFCofactorList`, and the counterpart of :class:`TFGene`. The
    table's four uniform columns are fields of their own — the **Gene id stem** it is
    keyed by, the symbol, the cofactor flag and which publisher listed the gene — and
    everything each publisher classified it with stays under that publisher's own
    namespaced name in :attr:`classifications`, because no two publishers carry the same
    columns and nothing here compares one's vocabulary with another's (ADR-0014).

    ``gene_ids`` is a tuple for the reason :class:`TFGene`'s is: one stem may name more
    than one gene id in one annotation and this never picks one. It is never empty — a
    stem the annotation carries no gene for is in :attr:`TFCofactorList.unresolved`
    instead of here.

    Attributes
    ----------
    gene_id_stem : str
        The **Gene id stem** the **Cofactor table** is keyed by.
    gene_ids : tuple of str
        Every gene id this annotation spells that stem with, ascending.
    symbol : str or None
        The gene symbol the table records, or ``None`` where it records none.
    is_cofactor : bool
        The table's own cofactor flag. ``True`` for every entry of a **TF cofactor
        list**, since no publisher shipping here releases a rejected set — a source that
        did would ship rejected rows, and they would be excluded rather than arriving
        here saying ``False``.
    source : str
        Which publisher listed the gene, in the table's own closed vocabulary —
        ``animaltfdb``, ``epifactors``, or ``both`` for a gene two of them listed. ``both``
        is agreement on **membership only** and never on how either classified it.
    classifications : mapping of str to (str or None)
        Every other column the table records for this gene, under the publisher's own
        namespaced snake_case name: ``animaltfdb_family`` and ``animaltfdb_category``
        for a gene AnimalTFDB listed, and the EpiFactors function, target, modification
        and complex name for one EpiFactors did. ``None`` is a cell that publisher left
        blank — which, for a gene the other publisher listed, is that publisher saying
        nothing rather than a value being lost.

    Examples
    --------
    >>> cofactor = TFCofactor(
    ...     gene_id_stem="ENSMUSG00000000085",
    ...     gene_ids=("ENSMUSG00000000085.16",),
    ...     symbol="Scmh1",
    ...     is_cofactor=True,
    ...     source="animaltfdb",
    ...     classifications={"animaltfdb_category": "Other Cofactors"},
    ... )
    >>> cofactor.classifications["animaltfdb_category"]
    'Other Cofactors'
    >>> cofactor.as_json()["gene_ids"]
    ['ENSMUSG00000000085.16']
    """

    gene_id_stem: str
    gene_ids: tuple[str, ...]
    symbol: str | None
    is_cofactor: bool
    source: str
    classifications: Mapping[str, str | None]

    def as_json(self) -> dict[str, Any]:
        """Return this cofactor as ``--json`` serializes it.

        Returns
        -------
        dict
            The fields above under their own names, with ``gene_ids`` as a list and
            ``classifications`` as a plain mapping under the publishers' own column names.
        """
        return {
            "gene_id_stem": self.gene_id_stem,
            "gene_ids": list(self.gene_ids),
            "symbol": self.symbol,
            "is_cofactor": self.is_cofactor,
            "source": self.source,
            "classifications": dict(self.classifications),
        }


@dataclass(frozen=True)
class TFCofactorList:
    """One **Assembly**'s **Transcription cofactor**s, in its annotation's own gene ids.

    :meth:`~genome.io.gtf.AnnotationRegistry.tf_cofactor_list`'s answer, and the
    counterpart of :class:`TFGeneList` in the same shape: the **Cofactor table**'s **Gene
    id stem**s resolved against one registered annotation, so the ids join to a counts
    matrix with nothing left to normalise.

    **Membership is this package's and classification is each publisher's.** A table built
    from two publishers is a union nobody else published, which is why :attr:`provenance`
    carries a record per publisher rather than one, and why a ``source`` of ``both`` on an
    entry says the two agreed the gene is a cofactor and nothing about how either
    classified it (ADR-0016).

    **What the table holds and this annotation does not is visible.** A stem no gene id
    here is of comes back in :attr:`unresolved` rather than being dropped.

    There is no empty one, for the reasons an absent table would give: an assembly whose
    species has no cofactor table, and one nothing names a species for, each raise a
    :class:`LookupError` of their own.

    Attributes
    ----------
    assembly : str
        The **Assembly** asked about.
    annotation : str
        The **Registered name** whose own gene ids these are.
    species : str
        The species the assembly's own metadata row names, which is what selected the
        table. Never passed in by a caller, so asking for one species' cofactors while
        holding another species' assembly is not expressible (ADR-0003).
    provenance : genome.tf.cofactor.table.CofactorProvenance
        Where the table came from: one record per publisher that contributed to it, plus
        the digest of the shipped bytes.
        :meth:`~genome.tf.cofactor.table.CofactorProvenance.attribution` renders the line
        to print beside anything it answered.
    cofactors : tuple of TFCofactor
        One entry per **Gene id stem** that named at least one gene id here, in the
        table's own row order.
    unresolved : tuple of str
        The stems this annotation carries no gene for, in table row order.

    Examples
    --------
    >>> from genome.tf.cofactor import cofactor_table
    >>> answer = TFCofactorList(
    ...     assembly="mm39",
    ...     annotation="gencode_vM39",
    ...     species="Mus musculus",
    ...     provenance=cofactor_table("Mus musculus").provenance,
    ...     cofactors=(
    ...         TFCofactor(
    ...             "ENSMUSG00000000085",
    ...             ("ENSMUSG00000000085.16",),
    ...             "Scmh1",
    ...             True,
    ...             "animaltfdb",
    ...             {},
    ...         ),
    ...     ),
    ...     unresolved=("ENSMUSG00000000275",),
    ... )
    >>> answer.gene_ids
    ['ENSMUSG00000000085.16']
    >>> answer.provenance.sources[0].publisher
    'AnimalTFDB'
    >>> answer.as_json()["unresolved"]
    ['ENSMUSG00000000275']
    """

    assembly: str
    annotation: str
    species: str
    provenance: CofactorProvenance
    cofactors: tuple[TFCofactor, ...]
    unresolved: tuple[str, ...]

    @property
    def gene_ids(self) -> list[str]:
        """Every gene id, cofactor order then id order — a fresh list each call.

        **Every** id, not one per gene, for the reason :attr:`TFGeneList.gene_ids` gives:
        flattening is where a reader would take the first id of a stem that names two and
        lose the other. :attr:`cofactors` is what says which gene an id came from, and
        what the publisher said about it.
        """
        return [gene_id for cofactor in self.cofactors for gene_id in cofactor.gene_ids]

    def as_json(self) -> dict[str, Any]:
        """Return this answer as ``--json`` serializes it.

        Returns
        -------
        dict
            ``assembly``, ``annotation``, ``species``, the table's ``provenance`` under
            its own field names with one entry per publisher under ``sources``,
            ``cofactors`` as a list of :meth:`TFCofactor.as_json` entries, the flattened
            ``gene_ids``, and ``unresolved`` as a list — the keys :class:`TFGeneList` uses,
            with the entries named for what they are.
        """
        provenance = asdict(self.provenance)
        # ``asdict`` leaves a tuple field a tuple and JSON has no tuple, so the ragged
        # per-publisher records — the one such field here — are written out as the list
        # they serialize to. A payload that did not survive its own round trip would be
        # one whose shape depended on whether anybody had serialized it yet.
        provenance["sources"] = list(provenance["sources"])
        return {
            "assembly": self.assembly,
            "annotation": self.annotation,
            "species": self.species,
            "provenance": provenance,
            "cofactors": [cofactor.as_json() for cofactor in self.cofactors],
            "gene_ids": self.gene_ids,
            "unresolved": list(self.unresolved),
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
