"""Filtering an **Xref source**'s rows by how it graded each assertion.

A pure module — it opens nothing, downloads nothing and knows no **Xref source** — and a
peer of :mod:`genome.xref.ids`: that one says how an identifier is spelled, this one says
what a caller may ask about the *evidence* behind a pair, and what happens when the answer
is nothing.

**Some publishers grade their cross-references and some do not**, so an evidence filter is
a capability of the source rather than of every set. A source whose file carries no such
column meets a filter with :class:`EvidenceNotRecordedError` rather than ignoring it, since
a filter silently dropped is a quality claim a caller believes and never made.

**A filter that keeps nothing raises.** This is not defensive tidying; it is the shape of a
real trap. Every human ``EntrezGene`` row in Ensembl release 116's TSV carries
``info_type=DEPENDENT`` and not one carries ``DIRECT`` — 552,633 rows, zero direct, and
mouse is the same at 358,853 — so the intuitive quality filter yields an empty set rather
than a smaller one. An empty **Xref set** answers every query with nothing, which is
indistinguishable from a gene list that matched none of it, so
:func:`check_evidence_filter` refuses to build one and names what the release actually
carries instead.

Examples
--------
>>> from genome.xref.evidence import check_evidence_filter, normalise_evidence
>>> normalise_evidence("dependent")
('DEPENDENT',)
>>> normalise_evidence(["DIRECT", "dependent", "DIRECT"])
('DEPENDENT', 'DIRECT')
>>> normalise_evidence(None)
()
>>> try:
...     check_evidence_filter(
...         wanted=("DIRECT",), seen={"DEPENDENT": 552633}, origin="release-116.tsv"
...     )
... except EmptyEvidenceFilterError as error:
...     print("DEPENDENT (552,633)" in str(error))
True
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping


class EvidenceNotRecordedError(LookupError):
    """The **Xref source**'s file grades nothing, so there is no evidence to filter on.

    A :class:`LookupError` for the same reason
    :class:`~genome.xref.xref.NamespaceNotCarriedError` is one: the caller named something
    this set does not carry, and the message names a source that does. Raised rather than
    ignored, because a filter that quietly does nothing leaves the caller believing their
    answer was graded when it was not.

    Examples
    --------
    >>> from genome.xref import ALLIANCE, XrefSet
    >>> XrefSet("Homo sapiens", ALLIANCE, evidence="DIRECT")        # doctest: +SKIP
    Traceback (most recent call last):
    EvidenceNotRecordedError: ...
    """


class EmptyEvidenceFilterError(LookupError):
    """The filter kept no rows, so the set would answer every query with nothing.

    A :class:`LookupError` and not a :class:`ValueError`: the evidence types named are ones
    this release resolves nothing under, which is the same kind of miss as an unknown
    **Namespace**, and the message answers it the same way — by naming what *is* there.

    Examples
    --------
    >>> from genome.xref import ENSEMBL_TSV, XrefSet
    >>> XrefSet("Homo sapiens", ENSEMBL_TSV, "116", evidence="DIRECT")  # doctest: +SKIP
    Traceback (most recent call last):
    EmptyEvidenceFilterError: ...
    """


def normalise_evidence(evidence: str | Iterable[str] | None) -> tuple[str, ...]:
    """Return the evidence types asked for, in one canonical spelling.

    Upper-cased because that is how every publisher surveyed writes them, stripped, emptied
    of blanks, deduplicated and **sorted** — the types are a set and their order carries no
    meaning, so sorting them is what lets the filter name a directory that two callers who
    asked in different orders both land on.

    Idempotent, and total: anything that is not asked for comes back as the empty tuple,
    which is what *no filter at all* is spelled as everywhere below.

    Parameters
    ----------
    evidence : str or iterable of str or None
        One evidence type, several, or ``None`` for no filter. A bare string is one type
        and never a sequence of one-character types.

    Returns
    -------
    tuple of str
        The types, upper-cased, unique and ascending. Empty when nothing was asked for.

    Examples
    --------
    >>> normalise_evidence("DIRECT")
    ('DIRECT',)
    >>> normalise_evidence(" dependent ")
    ('DEPENDENT',)
    >>> normalise_evidence(("DIRECT", "DEPENDENT"))
    ('DEPENDENT', 'DIRECT')
    >>> normalise_evidence(normalise_evidence(["direct", ""]))
    ('DIRECT',)
    """
    if evidence is None:
        return ()
    asked = [evidence] if isinstance(evidence, str) else list(evidence)
    return tuple(sorted({kind.strip().upper() for kind in asked if kind.strip()}))


def check_evidence_filter(*, wanted: tuple[str, ...], seen: Mapping[str, int], origin: str) -> None:
    """Refuse a filter that kept none of the rows, naming what the release carries.

    Called by a reader once it has counted, over the rows it would otherwise keep, how many
    carry each evidence type. It is a no-op when nothing was filtered, when something
    survived, or when the file held no rows at all — that last case is a file that is not
    the one the row pins, and it is answered further out by an error that names the URL.

    Parameters
    ----------
    wanted : tuple of str
        The filter, as :func:`normalise_evidence` spells it.
    seen : mapping of str to int
        Every evidence type the reader met, and how many rows carried it — counted over the
        rows the filter was applied to and never over rows the reader drops for other
        reasons, so the types offered back are ones a caller could actually have asked for.
    origin : str
        The file the counts came from; named in the message.

    Raises
    ------
    EmptyEvidenceFilterError
        If ``wanted`` is non-empty, ``seen`` is non-empty, and no type in ``seen`` was
        asked for.

    Examples
    --------
    >>> check_evidence_filter(wanted=(), seen={"DEPENDENT": 3}, origin="x") is None
    True
    >>> check_evidence_filter(wanted=("DEPENDENT",), seen={"DEPENDENT": 3}, origin="x") is None
    True
    >>> try:
    ...     check_evidence_filter(wanted=("DIRECT",), seen={"DEPENDENT": 3}, origin="x")
    ... except EmptyEvidenceFilterError as error:
    ...     print("without an evidence filter" in str(error))
    True
    """
    if not wanted or not seen or any(kind in wanted for kind in seen):
        return
    total = sum(seen.values())
    carried = ", ".join(f"{kind} ({count:,})" for kind, count in sorted(seen.items()))
    offered = " or ".join(sorted(seen))
    raise EmptyEvidenceFilterError(
        f"the evidence filter {'/'.join(wanted)} kept none of the {total:,} graded rows in "
        f"{origin}: they carry {carried} and nothing else. Filtering on evidence type here "
        f"empties the set rather than narrowing it, and a set that answers every query with "
        f"nothing is indistinguishable from a gene list that matched none of it — which is "
        f"why this refuses to build one. Ask for {offered}, or construct the set without an "
        f"evidence filter to read every assertion the release makes."
    )
