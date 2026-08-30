"""Where an assembly's FASTA comes from — the answer one name resolves to.

The **Source** as a value. Naming an assembly is the whole interface, so one name has to
answer *where do these bytes come from*, and there are three answers: a URL — the one the
assembly's **Assembly metadata** row pins, else the UCSC golden path derived from the name
— the local path or URL a caller handed over, or a recipe, *these components* (ADR-0008).
:func:`resolve_source` gives that answer once and the registration dispatches on which kind
came back, rather than each step asking the question again in its own words.

Deliberately the bottom of the stack: nothing here fetches, writes, or opens a
:class:`~genome.assembly.genome.Genome`. It reads the **Completion marker** already on disk, the
curated table, and the name — which is what lets :mod:`genome.assembly.download` import this at
the top of the file.

What a *finished* chimera wrote down about itself is :mod:`genome.assembly.components`, which
this module reads at one point only: the record here is believed before the name is
consulted, and telling a chimera's record from any other is what
:meth:`~genome.assembly.components.ChimeraDetails.from_record` is for. Nothing comes back the
other way.

Examples
--------
>>> from genome.assembly.registration import AssemblyDir
>>> nowhere = AssemblyDir.locate("hg38", "/tmp/definitely-not-a-build")
>>> resolve_source(nowhere, metadata=None, golden_path_url="https://example.org/hg38.fa.gz")
FetchedSource(url='https://example.org/hg38.fa.gz', derived=True)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from genome.assembly.chimera import ChimeraNamingError, derive_name, split_name
from genome.assembly.components import ChimeraDetails
from genome.assembly.metadata import AssemblyMetadata, lookup_assembly
from genome.assembly.registration import AssemblyDir, assembly_data_dir
from genome.store.completion import read_record


@dataclass(frozen=True)
class FetchedSource:
    """The assembly's FASTA is downloaded from a URL.

    One kind and not two, because a pinned URL and a derived one differ in exactly one
    consequence and it is carried here: validation is a property of the source (ADR-0003),
    so a URL the assembly's row named is fetched as it stands, while one derived from the
    assembly name has that name checked at UCSC first — a typo then fails on a ``HEAD``
    request rather than three minutes into a download.

    Attributes
    ----------
    url : str
        The URL the FASTA is fetched from.
    derived : bool
        Whether :attr:`url` was derived from the assembly name because no row pinned one.

    Examples
    --------
    >>> FetchedSource("https://example.org/hg38.fa.gz", derived=False).derived
    False
    """

    url: str
    derived: bool


@dataclass(frozen=True)
class SeededSource:
    """The caller handed the bytes over — a local path or a URL of their own.

    Never resolved, only declared: this is what ``Genome(path_or_url=...)`` and ``genome
    register <name> --source`` say, so nothing is inferred and UCSC is never consulted.
    The bytes are recorded rather than compared, and the assembly name degrades to a label
    for the directory they land in (ADR-0003).

    Attributes
    ----------
    location : str or pathlib.Path
        The local FASTA path or the http(s)/ftp/sftp URL that was named.

    Examples
    --------
    >>> SeededSource("/data/my_ref.fa").location
    '/data/my_ref.fa'
    """

    location: str | Path


@dataclass(frozen=True)
class ComponentSource:
    """A recipe rather than a location: this assembly's FASTA is *these components*.

    The third kind, and what makes a **Chimera** an assembly rather than a type of its own
    (ADR-0008). Nothing is fetched — the bytes are copied from components already prepared
    on this disk — so naming one is never a way to ask for its parts.

    Attributes
    ----------
    components : tuple of str
        The **Component** assembly names, in the sorted order the chimera's name spells
        them.

    Examples
    --------
    >>> ComponentSource(("ce11", "ecHT115")).components
    ('ce11', 'ecHT115')
    """

    components: tuple[str, ...]


#: An assembly's **Source**, in the three kinds one name can resolve to. The registration
#: dispatches on which of them came back from :func:`resolve_source`; a
#: :class:`SeededSource` reaches it from the caller instead, since a source somebody named
#: is not one anything has to work out.
Source = FetchedSource | SeededSource | ComponentSource


def fetched_source(metadata: AssemblyMetadata | None, golden_path_url: str) -> FetchedSource:
    """Return the URL this assembly's FASTA is downloaded from, and how that was decided.

    The metadata row's pinned source when it has one, and ``golden_path_url`` when it does
    not — the whole of *where does an ordinary assembly come from*, in one place, so that
    the URL a fetch uses and the question of whether the assembly name still has to be
    confirmed at UCSC cannot drift apart.

    Parameters
    ----------
    metadata : genome.assembly.metadata.AssemblyMetadata or None
        The assembly's curated row, or ``None`` for one the table does not list.
    golden_path_url : str
        The URL derived from the assembly name, used when nothing is pinned.

    Returns
    -------
    FetchedSource
        The URL to fetch, and whether it was derived rather than pinned.

    Examples
    --------
    >>> fetched_source(None, "https://example.org/hg38.fa.gz")
    FetchedSource(url='https://example.org/hg38.fa.gz', derived=True)
    """
    pinned = metadata.source_url if metadata is not None else None
    return FetchedSource(url=pinned if pinned else golden_path_url, derived=pinned is None)


def is_prepared(assembly: str) -> bool:
    """Return whether ``assembly`` is registered under the shared data root, by its record alone.

    By name and not by path, exactly as a chimera's recorded components are found again: a
    component is addressed by the key it was registered under. A directory holding files
    but no record is *not* prepared — nothing vouches for it.

    Parameters
    ----------
    assembly : str
        The assembly name to look for.

    Returns
    -------
    bool
        Whether a completion record is there.

    Examples
    --------
    >>> is_prepared("definitely-not-an-assembly")
    False
    """
    return read_record(assembly_data_dir(assembly)) is not None


def _could_be_a_component(assembly: str) -> bool:
    """Whether ``assembly`` names something a chimera could be built from.

    Prepared here, or listed in the curated table. The second clause is what separates a
    chimera's name from a free-form local key on a machine that holds neither: the shipped
    row for a chimera's components is why ``ce11_ecHT115`` reads as two assemblies while
    ``my_ref`` reads as one name somebody chose (ADR-0003).
    """
    return is_prepared(assembly) or lookup_assembly(assembly) is not None


def resolve_source(
    assembly_dir: AssemblyDir,
    *,
    metadata: AssemblyMetadata | None,
    golden_path_url: str,
) -> FetchedSource | ComponentSource:
    """Return where the assembly in ``assembly_dir`` gets its FASTA — a URL, or components.

    What a name means, decided once. Both ways in — :class:`~genome.assembly.genome.Genome` and
    :func:`~genome.assembly.download.register_assembly` — reach this through one registration
    step, which is what makes ``genome register <name>`` one command for every kind of
    **Source** and keeps the command line a thin client, since neither the resolution nor
    its refusals are written there. Four checks, in this order:

    1. **A record here.** An existing completion record already says what this assembly is,
       and is believed outright: a chimera's record rebuilds a chimera, and any other
       record — a plain ``hg38_mm10`` somebody seeded years ago — keeps the assembly
       whatever it was registered as. Only a record that is *lost* falls through to the
       name.
    2. **A source the caller named.** An explicitly seeded assembly is what the caller said
       it is, so it never reaches this function: ``--source`` is a :class:`SeededSource`
       handed straight to the registration. That route consults the record first as well —
       a registration already here is returned from it without the source being read again
       — so the first check stands whichever way in was taken, and the source answers only
       where there is nothing to believe.
    3. **The name.** It splits on ``_`` into two or more parts, and every one of them is
       either prepared here or listed in the curated table. That second clause is the whole
       separation between ``ce11_ecHT115`` and a free-form local key like ``my_ref`` on a
       machine holding neither (ADR-0003, ADR-0008).
    4. **Today's path**, for everything else — the FASTA is fetched from the pinned source
       or the derived golden-path URL, exactly as it always was.

    Rebuilding from scratch is not one of the checks and never overrules them: what an
    assembly *is* is not something a repair may change.

    Parameters
    ----------
    assembly_dir : genome.assembly.registration.AssemblyDir
        The **Assembly dir** in question. It carries the name being resolved and the
        directory whose record answers first.
    metadata : genome.assembly.metadata.AssemblyMetadata or None
        The assembly's curated row, or ``None`` for one the table does not list.
    golden_path_url : str
        The URL derived from the assembly name, used when nothing pins one.

    Returns
    -------
    FetchedSource or ComponentSource
        The components this assembly is concatenated from, or the URL its FASTA is
        downloaded from.

    Raises
    ------
    FileNotFoundError
        If the name spells a legal chimera's components in the wrong order. The message
        names the canonical spelling — a mis-ordered name being detectable is what leaves a
        ``--component`` flag with nothing left to buy.
    genome.store.completion.RegistrationError
        If the record here says it was built from components and cannot be read as one.
        The first check believes the record, so a record that cannot be believed is
        refused rather than falling through to the name — see
        :meth:`~genome.assembly.components.ChimeraDetails.from_details`.

    Examples
    --------
    >>> from genome.assembly.registration import AssemblyDir
    >>> nowhere = AssemblyDir.locate("my_ref", "/tmp/definitely-not-a-build")
    >>> resolve_source(nowhere, metadata=None, golden_path_url="https://example.org/x.fa.gz")
    FetchedSource(url='https://example.org/x.fa.gz', derived=True)
    """
    record = assembly_dir.read_record()
    if record is not None:
        details = ChimeraDetails.from_record(record)
        if details is None:
            return fetched_source(metadata, golden_path_url)
        return ComponentSource(tuple(details.components))
    try:
        candidates = split_name(assembly_dir.assembly)
    except ChimeraNamingError:
        return fetched_source(metadata, golden_path_url)
    if not all(_could_be_a_component(name) for name in candidates):
        return fetched_source(metadata, golden_path_url)
    canonical = derive_name(candidates)
    if canonical != assembly_dir.assembly:
        raise FileNotFoundError(
            f"nothing is registered as {assembly_dir.assembly!r} in {assembly_dir.path}, and "
            f"it is not how a chimera of {', '.join(candidates)} is spelled: a chimera's name "
            f"is its components sorted, so that one set of components means one "
            f"directory whatever order they are typed in (ADR-0008). Build it with "
            f"`genome register {canonical}`."
        )
    return ComponentSource(candidates)
