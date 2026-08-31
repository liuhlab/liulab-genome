"""What the curated table offers, set against what is prepared on this machine.

Two questions with one answer, and the second of them had nowhere to be asked before:
:func:`~genome.assembly.source.is_prepared` answers *is this one name prepared here* for a
name a caller already holds, and nothing enumerated the tree. "Which assemblies are here"
has no name to scope it, and its answer includes names no row lists — an assembly
registered from the UCSC golden path is a legitimate registration with no pinned row
behind it.

Reading only, and cheaply: the curated table, the directory names under the assembly tree,
and one **Completion marker** per directory that is there. Nothing is fetched, prepared,
built or created to answer, so a fresh install — where the tree does not exist at all — is
the ordinary case rather than the edge one.

**No integrity state.** A directory is registered when a record vouches for it and not
otherwise; whether the files that record claims are still what it claims is the question
``genome assembly verify`` exists to answer, deliberately, by re-reading the FASTA. Asking
it cheaply here would report *unchecked* in the words of *checked*, and asking it properly
would cost a whole genome's bytes on the one command a reader runs first.

The Annotation context's :class:`~genome.annotation.registry.AnnotationStatus` is the same
report one level down, and the shapes deliberately match — one row shape for every state,
columns that do not apply carried as ``None``, a ``state`` derived rather than stored.
Where they part company is that one: an annotation has no verify command, so ``broken``
folds into its listing because there is nowhere else for it to surface.

Examples
--------
>>> from genome.assembly.status import assembly_status
>>> status = assembly_status(root="/tmp/definitely-not-a-data-root")
>>> "hg38" in {row.assembly_name for row in status.assemblies}
True
>>> status.registered
()
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from genome.assembly.metadata import AssemblyMetadata, assembly_table
from genome.assembly.registration import AssemblyDir, assembly_root_dir


@dataclass(frozen=True)
class AssemblyStatusRow:
    """One assembly, in whichever of its states it is: offered, registered, or merely here.

    One shape for all of them, so a reader never has to ask which fields a row has — a
    name the table does not list carries the table's columns as ``None``, and one nothing
    is prepared for carries no :attr:`directory`.

    :attr:`registered` and :attr:`present` are not independent: a record lives *in* the
    directory it vouches for, so a registered assembly is always present, and the pair
    that needs saying out loud is the other one — **present and not registered**. That is
    the case enumerating the tree meets and asking after a single name never did: a
    directory that is neither a good registration nor absent. Reporting it as absent lies
    to somebody looking at a full disk, and reporting it as an assembly lies about what is
    trustworthy, so it is reported as what it is and :attr:`state` says so in three words.

    Attributes
    ----------
    assembly_name : str
        The assembly this row is about — the key its directory is addressed by.
    offered : bool
        Whether the curated table lists it.
    registered : bool
        Whether a record here vouches for it.
    present : bool
        Whether a directory for it is in the assembly tree, registered or not.
    directory : str or None
        Its **Assembly dir**, when one is there; ``None`` when nothing is prepared for it.
        A row that is not :attr:`present` carries none rather than the path one *would*
        take, which is the layout's answer and not a fact about this machine.
    species : str or None
        The species the table names; ``None`` for an assembly no row lists.
    ucsc_name : str or None
        Its name in UCSC's namespace, from the table; ``None`` for an unlisted one, and
        also for a listed reference UCSC has never carried.
    ncbi_name : str or None
        Its name in NCBI's namespace, from the table; ``None`` for an unlisted one.
    source_url : str or None
        Where the table says its FASTA is fetched from; ``None`` for an unlisted one, and
        for a **Chimera**, which is built rather than downloaded.
    sha256 : str or None
        The digest the table pins for the unpacked FASTA; ``None`` when it pins none.

    Examples
    --------
    >>> row = AssemblyStatusRow(
    ...     assembly_name="hg38",
    ...     offered=True,
    ...     registered=False,
    ...     present=False,
    ...     directory=None,
    ...     species="Homo sapiens",
    ...     ucsc_name="hg38",
    ...     ncbi_name="GRCh38",
    ...     source_url="https://example.org/hg38.fa.gz",
    ...     sha256=None,
    ... )
    >>> row.state
    'offered, not registered'
    >>> row.as_json()["offered"]
    True
    """

    assembly_name: str
    offered: bool
    registered: bool
    present: bool
    directory: str | None
    species: str | None
    ucsc_name: str | None
    ncbi_name: str | None
    source_url: str | None
    sha256: str | None

    @property
    def state(self) -> str:
        """Which of its states this row is in, in the words a surface prints.

        ``registered`` first, because it is the strongest thing that can be said and it
        settles the two weaker ones: a registered assembly is present, and whether the
        table also offers it is the only distinction left. ``here, not registered`` comes
        next for the same reason it exists at all — a directory is there, and answering
        with what the *table* says about the name would leave the reader to discover the
        directory some other way.

        Returns
        -------
        str
            One of ``registered``, ``registered, not offered``, ``here, not registered``,
            or ``offered, not registered``.
        """
        if self.registered:
            return "registered" if self.offered else "registered, not offered"
        if self.present:
            return "here, not registered"
        return "offered, not registered"

    def as_json(self) -> dict[str, Any]:
        """Return this row as ``--json`` serializes it: every attribute above, in order.

        :attr:`state` is not among them: it is read from :attr:`offered`,
        :attr:`registered` and :attr:`present`, which are all here, so writing it out too
        would be a second spelling of the same rule for a reader to disagree with.

        Returns
        -------
        dict
            The row's fields, under their own names.
        """
        return asdict(self)


@dataclass(frozen=True)
class AssemblyStatus:
    """What the curated table offers, set against what the assembly tree holds here.

    :func:`assembly_status`'s answer, and what ``genome assembly list`` prints. It is the
    first question a new user has — ``genome assembly register`` is the first command
    anyone runs and nothing else in the CLI says what may follow it — and the question
    somebody landing on a machine they did not set up has: what is already prepared on it.

    Attributes
    ----------
    directory : pathlib.Path
        The assembly tree's root, whether or not anything is there.
    assemblies : tuple of AssemblyStatusRow
        One row per name: the offered ones in table order, then whatever this tree holds
        that no row lists, in name order.

    Examples
    --------
    >>> status = assembly_status(root="/tmp/definitely-not-a-data-root")
    >>> status.registered
    ()
    >>> status.as_json()["directory"]
    '/tmp/definitely-not-a-data-root'
    """

    directory: Path
    assemblies: tuple[AssemblyStatusRow, ...]

    @property
    def registered(self) -> tuple[str, ...]:
        """The assemblies a record here vouches for, in the order they are reported.

        *What do I actually have on this machine* — the answer as data rather than as
        printed lines, for a caller who imports this instead of running the command.
        """
        return tuple(row.assembly_name for row in self.assemblies if row.registered)

    @property
    def summary(self) -> str:
        """The closing line: what is registered here, and the command that changes it.

        Two answers, and both name what to run next: nothing is registered, which is a
        fresh install's ordinary state and needs the command that prepares one; or
        something is, and then a reader who doubts one has a command for that too —
        ``verify`` is where integrity is settled, and is named here because this report
        deliberately does not settle it.

        Returns
        -------
        str
            One line, beginning ``registered here: ``.
        """
        registered = self.registered
        if not registered:
            return "registered here: (none) — prepare one with `genome assembly register <name>`"
        return (
            f"registered here: {len(registered)} — prepare another with "
            f"`genome assembly register <name>`, or re-check one with "
            f"`genome assembly verify <name>`"
        )

    @property
    def unregistered_note(self) -> str | None:
        """What a directory nothing vouches for means, when the tree holds one; else ``None``.

        Absent when there is nothing to explain, so the ordinary listing does not carry a
        sentence about a state nothing is in.

        Returns
        -------
        str or None
            One line about the ``here, not registered`` rows, or ``None`` when there are
            none.
        """
        if not any(row.present and not row.registered for row in self.assemblies):
            return None
        return (
            "`here, not registered` is a directory with no record beside it — nothing "
            "vouches for what is in it; `genome assembly register <name> --force` "
            "prepares it again from scratch"
        )

    def as_json(self) -> dict[str, Any]:
        """Return this report as ``--json`` serializes it.

        Returns
        -------
        dict
            The tree's ``directory`` as text, and ``assemblies`` as a list of
            :meth:`AssemblyStatusRow.as_json` rows.
        """
        return {
            "directory": str(self.directory),
            "assemblies": [row.as_json() for row in self.assemblies],
        }


def assembly_status(*, root: str | Path | None = None) -> AssemblyStatus:
    """Report what the curated table offers against what is prepared on this machine.

    What ``genome assembly list`` runs, and the whole of it: the command formats this and
    holds no rule of its own. Unlike :func:`~genome.annotation.registry.annotation_status`
    it takes no assembly to scope it — its scope is every assembly, offered or here.

    Nothing is fetched, prepared, built or created to answer. A tree that does not exist
    is a fresh install rather than a failure, and reports the table's rows with nothing
    registered.

    **What in the tree counts as an assembly** is a rule rather than a guess, because
    enumerating a directory means treating whatever is in it as the answer: the layout
    files one directory per assembly directly under the root and names it for the
    assembly, so an entry counts when it is a directory whose name does not begin with a
    dot. A file there belongs to no assembly, and no assembly is registered under a hidden
    name. Whether a directory that counts is *registered* is then the record's business
    and not the name's.

    Parameters
    ----------
    root : str or pathlib.Path, optional
        Override which directory the assembly tree is read from. Defaults to
        :func:`~genome.assembly.registration.assembly_root_dir`, the layout's answer. It
        is the tree itself, one level above the ``cache_dir`` that overrides a single
        **Assembly dir** — and it is not spelled ``cache_dir``, since the thing being
        pointed at is the lab's reference data rather than anything an eviction policy
        may delete.

    Returns
    -------
    AssemblyStatus
        The tree's root, and one :class:`AssemblyStatusRow` per name.

    Examples
    --------
    >>> status = assembly_status(root="/tmp/definitely-not-a-data-root")
    >>> next(row.state for row in status.assemblies if row.assembly_name == "hg38")
    'offered, not registered'
    """
    directory = Path(root).expanduser() if root is not None else assembly_root_dir()
    here = _assemblies_here(directory)
    offered = assembly_table()
    rows = [_status_row(record.assembly_name, table_row=record, here=here) for record in offered]
    listed = {record.assembly_name for record in offered}
    rows.extend(
        _status_row(name, table_row=None, here=here) for name in sorted(here.keys() - listed)
    )
    return AssemblyStatus(directory=directory, assemblies=tuple(rows))


def _assemblies_here(root: Path) -> dict[str, AssemblyDir]:
    """Return the assembly directories the tree holds, keyed by name, in name order.

    **What the tree is allowed to hold, written down rather than inferred from what
    happens to be there.** The layout files one directory per assembly directly under the
    root, named for the assembly, so an entry is an assembly's when it is a directory and
    its name does not begin with ``.``. A file there belongs to no assembly, and nothing
    registers an assembly under a hidden name, so a dotted entry is the platform's or a
    tool's rather than a registration's. A build's own working area is not among them
    either way: it lives *inside* the assembly it is building
    (:data:`~genome.store.completion.WORK_DIR_NAME`), a level below this one.

    Nothing is read to decide it: what is *in* each directory decides whether it is
    registered, and that is one record per directory, asked later and once.
    """
    if not root.is_dir():
        return {}
    return {
        entry.name: AssemblyDir(assembly=entry.name, path=entry)
        for entry in sorted(root.iterdir())
        if entry.is_dir() and not entry.name.startswith(".")
    }


def _status_row(
    name: str, *, table_row: AssemblyMetadata | None, here: dict[str, AssemblyDir]
) -> AssemblyStatusRow:
    """Return one :class:`AssemblyStatus` row, whichever of its states it is in."""
    directory = here.get(name)
    return AssemblyStatusRow(
        assembly_name=name,
        offered=table_row is not None,
        registered=directory is not None and directory.is_registered,
        present=directory is not None,
        directory=str(directory.path) if directory is not None else None,
        species=table_row.species if table_row is not None else None,
        ucsc_name=table_row.ucsc_name if table_row is not None else None,
        ncbi_name=table_row.ncbi_name if table_row is not None else None,
        source_url=table_row.source_url if table_row is not None else None,
        sha256=table_row.sha256 if table_row is not None else None,
    )
