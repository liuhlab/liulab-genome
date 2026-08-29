"""Where an assembly's files live, and the act of putting them there.

I/O boundary module, and the half of registering an assembly that has nothing to do
with where the bytes came from. Registering one is always the same four steps — stage
in the working area, place the FASTA under the assembly's own name, derive the
companions, write the **Completion marker** — and only the step that *produces* the
FASTA differs between one source and another. That common part is
:class:`AssemblyRegistration`; :class:`~genome.io.download.UCSCGenomeDownloader` is it
plus a fetch.

Nothing here reaches the network, and nothing here compares a digest against a pinned
one: an assembly that pins a checksum, and the fetching that makes a pin worth having,
belong to :mod:`genome.io.download`.

The **Assembly dir** layout lives here too — the data root, the per-assembly directory
and the subtrees other contexts own — because every one of those steps is expressed in
it and it cannot be read from anywhere else without a cycle. The names stay importable
from :mod:`genome.io.download`, which is where they used to live.

Most of the **Data dir** is that assembly tree, and not all of it: :func:`motif_data_dir`
and :func:`xref_data_dir` are filed *beside* it rather than inside it, because neither a
motif nor an identifier belongs to an assembly. They are here so that what lives under the
data root is legible in one file — what lives under ``motif/`` and ``xref/`` is each
context's own business and is spelled there.

:class:`AssemblyDir` is that layout as a **value**: where an assembly lives is settled
once, by :meth:`AssemblyDir.locate`, and then carried. A caller holding one — a
:class:`~genome.genome.Genome`, an :class:`~genome.aligner.aligner.Aligner` asking that
genome where its index goes — cannot resolve the directory differently from the caller
that opened it, which is what a free function reading the environment at each call site
could not promise.

Examples
--------
>>> import os
>>> from genome.io.registration import AssemblyRegistration
>>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
>>> AssemblyRegistration("hg38").cache_dir
PosixPath('/scratch/liulab/genome/hg38')
>>> del os.environ["LIULAB_DATA"]
"""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genome.io.completion import (
    CompletionRecord,
    build_record,
    check_registration,
    clear_work_dir,
    read_record,
    record_path,
    work_dir,
    write_record,
)
from genome.io.fasta import PREPARATION_TOOLS, GenomeFiles

#: Environment variable naming the lab data root directory.
LIULAB_DATA_ENV = "LIULAB_DATA"

#: Well-known lab data roots, tried in order when ``LIULAB_DATA`` is unset.
DEFAULT_LIULAB_DATA_PATHS = [
    "/share/lhqlab/liulab_data",
    "/large_storage/zhoulab/hanliu/liulab_data",
]

#: Subdirectory of the **Data dir** holding the assembly tree — one directory per
#: **Assembly** under it, and most of the data root.
ASSEMBLIES_SUBDIR = "genome"

#: Subdirectory of the **Data dir** holding motif data, and the first thing filed as a
#: *sibling* of the assembly tree rather than inside it: a **Motif** belongs to no
#: assembly, so there is no assembly directory it could go under.
MOTIF_SUBDIR = "motif"

#: Subdirectory of the **Data dir** holding **Xref set**s, a sibling of the assembly tree
#: beside ``motif/`` and for the same reason: an identifier is a name and not a place, so
#: it names no **Assembly** either. What lives *under* it — a directory per **Xref
#: source**, release and species — is the Xref context's business and is spelled there.
XREF_SUBDIR = "xref"

#: Subdirectory of an **Assembly dir** holding its annotations. The Assembly context
#: owns the layout, so the name lives here and the Annotation context reads it.
ANNOTATIONS_SUBDIR = "gtf"

#: Subdirectory of an **Assembly dir** holding its aligner indexes, likewise.
INDEXES_SUBDIR = "index"

#: The subtrees inside an **Assembly dir** that other contexts own. Each carries its own
#: completion record, so an assembly must not count them as files it failed to claim —
#: an annotation registered before its assembly is not an interrupted assembly.
_FOREIGN_SUBDIRS = frozenset({ANNOTATIONS_SUBDIR, INDEXES_SUBDIR})


def liulab_data_dir() -> Path:
    """Return the root directory for lab reference data.

    The location is read from the ``LIULAB_DATA`` environment variable. When that
    is unset (or empty), each entry in :data:`DEFAULT_LIULAB_DATA_PATHS` is checked
    in order and the first that exists is used as the root. If none exist, it falls
    back to ``~/liulab_data``. The path is expanded (``~`` resolved) but **not**
    created here — callers create the specific subdirectory they need on first write.

    Returns
    -------
    pathlib.Path
        The resolved lab data root.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> liulab_data_dir()
    PosixPath('/scratch/liulab')
    >>> del os.environ["LIULAB_DATA"]
    """
    env = os.environ.get(LIULAB_DATA_ENV)
    if env:
        return Path(env).expanduser()
    for candidate in DEFAULT_LIULAB_DATA_PATHS:
        path = Path(candidate).expanduser()
        if path.exists():
            return path
    return (Path.home() / "liulab_data").expanduser()


def assembly_data_dir(assembly: str) -> Path:
    """Return the directory holding all reference files for ``assembly``.

    Every file tied to a reference assembly (FASTA, indexes, annotations, …)
    lives under ``<liulab_data>/genome/<assembly>/`` so they stay co-located.

    Parameters
    ----------
    assembly : str
        Assembly name, e.g. ``"hg38"``.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/genome/<assembly>``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> assembly_data_dir("hg38")
    PosixPath('/scratch/liulab/genome/hg38')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / ASSEMBLIES_SUBDIR / assembly


def motif_data_dir() -> Path:
    """Return the directory holding motif data, which belongs to no assembly.

    ``<liulab_data>/motif/``, a **sibling** of the assembly tree rather than a tenant of
    it: a **Motif** is a pattern and not a place, so it names no **Assembly** and there is
    no per-assembly directory it could be filed under. Shared by every project on the
    machine, exactly as an assembly is. Nothing is created here — the caller that writes
    creates what it needs.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/motif``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> motif_data_dir()
    PosixPath('/scratch/liulab/motif')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / MOTIF_SUBDIR


def xref_data_dir() -> Path:
    """Return the directory holding **Xref set**s, which belong to no assembly.

    ``<liulab_data>/xref/``, a **sibling** of the assembly tree and of ``motif/``: an
    identifier is a name and not a place, so an **Xref set** names no **Assembly** and
    answers with no **Genome** open. Shared by every project on the machine. Nothing is
    created here — the caller that writes creates what it needs.

    Returns
    -------
    pathlib.Path
        ``<liulab_data>/xref``.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> xref_data_dir()
    PosixPath('/scratch/liulab/xref')
    >>> del os.environ["LIULAB_DATA"]
    """
    return liulab_data_dir() / XREF_SUBDIR


def assembly_repair_command(assembly: str, source: str | Path | None = None) -> str:
    """Return the command that registers ``assembly`` again from scratch.

    One spelling, wherever it is quoted: a broken **Assembly dir** names it, and so does a
    **Merged annotation** whose only repair is rebuilding the chimera that wrote it. A
    seeded assembly carries its own source into it — ``genome register tiny --force``
    would fetch from the golden path, which is not where such an assembly came from.

    Parameters
    ----------
    assembly : str
        The assembly to register again.
    source : str or pathlib.Path, optional
        Where its FASTA came from, for an assembly that was seeded rather than fetched.

    Returns
    -------
    str
        A command that runs as it stands.

    Examples
    --------
    >>> assembly_repair_command("hg38")
    'genome register hg38 --force'
    >>> assembly_repair_command("tiny", "/data/my ref.fa")
    "genome register tiny --force --source '/data/my ref.fa'"
    """
    base = f"genome register {assembly} --force"
    return base if source is None else f"{base} --source {shlex.quote(str(source))}"


@dataclass(frozen=True)
class AssemblyDir:
    """One **Assembly**'s directory, and everything the layout puts inside it.

    The **Assembly dir** as a value rather than a rule applied again at each call site.
    Where an assembly lives is decided once — by :meth:`locate`, which is the only
    implementation of *an explicit directory overrides the* **Data dir** *layout* — and
    the answer is then carried, so a caller holding one cannot resolve it differently
    from the caller that opened it. That is not hypothetical: an index derived from the
    data root rather than from the assembly it indexes lands beside a different
    assembly's files, and reads a completion record that is not there.

    Every subtree an assembly owns is spelled here, so the layout is legible in one
    place: the working area, the ``gtf/`` subtree the Annotation context files into, the
    ``index/`` subtree the Index context files into, and the four **Genome files**.

    Parameters
    ----------
    assembly : str
        The assembly name — the key this directory is addressed by, and the name its
        files carry.
    path : pathlib.Path
        The directory itself.

    Examples
    --------
    >>> import os
    >>> os.environ["LIULAB_DATA"] = "/scratch/liulab"
    >>> here = AssemblyDir.locate("hg38")
    >>> here.path
    PosixPath('/scratch/liulab/genome/hg38')
    >>> here.index_dir("chromap")
    PosixPath('/scratch/liulab/genome/hg38/index/chromap')
    >>> here.sibling("mm39").path
    PosixPath('/scratch/liulab/genome/mm39')
    >>> del os.environ["LIULAB_DATA"]
    """

    assembly: str
    path: Path

    @classmethod
    def locate(cls, assembly: str, cache_dir: str | Path | None = None) -> AssemblyDir:
        """Return where ``assembly`` lives: ``cache_dir`` when given, else the layout's answer.

        The one place the override rule is written. ``cache_dir`` is the directory
        itself and not a root to file under, which is what lets a caller put one
        assembly somewhere of its own without inventing a second **Data dir**.

        Parameters
        ----------
        assembly : str
            The assembly to locate.
        cache_dir : str or pathlib.Path, optional
            An explicit **Assembly dir**, overriding the layout.

        Returns
        -------
        AssemblyDir
            Where that assembly's files belong. Nothing is created and nothing is read.

        Examples
        --------
        >>> AssemblyDir.locate("hg38", "/tmp/elsewhere").path
        PosixPath('/tmp/elsewhere')
        """
        if cache_dir is not None:
            return cls(assembly=assembly, path=Path(cache_dir).expanduser())
        return cls(assembly=assembly, path=assembly_data_dir(assembly))

    def sibling(self, assembly: str) -> AssemblyDir:
        """Return another assembly's **Assembly dir**, found beside this one.

        How an assembly named in a record is found again — a **Chimera**'s components,
        say. A record carries names and never paths, which is what keeps a registered
        directory movable, so the name has to be resolved against something: this
        resolves it against *where the asking assembly is*, rather than against the
        **Data dir** the process happens to be pointed at. Under the ordinary layout the
        two agree, since every assembly is a sibling under ``<data dir>/genome/``. They
        part company exactly when one assembly was placed somewhere of its own, and then
        beside-me is the answer that finds anything at all.

        Parameters
        ----------
        assembly : str
            The other assembly's name.

        Returns
        -------
        AssemblyDir
            Its directory, beside this one. Nothing is created and nothing is read.

        Examples
        --------
        >>> AssemblyDir.locate("ce11_ecHT115", "/data/genome/ce11_ecHT115").sibling("ce11").path
        PosixPath('/data/genome/ce11')
        """
        return AssemblyDir(assembly=assembly, path=self.path.parent / assembly)

    @property
    def work_dir(self) -> Path:
        """The disposable working area a build stages in — see :func:`~genome.io.completion.work_dir`."""
        return work_dir(self.path)

    @property
    def record_path(self) -> Path:
        """Where this assembly's **Completion marker** is written, whether or not it exists."""
        return record_path(self.path)

    def read_record(self) -> CompletionRecord | None:
        """Return this assembly's completion record, or ``None`` when it has none."""
        return read_record(self.path)

    @property
    def annotations_root(self) -> Path:
        """The ``gtf/`` subtree, parent of every annotation directory."""
        return self.path / ANNOTATIONS_SUBDIR

    def annotation_dir(self, name: str) -> Path:
        """Return the directory the annotation registered as ``name`` is filed under."""
        return self.annotations_root / name

    @property
    def indexes_root(self) -> Path:
        """The ``index/`` subtree, parent of every **Index dir**."""
        return self.path / INDEXES_SUBDIR

    def index_dir(self, name: str) -> Path:
        """Return the **Index dir** of the index addressed as ``name``."""
        return self.indexes_root / name

    @property
    def genome_files(self) -> GenomeFiles:
        """The four **Genome files** this assembly's preparation produces, existing or not.

        Every way of producing the FASTA materializes it as ``<assembly>.fa`` and derives
        identically named companions, so one layout describes them all.
        """
        fasta = self.path / f"{self.assembly}.fa"
        return GenomeFiles(
            fasta=fasta,
            fai=fasta.with_name(fasta.name + ".fai"),
            twobit=self.path / f"{self.assembly}.2bit",
            chrom_sizes=self.path / f"{self.assembly}.chrom.sizes",
        )

    def completed_files(self, *, repair: str) -> GenomeFiles | None:
        """Return the prepared **Genome files** when the record vouches for them, else ``None``.

        *Is this assembly finished here?*, asked of the directory itself so that a caller
        holding one — a verification, say — needs no registration object to find out. The
        completion record is the only thing consulted: it must be there, and every file it
        claims must be present at the size it claims. That is one ``stat`` per file and no
        file contents, so reopening a prepared genome is instant.

        Parameters
        ----------
        repair : str
            The command quoted into the refusal a directory that cannot be trusted raises.

        Returns
        -------
        genome.io.fasta.GenomeFiles or None
            The four files, or ``None`` for an absent or empty directory — a fresh
            registration, which proceeds normally.

        Raises
        ------
        genome.io.completion.RegistrationError
            If the directory holds files with no record, or a record that disagrees with
            what is on disk (see :func:`~genome.io.completion.check_registration`).

        Examples
        --------
        >>> AssemblyDir.locate("hg38", "/tmp/definitely-not-a-build").completed_files(
        ...     repair="genome register hg38 --force"
        ... ) is None
        True
        """
        if check_registration(self.path, repair=repair, ignore=_FOREIGN_SUBDIRS) is None:
            return None
        return self.genome_files


class AssemblyRegistration:
    """One assembly's directory, and the steps that finish a registration in it.

    An assembly name plus the directory its files belong in, and everything a
    registration does against that pair whatever produced the FASTA: is it already
    registered, which paths does it own, where does the work happen, what repairs it,
    how does the FASTA get into place, and what does the record say afterwards.
    Subclass it and add the step that produces the FASTA.

    Every step here is a private one. The class is the seam, not a public surface: what
    a caller registers an assembly *with* is
    :func:`~genome.io.download.register_assembly` or :class:`~genome.genome.Genome`.

    Parameters
    ----------
    assembly : str
        The assembly name — the key its directory is addressed by, and the name its
        files and its record carry.
    cache_dir : str or pathlib.Path, optional
        Override the storage directory. Defaults to
        :func:`assembly_data_dir(assembly) <assembly_data_dir>`.

    Attributes
    ----------
    assembly : str
        The assembly name passed at construction.
    cache_dir : pathlib.Path
        The **Assembly dir** this registration fills.

    Examples
    --------
    >>> from pathlib import Path
    >>> registration = AssemblyRegistration("hg38", Path("/scratch/hg38"))
    >>> registration._repair_command()
    'genome register hg38 --force'
    """

    def __init__(self, assembly: str, cache_dir: str | Path | None = None) -> None:
        self.assembly = assembly
        self.dir: AssemblyDir = AssemblyDir.locate(assembly, cache_dir)

    @property
    def cache_dir(self) -> Path:
        """The **Assembly dir** this registration fills — :attr:`dir`'s path."""
        return self.dir.path

    @property
    def _work_dir(self) -> Path:
        """The disposable working area this assembly is built in.

        Inside :attr:`cache_dir` on purpose: same filesystem, so placing the finished
        FASTA is a rename rather than a copy, and it goes when the assembly does. See
        :func:`~genome.io.completion.work_dir`.
        """
        return self.dir.work_dir

    def _expected_genome_files(self) -> GenomeFiles:
        """Paths the FASTA pipeline produces for this assembly (whether or not they exist)."""
        return self.dir.genome_files

    def _repair_command(self, source: str | Path | None = None) -> str:
        """Return the command that re-registers this assembly from scratch.

        :func:`assembly_repair_command` for this registration's own assembly. Quoted
        verbatim into every error a broken directory raises, so it has to be a command
        that exists and does the job.
        """
        return assembly_repair_command(self.assembly, source)

    def _completed_genome(self, *, overwrite: bool, repair: str) -> GenomeFiles | None:
        """Return the prepared GenomeFiles when the record says so, else ``None``.

        :meth:`AssemblyDir.completed_files` for this registration's own directory, plus the
        one thing a *registration* adds to the question: ``overwrite`` skips it entirely,
        which is what makes it the repair.
        """
        if overwrite:
            return None
        return self.dir.completed_files(repair=repair)

    def _place_fasta(self, unpacked: Path) -> Path:
        """Move ``unpacked`` out of the working area to ``<assembly>.fa`` and return it.

        A rename within one filesystem, since the working area sits inside
        :attr:`cache_dir` — no second copy of a whole genome is ever made.
        """
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        fasta = self._expected_genome_files().fasta
        if unpacked.resolve() != fasta.resolve():
            unpacked.replace(fasta)
        return fasta

    def _record_completion(
        self,
        files: GenomeFiles,
        *,
        source_url: str | None,
        sha256: str | None,
        details: dict[str, Any] | None = None,
    ) -> None:
        """Write this assembly's completion record, then discard the working area.

        Called last, once every derived file exists — writing the record is what makes
        the registration finished, and the archive is only disposable after that. An
        interrupted run therefore leaves its work in place and repairs without
        producing a whole genome again.

        ``details`` is whatever is particular to how this FASTA was produced and cannot
        be read off the files themselves — a chimera's separator and components. A
        download has none: its source URL and digest are fields of their own.
        """
        record = build_record(
            self.cache_dir,
            kind="genome",
            name=self.assembly,
            files=[files.fasta, files.fai, files.twobit, files.chrom_sizes],
            source_url=source_url,
            sha256=sha256,
            tools=PREPARATION_TOOLS,
            details=details,
        )
        write_record(self.cache_dir, record)
        clear_work_dir(self.cache_dir)
