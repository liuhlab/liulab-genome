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
from pathlib import Path

from genome.io.completion import (
    build_record,
    check_registration,
    clear_work_dir,
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
    return liulab_data_dir() / "genome" / assembly


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
        self.cache_dir: Path = Path(assembly_data_dir(assembly) if cache_dir is None else cache_dir)

    @property
    def _work_dir(self) -> Path:
        """The disposable working area this assembly is built in.

        Inside :attr:`cache_dir` on purpose: same filesystem, so placing the finished
        FASTA is a rename rather than a copy, and it goes when the assembly does. See
        :func:`~genome.io.completion.work_dir`.
        """
        return work_dir(self.cache_dir)

    def _expected_genome_files(self) -> GenomeFiles:
        """Paths the FASTA pipeline produces for this assembly (whether or not they exist).

        Every way of producing the FASTA materializes it as ``<assembly>.fa`` and
        derives identically named companions, so a single layout describes them all.
        """
        fasta = self.cache_dir / f"{self.assembly}.fa"
        return GenomeFiles(
            fasta=fasta,
            fai=fasta.with_name(fasta.name + ".fai"),
            twobit=self.cache_dir / f"{self.assembly}.2bit",
            chrom_sizes=self.cache_dir / f"{self.assembly}.chrom.sizes",
        )

    def _repair_command(self, source: str | Path | None = None) -> str:
        """Return the command that re-registers this assembly from scratch.

        Quoted verbatim into every error a broken directory raises, so it has to be a
        command that exists and does the job. A seeded assembly carries its own source
        into it: ``genome register tiny --force`` would fetch from the golden path,
        which is not where such an assembly came from.
        """
        base = f"genome register {self.assembly} --force"
        return base if source is None else f"{base} --source {shlex.quote(str(source))}"

    def _completed_genome(self, *, overwrite: bool, repair: str) -> GenomeFiles | None:
        """Return the prepared GenomeFiles when the record says so, else ``None``.

        The completion record is the only thing consulted: it must be there, and every
        file it claims must be present at the size it claims. That is one ``stat`` per
        file and no file contents, so reopening a prepared genome is instant. An absent
        or empty directory answers ``None`` — a fresh registration, which proceeds
        normally — while a directory that cannot be trusted raises (see
        :func:`~genome.io.completion.check_registration`). ``overwrite`` skips the
        question entirely, which is what makes it the repair.
        """
        if overwrite:
            return None
        if check_registration(self.cache_dir, repair=repair, ignore=_FOREIGN_SUBDIRS) is None:
            return None
        return self._expected_genome_files()

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
        self, files: GenomeFiles, *, source_url: str | None, sha256: str | None
    ) -> None:
        """Write this assembly's completion record, then discard the working area.

        Called last, once every derived file exists — writing the record is what makes
        the registration finished, and the archive is only disposable after that. An
        interrupted run therefore leaves its work in place and repairs without
        producing a whole genome again.
        """
        record = build_record(
            self.cache_dir,
            kind="genome",
            name=self.assembly,
            files=[files.fasta, files.fai, files.twobit, files.chrom_sizes],
            source_url=source_url,
            sha256=sha256,
            tools=PREPARATION_TOOLS,
        )
        write_record(self.cache_dir, record)
        clear_work_dir(self.cache_dir)
