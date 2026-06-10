"""General aligner abstraction for building genome indexes.

An :class:`Aligner` wraps one external read-mapper (STAR, BWA, …) and knows how
to build that aligner's genome index for a :class:`~genome.genome.Genome`. The
base class owns the cross-aligner plumbing — installation checking, the on-disk
layout under ``<LIULAB_DATA>/genome/<assembly>/index/<name>/``, the ``.success``
flag, and the JSON parameter sidecar — while each concrete subclass supplies the
aligner-specific command and its exposed parameters via :meth:`Aligner.index`.

Only index construction is implemented here; mapping/alignment is out of scope.
"""

from __future__ import annotations

import json
import subprocess
import sys
from abc import ABC, abstractmethod
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

from genome.external import ToolNotFoundError, _resolve
from genome.io.download import assembly_data_dir

if TYPE_CHECKING:
    from genome.genome import Genome

#: Name of the marker file written into an index directory on success.
_SUCCESS_FLAG = ".success"


class Aligner(ABC):
    """Base class for an external aligner that can build a genome index.

    Subclasses set the class attributes :attr:`name` (the lowercase identifier
    used in the index path) and :attr:`binary` (the executable on ``PATH``), and
    implement :meth:`install_instructions`, :meth:`_detect_version`,
    :meth:`index`, and :attr:`_artifact`.

    Constructing an aligner assumes the tool is already installed: it resolves
    :attr:`binary` and queries its version, and on failure prints installation
    instructions and raises :class:`~genome.external.ToolNotFoundError`.

    Parameters
    ----------
    genome : genome.genome.Genome
        The genome whose reference FASTA will be indexed.
    """

    #: Lowercase identifier used in the index directory path (e.g. ``"star"``).
    name: str
    #: Executable name expected on ``PATH`` (e.g. ``"STAR"``).
    binary: str

    def __init__(self, genome: Genome) -> None:
        self._genome = genome
        try:
            self._executable: str = _resolve(self.binary)
            self._version: str = self._detect_version()
        except (ToolNotFoundError, OSError, subprocess.SubprocessError) as err:
            print(self.install_instructions(), file=sys.stderr)
            raise ToolNotFoundError(
                f"{self.binary!r} is required to build a {self.name} index but could not "
                f"be run. See the installation instructions above."
            ) from err

    # -- identity / layout ---------------------------------------------------

    @property
    def assembly(self) -> str:
        """The assembly name of the bound genome."""
        return self._genome.assembly

    @property
    def version(self) -> str:
        """The installed aligner version, detected at construction."""
        return self._version

    @property
    def index_dir(self) -> Path:
        """Directory holding this aligner's index for the assembly.

        ``<LIULAB_DATA>/genome/<assembly>/index/<name>/``.
        """
        return assembly_data_dir(self.assembly) / "index" / self.name

    @property
    def _flag_path(self) -> Path:
        """Path to the success marker written once an index completes."""
        return self.index_dir / _SUCCESS_FLAG

    @property
    def _metadata_path(self) -> Path:
        """Path to the JSON sidecar recording the index parameters."""
        return self.index_dir / f"{self.name}.index.json"

    @property
    def index_path(self) -> Path:
        """The built index file or prefix this aligner consumes.

        The exact flavour (a directory, a file, or a path prefix) is decided by
        the subclass via :attr:`_artifact`. Reading this property requires a
        completed build: it raises if the :attr:`_flag_path` success marker is
        absent.

        Raises
        ------
        RuntimeError
            If no successful index exists yet for this assembly.
        """
        if not self._flag_path.is_file():
            raise RuntimeError(
                f"No successful {self.name} index for {self.assembly!r} at "
                f"{self.index_dir} (missing success flag {_SUCCESS_FLAG!r}). "
                f"Build it first, e.g. Genome.build_{self.name}_index()."
            )
        return self._artifact

    # -- subclass contract ---------------------------------------------------

    @property
    @abstractmethod
    def _artifact(self) -> Path:
        """The index file or prefix this aligner consumes (flavour-specific)."""

    @abstractmethod
    def install_instructions(self) -> str:
        """Human-readable instructions for installing this aligner."""

    @abstractmethod
    def _detect_version(self) -> str:
        """Return the installed aligner version string."""

    @abstractmethod
    def index(self, *, overwrite: bool = False, **kwargs: Any) -> Path:
        """Build the genome index and return :attr:`index_path`."""

    # -- shared helpers for subclasses ---------------------------------------

    def _run(self, args: Sequence[str]) -> None:
        """Run the aligner binary with ``args`` inside :attr:`index_dir`.

        The index directory is used as the working directory so any log files
        the tool drops in the CWD stay co-located with the index. The tool's
        stdout/stderr are inherited, so its progress and any error messages
        stream live to the console rather than being captured.

        Raises
        ------
        RuntimeError
            If the tool exits non-zero. The tool's own output (already printed
            above) carries the diagnostic detail.
        """
        try:
            subprocess.run(
                [self._executable, *args],
                cwd=self.index_dir,
                check=True,
            )
        except subprocess.CalledProcessError as err:
            raise RuntimeError(
                f"{self.binary} failed (exit {err.returncode}); see its output above "
                f"for the error. Args: {list(args)!r}"
            ) from err

    def _write_metadata(self, *, command: Sequence[str], parameters: dict[str, Any]) -> None:
        """Write the JSON parameter sidecar next to the index files."""
        payload = {
            "aligner": self.name,
            "binary": self.binary,
            "version": self.version,
            "assembly": self.assembly,
            "fasta": str(self._genome.files.fasta),
            "command": list(command),
            "parameters": parameters,
        }
        self._metadata_path.write_text(json.dumps(payload, indent=2) + "\n")

    def _mark_success(self) -> None:
        """Write the :attr:`_flag_path` success marker."""
        self._flag_path.touch()
