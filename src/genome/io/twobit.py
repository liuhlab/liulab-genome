"""Read sequence from UCSC 2bit files via ``py2bit``.

I/O boundary module: it opens ``.2bit`` files (built by
:func:`genome.io.fasta.fasta_to_2bit`) through the ``py2bit`` extension rather
than parsing them in Python. ``py2bit``'s ``sequence(chrom, start, end)`` already
uses **0-based, half-open** coordinates — the package's canonical internal
convention — so no coordinate conversion happens here.

The :class:`TwoBit` class wraps an **open** file handle so repeated queries reuse
it (opening also loads the soft-mask index, so re-opening per query is wasteful).
Use it as a context manager, or call :meth:`TwoBit.close` when done.

Examples
--------
>>> from genome.io.twobit import TwoBit
>>> with TwoBit("sacCer3.2bit") as tb:           # doctest: +SKIP
...     tb.sequence("chrIV", 0, 10)
'ACACCACACC'
"""

from __future__ import annotations

from pathlib import Path
from types import TracebackType
from typing import Self

import py2bit


class TwoBit:
    """An open UCSC 2bit file, queried in 0-based half-open coordinates.

    The handle is opened at construction and held until :meth:`close` (or the
    context-manager exit). Reuse a single instance for many queries rather than
    re-opening per lookup.

    Parameters
    ----------
    path : str or pathlib.Path
        Path to a ``.2bit`` file.
    masked : bool, default True
        Preserve soft-masking: lower-case bases for repeat-masked regions. When
        ``False`` every base is upper-cased. (Soft-masking carries meaning, so it
        is kept by default; see :mod:`genome.seq`.)

    Attributes
    ----------
    path : pathlib.Path
        The 2bit file path.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    RuntimeError
        If ``path`` exists but cannot be opened as a 2bit file.

    Examples
    --------
    >>> tb = TwoBit("sacCer3.2bit")               # doctest: +SKIP
    >>> tb.sequence("chrIV", 0, 10)               # doctest: +SKIP
    'ACACCACACC'
    >>> tb.close()                                # doctest: +SKIP
    """

    def __init__(self, path: str | Path, *, masked: bool = True) -> None:
        self.path = Path(path)
        self._masked = masked
        if not self.path.is_file():
            raise FileNotFoundError(f"2bit file not found: {self.path}")
        try:
            handle = py2bit.open(str(self.path), masked)
        except RuntimeError as err:
            raise RuntimeError(f"could not open {self.path} as a 2bit file: {err}") from err
        self._handle = handle
        self._closed = False

    def chroms(self) -> dict[str, int]:
        """Return a ``{name: length}`` mapping of the sequences in the file."""
        self._ensure_open()
        return dict(self._handle.chroms())

    def sequence(self, chrom: str, start: int | None = None, end: int | None = None) -> str:
        """Return the sequence of ``chrom`` (or a sub-range) as a plain string.

        Coordinates are bounds-checked against the chromosome length: an ``end``
        past the end of the sequence raises rather than being silently clamped
        (``py2bit``'s default behavior).

        Parameters
        ----------
        chrom : str
            Sequence name as stored in the file.
        start : int, optional
            0-based start, inclusive. Defaults to the start of the chromosome.
        end : int, optional
            0-based end, exclusive (half-open). Defaults to the end of the
            chromosome.

        Returns
        -------
        str
            The bare sequence, case preserved when the file was opened with
            ``masked=True``.

        Raises
        ------
        ValueError
            If ``chrom`` is not in the file, ``start < 0``, ``end`` exceeds the
            chromosome length, or ``start > end``.

        Examples
        --------
        >>> tb.sequence("chrIV", 0, 10)           # doctest: +SKIP
        'ACACCACACC'
        """
        self._ensure_open()
        length = self._handle.chroms(chrom)
        if length is None:
            raise ValueError(f"unknown sequence {chrom!r} in 2bit file {self.path}.")
        resolved_start = 0 if start is None else start
        resolved_end = length if end is None else end
        if resolved_start < 0:
            raise ValueError(f"{chrom}: start must be >= 0 (0-based), got {resolved_start}.")
        if resolved_start > resolved_end:
            raise ValueError(f"{chrom}: start ({resolved_start}) is past end ({resolved_end}).")
        if resolved_end > length:
            raise ValueError(f"{chrom}: end ({resolved_end}) exceeds sequence length ({length}).")
        if resolved_start == resolved_end:
            return ""  # empty half-open interval; py2bit rejects start == end
        return self._handle.sequence(chrom, resolved_start, resolved_end)

    def nocheck_sequence(self, chrom: str, start: int | None = None, end: int | None = None) -> str:
        """Return the sequence of ``chrom`` (or a sub-range) without bounds-checking.

        Parameters
        ----------
        chrom : str
            Sequence name as stored in the file.
        start : int, optional
            0-based start, inclusive. Defaults to the start of the chromosome.
        end : int, optional
            0-based end, exclusive (half-open). Defaults to the end of the
            chromosome.

        Returns
        -------
        str
            The bare sequence, case preserved when the file was opened with
            ``masked=True``.

        Raises
        ------
        ValueError
            If ``chrom`` is not in the file.

        Examples
        --------
        >>> tb.nocheck_sequence("chrIV", 0, 10)    # doctest: +SKIP
        'ACACCACACC'
        """
        self._ensure_open()
        length = self._handle.chroms(chrom)
        if length is None:
            raise ValueError(f"unknown sequence {chrom!r} in 2bit file {self.path}.")
        resolved_start = 0 if start is None else start
        resolved_end = length if end is None else end
        if resolved_start == resolved_end:
            return ""  # empty half-open interval; py2bit rejects start == end
        return self._handle.sequence(chrom, resolved_start, resolved_end)

    def close(self) -> None:
        """Close the underlying file handle (idempotent)."""
        if getattr(self, "_closed", True):
            return
        self._handle.close()
        self._closed = True

    def _ensure_open(self) -> None:
        """Raise :class:`ValueError` if the handle has already been closed."""
        if getattr(self, "_closed", True):
            raise ValueError(f"operation on closed TwoBit file: {self.path}")

    def __enter__(self) -> Self:
        """Return ``self`` for use as a context manager."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        """Close the handle on context-manager exit."""
        self.close()

    def __del__(self) -> None:
        """Best-effort close when the object is garbage-collected."""
        self.close()

    def __repr__(self) -> str:
        """Return e.g. ``TwoBit('sacCer3.2bit', open)``."""
        state = "closed" if getattr(self, "_closed", True) else "open"
        return f"{type(self).__name__}({self.path.name!r}, {state})"
