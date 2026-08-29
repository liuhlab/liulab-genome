"""The **Motif** — what one transcription factor recognises, held as counts.

A motif is a leaf, in the sense a :class:`~genome.seq.DNA` is: it names no **Assembly**, no
**Region** and no **Strand**, and it holds no **Background**. Build one from a **Count
matrix** and it can answer everything about itself with nothing else loaded — no file, no
download, no network.

**The counts are the single source of truth.** Probabilities come from normalising a
column; log-odds come from a background and a pseudocount, which are *arguments and never
fields*. That is what keeps a background out of a motif's identity: one motif scored
against two backgrounds stays one motif, and two motifs sharing a **Motif id** can never
differ. Counts are float rather than int because JASPAR's own records carry fractional
values.

The matrix is **4 x L** — one row per base in :data:`BASES` order, one column per
position. Positions are columns everywhere in this package. logomaker is the one library
that wants the transpose, and :meth:`Motif.plot` is the only place that transpose is
taken.

Examples
--------
>>> import numpy as np
>>> from genome.tf.motif import Motif
>>> counts = np.array(
...     [
...         [0.0, 20.0, 0.0, 0.0, 20.0, 0.0],  # A
...         [0.0, 0.0, 0.0, 0.0, 0.0, 20.0],  # C
...         [20.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # G
...         [0.0, 0.0, 20.0, 20.0, 0.0, 0.0],  # T
...     ]
... )
>>> motif = Motif("MA9999.1", "Gattaca", counts)
>>> motif
Motif(motif_id='MA9999.1', motif_name='Gattaca', length=6, offset=0)
>>> motif.consensus
DNA('GATTAC')
>>> motif.information_content.round(3)
array([2., 2., 2., 2., 2., 2.])
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any

import numpy as np
import numpy.typing as npt

from genome.seq import DNA

if TYPE_CHECKING:  # pragma: no cover - typing only, so importing genome stays cheap
    from matplotlib.axes import Axes

#: The rows of a **Count matrix**, in order. The same order MOODS and logomaker use, so
#: nothing between here and the scan engine has to permute anything.
BASES: tuple[str, str, str, str] = ("A", "C", "G", "T")

#: The shortest motif this package will scan with, and the floor :meth:`Motif.trim` will
#: not go under. A 6-mer has only 4096 possible words, so its best possible match has
#: p = 1 / 4096 = 2.44e-4 and can never reach the default 1e-4 **Threshold** — a scan
#: engine asked for one anyway would clamp to the best attainable cutoff and over-call
#: silently. Trimming that produced a 6-column motif would therefore be trimming that
#: produced an unusable one.
MIN_MOTIF_LENGTH = 7

#: What :meth:`Motif.trim` calls uninformative, in bits. A flank whose most common base
#: holds barely half the observations carries about 0.21 bits and goes; one holding 55%
#: carries about 0.29 and stays.
_DEFAULT_TRIM_THRESHOLD = 0.25

#: The **Background** :meth:`Motif.log_odds` assumes when the caller names none, and the
#: reference :attr:`Motif.information_content` is always measured against.
_UNIFORM_BACKGROUND = np.full(4, 0.25)


@dataclass(frozen=True, eq=False, repr=False)
class Motif:
    """One motif: a **Count matrix**, the identity it is addressed by, and its annotation.

    Frozen, and frozen all the way down — the count matrix is copied at construction and
    marked read-only, so ``motif.counts[0, 0] = 1`` raises rather than quietly turning a
    motif into a different one behind its own hash.

    **Equality and hashing are split on purpose.** ``==`` compares every field, the matrix
    element by element (a numpy array cannot be compared with the tuple equality a
    dataclass generates, which is why ``eq=False``). The hash covers the identity only —
    id, name, offset and shape — never the several thousand floats behind them. Equal
    motifs always agree on those four, so the contract holds; and two motifs sharing a
    **Motif id** are meant to be the same motif, so collisions are theoretical.

    Parameters
    ----------
    motif_id : str
        The **Motif id** — JASPAR's versioned matrix accession, ``"MA0139.2"``. Never
        empty: it is what a motif is addressed by.
    motif_name : str
        The **Motif name** — the factor name, ``"CTCF"``, or ``"Ptf1a::Rbpj"`` for a
        dimer. A label and not a key; names collide.
    counts : array_like
        The **Count matrix**, 4 x L: one row per base in :data:`BASES` order, one column
        per position. Copied to ``float64``. Values must be finite and non-negative and
        every column must sum above zero.
    offset : int, default 0
        Where this motif's column zero sits in the frame of the motif with this id as the
        source published it, so ``position_in_this_frame + offset`` is the position in the
        full frame. Non-zero only on the result of :meth:`trim`.
    tax_group : str, default ``""``
        The **Tax group** the source filed this motif under — ``"vertebrates"``.
    tf_class : str, default ``""``
        The structural class of the factor, e.g. ``"C2H2 zinc finger factors"``.
    tf_family : str, default ``""``
        The family within that class, e.g. ``"More than 3 adjacent zinc fingers"``.
    uniprot_ids : iterable of str, default ``()``
        UniProt accessions for the factor — several for a dimer. Stored as a tuple.
    pubmed_ids : iterable of str, default ``()``
        PubMed ids for the experiment behind the matrix. Stored as a tuple.
    data_type : str, default ``""``
        How the matrix was measured, e.g. ``"ChIP-seq"``.

    Raises
    ------
    ValueError
        If ``motif_id`` is empty, ``offset`` is negative, or the count matrix is not a
        finite, non-negative 4 x L table with every column summing above zero.

    Notes
    -----
    An empty string or an empty tuple in one of the six annotations means the source
    stated nothing there — a de novo matrix out of a model carries none of them, and
    everything on this class still works.

    Examples
    --------
    >>> import numpy as np
    >>> motif = Motif(
    ...     "MA0139.2",
    ...     "CTCF",
    ...     np.array([[9.0, 1.0], [1.0, 1.0], [0.0, 7.0], [0.0, 1.0]]),
    ...     tax_group="vertebrates",
    ...     uniprot_ids=("P49711",),
    ... )
    >>> len(motif)
    2
    >>> motif.consensus
    DNA('AG')
    >>> motif.counts.flags.writeable
    False
    """

    motif_id: str
    motif_name: str
    counts: npt.NDArray[np.float64]
    offset: int = 0
    tax_group: str = ""
    tf_class: str = ""
    tf_family: str = ""
    uniprot_ids: tuple[str, ...] = ()
    pubmed_ids: tuple[str, ...] = ()
    data_type: str = ""

    def __post_init__(self) -> None:
        """Validate the identity and the matrix, then freeze an owned copy of the counts."""
        if not self.motif_id:
            raise ValueError(
                "motif_id must not be empty: a motif is addressed by its id, so give one "
                "(JASPAR's 'MA0139.2', or any id unique within the set you are building)."
            )
        if self.offset < 0:
            raise ValueError(
                f"offset must be >= 0, got {self.offset}: it is how far into the full "
                f"motif this one's column zero sits, and trimming only ever moves inward."
            )
        matrix = np.array(self.counts, dtype=np.float64)
        self._check_matrix(matrix)
        matrix.flags.writeable = False
        object.__setattr__(self, "counts", matrix)
        object.__setattr__(self, "uniprot_ids", _as_id_tuple(self.uniprot_ids, "uniprot_ids"))
        object.__setattr__(self, "pubmed_ids", _as_id_tuple(self.pubmed_ids, "pubmed_ids"))

    def _check_matrix(self, matrix: npt.NDArray[np.float64]) -> None:
        """Reject a count matrix that is the wrong shape, the wrong sign, or all zeros."""
        if matrix.ndim != 2 or matrix.shape[0] != len(BASES) or matrix.shape[1] < 1:
            raise ValueError(
                f"{self.motif_id}: a count matrix is 4 x L — one row per base in "
                f"{''.join(BASES)} order, one column per position — got shape "
                f"{matrix.shape}. Transpose it if your positions are rows."
            )
        if not np.isfinite(matrix).all():
            raise ValueError(
                f"{self.motif_id}: the count matrix holds non-finite values at "
                f"{_offending_columns(~np.isfinite(matrix).all(axis=0))}. Drop or repair "
                f"those positions before building a motif."
            )
        if (matrix < 0).any():
            raise ValueError(
                f"{self.motif_id}: the count matrix holds negative counts at "
                f"{_offending_columns((matrix < 0).any(axis=0))}. Counts are observations; "
                f"a weighted form belongs in log_odds(), not in the matrix."
            )
        sums = matrix.sum(axis=0)
        if not (sums > 0).all():
            raise ValueError(
                f"{self.motif_id}: the count matrix has no observations at "
                f"{_offending_columns(sums <= 0)}. Every position needs at least one "
                f"observed base, or its probabilities cannot be normalised."
            )

    # ----------------------------------------------------------------- identity

    def __len__(self) -> int:
        """Return L, the number of positions."""
        return int(self.counts.shape[1])

    @property
    def length(self) -> int:
        """Number of positions, ``L`` (alias of ``len(self)``).

        Examples
        --------
        >>> import numpy as np
        >>> Motif("MA9999.1", "x", np.ones((4, 11))).length
        11
        """
        return len(self)

    def __eq__(self, other: object) -> bool:
        """Compare every field, the count matrix element by element."""
        if not isinstance(other, Motif):
            return NotImplemented
        return (
            self.motif_id == other.motif_id
            and self.motif_name == other.motif_name
            and self.offset == other.offset
            and self.tax_group == other.tax_group
            and self.tf_class == other.tf_class
            and self.tf_family == other.tf_family
            and self.uniprot_ids == other.uniprot_ids
            and self.pubmed_ids == other.pubmed_ids
            and self.data_type == other.data_type
            and np.array_equal(self.counts, other.counts)
        )

    def __hash__(self) -> int:
        """Hash the identity — id, name, offset and shape — never the counts."""
        return hash((self.motif_id, self.motif_name, self.offset, self.counts.shape))

    def __repr__(self) -> str:
        """Return the identity and the shape, never the matrix itself."""
        return (
            f"Motif(motif_id={self.motif_id!r}, motif_name={self.motif_name!r}, "
            f"length={len(self)}, offset={self.offset})"
        )

    # ----------------------------------------------------- derived from counts

    @property
    def probabilities(self) -> npt.NDArray[np.float64]:
        """The counts with each column normalised to sum to 1, as a fresh 4 x L array.

        Derived on every call and never stored — the counts are the single source of
        truth, and a stored copy is a second one waiting to disagree.

        Returns
        -------
        numpy.ndarray
            Shape ``(4, L)``, ``float64``, each column summing to 1.

        Examples
        --------
        >>> import numpy as np
        >>> motif = Motif("MA9999.1", "x", np.array([[3.0], [1.0], [0.0], [0.0]]))
        >>> motif.probabilities.ravel()
        array([0.75, 0.25, 0.  , 0.  ])
        """
        return np.asarray(self.counts / self.counts.sum(axis=0), dtype=np.float64)

    def log_odds(
        self,
        background: Sequence[float] | npt.NDArray[np.float64] | None = None,
        pseudocount: float = 0.01,
    ) -> npt.NDArray[np.float64]:
        """Score each base at each position against a background, in bits.

        The background and the pseudocount are **arguments and never fields**, so the same
        motif can be scored against two backgrounds without becoming two motifs. Each
        entry is ``log2(p / q)``, where ``p`` is the pseudocounted column probability
        ``(count + pseudocount * q) / (column_sum + pseudocount)`` and ``q`` is the
        background frequency for that base — the same arithmetic the scan engine does,
        converted from natural log to bits.

        Parameters
        ----------
        background : sequence of float, optional
            Four frequencies over :data:`BASES`, each above zero and summing to 1.
            Uniform when omitted.
        pseudocount : float, default 0.01
            Added to each column, split by the background, so an unobserved base scores
            very low rather than ``-inf``. Must be above zero.

        Returns
        -------
        numpy.ndarray
            Shape ``(4, L)``, ``float64``, in bits.

        Raises
        ------
        ValueError
            If the background is not four positive frequencies summing to 1, or the
            pseudocount is not above zero.

        Examples
        --------
        >>> import numpy as np
        >>> motif = Motif("MA9999.1", "x", np.array([[100.0], [0.0], [0.0], [0.0]]))
        >>> motif.log_odds().round(2)[:, 0]
        array([  2.  , -13.29, -13.29, -13.29])
        >>> motif.log_odds([0.4, 0.1, 0.1, 0.4]).round(2)[0, 0]
        np.float64(1.32)
        """
        if pseudocount <= 0:
            raise ValueError(
                f"pseudocount must be > 0, got {pseudocount}: without one an unobserved "
                f"base scores -inf. Pass 0.01 unless you have a reason not to."
            )
        frequencies = _check_background(background)
        weighted = self.counts + pseudocount * frequencies[:, np.newaxis]
        probabilities = weighted / weighted.sum(axis=0)
        return np.asarray(np.log2(probabilities / frequencies[:, np.newaxis]), dtype=np.float64)

    @property
    def information_content(self) -> npt.NDArray[np.float64]:
        """How much each position says, in bits, in ``[0, 2]``.

        Measured against a uniform reference, which is what puts the ceiling at 2: a
        position fixed on one base says 2 bits, one that says nothing says 0. This is the
        y-axis :meth:`plot` draws and the quantity :meth:`trim` thresholds on, so the
        height you see is the number you set. Clipped to ``[0, 2]``, which only ever moves
        floating-point noise.

        Returns
        -------
        numpy.ndarray
            Shape ``(L,)``, ``float64``, one value per position.

        Examples
        --------
        >>> import numpy as np
        >>> fixed_then_flat = np.array([[8.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 2.0]])
        >>> Motif("MA9999.1", "x", fixed_then_flat).information_content
        array([2., 0.])
        """
        probabilities = self.probabilities
        # Masked rather than np.errstate'd: log2(0) is never evaluated, so no warning is
        # raised to be suppressed and no NaN is produced to be replaced.
        terms = np.zeros_like(probabilities)
        observed = probabilities > 0.0
        terms[observed] = probabilities[observed] * np.log2(probabilities[observed])
        return np.clip(2.0 + terms.sum(axis=0), 0.0, 2.0)

    @property
    def consensus(self) -> DNA:
        """The most common base at each position, as a typed :class:`~genome.seq.DNA`.

        A one-letter-per-position rendering of the counts, and nothing more: it is not the
        motif, and a **Motif hit** matched the matrix rather than this string. Positions
        saying nothing still get a letter, so read it beside
        :attr:`information_content`. Ties go to the first base in :data:`BASES` order, and
        no IUPAC ambiguity code is ever produced — this package's alphabet is ``ACGT``.

        Returns
        -------
        DNA
            One upper-case base per position, ``L`` long.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.array([[9.0, 1.0], [1.0, 1.0], [0.0, 7.0], [0.0, 1.0]])
        >>> Motif("MA9999.1", "x", counts).consensus
        DNA('AG')
        >>> Motif("MA9999.1", "x", np.ones((4, 3))).consensus     # every column tied
        DNA('AAA')
        """
        return DNA("".join(BASES[index] for index in self.counts.argmax(axis=0)))

    # -------------------------------------------------------------- reshaping

    def trim(
        self,
        threshold: float = _DEFAULT_TRIM_THRESHOLD,
        *,
        max_length: int | None = None,
        min_length: int = MIN_MOTIF_LENGTH,
    ) -> Motif:
        """Drop uninformative flanks, keeping the id, the name and a mapping back.

        **Only the ends move.** A position under ``threshold`` bits is dropped from the
        left while every position so far has been under it, and likewise from the right —
        so an uninformative spacer in the middle of a dimeric motif is kept and the motif
        can never be split in two.

        The result carries the same **Motif id** and **Motif name** and an ``offset`` such
        that ``position_in_trimmed_frame + offset`` is the position in the full motif's
        frame. Trimming a trimmed motif composes: the offsets add up, so a **Motif hit**
        found with either is readable in the full frame.

        Two bounds hold whatever the threshold says. It never returns a motif shorter than
        ``min_length`` — if the flank walk went too far, the window grows back one
        position at a time, always taking the more informative of the two neighbours, and
        the front one when they say the same, so the offset stays as small as it can. And
        a motif already shorter than ``min_length`` is returned untouched, since trimming
        cannot repair it and making it shorter would only make it worse.

        Parameters
        ----------
        threshold : float, default 0.25
            Bits below which a flanking position is uninformative.
        max_length : int, optional
            The most positions to keep. When the surviving window is longer, the less
            informative end is dropped one position at a time until it fits — the back
            one when the two ends say the same. Must be at least ``min_length``.
        min_length : int, default 7
            The fewest positions to keep, defaulting to :data:`MIN_MOTIF_LENGTH` — the
            shortest motif this package can scan with.

        Returns
        -------
        Motif
            The trimmed motif, or ``self`` when nothing was dropped.

        Raises
        ------
        ValueError
            If ``min_length`` is below 1, or ``max_length`` is below ``min_length``.

        Examples
        --------
        >>> import numpy as np
        >>> flat = np.full((4, 3), 5.0)                       # 0 bits per position
        >>> fixed = np.repeat(np.array([[9.0], [1.0], [0.0], [0.0]]), 8, axis=1)
        >>> motif = Motif("MA9999.1", "x", np.hstack([flat, fixed, flat]))
        >>> len(motif)
        14
        >>> trimmed = motif.trim()
        >>> trimmed
        Motif(motif_id='MA9999.1', motif_name='x', length=8, offset=3)
        >>> trimmed.trim().offset                             # trimming composes
        3
        """
        if min_length < 1:
            raise ValueError(f"min_length must be >= 1, got {min_length}.")
        if max_length is not None and max_length < min_length:
            raise ValueError(
                f"max_length ({max_length}) is below min_length ({min_length}), so no "
                f"length satisfies both. Raise max_length or lower min_length."
            )
        information = self.information_content
        length = len(self)
        # A motif already under the floor cannot be brought up to it, so the floor for
        # this call is whichever is smaller — trimming never makes a short motif shorter.
        floor = min(min_length, length)

        low, high = 0, length
        while low < high and information[low] < threshold:
            low += 1
        while high > low and information[high - 1] < threshold:
            high -= 1
        if low == high:
            # Every position was under the threshold. Anchor on the one that says most and
            # grow around it, rather than keeping whichever end the walks happened to meet.
            low = int(information.argmax())
            high = low + 1
        while high - low < floor:
            left = information[low - 1] if low > 0 else -np.inf
            right = information[high] if high < length else -np.inf
            if left >= right:
                low -= 1
            else:
                high += 1
        if max_length is not None:
            while high - low > max_length:
                # Ties go to the front, as they do when the window grows above: the
                # smaller offset is the one that needs less explaining.
                if information[low] < information[high - 1]:
                    low += 1
                else:
                    high -= 1

        if (low, high) == (0, length):
            return self
        return replace(self, counts=self.counts[:, low:high], offset=self.offset + low)

    # --------------------------------------------------------------- drawing

    def plot(self, ax: Axes | None = None, **kwargs: Any) -> Axes:
        """Draw the sequence logo, in bits, and return the axes it was drawn on.

        A new figure when ``ax`` is ``None``, the caller's axes when one is given — which
        is how a grid of motifs is built in one figure. The y-axis is
        :attr:`information_content`, so the height drawn is the same quantity
        :meth:`trim` thresholds on; the transform from counts to those heights is
        logomaker's own, rather than a second derivation of it here. The x-axis is this
        motif's own frame, ``0 .. L - 1``; add :attr:`offset` to read a trimmed motif's
        positions in the full frame.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw into. A new figure and axes are made when omitted.
        **kwargs
            Passed through to ``logomaker.Logo`` — ``color_scheme``, ``font_name``,
            ``shade_below`` and the rest of what it takes.

        Returns
        -------
        matplotlib.axes.Axes
            The axes drawn on, whether it was made here or handed in.

        Examples
        --------
        >>> import matplotlib
        >>> matplotlib.use("Agg")                              # no window, no display
        >>> import numpy as np
        >>> motif = Motif("MA9999.1", "x", np.array([[9.0, 1.0], [1.0, 1.0],
        ...                                          [0.0, 7.0], [0.0, 1.0]]))
        >>> axes = motif.plot()
        >>> axes.get_ylabel()
        'Information content (bits)'
        >>> low, high = axes.get_xlim()
        >>> float(high - low)                                  # one unit per position
        2.0
        """
        # Imported here, not at module scope: pyplot picks a backend on import and costs
        # about a second, and `import genome` must not pay for a plot nobody asked for.
        import logomaker
        import matplotlib.pyplot as plt
        import pandas as pd

        if ax is None:
            _, ax = plt.subplots()
        # logomaker wants positions as rows — the transpose of the 4 x L layout everything
        # else here uses, and the only place that transpose is taken.
        frame = pd.DataFrame(self.counts.T, columns=pd.Index(list(BASES)))
        heights = logomaker.transform_matrix(
            frame, from_type="counts", to_type="information", pseudocount=0
        )
        logomaker.Logo(heights, ax=ax, **kwargs)
        ax.set_xlabel("Position")
        ax.set_ylabel("Information content (bits)")
        # One unit per position and the full range of the quantity, both set rather than
        # left to autoscale: two logos in one grid are then read against one another.
        ax.set_xlim(-0.5, len(self) - 0.5)
        ax.set_ylim(0.0, 2.0)
        return ax


def _as_id_tuple(ids: Iterable[str], field: str) -> tuple[str, ...]:
    """Freeze an annotation's ids into a tuple, rejecting a bare string."""
    if isinstance(ids, str):
        raise ValueError(
            f"{field} must be an iterable of ids, not the string {ids!r}: pass "
            f"('{ids}',) for one, since a bare string would be read letter by letter."
        )
    return tuple(ids)


def _offending_columns(mask: npt.NDArray[np.bool_]) -> str:
    """Render the positions a mask flags, for an error message that names them."""
    positions = np.flatnonzero(mask).tolist()
    shown = ", ".join(str(position) for position in positions[:8])
    return f"position(s) {shown}" + (", ..." if len(positions) > 8 else "")


def _check_background(
    background: Sequence[float] | npt.NDArray[np.float64] | None,
) -> npt.NDArray[np.float64]:
    """Return a validated 4-frequency background array, uniform when none was given."""
    if background is None:
        return _UNIFORM_BACKGROUND
    frequencies = np.asarray(background, dtype=np.float64)
    if frequencies.shape != (len(BASES),):
        raise ValueError(
            f"background must be {len(BASES)} frequencies over {''.join(BASES)}, got "
            f"shape {frequencies.shape}. Higher-order backgrounds are not supported."
        )
    if not (frequencies > 0).all():
        raise ValueError(
            f"background frequencies must all be > 0, got {frequencies.tolist()}: a base "
            f"with frequency 0 makes every log-odds against it infinite."
        )
    if not np.isclose(frequencies.sum(), 1.0):
        raise ValueError(
            f"background frequencies must sum to 1, got {frequencies.sum()}. Normalise "
            f"them before passing them in."
        )
    return frequencies
