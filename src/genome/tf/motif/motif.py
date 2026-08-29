"""The **Motif** and the **Motif set** — what a factor recognises, and many held as one.

A motif is a leaf, in the sense a :class:`~genome.seq.DNA` is: it names no **Assembly**, no
**Region** and no **Strand**, and it holds no **Background**. Build one from a **Count
matrix** and it can answer everything about itself with nothing else loaded — no file, no
download, no network.

A :class:`MotifSet` is many of them held as one addressable group, and the container every
filter and — once they land — every scan and comparison is a method on. It is built from
*any* motifs, so a model's de novo matrices get the whole API and not a smaller one; a
**Release** read from disk is a motif set that also knows which release it is
(:class:`~genome.tf.motif.jaspar.JasparDatabase`). Nothing here reads a file or reaches the
network either: this module is still the pure core, and the loader is the edge.

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
>>> MotifSet([motif])["MA9999"]                      # the bare base id resolves too
Motif(motif_id='MA9999.1', motif_name='Gattaca', length=6, offset=0)
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, overload

import numpy as np
import numpy.typing as npt

from genome.seq import DNA
from genome.tf.motif.workers import DEFAULT_WORKERS

if TYPE_CHECKING:  # pragma: no cover - typing only, so importing genome stays cheap
    from pathlib import Path

    import pandas as pd
    from matplotlib.axes import Axes

    from genome.tf.motif.background import BackgroundArg
    from genome.tf.motif.compare import MotifComparison

#: The rows of a **Count matrix**, in order. The same order MOODS and logomaker use, so
#: nothing between here and the scan engine has to permute anything.
BASES: tuple[str, str, str, str] = ("A", "C", "G", "T")

#: The four annotations a **Motif** holds several of, in :class:`Motif`'s field order. The
#: source publishes one per half of a dimer and separates them with a semicolon, so each
#: is a tuple; ``tax_group`` and ``data_type`` are single strings and are not here.
PLURAL_ANNOTATIONS: tuple[str, str, str, str] = (
    "tf_class",
    "tf_family",
    "uniprot_ids",
    "pubmed_ids",
)

#: Every annotation :meth:`MotifSet.filter` takes as a keyword, mapped to whether its
#: values are prose (matched as a case-insensitive substring, since nobody remembers how
#: JASPAR spells a class) or ids (matched exactly, since a substring of an accession is
#: not a weaker accession — it is a different one).
_FILTERABLE: dict[str, bool] = {
    "tax_group": True,
    "tf_class": True,
    "tf_family": True,
    "data_type": True,
    "uniprot_ids": False,
    "pubmed_ids": False,
}

#: The shortest motif this package will scan with, and the floor :meth:`Motif.trim` will
#: not go under. A 6-mer has only 4096 possible words, so its best possible match has
#: p = 1 / 4096 = 2.44e-4 and can never reach the default 1e-4 **Threshold** — a scan
#: engine asked for one anyway would clamp to the best attainable cutoff and over-call
#: silently. Trimming that produced a 6-column motif would therefore be trimming that
#: produced an unusable one.
MIN_MOTIF_LENGTH = 7

#: The **Threshold** a scan is called at when the caller names none: one per-position
#: p-value, converted per **Motif** against the **Background** into the score that motif
#: must clear. Here rather than on the engine adapter because these signatures name it.
DEFAULT_THRESHOLD = 1e-4

#: What :meth:`MotifSet.scan` calls its sequence when the caller names none — a fixed
#: literal rather than a blank, so every row of every **Hit table** names its sequence.
DEFAULT_SEQUENCE_NAME = "sequence"

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
        The **Tax group** the source filed this motif under — ``"vertebrates"``. One
        value: a motif is filed under exactly one group.
    tf_class : iterable of str, default ``()``
        The structural classes of the factor, e.g. ``("C2H2 zinc finger factors",)``.
        Plural because a dimer has one per half — ``MA0119.1``'s ``NFIC::TLX1`` is a
        SMAD/NF-1 factor joined to a homeo domain factor — and the source separates them
        with a semicolon. Stored as a tuple.
    tf_family : iterable of str, default ``()``
        The families within those classes, e.g. ``("Nuclear factor 1", "NK")``. Plural
        for the same reason. Stored as a tuple.
    uniprot_ids : iterable of str, default ``()``
        UniProt accessions for the factor — several for a dimer. Stored as a tuple.
    pubmed_ids : iterable of str, default ``()``
        PubMed ids for the experiment behind the matrix. Stored as a tuple.
    data_type : str, default ``""``
        How the matrix was measured, e.g. ``"ChIP-seq"``. One value, and commas inside
        it are part of it: ``"PBM, CSA and/or DIP-chip"`` names one method.

    Raises
    ------
    ValueError
        If ``motif_id`` is empty, ``offset`` is negative, the count matrix is not a
        finite, non-negative 4 x L table with every column summing above zero, or one of
        the four plural annotations was given a bare string.

    Notes
    -----
    An empty string or an empty tuple in one of the six annotations means the source
    stated nothing there — a de novo matrix out of a model carries none of them, and
    everything on this class still works.

    **Four of the six are plural and two are not.** ``tf_class``, ``tf_family``,
    ``uniprot_ids`` and ``pubmed_ids`` hold a tuple, because the source publishes several
    of each for a dimer and separates them with a semicolon; ``tax_group`` and
    ``data_type`` hold one string each. A bare string handed to one of the plural four is
    refused rather than stored letter by letter.

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
    tf_class: tuple[str, ...] = ()
    tf_family: tuple[str, ...] = ()
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
        for field in PLURAL_ANNOTATIONS:
            object.__setattr__(self, field, _as_value_tuple(getattr(self, field), field))

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
            Uniform when omitted — a motif holds no input to derive one from, which is
            what a scan does instead.
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


class MotifNotFoundError(LookupError):
    """Nothing in this **Motif set** is addressed by that key.

    The plain absence: the key is not a **Motif id**, not a bare base id, and not a
    **Motif name** any motif in the set carries. Not an empty answer and never an empty
    collection — a set that does not hold CTCF says so rather than handing back nothing.

    A :class:`LookupError`, so it can be caught together with
    :class:`AmbiguousMotifNameError` and :class:`AmbiguousBaseIdError` and still be told
    apart from them: those two mean *too many*, and this one means *none*.

    Parameters
    ----------
    key : str
        The key that was asked for.
    size : int
        How many motifs the set holds, so an empty set explains itself.

    Attributes
    ----------
    key : str
        The key that was asked for.

    Examples
    --------
    >>> import numpy as np
    >>> empty = MotifSet([])
    >>> try:
    ...     empty["CTCF"]
    ... except LookupError as error:
    ...     print("holds no motifs" in str(error))
    True
    """

    def __init__(self, key: str, size: int) -> None:
        self.key = key
        held = (
            "this motif set holds no motifs at all"
            if size == 0
            else f"none of the {size} motifs in this motif set is addressed by it"
        )
        super().__init__(
            f"no motif is addressed by {key!r}: {held}. A motif is looked up by its full "
            f"motif id ('MA0139.2'), by the bare base id ('MA0139'), or by a motif name "
            f"that labels exactly one motif here ('CTCF'). Read `set.motif_ids` and "
            f"`set.motif_names` for what this set can answer to."
        )


class AmbiguousMotifNameError(LookupError):
    """A **Motif name** labels several motifs here, so it addresses none of them.

    A name is a label and not a key: 66 names collide in the 2024 **Release** and 71 in
    2026, so indexing a **Motif set** by one of them would have to pick one of four CTCFs
    and say nothing about the other three. It raises instead, naming every matching
    **Motif id** — which is what a caller addresses one by — and the call that hands back
    all of them.

    Parameters
    ----------
    name : str
        The **Motif name** that labels several motifs.
    motif_ids : iterable of str
        Every matching **Motif id**, in the set's own order.

    Attributes
    ----------
    name : str
        The name that was asked for.
    motif_ids : tuple of str
        Every matching **Motif id**.

    Examples
    --------
    >>> try:
    ...     raise AmbiguousMotifNameError("CTCF", ["MA0139.2", "MA1929.2"])
    ... except LookupError as error:
    ...     print("MA1929.2" in str(error))
    True
    """

    def __init__(self, name: str, motif_ids: Iterable[str]) -> None:
        self.name = name
        self.motif_ids: tuple[str, ...] = tuple(motif_ids)
        listed = ", ".join(self.motif_ids)
        super().__init__(
            f"the motif name {name!r} labels {len(self.motif_ids)} motifs here, so it "
            f"addresses none of them: {listed}. Index by the motif id of the one you mean, "
            f"or call set.by_name({name!r}) for all of them — a name is a label and only a "
            f"motif id is a key."
        )


class AmbiguousBaseIdError(LookupError):
    """A bare base id matches several **Motif id**s here, so it addresses none of them.

    Ordinarily impossible and checked anyway: a non-redundant **Release** ships exactly
    one version of each matrix, which is what makes ``MA0139`` mean ``MA0139.2`` without a
    caller having to remember the version. A set holding two versions of one matrix — a
    hand-built comparison of ``MA0139.1`` against ``MA0139.2``, or a release that is not
    the non-redundant one — breaks that, and it is said rather than resolved by whichever
    came first.

    Parameters
    ----------
    base_id : str
        The bare base id, versionless.
    motif_ids : iterable of str
        Every **Motif id** sharing it, in the set's own order.

    Attributes
    ----------
    base_id : str
        The base id that was asked for.
    motif_ids : tuple of str
        Every **Motif id** sharing it.

    Examples
    --------
    >>> try:
    ...     raise AmbiguousBaseIdError("MA0139", ["MA0139.1", "MA0139.2"])
    ... except LookupError as error:
    ...     print("MA0139.1" in str(error))
    True
    """

    def __init__(self, base_id: str, motif_ids: Iterable[str]) -> None:
        self.base_id = base_id
        self.motif_ids: tuple[str, ...] = tuple(motif_ids)
        listed = ", ".join(self.motif_ids)
        super().__init__(
            f"the base id {base_id!r} matches {len(self.motif_ids)} motifs here: {listed}. "
            f"Address the one you mean by its full motif id, version and all. A bare base "
            f"id resolves only where one version of a matrix is present, which is what a "
            f"non-redundant release guarantees and this set does not."
        )


class MotifSet:
    """**Motif**s held as one addressable group: indexed, filtered, and nothing else read.

    The container the whole motif API hangs off, and it is built from *any* motifs — a
    **Release** parsed off disk, a filtered part of one, or the de novo matrices a model
    found, which is what gives those everything a release can do rather than a smaller
    API. It reads no file and reaches no network; preparing a release does, and that is
    :class:`~genome.tf.motif.jaspar.JasparDatabase`.

    **Indexing always returns exactly one motif, never a union type.** ``set[key]``
    resolves a **Motif id**, then a bare base id, then a **Motif name** — in that order,
    and the first that matches wins. A name labelling several motifs raises rather than
    picking one of them; :meth:`by_name` is how all of them are had, and it hands back a
    tuple whether the name labels four motifs or one.

    Iterating yields :class:`Motif` objects rather than keys — the one place this reads
    unlike a dict, and it reads like the collection of motifs it is. Order is the order
    the motifs were given, everywhere: the tuple from :meth:`by_name`, the ids in an
    ambiguity error, and what a filtered set holds.

    Parameters
    ----------
    motifs : iterable of Motif
        The motifs to hold, in the order they should be held. May be empty — a filter
        that matched nothing is a real answer.

    Raises
    ------
    ValueError
        If two motifs share a **Motif id**. An id is what a motif is addressed by, so a
        set holding two of one id could not answer for either.

    Examples
    --------
    >>> import numpy as np
    >>> counts = np.array([[9.0, 1.0], [1.0, 1.0], [0.0, 7.0], [0.0, 1.0]])
    >>> ctcf = Motif("MA0139.2", "CTCF", counts, tf_class=("C2H2 zinc finger factors",))
    >>> prrx2 = Motif("MA0075.3", "PRRX2", counts, tf_class=("Homeo domain factors",))
    >>> motifs = MotifSet([ctcf, prrx2])
    >>> len(motifs)
    2
    >>> motifs["MA0139.2"] is motifs["MA0139"] is motifs["CTCF"]
    True
    >>> motifs.by_name("PRRX2")
    (Motif(motif_id='MA0075.3', motif_name='PRRX2', length=2, offset=0),)
    >>> motifs.filter(tf_class="zinc finger").motif_ids
    ('MA0139.2',)
    """

    def __init__(self, motifs: Iterable[Motif]) -> None:
        held = tuple(motifs)
        by_id: dict[str, Motif] = {}
        for motif in held:
            if motif.motif_id in by_id:
                raise ValueError(
                    f"two motifs in this set share the motif id {motif.motif_id!r}, so "
                    f"neither could be addressed by it. Give one of them an id of its own, "
                    f"or drop the duplicate — a motif name may be shared and an id may not."
                )
            by_id[motif.motif_id] = motif
        self._motifs = held
        self._by_id = by_id
        self._by_base = _group(held, base_id)
        self._by_name = _group(held, lambda motif: motif.motif_name)

    # ----------------------------------------------------------------- what is held

    @property
    def motifs(self) -> tuple[Motif, ...]:
        """Every **Motif** held, in the order it was given.

        Examples
        --------
        >>> import numpy as np
        >>> MotifSet([Motif("MA9999.1", "x", np.ones((4, 7)))]).motifs
        (Motif(motif_id='MA9999.1', motif_name='x', length=7, offset=0),)
        """
        return self._motifs

    @property
    def motif_ids(self) -> tuple[str, ...]:
        """Every **Motif id**, in the order the motifs were given. Unique, always.

        Examples
        --------
        >>> import numpy as np
        >>> MotifSet([Motif("MA9999.1", "x", np.ones((4, 7)))]).motif_ids
        ('MA9999.1',)
        """
        return tuple(self._by_id)

    @property
    def motif_names(self) -> tuple[str, ...]:
        """Every **Motif name**, one per motif and parallel to :attr:`motif_ids`.

        Duplicates are kept and never collapsed: a name labelling two motifs appears
        twice, because these are labels and there are as many of them as motifs.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.ones((4, 7))
        >>> pair = [Motif("MA0139.2", "CTCF", counts), Motif("MA1929.2", "CTCF", counts)]
        >>> MotifSet(pair).motif_names
        ('CTCF', 'CTCF')
        """
        return tuple(motif.motif_name for motif in self._motifs)

    def __len__(self) -> int:
        """Return how many motifs are held."""
        return len(self._motifs)

    def __iter__(self) -> Iterator[Motif]:
        """Iterate the **Motif**s themselves, in order — not their keys."""
        return iter(self._motifs)

    def __contains__(self, key: object) -> bool:
        """Answer whether a key matches anything here, or whether a **Motif** is one of ours.

        A key that matches *several* motifs is still in the set, though indexing by it
        raises: membership asks whether the set knows the key, and indexing asks it to
        name one motif.

        Examples
        --------
        >>> import numpy as np
        >>> motifs = MotifSet([Motif("MA0139.2", "CTCF", np.ones((4, 7)))])
        >>> "MA0139" in motifs, "CTCF" in motifs, "SOX2" in motifs
        (True, True, False)
        """
        if isinstance(key, Motif):
            return self._by_id.get(key.motif_id) == key
        if not isinstance(key, str):
            return False
        return key in self._by_id or key in self._by_base or key in self._by_name

    def __repr__(self) -> str:
        """Return how many motifs are held, never the motifs themselves."""
        return f"{type(self).__name__}(motifs={len(self._motifs)})"

    # --------------------------------------------------------------------- lookup

    def __getitem__(self, key: str) -> Motif:
        """Return the one **Motif** ``key`` addresses — id, base id, or unique name.

        Parameters
        ----------
        key : str
            A **Motif id** (``"MA0139.2"``), a bare base id (``"MA0139"``), or a **Motif
            name** labelling exactly one motif here (``"CTCF"``). Tried in that order.

        Returns
        -------
        Motif
            Exactly one motif, always — never a tuple and never a union type.

        Raises
        ------
        MotifNotFoundError
            If nothing here is addressed by ``key``.
        AmbiguousMotifNameError
            If ``key`` is a name labelling several motifs. It names every matching id.
        AmbiguousBaseIdError
            If ``key`` is a base id shared by several motifs, which a non-redundant
            **Release** cannot produce.

        Examples
        --------
        >>> import numpy as np
        >>> motifs = MotifSet([Motif("MA0139.2", "CTCF", np.ones((4, 7)))])
        >>> motifs["MA0139"].motif_name
        'CTCF'
        """
        found = self._by_id.get(key)
        if found is not None:
            return found
        for index, ambiguity in (
            (self._by_base, AmbiguousBaseIdError),
            (self._by_name, AmbiguousMotifNameError),
        ):
            matches = index.get(key)
            if matches is None:
                continue
            # Unpacked rather than indexed: a group always holds at least one motif, and
            # this says so in a way that needs no comment to prove.
            first, *others = matches
            if others:
                raise ambiguity(key, [motif.motif_id for motif in matches])
            return first
        raise MotifNotFoundError(key, len(self._motifs))

    def by_name(self, name: str) -> tuple[Motif, ...]:
        """Return every **Motif** labelled ``name``, always as a tuple.

        A tuple of one where the name is unique, so a caller writes one code path for the
        common case and the four-CTCFs case alike. Absence is not emptiness: a name
        nothing here carries raises rather than handing back ``()``.

        Parameters
        ----------
        name : str
            The **Motif name** to look up, matched exactly.

        Returns
        -------
        tuple of Motif
            Every motif with that name, in the set's own order. Never empty.

        Raises
        ------
        MotifNotFoundError
            If no motif here carries that name.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.ones((4, 7))
        >>> pair = [Motif("MA0139.2", "CTCF", counts), Motif("MA1929.2", "CTCF", counts)]
        >>> [motif.motif_id for motif in MotifSet(pair).by_name("CTCF")]
        ['MA0139.2', 'MA1929.2']
        """
        matches = self._by_name.get(name)
        if matches is None:
            raise MotifNotFoundError(name, len(self._motifs))
        return matches

    # ------------------------------------------------------------------ filtering

    def filter(
        self, predicate: Callable[[Motif], bool] | None = None, **annotations: str | Iterable[str]
    ) -> MotifSet:
        """Return a plain :class:`MotifSet` of the motifs that match.

        **A plain motif set, never a database**, even when called on one: a filtered
        **Release** is no longer that release, and calling it one would let a **Hit
        table** claim provenance it does not have. Everything a set does, the result
        still does.

        Every condition given must hold — the predicate *and* each annotation keyword.
        With nothing given, every motif matches.

        Parameters
        ----------
        predicate : callable, optional
            Called with each :class:`Motif`; kept where it returns true. For anything the
            keywords do not express — a length, a consensus, an information content.
        **annotations
            One or more of ``tax_group``, ``tf_class``, ``tf_family``, ``data_type``,
            ``uniprot_ids`` and ``pubmed_ids``. Each value is one string or an iterable of
            them, and a motif matches when *any* of them matches *any* value the motif
            holds for that annotation. The four prose annotations match on a
            case-insensitive substring, so ``tf_class="zinc finger"`` finds every spelling
            of one; the two id annotations match exactly, since a substring of an
            accession is a different accession rather than a looser one.

        Returns
        -------
        MotifSet
            The matching motifs, in this set's own order. Empty when nothing matched,
            which is a real answer and not an absence.

        Raises
        ------
        TypeError
            If a keyword is not one of the six annotations. The message lists them.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.ones((4, 9))
        >>> zinc = Motif("MA0139.2", "CTCF", counts, tf_class=("C2H2 zinc finger factors",))
        >>> homeo = Motif("MA0075.3", "PRRX2", counts, tf_class=("Homeo domain factors",))
        >>> motifs = MotifSet([zinc, homeo])
        >>> motifs.filter(tf_class="zinc finger").motif_ids
        ('MA0139.2',)
        >>> motifs.filter(lambda motif: len(motif) == 9).motif_ids
        ('MA0139.2', 'MA0075.3')
        >>> motifs.filter(tf_class=("zinc finger", "homeo")).motif_ids
        ('MA0139.2', 'MA0075.3')
        """
        unknown = sorted(set(annotations) - set(_FILTERABLE))
        if unknown:
            raise TypeError(
                f"filter() got no annotation named {', '.join(repr(key) for key in unknown)}. "
                f"The annotations it filters on are: {', '.join(_FILTERABLE)}. Anything else "
                f"— a length, a consensus — is what the predicate argument is for."
            )
        wanted = {key: _as_wanted(value) for key, value in annotations.items()}
        kept = [
            motif
            for motif in self._motifs
            if (predicate is None or predicate(motif))
            and all(_matches(motif, key, values) for key, values in wanted.items())
        ]
        return MotifSet(kept)

    # ------------------------------------------------------------------- scanning

    @overload
    def scan(
        self,
        sequence: str,
        name: str = ...,
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: None = ...,
        workers: int | None = ...,
    ) -> pd.DataFrame: ...

    @overload
    def scan(
        self,
        sequence: str,
        name: str = ...,
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: str | Path,
        workers: int | None = ...,
    ) -> Path: ...

    def scan(
        self,
        sequence: str,
        name: str = DEFAULT_SEQUENCE_NAME,
        *,
        threshold: float = DEFAULT_THRESHOLD,
        background: BackgroundArg = None,
        output: str | Path | None = None,
        workers: int | None = DEFAULT_WORKERS,
    ) -> pd.DataFrame | Path:
        """Scan one sequence with every motif here and return the **Hit table**.

        The quick case — one locus, checked once. It answers with exactly the table
        :meth:`scan_sequences` and :meth:`scan_fasta` answer with, down to the column
        order and the dtypes, so nothing downstream branches on how the scan was called.

        Both strands are scanned. Coordinates are **0-based half-open** and always in the
        forward frame, whichever strand matched, and **Strand** is ``+`` or ``-`` — never
        ``.``, because a scan knows which of the two it scored. The interval covers the
        bases the matrix scored, so a trimmed motif's hit is as long as the trimmed motif;
        :attr:`Motif.offset` maps a position *within* the motif and not the interval.

        **The sequence is upper-cased**, so a soft-masked one yields exactly the hits its
        upper-case equivalent does and there is no argument that would change that
        (ADR-0012).

        Parameters
        ----------
        sequence : str
            The bases to scan. A :class:`~genome.seq.DNA` or a plain string.
        name : str, default ``"sequence"``
            What the ``sequence_name`` column carries. A fixed literal by default, so
            every row names its sequence even when the caller did not.
        threshold : float, default 1e-4
            The **Threshold**: one per-position p-value, converted per motif against
            ``background`` into the score that motif must clear. In ``(0, 1)``.
        background : sequence of float or {"auto", "uniform", "derive"}, optional
            The **Background**: four frequencies over :data:`BASES`, above zero and
            summing to 1, or one of the three modes. **Automatic when omitted** — derived
            from ``sequence`` when it holds at least
            :data:`~genome.tf.motif.background.BACKGROUND_FLOOR` unambiguous bases,
            uniform below that. ``"uniform"`` pins it; ``"derive"`` derives whatever the
            input holds, floor or no floor. Whichever it is, it is recorded on the result.
        output : str or pathlib.Path, optional
            Where to stream the hits as Parquet instead of building a table. **A scan too
            large to hold goes to disk and hands back the path**; there is no row-count
            guard, because a genome-scale scan is the caller's decision. Read it back with
            :func:`~genome.tf.motif.parquet.read_hits`, which restores the dtypes and the
            provenance both — :func:`pandas.read_parquet` alone drops the provenance.
        workers : int, default 1
            How many processes to shard the scan across. **One by default**, so importing
            this package and calling a scan never starts a process unasked: under the spawn
            start method a pool re-imports the caller's script, and an unguarded one would
            re-execute itself. ``None`` resolves the count with
            :func:`~genome.tf.motif.workers.resolve_workers` — the Slurm allocation first,
            then process affinity, then the machine — which is what the command line
            passes. **More than one produces the identical table**, row for row; the choice
            is about wall time and nothing else.

        Returns
        -------
        pandas.DataFrame or pathlib.Path
            One row per **Motif hit**: ``motif_id``, ``motif_name``, ``sequence_name``,
            ``start``, ``end``, ``strand``, ``score`` — the score in bits, not a p-value.
            The scan's provenance is on ``frame.attrs``: the background, the threshold, the
            **Release** and **Tax group** where the set knows them, and which motifs were
            scanned and which were skipped. With ``output``, the path written.

        Raises
        ------
        ValueError
            If ``threshold`` is not in ``(0, 1)``, or ``background`` is neither one of the
            three modes nor four positive frequencies summing to 1.

        Notes
        -----
        Motifs shorter than :data:`MIN_MOTIF_LENGTH` are **not scanned** and are named in
        ``frame.attrs["motifs_skipped"]``. A 6-mer cannot reach the default threshold at
        all, and an engine asked for it anyway would fall back to that matrix's best
        attainable cutoff and over-call in silence.

        A **Release** answers with its release and tax group on the table; a set that
        :meth:`filter` returned answers with ``None`` for both, since a filtered release is
        no longer that release and must not claim to be.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.zeros((4, 8))
        >>> for column, base in enumerate("GATTACAG"):
        ...     counts["ACGT".index(base), column] = 100.0
        >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
        >>> hits = motifs.scan("TTTTTGATTACAGTTTTT")
        >>> hits[["motif_id", "sequence_name", "start", "end", "strand"]]
           motif_id sequence_name  start  end strand
        0  MA9999.1      sequence      5   13      +
        >>> hits.attrs["motifs_scanned"], hits.attrs["threshold"]
        (('MA9999.1',), 0.0001)
        """
        return self.scan_sequences(
            {name: sequence},
            threshold=threshold,
            background=background,
            output=output,
            workers=workers,
        )

    @overload
    def scan_sequences(
        self,
        sequences: Mapping[str, str],
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: None = ...,
        workers: int | None = ...,
    ) -> pd.DataFrame: ...

    @overload
    def scan_sequences(
        self,
        sequences: Mapping[str, str],
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: str | Path,
        workers: int | None = ...,
    ) -> Path: ...

    def scan_sequences(
        self,
        sequences: Mapping[str, str],
        *,
        threshold: float = DEFAULT_THRESHOLD,
        background: BackgroundArg = None,
        output: str | Path | None = None,
        workers: int | None = DEFAULT_WORKERS,
    ) -> pd.DataFrame | Path:
        """Scan named sequences — a peak set in one call — and return the **Hit table**.

        The same table :meth:`scan` returns, with ``sequence_name`` carrying the mapping's
        own keys and the rows in the mapping's own order. See :meth:`scan` for the
        coordinate convention, the strand rule, the upper-casing and the provenance.

        Parameters
        ----------
        sequences : mapping of str to str
            Name to bases. Scanned one at a time, so the peak cost is the longest sequence
            rather than all of them.
        threshold : float, default 1e-4
            The **Threshold**, as a per-position p-value.
        background : sequence of float or {"auto", "uniform", "derive"}, optional
            The **Background**. Automatic when omitted — see :meth:`scan`. Deciding reads
            a bounded prefix of ``sequences``, which is then scanned like the rest.
        output : str or pathlib.Path, optional
            Where to stream the hits as Parquet instead of building a table — see
            :meth:`scan`.
        workers : int, default 1
            How many processes to shard the scan across — see :meth:`scan`. Sequences are
            distributed whole, so a peak set parallelises without any of them being cut.

        Returns
        -------
        pandas.DataFrame or pathlib.Path
            The **Hit table**, empty of rows but not of schema when nothing matched. With
            ``output``, the path written.

        Raises
        ------
        ValueError
            If ``threshold`` is not in ``(0, 1)``, or the background is neither a mode nor
            four positive frequencies summing to 1.

        Examples
        --------
        >>> import numpy as np
        >>> counts = np.zeros((4, 8))
        >>> for column, base in enumerate("GATTACAG"):
        ...     counts["ACGT".index(base), column] = 100.0
        >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
        >>> peaks = {"peak1": "TTTTT" + "GATTACAG" + "TTTTT",
        ...          "peak2": "CCCCC" + "CTGTAATC" + "CCCCC"}   # the same site, flipped
        >>> hits = motifs.scan_sequences(peaks)
        >>> list(zip(hits["sequence_name"], hits["strand"], hits["start"], hits["end"]))
        [('peak1', '+', 5, 13), ('peak2', '-', 5, 13)]
        """
        # Imported here rather than at module scope: the engine adapter imports this
        # module, so naming it up there would be a cycle — and this stays the pure core.
        from genome.tf.motif.scan import scan_stream

        return scan_stream(
            self,
            sequences.items(),
            threshold=threshold,
            background=background,
            output=output,
            workers=workers,
        )

    @overload
    def scan_fasta(
        self,
        path: str | Path,
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: None = ...,
        workers: int | None = ...,
    ) -> pd.DataFrame: ...

    @overload
    def scan_fasta(
        self,
        path: str | Path,
        *,
        threshold: float = ...,
        background: BackgroundArg = ...,
        output: str | Path,
        workers: int | None = ...,
    ) -> Path: ...

    def scan_fasta(
        self,
        path: str | Path,
        *,
        threshold: float = DEFAULT_THRESHOLD,
        background: BackgroundArg = None,
        output: str | Path | None = None,
        workers: int | None = DEFAULT_WORKERS,
    ) -> pd.DataFrame | Path:
        r"""Scan every record of a FASTA and return the **Hit table**.

        The same table :meth:`scan` returns. Records are read and scanned one at a time,
        so the file is never held whole; plain or gzipped.

        **A record's name is its header up to the first whitespace**, which is what STAR
        and chromap write into an alignment produced from the same file — so this table
        joins against that alignment with nobody renaming anything.

        Parameters
        ----------
        path : str or pathlib.Path
            The FASTA to scan, ``.fa`` or ``.fa.gz``.
        threshold : float, default 1e-4
            The **Threshold**, as a per-position p-value.
        background : sequence of float or {"auto", "uniform", "derive"}, optional
            The **Background**. Automatic when omitted — see :meth:`scan`. Deciding reads
            records off the front of the file, which are then scanned like the rest, so
            the file is still read exactly once.
        output : str or pathlib.Path, optional
            Where to stream the hits as Parquet instead of building a table — see
            :meth:`scan`. The whole-genome case this exists for: a FASTA in, a Parquet
            out, and nothing the size of either in memory.
        workers : int, default 1
            How many processes to shard the scan across — see :meth:`scan`. A record long
            enough is cut into pieces with an overlap, so one chromosome still uses the
            whole allocation.

        Returns
        -------
        pandas.DataFrame or pathlib.Path
            The **Hit table**, with ``sequence_name`` carrying the truncated record names.
            With ``output``, the path written.

        Raises
        ------
        FileNotFoundError
            If ``path`` does not exist.
        ValueError
            If the file is not FASTA, a record carries no name, ``threshold`` is not in
            ``(0, 1)``, or the background is neither a mode nor four positive frequencies
            summing to 1.

        Examples
        --------
        >>> import tempfile
        >>> from pathlib import Path
        >>> import numpy as np
        >>> counts = np.zeros((4, 8))
        >>> for column, base in enumerate("GATTACAG"):
        ...     counts["ACGT".index(base), column] = 100.0
        >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     fasta = Path(directory) / "peaks.fa"
        ...     _ = fasta.write_text(">peak1 chrI:100-118 of nowhere\nTTTTTGATTACAGTTTTT\n")
        ...     hits = motifs.scan_fasta(fasta)
        >>> list(zip(hits["sequence_name"], hits["start"], hits["strand"]))
        [('peak1', 5, '+')]
        """
        from genome.tf.motif.scan import read_fasta, scan_stream

        return scan_stream(
            self,
            read_fasta(path),
            threshold=threshold,
            background=background,
            output=output,
            workers=workers,
        )

    # ------------------------------------------------------------------ comparison

    def compare(
        self, queries: Motif | Iterable[Motif], *, top: int | None = None
    ) -> MotifComparison:
        """Ask what one or more motifs look like, against the motifs held here.

        **The motifs held here are the targets, and the argument is the queries** — read
        it as *compare these against this release*. The use case is naming: a chromBPNet
        or TF-MoDISco run hands back matrices with no names on them, and
        ``release.compare(de_novo)`` says which published motif each one most resembles.

        The comparison is tomtom's, from `memelite`, handed the same 4 x L probability
        matrices :attr:`Motif.probabilities` produces. What comes back is a labelled
        array indexed by **Motif id** on both axes; see :class:`MotifComparison` for its
        two shapes and :meth:`MotifComparison.to_frame` for the flat table.

        Parameters
        ----------
        queries : Motif or iterable of Motif
            One motif, several, or a whole :class:`MotifSet`. Two queries sharing a
            **Motif id** are refused: the array's query axis is labelled with them, so a
            repeated label could answer for neither.
        top : int, optional
            Keep only this many targets per query, best first. It is *not* a convenience
            over the complete answer — it sends the work down tomtom's faster
            nearest-neighbour path, which never scores the targets that lose. The result
            is then **ragged**: its target axis is per query, and it **cannot be widened
            without recomputing**, which is accepted rather than a defect. Omit it for
            the complete query x target array.

        Returns
        -------
        MotifComparison
            The labelled array and the methods that read it.

        Raises
        ------
        ValueError
            If ``queries`` is empty, two queries share a **Motif id**, this set holds no
            motifs, or ``top`` is below 1 or above the number of motifs held here.

        Notes
        -----
        **A motif compared against itself aligns to itself perfectly** — offset 0, the
        whole length overlapping, on the ``+`` strand, and no target scores higher. It is
        usually ranked first too, and the exception is worth knowing: TOMTOM's p-value
        rewards a short dense alignment, so a long motif that embeds a shorter one can
        rank the shorter one above itself. Both of the fixture's 31- and 33-column CTCF
        matrices do exactly that with the 15-column ``MA0139.2`` they contain. That is a
        property of the statistic, not of this wrapper, and it is what a caller naming a
        de novo pattern should expect to see from a family of nested matrices.

        Examples
        --------
        >>> import numpy as np
        >>> def spelled(motif_id, bases):
        ...     columns = [[19.0 if b == l else 1.0 for l in bases] for b in BASES]
        ...     return Motif(motif_id, "", np.array(columns))
        >>> published = MotifSet([spelled("MA0001.1", "ACGTACGTA"),
        ...                       spelled("MA0002.1", "TTTTTTTTT")])
        >>> published.compare(spelled("pattern_0", "ACGTACGTA")).to_frame()["target"]
        0    MA0001.1
        Name: target, dtype: object
        >>> published.compare(published, top=1).is_ragged
        True
        """
        # Imported here, not at module scope: the comparison engine pulls in numba, and
        # `import genome` must not pay for it — the same reason plot() defers logomaker.
        from genome.tf.motif.compare import _compare

        return _compare(queries, self, top=top)


def base_id(motif: Motif) -> str:
    """Return a **Motif id** without its version — ``MA0139.2`` becomes ``MA0139``.

    Examples
    --------
    >>> import numpy as np
    >>> base_id(Motif("MA0139.2", "CTCF", np.ones((4, 7))))
    'MA0139'
    """
    return motif.motif_id.split(".", 1)[0]


def _group(motifs: tuple[Motif, ...], key: Callable[[Motif], str]) -> dict[str, tuple[Motif, ...]]:
    """Index motifs by a key several of them may share, keeping the order they came in."""
    grouped: dict[str, list[Motif]] = {}
    for motif in motifs:
        grouped.setdefault(key(motif), []).append(motif)
    return {key_: tuple(group) for key_, group in grouped.items()}


def _as_wanted(value: str | Iterable[str]) -> tuple[str, ...]:
    """Read one filter keyword's value as the alternatives it offers."""
    return (value,) if isinstance(value, str) else tuple(value)


def _matches(motif: Motif, annotation: str, wanted: tuple[str, ...]) -> bool:
    """Answer whether one motif's annotation matches any of the alternatives asked for."""
    held = getattr(motif, annotation)
    values = (held,) if isinstance(held, str) else held
    if _FILTERABLE[annotation]:
        return any(want.lower() in value.lower() for value in values for want in wanted)
    return any(value in wanted for value in values)


def _as_value_tuple(values: Iterable[str], field: str) -> tuple[str, ...]:
    """Freeze one plural annotation into a tuple, rejecting a bare string."""
    if isinstance(values, str):
        raise ValueError(
            f"{field} must be an iterable of values, not the string {values!r}: pass "
            f"('{values}',) for one, since a bare string would be read letter by letter. "
            f"{field} is plural because the source publishes several for a dimer."
        )
    return tuple(values)


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
