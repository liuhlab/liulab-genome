"""Comparing motifs — asking what a **Motif** looks like, in a **Release**'s own words.

The question this answers is *naming*: a chromBPNet or TF-MoDISco run hands back matrices
with no names on them, and the way to name one is to ask which published motif it most
resembles. :meth:`~genome.tf.motif.motif.MotifSet.compare` takes one **Motif**, several,
or a whole **Motif set** as the queries and compares them against the motifs of the set it
is called on.

**The comparison itself is tomtom's**, from `memelite`, which takes probability matrices
in the same 4 x L layout :attr:`~genome.tf.motif.motif.Motif.probabilities` already
produces — so nothing here permutes, transposes or rescales anything on the way in. This
is the *only* place `memelite` is used: its scanner was measured against MOODS and
rejected, and that measurement is in ``docs/research/``.

**What comes back is a labelled array**, not a table — :class:`MotifComparison` wraps an
:class:`xarray.Dataset` indexed by **Motif id** on both axes, which is what makes
``data.sel(query="pattern_0", target="MA0139.2")`` the natural way to ask about one pair
and ``data["neg_log10_p"]`` a similarity matrix ready to cluster. :meth:`
MotifComparison.to_frame` flattens it to one row per pair for the common question.

**Negative log10 p is stored, never raw p.** The array holds half precision, whose
smallest normal value is 6.1e-5 — a p of 1e-20 stored raw would flush to zero and take
every motif's best match down with it. Stored as 20.0 it is an ordinary small number, and
that inversion is the whole reason the column is named for the transform rather than for
the p-value.

Examples
--------
>>> import numpy as np
>>> from genome.tf.motif import Motif, MotifSet
>>> def spelled(motif_id, bases):
...     columns = [[19.0 if base == letter else 1.0 for letter in bases] for base in "ACGT"]
...     return Motif(motif_id, "", np.array(columns))
>>> published = MotifSet([spelled("MA0001.1", "ACGTACGTA"), spelled("MA0002.1", "TTTTTTTTT")])
>>> de_novo = spelled("pattern_0", "ACGTACGTA")
>>> named = published.compare(de_novo).to_frame()
>>> named.loc[0, ["query", "target", "strand"]].tolist()
['pattern_0', 'MA0001.1', '+']
>>> int(named.loc[0, "overlap"])                  # the whole 9 columns aligned
9
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:  # pragma: no cover - typing only, so importing genome stays cheap
    import pandas as pd
    import xarray as xr

    from genome.tf.motif.motif import Motif, MotifSet

#: The five things a comparison says about one query-target pair, and the data variables
#: of the labelled array, in order. Fixed: nothing downstream should have to branch on
#: which of them a particular comparison happened to carry.
COMPARISON_VARIABLES: tuple[str, ...] = (
    "neg_log10_p",
    "score",
    "offset",
    "overlap",
    "strand",
)

#: The columns :meth:`MotifComparison.to_frame` produces, in order. ``query`` and
#: ``target`` are **Motif id**s and ``rank`` is 0 for the best target of each query.
FRAME_COLUMNS: tuple[str, ...] = ("query", "target", "rank", *COMPARISON_VARIABLES)

#: How each variable is stored. Half precision for the p-value is what makes storing its
#: negative log10 necessary rather than merely tidy; the score keeps single precision
#: because it is the tiebreaker when two targets round to the same p, and a tiebreaker
#: that ties is not one. Offsets and overlaps are bounded by motif lengths, which are
#: tens of columns, so a 16-bit integer is not a squeeze.
_DTYPES: dict[str, npt.DTypeLike] = {
    "neg_log10_p": np.float16,
    "score": np.float32,
    "offset": np.int16,
    "overlap": np.int16,
}

#: tomtom reports the strand as 0 for the target as given and 1 for its reverse
#: complement. A comparison always knows which of the two it aligned, so the **Strand**
#: here is always ``+`` or ``-`` and never ``.``.
_STRAND_SYMBOLS = np.array(["+", "-"], dtype="<U1")


class RaggedComparisonError(ValueError):
    """A limited comparison was asked about targets it never scored.

    Raised by :attr:`MotifComparison.target_ids` and by
    :meth:`MotifComparison.to_frame` when a comparison run with a ``top`` limit is asked
    for a shared target axis or for more targets per query than it kept.

    **A limited comparison cannot be widened without recomputing, and that is accepted
    rather than a defect.** Passing ``top`` sends the work down tomtom's nearest-neighbour
    path, which never scores the targets that lose — so the missing cells were not
    discarded, they were never computed, and no amount of rearranging the array will
    produce them. Recompute with a larger ``top``, or with none at all for the complete
    query x target array.

    Examples
    --------
    >>> raise RaggedComparisonError("kept 3 targets per query")
    Traceback (most recent call last):
        ...
    genome.tf.motif.compare.RaggedComparisonError: kept 3 targets per query
    """


class MotifComparison:
    """How a set of query **Motif**s compares against a set of target motifs.

    A thin wrapper over an :class:`xarray.Dataset` — the array is the answer, and this
    class exists to say which of its two shapes it is in, to keep the ranking rule in one
    place, and to refuse the one question a limited comparison cannot answer. Reach for
    :attr:`data` whenever you want the array itself.

    **Two shapes, and every method here handles both.**

    - *Complete*, from ``compare(queries)``: dimensions ``(query, target)``, both labelled
      with **Motif id**s, every pair scored. :attr:`target_ids` answers, and
      ``data.sel(query=..., target=...)`` reaches any cell.
    - *Limited*, from ``compare(queries, top=n)``: dimensions ``(query, rank)``, where
      rank 0 is each query's best target. **The target axis is per query** — rank 0 names
      a different motif for each row — so the target ids ride along as a ``target`` data
      variable rather than as a coordinate, and there is no shared axis to index or to
      widen.

    Parameters
    ----------
    data : xarray.Dataset
        The labelled array, in either shape. Must carry a ``query`` dimension labelled
        with **Motif id**s, every variable in :data:`COMPARISON_VARIABLES`, and either a
        ``target`` coordinate (complete) or a ``rank`` dimension with a ``target``
        variable over ``(query, rank)`` (limited). ``attrs["targets_compared"]`` records
        how many targets the comparison ran against, which is the only thing a limited
        array cannot say for itself.

    Raises
    ------
    ValueError
        If ``data`` is in neither shape. The message names what is missing.

    Examples
    --------
    >>> import numpy as np
    >>> from genome.tf.motif import Motif, MotifSet
    >>> counts = np.array([[19.0, 1.0, 1.0], [1.0, 19.0, 1.0],
    ...                    [1.0, 1.0, 19.0], [1.0, 1.0, 1.0]])
    >>> published = MotifSet([Motif("MA0001.1", "", np.tile(counts, 3))])
    >>> comparison = published.compare(Motif("pattern_0", "", np.tile(counts, 3)))
    >>> comparison
    MotifComparison(queries=1, targets=1, top=None)
    >>> comparison.is_ragged
    False
    >>> float(comparison.data["neg_log10_p"].sel(query="pattern_0", target="MA0001.1")) > 0
    True
    """

    def __init__(self, data: xr.Dataset) -> None:
        self._data = _checked(data)

    # ----------------------------------------------------------------- what is held

    @property
    def data(self) -> xr.Dataset:
        """The labelled array itself, indexable by **Motif id**.

        Examples
        --------
        >>> import numpy as np
        >>> from genome.tf.motif import Motif, MotifSet
        >>> published = MotifSet([Motif("MA0001.1", "", np.eye(4)[:, [0, 1, 2, 3, 0, 1, 2]] + 1)])
        >>> sorted(published.compare(published).data.data_vars)
        ['neg_log10_p', 'offset', 'overlap', 'score', 'strand']
        """
        return self._data

    @property
    def is_ragged(self) -> bool:
        """Whether the target axis is per query, which a ``top`` limit makes it.

        Examples
        --------
        >>> import numpy as np
        >>> from genome.tf.motif import Motif, MotifSet
        >>> published = MotifSet([Motif("MA0001.1", "", np.ones((4, 8)) + np.eye(4, 8))])
        >>> published.compare(published).is_ragged
        False
        >>> published.compare(published, top=1).is_ragged
        True
        """
        return "rank" in self._data.dims

    @property
    def top(self) -> int | None:
        """The limit this comparison was run with — ``None`` when every pair was scored."""
        return int(self._data.sizes["rank"]) if self.is_ragged else None

    @property
    def query_ids(self) -> tuple[str, ...]:
        """Every query **Motif id**, in the order the queries were given."""
        return tuple(str(value) for value in self._data["query"].to_numpy())

    @property
    def target_ids(self) -> tuple[str, ...]:
        """Every target **Motif id**, in the order the target set held them.

        Returns
        -------
        tuple of str
            The shared target axis.

        Raises
        ------
        RaggedComparisonError
            If this comparison was limited. A limited comparison has no shared target
            axis: rank 0 names a different motif for each query.
        """
        if self.is_ragged:
            raise RaggedComparisonError(
                f"this comparison was limited to top={self.top}, so its target axis is per "
                f"query and there is no shared one to return: rank 0 names a different motif "
                f"for each of the {len(self.query_ids)} queries. Read "
                f"comparison.data['target'] for the ids it kept, or recompute with "
                f"compare(..., top=None) for a complete query x target array."
            )
        return tuple(str(value) for value in self._data["target"].to_numpy())

    @property
    def targets_compared(self) -> int:
        """How many target motifs the comparison ran against, whichever shape it is in."""
        return int(self._data.attrs["targets_compared"])

    def __repr__(self) -> str:
        """Return the two sizes and the limit, never the array itself."""
        return (
            f"{type(self).__name__}(queries={len(self.query_ids)}, "
            f"targets={self.targets_compared}, top={self.top})"
        )

    # -------------------------------------------------------------------- flattening

    def to_frame(self, top: int | None = 1) -> pd.DataFrame:
        """Flatten to one row per query-target pair, best first within each query.

        The default is the one question most callers have — *what does each of my motifs
        look like* — so it hands back one row per query. Ranking is by
        ``neg_log10_p`` descending, ties broken by ``score`` descending and then by the
        target set's own order, which is also the order tomtom's nearest-neighbour path
        returns, so the two shapes agree on which target is best.

        Parameters
        ----------
        top : int or None, default 1
            How many targets to keep per query, best first. ``None`` keeps every pair the
            comparison holds. A limit above the number of targets is not an error for a
            complete comparison — nothing is missing from one, so it simply returns every
            pair.

        Returns
        -------
        pandas.DataFrame
            Columns :data:`FRAME_COLUMNS`, with a fresh range index. Numeric dtypes are
            the array's own, unchanged: the frame is the array flattened, not a copy of
            it at a different precision.

        Raises
        ------
        ValueError
            If ``top`` is below 1.
        RaggedComparisonError
            If this comparison was limited and ``top`` asks for more targets per query
            than it kept. Those pairs were never scored, so widening means recomputing.

        Examples
        --------
        >>> import numpy as np
        >>> from genome.tf.motif import Motif, MotifSet
        >>> def spelled(motif_id, bases):
        ...     columns = [[19.0 if b == l else 1.0 for l in bases] for b in "ACGT"]
        ...     return Motif(motif_id, "", np.array(columns))
        >>> published = MotifSet([spelled("MA0001.1", "ACGTACGTA"),
        ...                       spelled("MA0002.1", "TTTTTTTTT")])
        >>> published.compare(spelled("pattern_0", "ACGTACGTA")).to_frame()["target"]
        0    MA0001.1
        Name: target, dtype: object
        """
        import pandas as pd

        if top is not None and top < 1:
            raise ValueError(
                f"top must be >= 1 or None, got {top}: it is how many targets to keep per "
                f"query. Pass None for every pair."
            )
        available = int(self._data.sizes["rank" if self.is_ragged else "target"])
        if top is not None and top > available and available < self.targets_compared:
            raise RaggedComparisonError(
                f"to_frame(top={top}) asks for {top} targets per query and this comparison "
                f"kept {available}: the other "
                f"{self.targets_compared - available} of its {self.targets_compared} targets "
                f"were never scored, because a top limit takes tomtom's nearest-neighbour "
                f"path. A limited comparison cannot be widened without recomputing. Ask for "
                f"at most {available}, or recompute with compare(..., top={top})."
            )
        keep = available if top is None else min(top, available)

        order = self._order()[:, :keep]
        targets = np.take_along_axis(self._targets(), order, axis=1)
        columns: dict[str, npt.NDArray[np.generic]] = {
            "query": np.repeat(np.asarray(self.query_ids), keep),
            "target": targets.ravel(),
            "rank": np.tile(np.arange(keep, dtype=np.int16), len(self.query_ids)),
        }
        for name in COMPARISON_VARIABLES:
            values = self._data[name].to_numpy()
            columns[name] = np.take_along_axis(values, order, axis=1).ravel()
        # Inserted in FRAME_COLUMNS order above, so the frame is already in it — asserted
        # by a test rather than restated here, where it would only be a second source.
        return pd.DataFrame(columns)

    # ----------------------------------------------------------------------- internals

    def _order(self) -> npt.NDArray[np.intp]:
        """Return, per query, the column order that puts the best target first."""
        if self.is_ragged:
            # tomtom's nearest-neighbour path already sorted these by p-value, and rank
            # is what that order was stored as.
            ranks = np.arange(self._data.sizes["rank"], dtype=np.intp)
            return np.tile(ranks, (self._data.sizes["query"], 1))
        strength = self._data["neg_log10_p"].to_numpy()
        score = self._data["score"].to_numpy()
        # Ties fall back to the target set's own order, so a comparison of one set
        # against another is reproducible rather than left to the sort's internals.
        given = np.tile(np.arange(strength.shape[1]), (strength.shape[0], 1))
        return np.lexsort((given, -score, -strength.astype(np.float32)), axis=1)

    def _targets(self) -> npt.NDArray[np.str_]:
        """Return the target **Motif id** of every cell, whichever shape this is in."""
        if self.is_ragged:
            return np.asarray(self._data["target"].to_numpy(), dtype=np.str_)
        return np.tile(np.asarray(self.target_ids), (self._data.sizes["query"], 1))


def _checked(data: xr.Dataset) -> xr.Dataset:
    """Refuse a dataset that is in neither of the two shapes, naming what is missing."""
    if "query" not in data.dims:
        raise ValueError(
            "a comparison is a labelled array with a 'query' dimension holding motif ids, "
            f"and this one has {tuple(data.dims)}. Build it with MotifSet.compare() rather "
            "than by hand unless you are reproducing its layout exactly."
        )
    missing = [name for name in COMPARISON_VARIABLES if name not in data.data_vars]
    if missing:
        raise ValueError(
            f"a comparison carries {', '.join(COMPARISON_VARIABLES)} and this one is missing "
            f"{', '.join(missing)}. Build it with MotifSet.compare()."
        )
    complete = "target" in data.coords
    limited = "rank" in data.dims and "target" in data.data_vars
    if not (complete or limited):
        raise ValueError(
            "a comparison is either complete — dimensions (query, target), both labelled "
            "with motif ids — or limited — dimensions (query, rank) with a 'target' variable "
            f"naming a motif per cell — and this one is neither, with {tuple(data.dims)}. "
            "Build it with MotifSet.compare()."
        )
    if "targets_compared" not in data.attrs:
        raise ValueError(
            "a comparison records how many targets it ran against in "
            "attrs['targets_compared'], and this one does not. A limited comparison cannot "
            "say it any other way, which is how to_frame() knows a widening from a request "
            "for everything."
        )
    return data


def _neg_log10(p_values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Return ``-log10(p)``, infinite where a p-value underflowed to zero."""
    # Masked rather than errstate'd: log10(0) is never evaluated, so no warning is raised
    # to be suppressed — and this suite turns warnings into errors.
    strength = np.full(p_values.shape, np.inf, dtype=np.float64)
    positive = p_values > 0.0
    strength[positive] = -np.log10(p_values[positive])
    return strength


def _compare(
    queries: Motif | Iterable[Motif], targets: MotifSet, *, top: int | None
) -> MotifComparison:
    """Run tomtom over the queries against the targets and label what it returns."""
    # Imported here, not at module scope: memelite pulls in numba, and `import genome`
    # must not pay half a second for a comparison nobody asked for.
    import xarray as xr
    from memelite import tomtom

    from genome.tf.motif.motif import Motif, MotifSet

    given = (queries,) if isinstance(queries, Motif) else tuple(queries)
    if not given:
        raise ValueError(
            "no motifs to compare: compare() takes one Motif, several, or a MotifSet as its "
            "queries, and was given none. A filter that matched nothing is a real answer — "
            "check what it kept before comparing."
        )
    # Built rather than iterated, because a MotifSet is what refuses two queries sharing an
    # id — and a labelled array whose query axis repeats a label could answer for neither.
    query_set = MotifSet(given)
    if not len(targets):
        raise ValueError(
            "this motif set holds no motifs, so there is nothing to compare against. Compare "
            "against a JasparDatabase, or against a set that kept something."
        )
    n_targets = len(targets)
    if top is not None:
        if top < 1:
            raise ValueError(
                f"top must be >= 1 or None, got {top}: it is how many targets to keep per "
                f"query. Pass None to score every pair."
            )
        if top > n_targets:
            raise ValueError(
                f"top={top} asks for more targets per query than there are: this set holds "
                f"only {n_targets} targets. Ask for at most {n_targets}, or pass top=None to "
                f"score every pair."
            )

    query_ids = list(query_set.motif_ids)
    target_ids = list(targets.motif_ids)
    returned = tomtom(
        [motif.probabilities for motif in query_set],
        [motif.probabilities for motif in targets],
        n_nearest=top,
    )
    p_values, score, offset, overlap, strand = returned[:5]

    dims = ("query", "rank") if top is not None else ("query", "target")
    variables: dict[str, tuple[tuple[str, str], npt.NDArray[np.generic]]] = {
        "neg_log10_p": (dims, _neg_log10(np.asarray(p_values)).astype(_DTYPES["neg_log10_p"])),
        "score": (dims, np.asarray(score).astype(_DTYPES["score"])),
        "offset": (dims, np.rint(offset).astype(_DTYPES["offset"])),
        "overlap": (dims, np.rint(overlap).astype(_DTYPES["overlap"])),
        "strand": (dims, _STRAND_SYMBOLS[np.rint(strand).astype(np.intp)]),
    }
    coords: dict[str, object] = {"query": query_ids}
    if top is None:
        coords["target"] = target_ids
    else:
        # The per-query target axis. tomtom hands back indices into the target order,
        # sorted by p-value, so rank 0 is each query's best and the names ride along as a
        # variable — there is no shared axis to make a coordinate out of.
        chosen = np.rint(returned[5]).astype(np.intp)
        coords["rank"] = np.arange(top, dtype=np.int16)
        variables["target"] = (dims, np.asarray(target_ids)[chosen])
    return MotifComparison(
        xr.Dataset(variables, coords=coords, attrs={"targets_compared": n_targets})
    )
