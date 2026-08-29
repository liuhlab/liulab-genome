"""The Parquet sink a **Hit table** too large to hold streams to, and the reader for it.

Scanning hg38 with a full vertebrate release is roughly 550 million **Motif hit**s; at the
19 bytes a row the fixed dtypes cost, that is 10.5 GB, which is not a DataFrame on any
machine in the lab. So a scan handed an output path writes there instead and answers with
the path. **There is no row-count guard and no refusal** — a genome-scale scan is the
caller's decision, and a package that second-guessed it would only be guessing.

**Batches are written as they are produced.** :func:`write_hits` drains the same
per-sequence iterator :func:`~genome.tf.motif.scan._collect` drains, one row group per
batch, so what is held at once is one sequence's hits and never the whole answer. That is
the seam the scan loop was built for: the sink replaces the collector and nothing else
moves.

**The provenance travels in the file, because ``frame.attrs`` does not.** pandas drops
``attrs`` on the way through pyarrow, and a **Hit table** without its **Background**,
**Threshold**, **Release** and **Tax group** is the same table with its meaning removed. So
it is written into the Parquet file's own key-value metadata under
:data:`HIT_PROVENANCE_KEY` and put back on ``frame.attrs`` by :func:`read_hits`.

**The dtypes are written explicitly rather than inferred per batch.** A categorical
column's index width follows its cardinality, so two batches would otherwise disagree on
the schema — ``int8`` for the batch with four sequence names and ``int16`` for the one with
four hundred. The schema here pins ``int32`` for all of them, and :func:`read_hits` sorts
each column's categories back into the order ``astype("category")`` would have produced, so
what comes off the disk equals the in-memory table down to its category order.

Examples
--------
>>> import tempfile
>>> from pathlib import Path
>>> import numpy as np
>>> from genome.tf.motif import Motif, MotifSet
>>> counts = np.zeros((4, 8))
>>> for column, base in enumerate("GATTACAG"):
...     counts["ACGT".index(base), column] = 100.0
>>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
>>> with tempfile.TemporaryDirectory() as directory:
...     written = motifs.scan("TTTTTGATTACAGTTTTT", "chrTest",
...                           output=Path(directory) / "hits.parquet")
...     hits = read_hits(written)
>>> hits[["motif_id", "sequence_name", "start", "end", "strand"]]
   motif_id sequence_name  start  end strand
0  MA9999.1       chrTest      5   13      +
>>> hits.attrs["background"], hits.attrs["motifs_scanned"]
((0.25, 0.25, 0.25, 0.25), ('MA9999.1',))
"""

from __future__ import annotations

import json
from collections.abc import Hashable, Iterable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from genome.tf.motif.scan import HIT_DTYPES, empty_hits

#: Where the scan's provenance lives in the Parquet file's key-value metadata — the six
#: things :data:`~genome.tf.motif.scan.HIT_PROVENANCE` names, as JSON. Its own key rather
#: than pandas' ``attributes`` block, which pandas neither writes nor reads back.
HIT_PROVENANCE_KEY = b"genome.hit_provenance"

#: How the **Hit table**'s dtypes are spelled in Arrow. Built from
#: :func:`~genome.tf.motif.scan.empty_hits` so that pandas' own column metadata comes along
#: — that is what makes a category read back as a category — then widened to ``int32``
#: dictionary indices, which an empty frame would otherwise leave at ``int8`` and a large
#: batch would overflow.
_SCHEMA: pa.Schema = pa.schema(
    [
        field.with_type(pa.dictionary(pa.int32(), pa.string()))
        if pa.types.is_dictionary(field.type)
        else field
        for field in pa.Table.from_pandas(empty_hits(), preserve_index=False).schema
    ],
    metadata=pa.Table.from_pandas(empty_hits(), preserve_index=False).schema.metadata,
)


def write_hits(
    batches: Iterable[pd.DataFrame],
    output: str | Path,
    provenance: Mapping[Hashable, Any],
) -> Path:
    """Stream **Hit table** batches to Parquet and return where they went.

    Each batch is written as it arrives, so the whole result is never held. An empty batch
    contributes no row group; a scan that found nothing still writes a file, with the
    schema and the provenance and no rows, because that is a real answer.

    Parameters
    ----------
    batches : iterable of pandas.DataFrame
        One frame per named sequence, each with :data:`~genome.tf.motif.scan.HIT_COLUMNS`
        and :data:`~genome.tf.motif.scan.HIT_DTYPES`.
    output : str or pathlib.Path
        Where to write. The parent directory is created if it does not exist.
    provenance : mapping
        What :data:`~genome.tf.motif.scan.HIT_PROVENANCE` names, as
        :func:`read_hits` will put it back on ``frame.attrs``.

    Returns
    -------
    pathlib.Path
        ``output``, as a path.

    Raises
    ------
    IsADirectoryError
        If ``output`` names a directory that already exists.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> from genome.tf.motif.scan import empty_hits
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     written = write_hits([empty_hits()], Path(directory) / "none.parquet",
    ...                          {"threshold": 0.0001})
    ...     read_hits(written).attrs
    {'threshold': 0.0001}
    """
    path = Path(output)
    if path.is_dir():
        raise IsADirectoryError(
            f"cannot write hits to {path}: it is a directory. Name the file itself, "
            f"e.g. {path / 'hits.parquet'}."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(_SCHEMA.metadata or {})
    metadata[HIT_PROVENANCE_KEY] = json.dumps(dict(provenance)).encode("utf-8")
    with pq.ParquetWriter(path, _SCHEMA.with_metadata(metadata)) as writer:
        for batch in batches:
            if batch.empty:
                continue
            writer.write_table(pa.Table.from_pandas(batch, schema=_SCHEMA, preserve_index=False))
    return path


def read_hits(path: str | Path) -> pd.DataFrame:
    """Read a **Hit table** back from Parquet, provenance and dtypes included.

    The counterpart of :func:`write_hits`, and the only reader that restores
    ``frame.attrs``: :func:`pandas.read_parquet` alone gives the rows and drops what the
    scan was, which is the same table with its meaning removed.

    Parameters
    ----------
    path : str or pathlib.Path
        A Parquet file written by :func:`write_hits`.

    Returns
    -------
    pandas.DataFrame
        The **Hit table**, equal to the in-memory form of the same scan — column order,
        dtypes and category order included — with
        :data:`~genome.tf.motif.scan.HIT_PROVENANCE` on ``frame.attrs``. A file written by
        something other than :func:`write_hits` simply carries no provenance, and
        ``attrs`` is then empty.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import numpy as np
    >>> from genome.tf.motif import Motif, MotifSet
    >>> counts = np.zeros((4, 8))
    >>> for column, base in enumerate("GATTACAG"):
    ...     counts["ACGT".index(base), column] = 100.0
    >>> motifs = MotifSet([Motif("MA9999.1", "Gattacag", counts)])
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     written = motifs.scan("TTTTTGATTACAGTTTTT", output=Path(directory) / "h.parquet")
    ...     hits = read_hits(written)
    >>> {name: str(dtype) for name, dtype in hits.dtypes.items()}["score"]
    'float16'
    >>> hits.attrs["threshold"]
    0.0001
    """
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(
            f"hit table not found: {source}. Scan with output={source} to write one."
        )
    frame = pd.read_parquet(source)
    for name, dtype in HIT_DTYPES.items():
        if dtype == "category" and name in frame.columns:
            # Arrow keeps the categories in the order the row groups introduced them;
            # astype("category") sorts them, and the in-memory table is what this must equal.
            column = frame[name]
            frame[name] = column.cat.set_categories(sorted(column.cat.categories))
    frame.attrs = _provenance_of(source)
    return frame


def _provenance_of(path: Path) -> dict[Hashable, Any]:
    """Read the scan's provenance out of the file's metadata, or ``{}`` if it carries none.

    Every JSON array becomes a tuple again, which is what a **Hit table**'s background and
    its two motif lists are in memory.
    """
    metadata = pq.read_schema(path).metadata or {}
    raw = metadata.get(HIT_PROVENANCE_KEY)
    if raw is None:
        return {}
    recorded = json.loads(raw.decode("utf-8"))
    return {
        key: tuple(value) if isinstance(value, list) else value for key, value in recorded.items()
    }
