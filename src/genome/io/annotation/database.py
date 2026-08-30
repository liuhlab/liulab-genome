"""The **Annotation database** — building one, and reading its gene rows back.

The `gffutils` adapter, and the only module in the package that imports it. Two things are
asked of the library and they are both here: building the SQLite database beside a placed
**GTF**, and reading every gene feature's id out of one that was built. Keeping them
together is what makes the dependency one import in one small module rather than two lines
buried in the largest one — the API surface never says "gffutils", and swapping what backs
an **Annotation database** would be a change to this file.

**Feature inference** is off by default in both directions, deliberately: GENCODE, Ensembl
and RefSeq all declare ``gene`` and ``transcript`` rows already, and reconstructing them
from exon lines is the library's slow path.

Examples
--------
>>> from genome.io.annotation.database import LIBRARY_VERSION_KEY
>>> LIBRARY_VERSION_KEY
'gffutils_version'
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import gffutils

#: What a **Completion marker**'s ``details`` calls the library version that built the
#: database. Recorded there rather than in ``tool_versions``, which is for **External
#: tool**s: a tool resolved on ``PATH``, version-detected by running it, and installable by
#: a command an error can name. `gffutils` is an installed Python library and none of that
#: applies to it, so recording it there would blur the one word the package keeps sharp for
#: binaries it shells out to.
LIBRARY_VERSION_KEY = "gffutils_version"

#: What is asked of an **Annotation database** to resolve **Gene id stem**s: the id of
#: every gene feature, ascending, and nothing else. That id is the table's primary key —
#: so the ordering is an index walk — and for a GTF it is the ``gene_id`` gffutils keys
#: gene features by, which is the annotation's own spelling of its gene ids. It goes
#: through gffutils' own ``execute`` and is read a row at a time off the cursor: building
#: a feature object per row would parse an attribute blob per gene for a value already in
#: hand, and a GENCODE annotation has some 78,000 of them.
_GENE_IDS_QUERY = "SELECT id FROM features WHERE featuretype = 'gene' ORDER BY id"


def build_database(
    gtf: Path,
    database: Path,
    *,
    disable_infer_genes: bool = True,
    disable_infer_transcripts: bool = True,
) -> str:
    """Build the **Annotation database** at ``database`` from the **GTF** at ``gtf``.

    The one build in the package. An older database left by an interrupted or forced
    re-registration is replaced rather than refused, since this is reached only when the
    annotation is being built; the SQLite connection is closed behind us, because what the
    caller needs is the file and not a handle.

    Parameters
    ----------
    gtf : pathlib.Path
        The placed GTF to build from, read by the library rather than by this package.
    database : pathlib.Path
        Where to write the SQLite database, beside the GTF as ``<name>.db``.
    disable_infer_genes : bool, default True
        Do not reconstruct ``gene`` features from exon lines — **Feature inference** off.
    disable_infer_transcripts : bool, default True
        Do not reconstruct ``transcript`` features from exon lines.

    Returns
    -------
    str
        The version of the library that built it, for a **Completion marker** to record
        under :data:`LIBRARY_VERSION_KEY`.

    Examples
    --------
    >>> from pathlib import Path
    >>> build_database(Path("mine.gtf"), Path("mine.db"))     # doctest: +SKIP
    '0.13'
    """
    built = gffutils.create_db(
        str(gtf),
        str(database),
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts,
    )
    # The on-disk database is now fully written; release the SQLite connection
    # so we don't leak an open file handle (the build is the only thing we need).
    built.conn.close()
    return gffutils.__version__


def gene_ids(database: Path) -> Iterator[str]:
    """Yield the id of every gene feature in ``database``, ascending, one row at a time.

    The database is *queried* rather than read: a cursor over the gene features alone,
    walked in primary-key order, so a GENCODE-sized annotation costs an index walk of its
    gene rows and nothing is ever held in memory. A caller that stops early closes the
    generator, which closes the connection with it — as does exhausting it, and as does an
    exception on the way out.

    Parameters
    ----------
    database : pathlib.Path
        The **Annotation database** to walk. It is opened read-only for the walk and
        closed when the walk ends.

    Yields
    ------
    str
        Each gene feature's id, in ascending order, which for a GTF is the annotation's
        own spelling of that gene's id.

    Examples
    --------
    >>> from pathlib import Path
    >>> next(gene_ids(Path("gencode_v50.db")))                # doctest: +SKIP
    'ENSG00000000003.16'
    """
    handle = gffutils.FeatureDB(str(database))
    try:
        for row in handle.execute(_GENE_IDS_QUERY):
            yield str(row["id"])
    finally:
        # The registry hands back an answer and not a handle, so the SQLite connection
        # this opened closes here — including on the way out of an exception.
        handle.conn.close()
