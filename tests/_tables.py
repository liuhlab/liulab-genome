"""Shared fixtures for the shipped-table tests: ``test_link.py``, ``test_tf_link.py``,
``test_tf_census.py`` and ``test_tf_cofactor.py``.

Two things live here rather than being reinvented once per file. First, a helper that
builds a fake shipped table's text and hands it to a parse entry point, rather than
writing a file to disk, which is what each entry point is public (or, for the census,
package-private but already that module's own testing seam) for. Not a fixture and not a
conftest addition — conftest.py is shared by every test module in the suite, and this is
one specific convenience four of them happen to want the same shape of. The one real
difference between the shapes is whether the header can be overridden: the two
uniform-four tables (census, cofactor) each exercise a "the header does not lead with the
required columns" refusal that needs a header of its own; the link table's
malformed-header cases build their raw text directly instead and never call this helper
for that case, so it never exercises the override. Kept as one optional parameter rather
than two helpers, since every other line is the same shape.

Second, the representative sample used everywhere a per-table-and-release or a
per-species check is run against a sample rather than every shipped file — one human
table and one mouse table, on different releases, since those are the two boundaries
(species, release) most of these checks are about. Named once so that widening or
narrowing the sample is one edit rather than three.
"""

from __future__ import annotations

from collections.abc import Sequence

#: One human table and one mouse table, on different releases: enough to cross the
#: species boundary and the release boundary at once for a row- or table-level rule,
#: without re-running it against every shipped file — the pinned-count tests in
#: ``test_tf_link.py`` and ``test_tf_cofactor.py`` already guard the full population.
REPRESENTATIVE_TABLES = (("homo_sapiens", "2026"), ("mus_musculus", "2024"))


def table_text(
    columns: Sequence[str],
    *rows: str | Sequence[str],
    header: Sequence[str] | str | None = None,
) -> str:
    """Return a table's text: ``header`` (or ``columns``) then one line per row.

    A row is either a pre-joined line (``"g1\\tA\\tyes\\tX"``) or a sequence of cells to
    join with a tab, so a caller that builds rows cell-by-cell and one that pastes a
    literal line both work — without pretending the two conventions are one.

    >>> table_text(("a", "b"), "1\\t2", ("3", "4"))
    'a\\tb\\n1\\t2\\n3\\t4\\n'
    """
    if header is not None:
        head = header if isinstance(header, str) else "\t".join(header)
    else:
        head = "\t".join(columns)
    lines = [head, *(row if isinstance(row, str) else "\t".join(row) for row in rows)]
    return "\n".join(lines) + "\n"
