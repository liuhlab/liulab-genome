"""The motif libraries answer in every environment, before any code depends on them.

Nothing under ``src/`` imports these yet — the motif subpackage is not written. They sit
in the core dependency table rather than behind a feature so that no lane of the suite
can be without them, which is the same reasoning that puts the GTF library there. A
library is not an **External tool**: it is not resolved on ``PATH`` and not version-
checked at call time, so there is nothing here to skip on and these tests belong in the
unit lane, unmarked.

Each test does a little work rather than only importing, because three of the four are
compiled or draw pixels, and a name that resolves is not yet a library that runs. The
scan engine among them was chosen by measurement; the numbers and the method are in
``docs/research/motif-scan-engine-2026-08-28.md``.
"""

from __future__ import annotations

from pathlib import Path

import logomaker
import matplotlib
import matplotlib.pyplot as plt
import memelite
import MOODS.parsers
import MOODS.scan
import MOODS.tools
import numpy as np
import pandas as pd
import xarray as xr

matplotlib.use("Agg")  # headless: no test may open a window

#: A 4 x 6 count matrix over A, C, G, T reading GATTAC with no ambiguity anywhere.
_COUNTS = [
    [0.0, 20.0, 0.0, 0.0, 20.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 20.0],
    [20.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 20.0, 20.0, 0.0, 0.0],
]


def test_moods_finds_a_planted_site_at_the_position_it_was_planted() -> None:
    # The compiled extension, its threshold arithmetic and its scanner, in one pass.
    background = MOODS.tools.flat_bg(4)
    matrix = MOODS.tools.log_odds(_COUNTS, background, 0.01)
    threshold = MOODS.tools.threshold_from_p(matrix, background, 1e-4)
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs([matrix], background, [threshold])

    sequence = "TTTTTTTTTTGATTACTTTTTTTTTT"
    (matches,) = scanner.scan(sequence)

    assert [match.pos for match in matches] == [sequence.index("GATTAC")]


def test_moods_parsers_reads_a_count_matrix_from_disk(tmp_path: Path) -> None:
    # MOODS.parsers is the third submodule the adapter will reach for.
    path = tmp_path / "gattac.pfm"
    path.write_text("\n".join(" ".join(str(x) for x in row) for row in _COUNTS) + "\n")

    matrix = MOODS.parsers.pfm(str(path))

    assert len(matrix) == 4
    assert all(len(row) == 6 for row in matrix)


def test_logomaker_draws_a_logo_as_wide_as_the_motif() -> None:
    # logomaker wants positions as rows, the transpose of the layout everything else uses.
    frame = pd.DataFrame(np.asarray(_COUNTS).T, columns=pd.Index(["A", "C", "G", "T"]))
    information = logomaker.transform_matrix(frame, from_type="counts", to_type="information")

    figure, axes = plt.subplots()
    try:
        logomaker.Logo(information, ax=axes)
        assert len(axes.get_xticks()) > 0
        assert information.shape == (6, 4)
    finally:
        plt.close(figure)


def test_xarray_labels_an_array_by_motif_id() -> None:
    # The shape a motif-versus-motif comparison result is carried in.
    scores = xr.DataArray(
        np.eye(2, dtype="float32"),
        dims=("query", "target"),
        coords={"query": ["MA0004.1", "MA0006.1"], "target": ["MA0004.1", "MA0006.1"]},
    )

    assert scores.sel(query="MA0004.1", target="MA0004.1") == 1.0


def test_memelite_exposes_the_comparison_entry_point() -> None:
    # Only tomtom is wanted from memelite; its scanner lost the measurement.
    assert callable(memelite.tomtom)
    assert callable(memelite.fimo)
