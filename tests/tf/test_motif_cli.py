"""Tests for the ``genome motif`` sub-app."""

from __future__ import annotations

import json as _json
from pathlib import Path
from types import MappingProxyType

import pytest

from genome.cli import app
from genome.tf.motif import MIN_MOTIF_LENGTH, hit_count, provenance_of, read_hits
from genome.tf.motif import jaspar as jaspar_mod
from genome.tf.motif.jaspar import MOTIF_COUNTS

from .._cli import output as _output
from .._cli import runner
from ..conftest import FakeFetch
from .test_jaspar import FIXTURE as _MOTIF_FIXTURE
from .test_jaspar import FIXTURE_COUNT as _MOTIF_COUNT
from .test_jaspar import FIXTURE_MOTIFS as _MOTIF_RECORDS
from .test_scan import FIXTURE as _PLANTED_FASTA

#: Which of the committed motifs a scan leaves out, and how many are left to scan with.
#: Read off the fixture table rather than written down again, so a changed fixture moves the
#: expected summary with it: a motif under the minimum length cannot reach the default
#: **Threshold** at all, and is named among the skipped rather than called at something looser.
_MOTIFS_TOO_SHORT = [
    motif_id for motif_id, _name, length, _tax in _MOTIF_RECORDS if length < MIN_MOTIF_LENGTH
]
_MOTIFS_LONG_ENOUGH = _MOTIF_COUNT - len(_MOTIFS_TOO_SHORT)


@pytest.fixture
def motif_release(fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> FakeFetch:
    """Serve the committed transfac records as whichever **Release** is asked for.

    The arrangement ``tests/tf/test_jaspar.py`` uses: the count check that stands where a
    **Completion marker** stands elsewhere is never switched off, only pointed at what the
    fake fetch actually serves.
    """
    monkeypatch.setattr(
        jaspar_mod, "MOTIF_COUNTS", MappingProxyType(dict.fromkeys(MOTIF_COUNTS, _MOTIF_COUNT))
    )
    fake_fetch.serve(_MOTIF_FIXTURE)
    return fake_fetch


@pytest.fixture
def planted_fasta(data_dir: Path) -> Path:
    """The committed FASTA with motifs planted at positions ``tests/data/README.md`` lists."""
    return data_dir / _PLANTED_FASTA


@pytest.mark.xdist_group("spawns_mixin")
class TestMotifScan:
    """The batch scan: what the summary says, where each half of the answer goes.

    Offline like every other test here — the ``fake_fetch`` fixture serves the committed
    transfac records as whatever release is asked for, with the count check pointed at
    them — and unmarked, since a process pool is not a binary this package ships.
    """

    def _summary(self, result: object) -> dict[str, object]:
        """Parse the JSON summary, which is the whole of what the command puts on stdout."""
        payload = _json.loads(getattr(result, "stdout", ""))
        assert isinstance(payload, dict)
        return payload

    def test_the_json_summary_carries_the_run(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "hits.parquet"

        result = runner.invoke(
            app,
            [
                "motif",
                "scan",
                str(planted_fasta),
                str(output),
                "--release",
                "2024",
                "--tax-group",
                "all",
                "--workers",
                "1",
                "--json",
            ],
        )

        assert result.exit_code == 0, _output(result)
        summary = self._summary(result)
        assert summary == {
            "release": "2024",
            "tax_group": "all",
            "motifs_scanned": _MOTIFS_LONG_ENOUGH,
            "motifs_skipped": _MOTIFS_TOO_SHORT,
            # Two 600-base records are far under the derivation floor, so 'auto' is uniform
            # — and says so rather than leaving it to be assumed.
            "background": [0.25, 0.25, 0.25, 0.25],
            "threshold": 1e-4,
            "sequences_scanned": 2,
            "hits_written": hit_count(output),
            "workers": 1,
            "output": str(output),
        }

    def test_the_summary_goes_to_stdout_hits_to_the_file_and_the_human_form_agrees(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        # The whole reason the hits are written rather than printed: stdout is one JSON
        # document and nothing else, whatever the scan found.
        output = tmp_path / "hits.parquet"
        result = runner.invoke(
            app, ["motif", "scan", str(planted_fasta), str(output), "--workers", "1", "--json"]
        )
        summary = self._summary(result)
        written = read_hits(output)
        assert len(written) > 0
        assert summary["hits_written"] == len(written)
        assert "plantedI" in set(written["sequence_name"])
        # The table's own provenance is what the summary was read off, so the two agree.
        assert written.attrs["release"] == summary["release"]

        human_output = tmp_path / "human.parquet"
        human = runner.invoke(
            app, ["motif", "scan", str(planted_fasta), str(human_output), "--workers", "1"]
        )
        assert human.exit_code == 0, _output(human)
        assert f"scanned 2 sequences with {_MOTIFS_LONG_ENOUGH} motifs" in human.stdout
        assert str(human_output) in human.stdout
        # Which motifs were left out is printed, not silently dropped: an absent factor
        # is explainable only if the scan says it never scanned for it.
        for motif_id in _MOTIFS_TOO_SHORT:
            assert motif_id in human.stdout

    def test_a_missing_fasta_a_bad_argument_is_refused_and_an_empty_scan_still_writes(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        missing = tmp_path / "nowhere.fa"
        missing_output = tmp_path / "missing.parquet"
        missing_result = runner.invoke(
            app, ["motif", "scan", str(missing), str(missing_output), "--workers", "1", "--json"]
        )
        assert missing_result.exit_code == 1
        assert "not found" in _output(missing_result)
        assert str(missing) in _output(missing_result)
        # Nothing half-written to be mistaken for an answer.
        assert not missing_output.exists()

        empty = tmp_path / "unreadable.fa"
        empty.write_text(">nothing\n" + "N" * 400 + "\n")
        empty_output = tmp_path / "empty.parquet"
        empty_result = runner.invoke(
            app, ["motif", "scan", str(empty), str(empty_output), "--workers", "1", "--json"]
        )
        assert empty_result.exit_code == 0, _output(empty_result)
        assert self._summary(empty_result)["hits_written"] == 0
        assert empty_output.is_file()

        bad_release = runner.invoke(
            app,
            [
                "motif",
                "scan",
                str(planted_fasta),
                str(tmp_path / "hits.parquet"),
                "--release",
                "2019",
            ],
        )
        assert bad_release.exit_code == 1
        assert "2024, 2026" in _output(bad_release)

        bad_threshold = runner.invoke(
            app,
            [
                "motif",
                "scan",
                str(planted_fasta),
                str(tmp_path / "hits.parquet"),
                "--threshold",
                "5",
            ],
        )
        assert bad_threshold.exit_code == 1
        assert "p-value" in _output(bad_threshold)

        output = tmp_path / "hits.parquet"
        zero_workers = runner.invoke(
            app, ["motif", "scan", str(planted_fasta), str(output), "--workers", "0"]
        )
        assert zero_workers.exit_code == 1
        assert "at least 1" in _output(zero_workers)
        assert not output.exists()

    def test_a_background_mode_reaches_the_scan_or_is_refused_before_anything_runs(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        # An invalid mode is refused first, on a set that has fetched nothing yet: not
        # even the release is fetched before an argument this basic is checked.
        bad_output = tmp_path / "bad.parquet"
        bad = runner.invoke(
            app, ["motif", "scan", str(planted_fasta), str(bad_output), "--background", "gc"]
        )
        assert bad.exit_code == 2
        assert not bad_output.exists()
        assert not motif_release.calls

        # The parameter that decides the answer more than any other, and the summary
        # reports the one actually used rather than the one asked for.
        output = tmp_path / "hits.parquet"
        result = runner.invoke(
            app,
            [
                "motif",
                "scan",
                str(planted_fasta),
                str(output),
                "--background",
                "derive",
                "--workers",
                "1",
                "--json",
            ],
        )
        background = self._summary(result)["background"]
        assert background != [0.25, 0.25, 0.25, 0.25]
        assert list(provenance_of(output)["background"]) == background

    def test_the_worker_count_defaults_to_the_allocation(
        self,
        motif_release: FakeFetch,
        planted_fasta: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        # The one place the command deliberately differs from the library, which defaults
        # to serial: a console script is a proper entry point, so it takes what it was
        # given — the allocation first, and never the machine's cores over it.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
        shared, serial = tmp_path / "shared.parquet", tmp_path / "serial.parquet"

        allocated = runner.invoke(app, ["motif", "scan", str(planted_fasta), str(shared), "--json"])
        alone = runner.invoke(
            app, ["motif", "scan", str(planted_fasta), str(serial), "--workers", "1", "--json"]
        )

        assert self._summary(allocated)["workers"] == 2
        assert self._summary(alone)["workers"] == 1
        # And the choice is about wall time and nothing else.
        assert read_hits(shared).equals(read_hits(serial))

    def test_the_progress_display_is_suppressed_under_the_json_flag(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        runner.invoke(
            app,
            [
                "motif",
                "scan",
                str(planted_fasta),
                str(tmp_path / "hits.parquet"),
                "--workers",
                "1",
                "--json",
            ],
        )

        assert motif_release.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(
        self, motif_release: FakeFetch, planted_fasta: Path, tmp_path: Path
    ) -> None:
        runner.invoke(
            app,
            ["motif", "scan", str(planted_fasta), str(tmp_path / "hits.parquet"), "--workers", "1"],
        )

        assert motif_release.last.progressbar is True
