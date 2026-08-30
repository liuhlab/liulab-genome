"""Tests for genome.tf.motif.scan — the MOODS adapter and the **Hit table** it answers with.

**MOODS is called for real here.** It is a core dependency, deterministic, and fast enough
on a 1.2 kb fixture that a fake would only test the adapter against our own assumptions
about the engine — which is the class of bug the benchmark behind the engine choice found.
So the one claim this whole feature rests on is checked against the engine itself: that a
position reported for a reverse-complement matrix is already a forward-frame start. The
fixture plants one site and its own reverse complement and asserts both land on the bases
they were written over.

The committed bytes of ``tests/data/planted_motifs.fa`` are asserted here rather than
trusted from the README — the backbone against ``tiny.fa``'s own bases, the planted words
against their offsets, and the soft-masked window against its bounds.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from types import MappingProxyType
from typing import Any

import MOODS.tools
import numpy as np
import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from genome.seq import DNA
from genome.tf.motif import (
    BACKGROUND_FLOOR,
    DEFAULT_SEQUENCE_NAME,
    DEFAULT_THRESHOLD,
    HIT_COLUMNS,
    HIT_DTYPES,
    HIT_PROVENANCE,
    MIN_MOTIF_LENGTH,
    UNIFORM_BACKGROUND,
    FastaFormatError,
    JasparDatabase,
    Motif,
    MotifSet,
    parse_transfac,
)
from genome.tf.motif import jaspar as jaspar_mod
from genome.tf.motif.jaspar import MOTIF_COUNTS
from genome.tf.motif.scan import empty_hits, read_fasta, scan_stream

from .conftest import FakeFetch

#: The committed motif records every scan here is run with — the same ten #117 uses.
MOTIF_FIXTURE = "tiny_jaspar_transfac.txt"

#: The committed FASTA. See ``tests/data/README.md``.
FIXTURE = "planted_motifs.fa"

#: What each record's backbone is a copy of: a record of ``tiny.fa`` and the 0-based
#: half-open window of it taken, before anything was planted over it.
BACKBONE: tuple[tuple[str, str, int, int], ...] = (
    ("plantedI", "chrI", 0, 600),
    ("plantedII", "chrII", 0, 600),
)

#: What was planted, where, and on which strand — the table ``tests/data/README.md``
#: carries in prose. Each row is the record, the 0-based half-open interval the word
#: occupies, the **Strand** a scan must report it on, the **Motif id** whose consensus it
#: is, and the bases as the file holds them, case and all.
PLANTED: tuple[tuple[str, int, int, str, str, str], ...] = (
    ("plantedI", 100, 115, "+", "MA0139.2", "GCCACCAGGGGGCGC"),
    ("plantedI", 300, 315, "-", "MA0139.2", "GCGCCCCCTGGTGGC"),
    ("plantedII", 200, 209, "+", "MA0789.1", "tatgcaaat"),
)

#: The one soft-masked stretch, 0-based half-open — one wrapped line of ``plantedII``,
#: holding the third planted site.
MASKED = (180, 240)

#: How the committed FASTA is wrapped.
WRAP = 60

#: The two records of the motif fixture below :data:`MIN_MOTIF_LENGTH`, so no scan may
#: call them and every scan must name them.
TOO_SHORT: tuple[str, ...] = ("MA2355.1", "MA0261.1")


def read_records(text: str) -> dict[str, str]:
    """Read a FASTA into ``{name: bases}`` — the test's own reader, not the package's."""
    records: dict[str, list[str]] = {}
    name = ""
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = []
        else:
            records[name].append(line)
    return {name: "".join(lines) for name, lines in records.items()}


def revcomp(bases: str) -> str:
    """Reverse-complement, case preserved."""
    return str(DNA(bases).reverse_complement())


def word_motif(bases: str, motif_id: str = "MA9999.1", name: str = "Testin") -> Motif:
    """A motif fixed on every base of ``bases`` — 2 bits a position."""
    counts = np.zeros((4, len(bases)))
    for column, base in enumerate(bases.upper()):
        counts["ACGT".index(base), column] = 100.0
    return Motif(motif_id, name, counts)


@pytest.fixture
def motifs(data_dir: Path) -> MotifSet:
    """The ten committed JASPAR records as a plain **Motif set**."""
    return MotifSet(parse_transfac((data_dir / MOTIF_FIXTURE).read_text(encoding="utf-8")))


@pytest.fixture
def planted(data_dir: Path) -> Path:
    """The committed FASTA with motifs planted in it."""
    return data_dir / FIXTURE


@pytest.fixture
def planted_records(planted: Path) -> dict[str, str]:
    """The committed FASTA read into ``{name: bases}``, case as committed."""
    return read_records(planted.read_text(encoding="utf-8"))


@pytest.fixture
def hits(motifs: MotifSet, planted: Path) -> pd.DataFrame:
    """The committed FASTA scanned with the committed motifs, at every default."""
    return motifs.scan_fasta(planted)


@pytest.fixture
def release(fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch) -> JasparDatabase:
    """The same ten motifs, but arriving as a **Release** that knows what it is."""
    monkeypatch.setattr(
        jaspar_mod, "MOTIF_COUNTS", MappingProxyType(dict.fromkeys(MOTIF_COUNTS, 10))
    )
    fake_fetch.serve(MOTIF_FIXTURE)
    return JasparDatabase("2024", "all")


def wide(bases: str, copies: int = 20) -> dict[str, str]:
    """``copies`` named sequences of the same bases — an input over the derivation floor.

    Every record identical, so the composition of any prefix is the composition of all of
    them and the expected background does not depend on how many were read while deciding.
    """
    return {f"peak{index}": bases for index in range(copies)}


def composition(bases: str) -> list[float]:
    """The four base frequencies of one sequence, ambiguous bases ignored, in ACGT order."""
    counts = [bases.upper().count(base) for base in "ACGT"]
    return [count / sum(counts) for count in counts]


def rows(frame: pd.DataFrame) -> list[tuple[Any, ...]]:
    """Every row as a plain tuple, so two tables compare without their category dtypes."""
    return [tuple(row) for row in frame.itertuples(index=False, name=None)]


def sites(frame: pd.DataFrame) -> set[tuple[Any, ...]]:
    """Every hit as ``(motif_id, sequence_name, start, end, strand)`` — the score dropped."""
    return {
        (motif_id, sequence_name, int(start), int(end), strand)
        for motif_id, _name, sequence_name, start, end, strand, _score in rows(frame)
    }


def of(frame: pd.DataFrame, **equals: object) -> list[tuple[Any, ...]]:
    """The rows whose named columns hold the given values, as plain tuples."""
    wanted = [(list(HIT_COLUMNS).index(column), value) for column, value in equals.items()]
    return [row for row in rows(frame) if all(row[index] == value for index, value in wanted)]


# ---------------------------------------------------------------------------
# The committed bytes: everything the README says about the fixture
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    def test_the_fixture_shape_headers_and_format(
        self, planted: Path, planted_records: dict[str, str]
    ) -> None:
        assert list(planted_records) == ["plantedI", "plantedII"]
        assert [len(bases) for bases in planted_records.values()] == [600, 600]
        # What the whitespace truncation is here to be tested against.
        headers = [
            line for line in planted.read_text(encoding="utf-8").splitlines() if line[:1] == ">"
        ]
        assert headers == [
            ">plantedI",
            ">plantedII  sacCer3 chrII:1-600, bases 180-240 soft-masked",
        ]
        for bases in planted_records.values():
            assert DNA.outside_alphabet(bases) == []
        widths = {
            len(line)
            for line in planted.read_text(encoding="utf-8").splitlines()
            if not line.startswith(">")
        }
        assert widths == {WRAP}

    def test_the_backbone_is_sacCer3_and_only_the_planted_words_are_not(  # noqa: N802
        self, data_dir: Path, planted_records: dict[str, str]
    ) -> None:
        # The README promises real bases everywhere except the planted words; asserted by
        # putting the source bases back and getting tiny.fa's own window out.
        tiny = read_records((data_dir / "tiny.fa").read_text(encoding="utf-8"))
        for name, source, start, end in BACKBONE:
            restored = list(planted_records[name].upper())
            for record, low, high, _strand, _motif_id, _word in PLANTED:
                if record == name:
                    restored[low:high] = tiny[source][start + low : start + high]
            assert "".join(restored) == tiny[source][start:end]

    def test_each_planted_word_sits_at_its_offset_and_the_reverse_is_the_forward_flipped(
        self, planted_records: dict[str, str], motifs: MotifSet
    ) -> None:
        for record, start, end, strand, motif_id, word in PLANTED:
            assert planted_records[record][start:end] == word
            expected = str(motifs[motif_id].consensus)
            assert word.upper() == (expected if strand == "+" else revcomp(expected))
        forward, reverse = PLANTED[0], PLANTED[1]
        bases = planted_records["plantedI"]
        assert revcomp(bases[reverse[1] : reverse[2]]) == bases[forward[1] : forward[2]]

    def test_the_masked_window_bounds_and_holds_the_third_planted_site(
        self, planted_records: dict[str, str]
    ) -> None:
        low, high = MASKED
        bases = planted_records["plantedII"]
        assert bases[low:high].islower()
        assert bases[:low].isupper()
        assert bases[high:].isupper()
        assert planted_records["plantedI"].isupper()
        _record, start, end, _strand, _motif_id, _word = PLANTED[2]
        assert low <= start < end <= high


# ---------------------------------------------------------------------------
# Reading a FASTA
# ---------------------------------------------------------------------------


class TestReadFasta:
    def test_reading_the_fixture_gives_names_and_bases_as_written_gzipped_or_not(
        self, planted: Path, planted_records: dict[str, str], tmp_path: Path
    ) -> None:
        # What STAR and chromap write into an alignment from the same file, so a hit table
        # joins against it without anyone renaming anything.
        assert [name for name, _ in read_fasta(planted)] == ["plantedI", "plantedII"]
        assert dict(read_fasta(planted)) == planted_records

        zipped = tmp_path / "planted.fa.gz"
        with gzip.open(zipped, "wt", encoding="utf-8") as handle:
            handle.write(planted.read_text(encoding="utf-8"))
        assert dict(read_fasta(zipped)) == dict(read_fasta(planted))

    def test_malformed_or_missing_files_are_refused_but_an_empty_one_yields_no_records(
        self, tmp_path: Path
    ) -> None:
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            list(read_fasta(tmp_path / "nope.fa"))

        bare = tmp_path / "bare.fa"
        bare.write_text("ACGTACGT\n")
        with pytest.raises(FastaFormatError, match="line 1"):
            list(read_fasta(bare))

        nameless = tmp_path / "nameless.fa"
        nameless.write_text(">chrI\nACGT\n>\nACGT\n")
        with pytest.raises(FastaFormatError, match="no name"):
            list(read_fasta(nameless))

        empty = tmp_path / "empty.fa"
        empty.write_text("")
        assert list(read_fasta(empty)) == []


# ---------------------------------------------------------------------------
# The schema, which is the contract
# ---------------------------------------------------------------------------


class TestSchema:
    def test_the_schema_and_index_are_the_contract_hits_or_not(
        self, hits: pd.DataFrame, motifs: MotifSet, planted: Path
    ) -> None:
        assert list(hits.columns) == list(HIT_COLUMNS)
        assert list(HIT_COLUMNS) == [
            "motif_id",
            "motif_name",
            "sequence_name",
            "start",
            "end",
            "strand",
            "score",
        ]
        assert {name: str(dtype) for name, dtype in hits.dtypes.items()} == dict(HIT_DTYPES)
        assert list(hits.index) == list(range(len(hits)))

        # And a scan that finds nothing, or has nothing to scan with, keeps the same
        # schema — the empty table built by hand included.
        found = motifs.scan("N" * 300)
        assert len(found) == 0
        assert list(found.columns) == list(HIT_COLUMNS)
        assert {name: str(dtype) for name, dtype in found.dtypes.items()} == dict(HIT_DTYPES)

        no_motifs = MotifSet([]).scan_fasta(planted)
        assert len(no_motifs) == 0
        assert no_motifs.attrs["motifs_scanned"] == ()

        assert list(empty_hits().columns) == list(HIT_COLUMNS)
        assert {name: str(dtype) for name, dtype in empty_hits().dtypes.items()} == dict(HIT_DTYPES)


# ---------------------------------------------------------------------------
# The planted sites, and the frame every coordinate is in
# ---------------------------------------------------------------------------


class TestPlantedSites:
    def test_every_planted_site_is_found_where_it_was_planted(self, hits: pd.DataFrame) -> None:
        for record, start, end, strand, motif_id, _word in PLANTED:
            assert (motif_id, record, start, end, strand) in sites(hits)

    def test_a_reverse_site_covers_the_same_bases_as_its_forward_equivalent(
        self, hits: pd.DataFrame, planted_records: dict[str, str]
    ) -> None:
        # The whole claim the adapter rests on: a reverse-complement matrix reports a
        # forward-frame start, so the two intervals below are read off the same strand of
        # the same record and one is the other flipped.
        found = of(hits, motif_id="MA0139.2", sequence_name="plantedI")
        by_strand = {row[5]: (int(row[3]), int(row[4])) for row in found}
        assert set(by_strand) == {"+", "-"}
        bases = planted_records["plantedI"]
        forward = bases[slice(*by_strand["+"])]
        reverse = bases[slice(*by_strand["-"])]
        assert revcomp(reverse) == forward
        assert by_strand["+"][1] - by_strand["+"][0] == by_strand["-"][1] - by_strand["-"][0]

    def test_every_hit_obeys_the_interval_and_provenance_invariants(
        self, hits: pd.DataFrame, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        for motif_id, _name, sequence, start, end, strand, _score in rows(hits):
            assert int(end) - int(start) == len(motifs[str(motif_id)])
            assert 0 <= int(start) < int(end) <= len(planted_records[str(sequence)])
            assert strand in {"+", "-"}
            assert motif_id in hits.attrs["motifs_scanned"]
        assert {"+", "-"} <= set(hits["strand"])

    def test_the_score_is_log_odds_in_bits_and_not_a_p_value(
        self, hits: pd.DataFrame, motifs: MotifSet
    ) -> None:
        # The planted word is the consensus, so its score is the matrix's best possible
        # one — in bits, which is what makes that number checkable at all.
        best = float(motifs["MA0139.2"].log_odds().max(axis=0).sum())
        scored = [float(row[6]) for row in of(hits, motif_id="MA0139.2")]
        assert max(scored) == pytest.approx(best, abs=0.02)
        assert best > 1.0  # a p-value could not be


# ---------------------------------------------------------------------------
# What the engine is trusted for, checked against the engine
# ---------------------------------------------------------------------------


class TestForwardFrame:
    def test_a_word_planted_at_a_known_offset_is_found_there_on_either_strand(self) -> None:
        word = "GATTACAGTC"
        filler = "AAACCCTTT" * 5
        for strand, planted in [("+", word), ("-", revcomp(word))]:
            found = MotifSet([word_motif(word)]).scan(filler + planted + filler)
            assert (
                "MA9999.1",
                DEFAULT_SEQUENCE_NAME,
                len(filler),
                len(filler) + len(word),
                strand,
            ) in sites(found)

    def test_the_boundaries_of_the_sequence_are_handled(self) -> None:
        word = "GATTACAGTC"
        filler = "AAACCCTTT" * 5
        # The engine was chosen partly because it scans the last window; a site at the
        # final position must not need a base appended after it to be found.
        at_the_end = MotifSet([word_motif(word)]).scan(filler + word)
        assert (int(at_the_end["end"].max())) == len(filler) + len(word)
        # And a sequence shorter than the motif yields no hit rather than raising.
        too_short = MotifSet([word_motif(word)]).scan("ACGT")
        assert len(too_short) == 0

    @settings(max_examples=30)
    @given(
        word=st.text(alphabet="ACGT", min_size=8, max_size=14),
        left=st.text(alphabet="ACGT", min_size=0, max_size=30),
        right=st.text(alphabet="ACGT", min_size=0, max_size=30),
    )
    def test_a_reverse_complement_hit_is_a_forward_frame_start(
        self, word: str, left: str, right: str
    ) -> None:
        # Whatever the word, planting its reverse complement puts a '-' hit on exactly the
        # bases it was written over — nothing is subtracted from what the engine reports.
        found = MotifSet([word_motif(word)]).scan(left + revcomp(word) + right)
        assert (
            "MA9999.1",
            DEFAULT_SEQUENCE_NAME,
            len(left),
            len(left) + len(word),
            "-",
        ) in sites(found)

    @settings(max_examples=30)
    @given(sequence=st.text(alphabet="ACGTN", max_size=200))
    def test_every_hit_lands_inside_the_sequence_it_was_found_in(self, sequence: str) -> None:
        found = MotifSet([word_motif("GATTACAGTC")]).scan(sequence)
        for _id, _name, _sequence, start, end, _strand, _score in rows(found):
            assert int(start) >= 0
            assert int(end) - int(start) == 10
            assert int(end) <= len(sequence)


# ---------------------------------------------------------------------------
# Soft-masking, which a scan discards without being asked (ADR-0012)
# ---------------------------------------------------------------------------


class TestSoftMasking:
    def test_the_fixture_softmask_makes_no_difference_to_the_scan(
        self, motifs: MotifSet, hits: pd.DataFrame, planted_records: dict[str, str]
    ) -> None:
        masked = planted_records["plantedII"]
        assert not masked.isupper()  # the fixture really does carry masking
        assert rows(motifs.scan(masked, "x")) == rows(motifs.scan(masked.upper(), "x"))
        # And the real planted site under that mask is still found despite its case.
        record, start, end, strand, motif_id, word = PLANTED[2]
        assert word.islower()
        assert (motif_id, record, start, end, strand) in sites(hits)

    def test_a_wholly_lower_case_file_scans_as_its_upper_case_one(
        self, motifs: MotifSet, planted: Path, tmp_path: Path
    ) -> None:
        text = planted.read_text(encoding="utf-8")
        lowered = "\n".join(
            line if line.startswith(">") else line.lower() for line in text.splitlines()
        )
        path = tmp_path / "lowered.fa"
        path.write_text(lowered + "\n")
        assert rows(motifs.scan_fasta(path)) == rows(motifs.scan_fasta(planted))

    @settings(max_examples=30)
    @given(sequence=st.text(alphabet="ACGT", max_size=150))
    def test_masking_changes_no_hit_whatever_the_sequence(self, sequence: str) -> None:
        motifs = MotifSet([word_motif("GATTACAGTC")])
        assert rows(motifs.scan(sequence)) == rows(motifs.scan(sequence.lower()))


# ---------------------------------------------------------------------------
# The minimum length, and the motifs it leaves out
# ---------------------------------------------------------------------------


class TestSkippedMotifs:
    def test_short_motifs_are_never_scanned_and_are_always_named(
        self, hits: pd.DataFrame, motifs: MotifSet
    ) -> None:
        assert set(hits["motif_id"]).isdisjoint(TOO_SHORT)
        assert hits.attrs["motifs_skipped"] == TOO_SHORT
        assert hits.attrs["motifs_scanned"] == tuple(
            motif.motif_id for motif in motifs if motif.motif_id not in TOO_SHORT
        )
        assert set(hits.attrs["motifs_scanned"]) | set(hits.attrs["motifs_skipped"]) == set(
            motifs.motif_ids
        )

    def test_the_boundary_is_seven_positions_and_an_all_short_set_scans_empty(
        self, planted: Path
    ) -> None:
        short = word_motif("GATTAC", "MA0001.1")  # six, and unreachable at any threshold
        long = word_motif("GATTACA", "MA0002.1")  # seven, the shortest that is scannable
        assert (len(short), len(long)) == (MIN_MOTIF_LENGTH - 1, MIN_MOTIF_LENGTH)
        found = MotifSet([short, long]).scan_fasta(planted)
        assert found.attrs["motifs_skipped"] == ("MA0001.1",)
        assert found.attrs["motifs_scanned"] == ("MA0002.1",)

        only_short = MotifSet([short]).scan_fasta(planted)
        assert len(only_short) == 0
        assert only_short.attrs["motifs_skipped"] == ("MA0001.1",)


# ---------------------------------------------------------------------------
# The threshold, and what the score column is not
# ---------------------------------------------------------------------------


class TestThreshold:
    def test_stricter_is_a_subset_looser_is_a_superset_and_out_of_range_is_refused(
        self, motifs: MotifSet, planted: Path
    ) -> None:
        assert DEFAULT_THRESHOLD == 1e-4
        loose = sites(motifs.scan_fasta(planted, threshold=1e-4))
        strict = sites(motifs.scan_fasta(planted, threshold=1e-6))
        assert strict < loose
        assert loose < sites(motifs.scan_fasta(planted, threshold=1e-3))

        for threshold in (0.0, 1.0):
            with pytest.raises(ValueError, match="per-position p-value"):
                motifs.scan("ACGT" * 20, threshold=threshold)

    def test_a_thresholds_conversion_is_cached_per_value_and_not_reused_across_values(
        self, motifs: MotifSet, planted: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # The engine's one slow step — seconds for a full vertebrate release — and a pure
        # function of (matrices, background, p). The engine still runs; it is counted.
        calls: list[float] = []
        real = MOODS.tools.threshold_from_p

        def counted(matrix: object, background: object, p: float, *args: object) -> float:
            calls.append(p)
            return real(matrix, background, p, *args)

        monkeypatch.setattr(MOODS.tools, "threshold_from_p", counted)
        first = motifs.scan_fasta(planted)
        converted = len(calls)
        assert converted > 0
        second = motifs.scan_fasta(planted)
        assert len(calls) == converted
        assert rows(second) == rows(first)

        motifs.scan_fasta(planted, threshold=1e-6)
        assert len(calls) == 2 * converted


# ---------------------------------------------------------------------------
# The background, and the provenance the table carries
# ---------------------------------------------------------------------------


class TestBackgroundAndProvenance:
    def test_provenance_keys_are_always_present_hit_or_not(
        self, hits: pd.DataFrame, motifs: MotifSet
    ) -> None:
        assert set(hits.attrs) == set(HIT_PROVENANCE)
        empty = motifs.scan("N" * 200)
        assert len(empty) == 0
        assert set(empty.attrs) == set(HIT_PROVENANCE)
        assert empty.attrs["motifs_skipped"] == TOO_SHORT

    def test_background_mode_auto_uniform_and_derive_pick_correctly_relative_to_the_floor(
        self,
        motifs: MotifSet,
        hits: pd.DataFrame,
        planted: Path,
        planted_records: dict[str, str],
    ) -> None:
        # 1200 bases: a composition estimated from that few would distort its own cutoffs,
        # so the actual fixture scan — under the floor — stays uniform.
        assert sum(len(bases) for bases in planted_records.values()) < BACKGROUND_FLOOR
        assert hits.attrs["background"] == UNIFORM_BACKGROUND

        bases = planted_records["plantedI"]
        peaks = wide(bases)
        assert sum(len(sequence) for sequence in peaks.values()) > BACKGROUND_FLOOR

        # Auto derives once the floor is crossed, from the very bases it was given.
        over_floor = motifs.scan_sequences(peaks)
        assert over_floor.attrs["background"] != UNIFORM_BACKGROUND
        for recorded, wanted in zip(
            over_floor.attrs["background"], composition(bases), strict=True
        ):
            assert recorded == pytest.approx(wanted, abs=0.002)

        # And either mode can be asked for explicitly regardless of which side of the
        # floor the input falls on: uniform over the floor, and derive under it — from
        # the whole fixture this time, so the composition it derives is checked too.
        forced_uniform = motifs.scan_sequences(peaks, background="uniform")
        assert forced_uniform.attrs["background"] == UNIFORM_BACKGROUND
        forced_derive = motifs.scan_fasta(planted, background="derive")
        assert forced_derive.attrs["background"] != UNIFORM_BACKGROUND
        whole = "".join(planted_records.values())
        for recorded, wanted in zip(
            forced_derive.attrs["background"], composition(whole), strict=True
        ):
            assert recorded == pytest.approx(wanted, abs=0.002)

    def test_an_explicit_background_is_recorded_and_wins_over_the_input(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str]
    ) -> None:
        given_background = motifs.scan_fasta(planted, background=[0.3, 0.2, 0.2, 0.3])
        assert given_background.attrs["background"] == (0.3, 0.2, 0.2, 0.3)

        over_floor = motifs.scan_sequences(
            wide(planted_records["plantedI"]), background=[0.3, 0.2, 0.2, 0.3]
        )
        assert over_floor.attrs["background"] == (0.3, 0.2, 0.2, 0.3)

    def test_the_recorded_background_reproduces_the_scan_exactly(
        self, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        # The whole point of recording it: handing the recorded value back must reproduce
        # the scan exactly, or two runs could not be reconciled from what they carry.
        peaks = wide(planted_records["plantedI"])
        derived = motifs.scan_sequences(peaks)
        replayed = motifs.scan_sequences(peaks, background=list(derived.attrs["background"]))
        pd.testing.assert_frame_equal(replayed, derived)
        assert replayed.attrs["background"] == derived.attrs["background"]

    def test_the_background_choice_changes_which_hits_are_found(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str]
    ) -> None:
        # Why this is automatic rather than an option nobody remembers to pass, and why it
        # is recorded rather than assumed: an AT-rich null and a uniform one do not agree.
        peaks = wide(planted_records["plantedI"])
        assert sites(motifs.scan_sequences(peaks)) != sites(
            motifs.scan_sequences(peaks, background="uniform")
        )
        uniform = sites(motifs.scan_fasta(planted))
        at_rich = sites(motifs.scan_fasta(planted, background=[0.35, 0.15, 0.15, 0.35]))
        assert uniform != at_rich

    @pytest.mark.parametrize(
        ("background", "message"),
        [
            ([0.25, 0.25, 0.25], "4 frequencies"),
            ([0.0, 0.25, 0.25, 0.5], "must all be > 0"),
            ([0.3, 0.3, 0.3, 0.3], "must sum to 1"),
        ],
    )
    def test_a_background_that_is_not_four_frequencies_is_refused(
        self, motifs: MotifSet, background: list[float], message: str
    ) -> None:
        with pytest.raises(ValueError, match=message):
            motifs.scan("ACGT" * 20, background=background)

    def test_an_unknown_background_mode_names_the_ones_there_are(self, motifs: MotifSet) -> None:
        with pytest.raises(ValueError, match="auto, uniform, derive"):
            motifs.scan("ACGT" * 20, background="gc")  # type: ignore[arg-type]

    def test_deciding_the_background_does_not_eat_the_first_records(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str]
    ) -> None:
        # A FASTA is read once. Records pulled off it while the background was being
        # decided must still be scanned, or a scan would silently skip its first records.
        found = motifs.scan_fasta(planted, background="derive")
        assert set(found["sequence_name"]) == set(planted_records)

    def test_release_identity_is_recorded_only_for_an_unfiltered_release(
        self, release: JasparDatabase, motifs: MotifSet, planted: Path
    ) -> None:
        from_release = release.scan_fasta(planted)
        assert (from_release.attrs["release"], from_release.attrs["tax_group"]) == (
            "2024",
            "all",
        )
        # filter() hands back a plain motif set, so a filtered release is no longer that
        # release and its table must not claim to be — and neither does a de novo set.
        filtered = release.filter(tax_group="vertebrates").scan_fasta(planted)
        assert (filtered.attrs["release"], filtered.attrs["tax_group"]) == (None, None)
        de_novo = motifs.scan_fasta(planted)
        assert (de_novo.attrs["release"], de_novo.attrs["tax_group"]) == (None, None)


# ---------------------------------------------------------------------------
# One table, whatever the scan was handed
# ---------------------------------------------------------------------------


class TestEntryPointsAgree:
    def test_the_single_sequence_form_names_its_sequence_by_default_or_by_request(
        self, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        found = motifs.scan("ACGT" * 40)
        assert DEFAULT_SEQUENCE_NAME == "sequence"
        assert set(found["sequence_name"]) <= {DEFAULT_SEQUENCE_NAME}

        named = motifs.scan(planted_records["plantedI"], "chosen")
        assert set(named["sequence_name"]) == {"chosen"}

    def test_the_mapping_fasta_stream_and_one_sequence_forms_agree_and_batches_stream_in_order(
        self,
        motifs: MotifSet,
        planted: Path,
        hits: pd.DataFrame,
        planted_records: dict[str, str],
    ) -> None:
        pd.testing.assert_frame_equal(
            motifs.scan_fasta(planted), motifs.scan_sequences(planted_records)
        )
        pd.testing.assert_frame_equal(
            scan_stream(motifs, planted_records.items()), motifs.scan_fasta(planted)
        )
        alone = motifs.scan(planted_records["plantedI"], "plantedI")
        assert rows(alone) == of(hits, sequence_name="plantedI")
        # And a DNA scans as the string it is.
        bases = planted_records["plantedI"]
        assert rows(motifs.scan(DNA(bases), "x")) == rows(motifs.scan(bases, "x"))

        # The loop a Parquet sink and a parallel source attach to: one batch per named
        # sequence, concatenated in the order the source yielded them.
        pieces = [("second", "ACGT" * 30), ("first", "TTTTGATTACAGTTTT")]
        found = scan_stream(MotifSet([word_motif("GATTACAG")]), pieces)
        assert list(dict.fromkeys(found["sequence_name"])) == ["first"]
        both = scan_stream(MotifSet([word_motif("GATTACAG")]), [pieces[1], pieces[1]])
        assert len(both) == 2 * len(found)

        # A caller sharding one sequence wants the pieces to share a name.
        sequence = "TTTTGATTACAGTTTT"
        shared_name = scan_stream(
            MotifSet([word_motif("GATTACAG")]), [("shard", sequence), ("shard", sequence)]
        )
        assert set(shared_name["sequence_name"]) == {"shard"}
        assert len(shared_name) == 2
