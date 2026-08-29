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

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.seq import DNA
from genome.tf.motif import (
    DEFAULT_SEQUENCE_NAME,
    DEFAULT_THRESHOLD,
    HIT_COLUMNS,
    HIT_DTYPES,
    HIT_PROVENANCE,
    MIN_MOTIF_LENGTH,
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
    def test_the_fixture_holds_two_records_of_six_hundred_bases(
        self, planted_records: dict[str, str]
    ) -> None:
        assert list(planted_records) == ["plantedI", "plantedII"]
        assert [len(bases) for bases in planted_records.values()] == [600, 600]

    def test_one_header_carries_a_description_after_its_name(self, planted: Path) -> None:
        # What the whitespace truncation is here to be tested against.
        headers = [
            line for line in planted.read_text(encoding="utf-8").splitlines() if line[:1] == ">"
        ]
        assert headers == [
            ">plantedI",
            ">plantedII  sacCer3 chrII:1-600, bases 180-240 soft-masked",
        ]

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

    @pytest.mark.parametrize(("record", "start", "end", "strand", "motif_id", "word"), PLANTED)
    def test_each_planted_word_sits_at_the_offset_the_readme_names(
        self,
        planted_records: dict[str, str],
        motifs: MotifSet,
        record: str,
        start: int,
        end: int,
        strand: str,
        motif_id: str,
        word: str,
    ) -> None:
        assert planted_records[record][start:end] == word
        expected = str(motifs[motif_id].consensus)
        assert word.upper() == (expected if strand == "+" else revcomp(expected))

    def test_the_reverse_site_is_the_forward_one_flipped(
        self, planted_records: dict[str, str]
    ) -> None:
        forward, reverse = PLANTED[0], PLANTED[1]
        bases = planted_records["plantedI"]
        assert revcomp(bases[reverse[1] : reverse[2]]) == bases[forward[1] : forward[2]]

    def test_one_window_of_plantedII_is_soft_masked_and_nothing_else_is(  # noqa: N802
        self, planted_records: dict[str, str]
    ) -> None:
        low, high = MASKED
        bases = planted_records["plantedII"]
        assert bases[low:high].islower()
        assert bases[:low].isupper()
        assert bases[high:].isupper()
        assert planted_records["plantedI"].isupper()

    def test_the_masked_window_holds_the_third_planted_site(
        self, planted_records: dict[str, str]
    ) -> None:
        _record, start, end, _strand, _motif_id, _word = PLANTED[2]
        assert MASKED[0] <= start < end <= MASKED[1]

    def test_the_fixture_holds_only_the_four_bases(self, planted_records: dict[str, str]) -> None:
        for bases in planted_records.values():
            assert DNA.outside_alphabet(bases) == []

    def test_the_fixture_is_wrapped_at_sixty(self, planted: Path) -> None:
        widths = {
            len(line)
            for line in planted.read_text(encoding="utf-8").splitlines()
            if not line.startswith(">")
        }
        assert widths == {WRAP}


# ---------------------------------------------------------------------------
# Reading a FASTA
# ---------------------------------------------------------------------------


class TestReadFasta:
    def test_a_record_name_stops_at_the_first_whitespace(self, planted: Path) -> None:
        # What STAR and chromap write into an alignment from the same file, so a hit table
        # joins against it without anyone renaming anything.
        assert [name for name, _ in read_fasta(planted)] == ["plantedI", "plantedII"]

    def test_the_bases_come_back_joined_and_in_the_case_they_were_written(
        self, planted: Path, planted_records: dict[str, str]
    ) -> None:
        assert dict(read_fasta(planted)) == planted_records

    def test_a_gzipped_fasta_reads_the_same(self, planted: Path, tmp_path: Path) -> None:
        zipped = tmp_path / "planted.fa.gz"
        with gzip.open(zipped, "wt", encoding="utf-8") as handle:
            handle.write(planted.read_text(encoding="utf-8"))
        assert dict(read_fasta(zipped)) == dict(read_fasta(planted))

    def test_a_missing_file_names_itself(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match="FASTA file not found"):
            list(read_fasta(tmp_path / "nope.fa"))

    def test_bases_before_any_header_are_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "bare.fa"
        path.write_text("ACGTACGT\n")
        with pytest.raises(FastaFormatError, match="line 1"):
            list(read_fasta(path))

    def test_a_header_with_no_name_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "nameless.fa"
        path.write_text(">chrI\nACGT\n>\nACGT\n")
        with pytest.raises(FastaFormatError, match="no name"):
            list(read_fasta(path))

    def test_an_empty_file_yields_no_records(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.fa"
        path.write_text("")
        assert list(read_fasta(path)) == []


# ---------------------------------------------------------------------------
# The schema, which is the contract
# ---------------------------------------------------------------------------


class TestSchema:
    def test_the_columns_are_the_seven_in_that_order(self, hits: pd.DataFrame) -> None:
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

    def test_the_dtypes_are_the_compact_ones(self, hits: pd.DataFrame) -> None:
        assert {name: str(dtype) for name, dtype in hits.dtypes.items()} == dict(HIT_DTYPES)

    def test_a_scan_that_found_nothing_has_the_same_schema(self, motifs: MotifSet) -> None:
        found = motifs.scan("N" * 300)
        assert len(found) == 0
        assert list(found.columns) == list(HIT_COLUMNS)
        assert {name: str(dtype) for name, dtype in found.dtypes.items()} == dict(HIT_DTYPES)

    def test_an_empty_motif_set_scans_to_an_empty_table(self, planted: Path) -> None:
        found = MotifSet([]).scan_fasta(planted)
        assert len(found) == 0
        assert found.attrs["motifs_scanned"] == ()

    def test_the_empty_table_is_the_schema_with_no_rows(self) -> None:
        assert list(empty_hits().columns) == list(HIT_COLUMNS)
        assert {name: str(dtype) for name, dtype in empty_hits().dtypes.items()} == dict(HIT_DTYPES)

    def test_the_index_is_a_fresh_range_across_records(self, hits: pd.DataFrame) -> None:
        assert list(hits.index) == list(range(len(hits)))


# ---------------------------------------------------------------------------
# The planted sites, and the frame every coordinate is in
# ---------------------------------------------------------------------------


class TestPlantedSites:
    @pytest.mark.parametrize(("record", "start", "end", "strand", "motif_id", "word"), PLANTED)
    def test_the_planted_site_is_found_where_it_was_planted(
        self,
        hits: pd.DataFrame,
        record: str,
        start: int,
        end: int,
        strand: str,
        motif_id: str,
        word: str,
    ) -> None:
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

    def test_every_interval_is_as_long_as_the_motif_that_made_it(
        self, hits: pd.DataFrame, motifs: MotifSet
    ) -> None:
        for motif_id, _name, _sequence, start, end, _strand, _score in rows(hits):
            assert int(end) - int(start) == len(motifs[str(motif_id)])

    def test_every_interval_lies_inside_its_record(
        self, hits: pd.DataFrame, planted_records: dict[str, str]
    ) -> None:
        for _id, _name, sequence, start, end, _strand, _score in rows(hits):
            assert 0 <= int(start) < int(end) <= len(planted_records[str(sequence)])

    def test_every_strand_is_plus_or_minus_and_never_unknown(self, hits: pd.DataFrame) -> None:
        assert set(hits["strand"]) <= {"+", "-"}
        assert "." not in set(hits["strand"])
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

    def test_a_hit_is_never_reported_for_a_motif_that_was_not_scanned(
        self, hits: pd.DataFrame
    ) -> None:
        assert set(hits["motif_id"]) <= set(hits.attrs["motifs_scanned"])


# ---------------------------------------------------------------------------
# What the engine is trusted for, checked against the engine
# ---------------------------------------------------------------------------


class TestForwardFrame:
    @pytest.mark.parametrize("strand", ["+", "-"])
    def test_a_word_planted_at_a_known_offset_is_found_at_that_offset(self, strand: str) -> None:
        word = "GATTACAGTC"
        filler = "AAACCCTTT" * 5
        planted = word if strand == "+" else revcomp(word)
        found = MotifSet([word_motif(word)]).scan(filler + planted + filler)
        assert (
            "MA9999.1",
            DEFAULT_SEQUENCE_NAME,
            len(filler),
            len(filler) + len(word),
            strand,
        ) in sites(found)

    def test_the_last_window_of_a_sequence_is_scanned(self) -> None:
        # The engine was chosen partly because it scans it; a site at the final position
        # must not need a base appended after it to be found.
        word = "GATTACAGTC"
        filler = "AAACCCTTT" * 5
        found = MotifSet([word_motif(word)]).scan(filler + word)
        assert (int(found["end"].max())) == len(filler) + len(word)

    def test_a_sequence_shorter_than_the_motif_yields_no_hit(self) -> None:
        found = MotifSet([word_motif("GATTACAGTC")]).scan("ACGT")
        assert len(found) == 0

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
    def test_the_committed_masked_record_scans_as_its_upper_case_equivalent(
        self, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        masked = planted_records["plantedII"]
        assert not masked.isupper()  # the fixture really does carry masking
        assert rows(motifs.scan(masked, "x")) == rows(motifs.scan(masked.upper(), "x"))

    def test_the_masked_planted_site_is_found_despite_its_case(
        self, hits: pd.DataFrame, planted_records: dict[str, str]
    ) -> None:
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

    @given(sequence=st.text(alphabet="ACGT", max_size=150))
    def test_masking_changes_no_hit_whatever_the_sequence(self, sequence: str) -> None:
        motifs = MotifSet([word_motif("GATTACAGTC")])
        assert rows(motifs.scan(sequence)) == rows(motifs.scan(sequence.lower()))


# ---------------------------------------------------------------------------
# The minimum length, and the motifs it leaves out
# ---------------------------------------------------------------------------


class TestSkippedMotifs:
    def test_the_two_short_records_are_never_scanned(self, hits: pd.DataFrame) -> None:
        assert set(hits["motif_id"]).isdisjoint(TOO_SHORT)

    def test_they_are_named_on_the_result(self, hits: pd.DataFrame) -> None:
        assert hits.attrs["motifs_skipped"] == TOO_SHORT

    def test_the_other_eight_were_scanned(self, hits: pd.DataFrame, motifs: MotifSet) -> None:
        assert hits.attrs["motifs_scanned"] == tuple(
            motif.motif_id for motif in motifs if motif.motif_id not in TOO_SHORT
        )

    def test_scanned_and_skipped_together_are_the_whole_set(
        self, hits: pd.DataFrame, motifs: MotifSet
    ) -> None:
        assert hits.attrs["motifs_scanned"] + hits.attrs["motifs_skipped"] != ()
        assert set(hits.attrs["motifs_scanned"]) | set(hits.attrs["motifs_skipped"]) == set(
            motifs.motif_ids
        )

    def test_the_boundary_is_seven_positions(self, planted: Path) -> None:
        short = word_motif("GATTAC", "MA0001.1")  # six, and unreachable at any threshold
        long = word_motif("GATTACA", "MA0002.1")  # seven, the shortest that is scannable
        assert (len(short), len(long)) == (MIN_MOTIF_LENGTH - 1, MIN_MOTIF_LENGTH)
        found = MotifSet([short, long]).scan_fasta(planted)
        assert found.attrs["motifs_skipped"] == ("MA0001.1",)
        assert found.attrs["motifs_scanned"] == ("MA0002.1",)

    def test_a_set_of_only_short_motifs_scans_to_an_empty_table(self, planted: Path) -> None:
        found = MotifSet([word_motif("GATTAC", "MA0001.1")]).scan_fasta(planted)
        assert len(found) == 0
        assert found.attrs["motifs_skipped"] == ("MA0001.1",)


# ---------------------------------------------------------------------------
# The threshold, and what the score column is not
# ---------------------------------------------------------------------------


class TestThreshold:
    def test_the_default_is_a_per_position_p_value(self, hits: pd.DataFrame) -> None:
        assert DEFAULT_THRESHOLD == 1e-4
        assert hits.attrs["threshold"] == DEFAULT_THRESHOLD

    def test_a_stricter_threshold_keeps_a_subset(self, motifs: MotifSet, planted: Path) -> None:
        loose = sites(motifs.scan_fasta(planted, threshold=1e-4))
        strict = sites(motifs.scan_fasta(planted, threshold=1e-6))
        assert strict < loose

    def test_a_looser_threshold_keeps_a_superset(self, motifs: MotifSet, planted: Path) -> None:
        assert sites(motifs.scan_fasta(planted)) < sites(motifs.scan_fasta(planted, threshold=1e-3))

    @pytest.mark.parametrize("threshold", [0.0, 1.0, -1e-4, 2.0])
    def test_a_threshold_outside_zero_to_one_says_it_is_a_p_value(
        self, motifs: MotifSet, threshold: float
    ) -> None:
        with pytest.raises(ValueError, match="per-position p-value"):
            motifs.scan("ACGT" * 20, threshold=threshold)


# ---------------------------------------------------------------------------
# The background, and the provenance the table carries
# ---------------------------------------------------------------------------


class TestBackgroundAndProvenance:
    def test_every_provenance_key_is_present(self, hits: pd.DataFrame) -> None:
        assert set(hits.attrs) == set(HIT_PROVENANCE)

    def test_the_background_defaults_to_uniform_and_is_recorded(self, hits: pd.DataFrame) -> None:
        assert hits.attrs["background"] == (0.25, 0.25, 0.25, 0.25)

    def test_a_background_given_is_the_one_recorded(self, motifs: MotifSet, planted: Path) -> None:
        found = motifs.scan_fasta(planted, background=[0.3, 0.2, 0.2, 0.3])
        assert found.attrs["background"] == (0.3, 0.2, 0.2, 0.3)

    def test_the_background_decides_the_answer(self, motifs: MotifSet, planted: Path) -> None:
        # Recorded rather than assumed because it moves the hits: an AT-rich null and a
        # uniform one do not agree on this fixture.
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

    def test_a_release_records_which_release_it_is(
        self, release: JasparDatabase, planted: Path
    ) -> None:
        found = release.scan_fasta(planted)
        assert (found.attrs["release"], found.attrs["tax_group"]) == ("2024", "all")

    def test_a_filtered_release_records_neither(
        self, release: JasparDatabase, planted: Path
    ) -> None:
        # filter() hands back a plain motif set, so a filtered release is no longer that
        # release and its table must not claim to be.
        found = release.filter(tax_group="vertebrates").scan_fasta(planted)
        assert (found.attrs["release"], found.attrs["tax_group"]) == (None, None)

    def test_a_de_novo_set_records_neither(self, motifs: MotifSet, planted: Path) -> None:
        found = motifs.scan_fasta(planted)
        assert (found.attrs["release"], found.attrs["tax_group"]) == (None, None)

    def test_an_empty_result_still_carries_its_provenance(self, motifs: MotifSet) -> None:
        found = motifs.scan("N" * 200)
        assert len(found) == 0
        assert set(found.attrs) == set(HIT_PROVENANCE)
        assert found.attrs["motifs_skipped"] == TOO_SHORT


# ---------------------------------------------------------------------------
# One table, whatever the scan was handed
# ---------------------------------------------------------------------------


class TestEntryPointsAgree:
    def test_the_single_sequence_form_names_its_sequence(self, motifs: MotifSet) -> None:
        found = motifs.scan("ACGT" * 40)
        assert DEFAULT_SEQUENCE_NAME == "sequence"
        assert set(found["sequence_name"]) <= {DEFAULT_SEQUENCE_NAME}

    def test_the_single_sequence_form_takes_a_name(
        self, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        found = motifs.scan(planted_records["plantedI"], "chosen")
        assert set(found["sequence_name"]) == {"chosen"}

    def test_a_mapping_and_a_fasta_agree(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str]
    ) -> None:
        pd.testing.assert_frame_equal(
            motifs.scan_fasta(planted), motifs.scan_sequences(planted_records)
        )

    def test_one_sequence_agrees_with_the_file_it_came_from(
        self, motifs: MotifSet, hits: pd.DataFrame, planted_records: dict[str, str]
    ) -> None:
        alone = motifs.scan(planted_records["plantedI"], "plantedI")
        assert rows(alone) == of(hits, sequence_name="plantedI")

    def test_a_dna_scans_as_the_string_it_is(
        self, motifs: MotifSet, planted_records: dict[str, str]
    ) -> None:
        bases = planted_records["plantedI"]
        assert rows(motifs.scan(DNA(bases), "x")) == rows(motifs.scan(bases, "x"))

    def test_the_stream_form_is_the_one_the_three_share(
        self, motifs: MotifSet, planted: Path, planted_records: dict[str, str]
    ) -> None:
        pd.testing.assert_frame_equal(
            scan_stream(motifs, planted_records.items()), motifs.scan_fasta(planted)
        )

    def test_the_batches_are_drained_in_the_order_they_arrive(self, motifs: MotifSet) -> None:
        # The loop a Parquet sink and a parallel source attach to: one batch per named
        # sequence, concatenated in the order the source yielded them.
        pieces = [("second", "ACGT" * 30), ("first", "TTTTGATTACAGTTTT")]
        found = scan_stream(MotifSet([word_motif("GATTACAG")]), pieces)
        assert list(dict.fromkeys(found["sequence_name"])) == ["first"]
        both = scan_stream(MotifSet([word_motif("GATTACAG")]), [pieces[1], pieces[1]])
        assert len(both) == 2 * len(found)

    def test_the_same_name_twice_is_the_callers_business(self, motifs: MotifSet) -> None:
        # A caller sharding one sequence wants the pieces to share a name.
        sequence = "TTTTGATTACAGTTTT"
        found = scan_stream(
            MotifSet([word_motif("GATTACAG")]), [("shard", sequence), ("shard", sequence)]
        )
        assert set(found["sequence_name"]) == {"shard"}
        assert len(found) == 2
