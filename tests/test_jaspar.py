"""Tests for genome.tf.motif.jaspar — the transfac parser and the JASPAR database.

Offline throughout: the ``fake_fetch`` fixture stands in for the package's one fetch step
and serves ``tests/data/tiny_jaspar_transfac.txt``, and the URL asserted here is the one
the package *built*, read back off the recorded call rather than off the network. The
autouse data-root fixture puts the cache under the test's own directory, so the layout is
exercised for real.

The fixture is ten real records, and everything ``tests/data/README.md`` says about them
is asserted here rather than trusted — the traps especially: values separated by a
semicolon and never by a comma, commas that live *inside* one value, an annotation the
source left empty, and counts that are not integers.

There is no **Completion marker** for a motif download, deliberately, so the two things
standing in its place get tests of their own: the rename into place, which is why an
interrupted download never occupies the final name, and the count check, which is why a
short file raises instead of scanning with half a release.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import itertools
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.io import fetch as fetch_mod
from genome.tf.motif import MIN_MOTIF_LENGTH, AmbiguousMotifNameError, Motif, MotifSet
from genome.tf.motif import jaspar as jaspar_mod
from genome.tf.motif.jaspar import (
    DEFAULT_RELEASE,
    DEFAULT_TAX_GROUP,
    JASPAR_RELEASES,
    JASPAR_TAX_GROUPS,
    MOTIF_COUNTS,
    JasparDatabase,
    JasparReleaseError,
    TransfacError,
    jaspar_data_dir,
    jaspar_filename,
    jaspar_url,
    parse_transfac,
)

from .conftest import FakeFetch

#: The committed fixture, cut from JASPAR 2024's `all` union file. See tests/data/README.md.
FIXTURE = "tiny_jaspar_transfac.txt"

#: What the fixture is: every record, in file order, with the four facts a reader of the
#: README needs to be able to check. The rules each one exists to break are asserted one
#: by one below, since that is what it is in the set for.
FIXTURE_MOTIFS: tuple[tuple[str, str, int, str], ...] = (
    ("MA0119.1", "NFIC::TLX1", 14, "vertebrates"),
    ("MA0789.1", "POU3F4", 9, "vertebrates"),
    ("MA0079.5", "SP1", 9, "vertebrates"),
    ("MA0139.2", "CTCF", 15, "vertebrates"),
    ("MA1929.2", "CTCF", 31, "vertebrates"),
    ("MA1930.2", "CTCF", 33, "vertebrates"),
    ("MA2355.1", "PK06791.1", 6, "plants"),
    ("MA0261.1", "lin-14", 6, "nematodes"),
    ("MA0283.1", "CHA4", 8, "fungi"),
    ("MA1407.2", "bZIP14", 8, "diatoms"),
)

#: How many motifs the fixture holds — the count every database test is held to, in place
#: of the release's real one.
FIXTURE_COUNT = len(FIXTURE_MOTIFS)


@pytest.fixture
def fixture_text(data_dir: Path) -> str:
    """The committed transfac fixture, as text."""
    return (data_dir / FIXTURE).read_text(encoding="utf-8")


@pytest.fixture
def fixture_motifs(fixture_text: str) -> tuple[Motif, ...]:
    """The committed fixture, parsed."""
    return parse_transfac(fixture_text)


@pytest.fixture
def jaspar_counts(monkeypatch: pytest.MonkeyPatch) -> None:
    """Expect the fixture's ten motifs from every release, so any of them can be served.

    The count check is what stands where a **Completion marker** stands elsewhere, so it
    is never switched off — only pointed at what the fake fetch actually serves. Every
    key of the real table survives, so a test may still name any release and tax group.
    """
    monkeypatch.setattr(
        jaspar_mod, "MOTIF_COUNTS", MappingProxyType(dict.fromkeys(MOTIF_COUNTS, FIXTURE_COUNT))
    )


@pytest.fixture
def served(fake_fetch: FakeFetch, jaspar_counts: None) -> FakeFetch:
    """A fetch step serving the transfac fixture, with the count check pointed at it."""
    fake_fetch.serve(FIXTURE)
    return fake_fetch


# ---------------------------------------------------------------------------
# The committed bytes: everything the README says about the fixture
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    def test_the_fixture_holds_the_records_the_readme_lists_in_that_order(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        read = tuple((m.motif_id, m.motif_name, len(m), m.tax_group) for m in fixture_motifs)
        assert read == FIXTURE_MOTIFS

    def test_every_record_carries_all_six_annotations_or_says_nothing(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        for motif in fixture_motifs:
            assert motif.tax_group  # the one annotation no record leaves empty
            assert motif.pubmed_ids
            assert motif.data_type

    def test_three_records_share_the_name_ctcf(self, fixture_motifs: tuple[Motif, ...]) -> None:
        names = [m.motif_name for m in fixture_motifs]
        assert names.count("CTCF") == 3

    def test_the_dimer_carries_two_classes_two_families_and_two_accessions(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        # MA0119.1 NFIC::TLX1 — the semicolon-separated case, and why tf_class is plural.
        dimer = MotifSet(fixture_motifs)["MA0119.1"]
        assert dimer.motif_name == "NFIC::TLX1"
        assert dimer.tf_class == (
            "SMAD/NF-1 DNA-binding domain factors",
            "Homeo domain factors",
        )
        assert dimer.tf_family == ("Nuclear factor 1", "NK")
        assert dimer.uniprot_ids == ("P08651", "P31314")

    def test_one_record_carries_two_pubmed_ids(self, fixture_motifs: tuple[Motif, ...]) -> None:
        assert MotifSet(fixture_motifs)["MA0789.1"].pubmed_ids == ("8876240", "2350782")

    def test_one_record_carries_fractional_counts(self, fixture_motifs: tuple[Motif, ...]) -> None:
        # What the transfac serialization keeps and `.jaspar` rounds away.
        sp1 = MotifSet(fixture_motifs)["MA0079.5"]
        assert sp1.counts[0, 0] == 1.05485
        assert not (sp1.counts == sp1.counts.astype(int)).all()

    def test_the_other_nine_records_hold_whole_counts(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        fractional = [
            m.motif_id for m in fixture_motifs if not (m.counts == m.counts.astype(int)).all()
        ]
        assert fractional == ["MA0079.5"]

    def test_a_comma_inside_one_class_is_part_of_it_and_never_a_separator(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        # Splitting on the comma would make two classes out of one, silently.
        assert MotifSet(fixture_motifs)["MA2355.1"].tf_class == (
            "C3H(C),C2HC zinc-fingers like factors",
        )

    def test_a_comma_inside_a_data_type_is_part_of_it_too(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        assert MotifSet(fixture_motifs)["MA0283.1"].data_type == "PBM, CSA and/or DIP-chip"

    @pytest.mark.parametrize(
        ("motif_id", "empty"),
        [
            ("MA2355.1", ("uniprot_ids",)),
            ("MA0261.1", ("tf_class", "tf_family")),
            ("MA0283.1", ("tf_family",)),
        ],
    )
    def test_an_annotation_the_source_left_blank_is_an_empty_tuple(
        self, fixture_motifs: tuple[Motif, ...], motif_id: str, empty: tuple[str, ...]
    ) -> None:
        motif = MotifSet(fixture_motifs)[motif_id]
        for field in empty:
            assert getattr(motif, field) == ()

    def test_two_records_are_below_the_minimum_scannable_length(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        short = [m.motif_id for m in fixture_motifs if len(m) < MIN_MOTIF_LENGTH]
        assert short == ["MA2355.1", "MA0261.1"]

    def test_the_long_record_has_the_least_informative_flanks(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        long = MotifSet(fixture_motifs)["MA1930.2"]
        assert len(long) == 33 == max(len(m) for m in fixture_motifs)
        bits = long.information_content
        assert bits[0] == pytest.approx(0.36, abs=0.005)
        assert bits[-1] == pytest.approx(0.31, abs=0.005)

    def test_trimming_the_long_record_takes_its_flanks_and_keeps_its_middle(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        # Real data for the rule that trimming acts only on the ends: this matrix carries
        # twelve interior positions under 0.25 bits, and every one of them survives.
        long = MotifSet(fixture_motifs)["MA1930.2"]
        trimmed = long.trim(0.4)
        assert (len(trimmed), trimmed.offset) == (30, 1)
        assert (trimmed.information_content < 0.25).sum() == 12

    def test_only_one_record_is_a_diatom(self, fixture_motifs: tuple[Motif, ...]) -> None:
        # The degenerate tax group: JASPAR really does publish one diatom matrix.
        diatoms = MotifSet(fixture_motifs).filter(tax_group="diatoms")
        assert diatoms.motif_ids == ("MA1407.2",)

    def test_five_tax_groups_are_represented(self, fixture_motifs: tuple[Motif, ...]) -> None:
        assert {m.tax_group for m in fixture_motifs} == {
            "vertebrates",
            "plants",
            "nematodes",
            "fungi",
            "diatoms",
        }


# ---------------------------------------------------------------------------
# The parser
# ---------------------------------------------------------------------------

ONE_RECORD = """AC MA0260.1
XX
ID che-1
XX
DE MA0260.1 che-1 ; From JASPAR
PO\tA\tC\tG\tT
01\t0.0\t0.0\t37.0\t0.0
02\t37.0\t0.0\t0.0\t0.0
XX
CC tax_group:nematodes
CC tf_family:More than 3 adjacent zinc fingers
CC tf_class:C2H2 zinc finger factors
CC pubmed_ids:17606643
CC uniprot_ids:Q966L8
CC data_type:COMPILED
XX
//
"""


class TestParseTransfac:
    def test_one_record_becomes_one_motif(self) -> None:
        (motif,) = parse_transfac(ONE_RECORD)
        assert (motif.motif_id, motif.motif_name) == ("MA0260.1", "che-1")
        assert motif.counts.shape == (4, 2)

    def test_the_matrix_is_four_rows_by_l_columns(self) -> None:
        (motif,) = parse_transfac(ONE_RECORD)
        # Column zero observed G 37 times; positions are columns, bases are rows.
        assert motif.counts[:, 0].tolist() == [0.0, 0.0, 37.0, 0.0]
        assert motif.consensus == "GA"

    def test_all_six_annotations_are_read(self) -> None:
        (motif,) = parse_transfac(ONE_RECORD)
        assert motif.tax_group == "nematodes"
        assert motif.tf_class == ("C2H2 zinc finger factors",)
        assert motif.tf_family == ("More than 3 adjacent zinc fingers",)
        assert motif.uniprot_ids == ("Q966L8",)
        assert motif.pubmed_ids == ("17606643",)
        assert motif.data_type == "COMPILED"

    def test_values_split_on_a_semicolon(self) -> None:
        text = ONE_RECORD.replace("CC tf_class:C2H2", "CC tf_class:Alpha; Beta ;C2H2")
        (motif,) = parse_transfac(text)
        assert motif.tf_class == ("Alpha", "Beta", "C2H2 zinc finger factors")

    def test_values_never_split_on_a_comma(self) -> None:
        text = ONE_RECORD.replace("CC tf_class:C2H2", "CC tf_class:Zinc finger, BED-type; C2H2")
        (motif,) = parse_transfac(text)
        assert motif.tf_class == ("Zinc finger, BED-type", "C2H2 zinc finger factors")

    def test_an_empty_value_means_the_source_stated_nothing(self) -> None:
        text = ONE_RECORD.replace("CC tf_family:More than 3 adjacent zinc fingers", "CC tf_family:")
        (motif,) = parse_transfac(text)
        assert motif.tf_family == ()

    def test_a_trailing_separator_yields_no_empty_value(self) -> None:
        text = ONE_RECORD.replace("CC uniprot_ids:Q966L8", "CC uniprot_ids:Q966L8; ")
        (motif,) = parse_transfac(text)
        assert motif.uniprot_ids == ("Q966L8",)

    def test_empty_text_yields_no_motifs(self) -> None:
        assert parse_transfac("") == ()
        assert parse_transfac("\n\n") == ()

    def test_it_reads_nothing_off_disk(self, fixture_text: str) -> None:
        # A pure function from text to motifs: no path, no release, no network.
        assert len(parse_transfac(fixture_text)) == FIXTURE_COUNT

    def test_a_record_without_an_accession_is_refused(self) -> None:
        text = ONE_RECORD.replace("AC MA0260.1\n", "")
        with pytest.raises(TransfacError, match="record 1"):
            parse_transfac(text)

    def test_a_record_without_a_matrix_is_refused_and_named(self) -> None:
        text = "AC MA0001.1\nXX\nID x\nXX\nCC tax_group:fungi\nXX\n//\n"
        with pytest.raises(TransfacError, match=r"MA0001\.1"):
            parse_transfac(text)

    def test_a_header_in_another_base_order_is_refused_rather_than_transposed(self) -> None:
        text = ONE_RECORD.replace("PO\tA\tC\tG\tT", "PO\tA\tG\tC\tT")
        with pytest.raises(TransfacError, match="another order"):
            parse_transfac(text)

    def test_a_count_row_that_is_not_numbers_is_refused(self) -> None:
        text = ONE_RECORD.replace("01\t0.0\t0.0\t37.0\t0.0", "01\t0.0\tx\t37.0\t0.0")
        with pytest.raises(TransfacError, match="four numbers"):
            parse_transfac(text)

    def test_a_final_record_missing_its_terminator_is_still_read(self) -> None:
        # Losing it silently is the truncation the count check exists to catch.
        assert len(parse_transfac(ONE_RECORD.replace("//\n", ""))) == 1


# ---------------------------------------------------------------------------
# Parser properties
# ---------------------------------------------------------------------------


def render(motif: Motif) -> str:
    """Render one motif back into a transfac record, as JASPAR serializes it."""
    rows = "\n".join(
        "\t".join([f"{position + 1:02d}", *(repr(float(count)) for count in column)])
        for position, column in enumerate(motif.counts.T)
    )
    annotations = "\n".join(
        [
            f"CC tax_group:{motif.tax_group}",
            f"CC tf_class:{'; '.join(motif.tf_class)}",
            f"CC tf_family:{'; '.join(motif.tf_family)}",
            f"CC pubmed_ids:{'; '.join(motif.pubmed_ids)}",
            f"CC uniprot_ids:{'; '.join(motif.uniprot_ids)}",
            f"CC data_type:{motif.data_type}",
        ]
    )
    return (
        f"AC {motif.motif_id}\nXX\nID {motif.motif_name}\nXX\n"
        f"DE {motif.motif_id} {motif.motif_name} ; From JASPAR\n"
        f"PO\tA\tC\tG\tT\n{rows}\nXX\n{annotations}\nXX\n//\n"
    )


#: Annotation text a transfac file could hold: printable, no separator, and nothing that
#: is only whitespace, since the line-based format cannot express any of those. Commas are
#: deliberately in the alphabet — a comma inside a value is exactly what must survive.
_value = (
    st.text(
        alphabet=st.characters(min_codepoint=32, max_codepoint=126, blacklist_characters=";"),
        min_size=1,
        max_size=20,
    )
    .map(str.strip)
    .filter(bool)
)

_counts = st.lists(
    st.lists(
        st.floats(min_value=0.0, max_value=1e5, allow_nan=False, allow_infinity=False),
        min_size=4,
        max_size=4,
    ).filter(lambda column: sum(column) > 0),
    min_size=1,
    max_size=8,
)


@st.composite
def transfac_motifs(draw: st.DrawFn) -> Motif:
    """A motif shaped like something a transfac file could carry."""
    counts = np.array(draw(_counts), dtype=float).T
    return Motif(
        draw(_value),
        draw(_value),
        counts,
        tax_group=draw(_value),
        tf_class=tuple(draw(st.lists(_value, max_size=3))),
        tf_family=tuple(draw(st.lists(_value, max_size=3))),
        uniprot_ids=tuple(draw(st.lists(_value, max_size=3))),
        pubmed_ids=tuple(draw(st.lists(_value, max_size=3))),
        data_type=draw(_value),
    )


class TestParserProperties:
    @given(motifs=st.lists(transfac_motifs(), min_size=1, max_size=4))
    def test_rendering_and_parsing_back_is_identity(self, motifs: list[Motif]) -> None:
        assert list(parse_transfac("".join(render(motif) for motif in motifs))) == motifs

    def test_parsing_is_total_over_the_fixture(self, fixture_motifs: tuple[Motif, ...]) -> None:
        assert len(fixture_motifs) == FIXTURE_COUNT

    @pytest.mark.parametrize("motif_id", [entry[0] for entry in FIXTURE_MOTIFS])
    def test_every_record_yields_four_rows_of_equal_length(
        self, fixture_motifs: tuple[Motif, ...], motif_id: str
    ) -> None:
        motif = MotifSet(fixture_motifs)[motif_id]
        assert motif.counts.ndim == 2
        assert motif.counts.shape[0] == 4
        assert len({len(row) for row in motif.counts}) == 1

    @pytest.mark.parametrize("motif_id", [entry[0] for entry in FIXTURE_MOTIFS])
    def test_every_column_sum_is_positive(
        self, fixture_motifs: tuple[Motif, ...], motif_id: str
    ) -> None:
        assert (MotifSet(fixture_motifs)[motif_id].counts.sum(axis=0) > 0).all()

    @pytest.mark.parametrize("motif_id", [entry[0] for entry in FIXTURE_MOTIFS])
    def test_probabilities_and_back_preserve_proportions(
        self, fixture_motifs: tuple[Motif, ...], motif_id: str
    ) -> None:
        motif = MotifSet(fixture_motifs)[motif_id]
        restored = motif.probabilities * motif.counts.sum(axis=0)
        assert restored == pytest.approx(motif.counts, rel=1e-9)


# ---------------------------------------------------------------------------
# URLs, file names and the arguments they take
# ---------------------------------------------------------------------------


class TestReleaseVocabulary:
    def test_the_eight_tax_groups(self) -> None:
        assert JASPAR_TAX_GROUPS == (
            "vertebrates",
            "plants",
            "insects",
            "nematodes",
            "fungi",
            "urochordates",
            "diatoms",
            "all",
        )

    def test_two_releases_and_the_newer_is_the_default(self) -> None:
        assert JASPAR_RELEASES == ("2024", "2026")
        assert DEFAULT_RELEASE == "2026"

    def test_vertebrates_is_the_default_tax_group(self) -> None:
        assert DEFAULT_TAX_GROUP == "vertebrates"

    def test_a_count_is_held_for_every_release_and_tax_group(self) -> None:
        assert set(MOTIF_COUNTS) == set(itertools.product(JASPAR_RELEASES, JASPAR_TAX_GROUPS))

    @pytest.mark.parametrize("release", JASPAR_RELEASES)
    def test_all_holds_exactly_the_union_of_the_other_seven(self, release: str) -> None:
        named = [group for group in JASPAR_TAX_GROUPS if group != "all"]
        assert MOTIF_COUNTS[release, "all"] == sum(MOTIF_COUNTS[release, g] for g in named)

    def test_diatoms_really_does_hold_one_motif(self) -> None:
        assert MOTIF_COUNTS["2024", "diatoms"] == MOTIF_COUNTS["2026", "diatoms"] == 1

    def test_an_unknown_release_names_the_ones_there_are(self) -> None:
        with pytest.raises(ValueError, match="no JASPAR release '2020'"):
            jaspar_url("2020", "vertebrates")

    def test_an_unknown_tax_group_names_the_ones_there_are(self) -> None:
        with pytest.raises(ValueError, match="no JASPAR tax group 'mammals'"):
            jaspar_url("2024", "mammals")


class TestUrlAndFilename:
    def test_a_named_taxon(self) -> None:
        assert jaspar_url("2024", "vertebrates") == (
            "https://jaspar.elixir.no/download/data/2024/CORE/"
            "JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt"
        )

    def test_the_union_file_drops_the_taxon_segment(self) -> None:
        assert jaspar_url("2026", "all") == (
            "https://jaspar.elixir.no/download/data/2026/CORE/"
            "JASPAR2026_CORE_non-redundant_pfms_transfac.txt"
        )

    def test_the_cached_name_carries_the_release_and_the_tax_group(self) -> None:
        assert jaspar_filename("2024", "vertebrates") == (
            "JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt"
        )
        # Unlike the published union file, whose own name says neither.
        assert jaspar_filename("2026", "all") == (
            "JASPAR2026_CORE_all_non-redundant_pfms_transfac.txt"
        )

    def test_every_cached_name_is_distinct_so_the_flat_cache_works(self) -> None:
        names = {
            jaspar_filename(release, tax_group)
            for release, tax_group in itertools.product(JASPAR_RELEASES, JASPAR_TAX_GROUPS)
        }
        assert len(names) == len(JASPAR_RELEASES) * len(JASPAR_TAX_GROUPS)


# ---------------------------------------------------------------------------
# The database: downloading, caching, and what it knows about itself
# ---------------------------------------------------------------------------


class TestJasparDatabaseDownload:
    @pytest.mark.parametrize("tax_group", JASPAR_TAX_GROUPS)
    @pytest.mark.parametrize("release", JASPAR_RELEASES)
    def test_the_url_it_asked_for(self, served: FakeFetch, release: str, tax_group: str) -> None:
        # Read off the recorded fetch call, which is the only thing the package built.
        JasparDatabase(release, tax_group)
        segment = "" if tax_group == "all" else f"{tax_group}_"
        assert served.last.url == (
            f"https://jaspar.elixir.no/download/data/{release}/CORE/"
            f"JASPAR{release}_CORE_{segment}non-redundant_pfms_transfac.txt"
        )

    def test_the_file_lands_where_the_layout_puts_it(
        self, served: FakeFetch, liulab_data: Path
    ) -> None:
        database = JasparDatabase("2024", "nematodes")
        expected = (
            liulab_data
            / "motif"
            / "jaspar"
            / "JASPAR2024_CORE_nematodes_non-redundant_pfms_transfac.txt"
        )
        assert database.path == expected
        assert expected.is_file()
        assert jaspar_data_dir() == expected.parent

    def test_motif_data_is_a_sibling_of_the_assembly_tree(
        self, served: FakeFetch, liulab_data: Path
    ) -> None:
        JasparDatabase("2024", "nematodes")
        assert (liulab_data / "motif").is_dir()
        assert not (liulab_data / "genome").exists()

    def test_the_cache_is_flat(self, served: FakeFetch) -> None:
        JasparDatabase("2024", "nematodes")
        JasparDatabase("2026", "all")
        cached = sorted(path.name for path in jaspar_data_dir().iterdir())
        assert cached == [
            "JASPAR2024_CORE_nematodes_non-redundant_pfms_transfac.txt",
            "JASPAR2026_CORE_all_non-redundant_pfms_transfac.txt",
        ]

    def test_a_second_construction_fetches_nothing(self, served: FakeFetch) -> None:
        first = JasparDatabase("2024", "vertebrates")
        second = JasparDatabase("2024", "vertebrates")
        assert len(served.calls) == 1
        assert first.path == second.path
        assert first.motif_ids == second.motif_ids

    def test_another_tax_group_is_another_file(self, served: FakeFetch) -> None:
        JasparDatabase("2024", "vertebrates")
        JasparDatabase("2024", "plants")
        assert len(served.calls) == 2

    def test_an_explicit_cache_dir_overrides_the_layout(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        elsewhere = tmp_path / "somewhere-else"
        database = JasparDatabase("2024", "fungi", cache_dir=elsewhere)
        assert database.path.parent == elsewhere

    def test_an_interrupted_download_never_occupies_the_final_name(
        self, monkeypatch: pytest.MonkeyPatch, jaspar_counts: None, data_dir: Path
    ) -> None:
        # Half a file arrives and the fetch then dies. The final name must be untouched:
        # this is what stands in for a completion record.
        half = (data_dir / FIXTURE).read_text(encoding="utf-8")[:2000]

        def die(url: str, dest_dir: Path, **kwargs: Any) -> Path:
            dest_dir.mkdir(parents=True, exist_ok=True)
            (dest_dir / str(kwargs["fname"])).write_text(half)
            raise ConnectionError("the network went away")

        monkeypatch.setattr(fetch_mod, "fetch_url", die)
        with pytest.raises(ConnectionError):
            JasparDatabase("2024", "vertebrates")
        final = jaspar_data_dir() / jaspar_filename("2024", "vertebrates")
        assert not final.exists()

    def test_what_an_interrupted_download_left_behind_is_never_adopted(
        self, monkeypatch: pytest.MonkeyPatch, jaspar_counts: None, data_dir: Path
    ) -> None:
        # The half file is still in the working area when the next construction starts;
        # it must be swept rather than picked up as though it had finished.
        half = (data_dir / FIXTURE).read_text(encoding="utf-8")[:2000]
        work = jaspar_data_dir() / ".work"
        work.mkdir(parents=True)
        part = work / f"{jaspar_filename('2024', 'vertebrates')}.part"
        part.write_text(half)

        fake = FakeFetch(FIXTURE)
        monkeypatch.setattr(fetch_mod, "fetch_url", fake)
        database = JasparDatabase("2024", "vertebrates")
        assert len(database) == FIXTURE_COUNT
        assert not part.exists()

    def test_the_working_area_is_gone_once_the_file_is_in_place(self, served: FakeFetch) -> None:
        JasparDatabase("2024", "vertebrates")
        assert not (jaspar_data_dir() / ".work").exists()

    def test_the_progress_bar_is_the_callers_choice(self, served: FakeFetch) -> None:
        JasparDatabase("2024", "vertebrates", progressbar=False)
        assert served.last.progressbar is False


class TestJasparDatabaseIdentity:
    def test_it_knows_which_release_it_is(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "nematodes")
        assert (database.release, database.tax_group) == ("2024", "nematodes")
        assert database.source_url == jaspar_url("2024", "nematodes")

    def test_the_defaults_are_the_newest_release_and_vertebrates(self, served: FakeFetch) -> None:
        database = JasparDatabase()
        assert (database.release, database.tax_group) == (DEFAULT_RELEASE, DEFAULT_TAX_GROUP)

    def test_repr_names_the_release_and_the_size(self, served: FakeFetch) -> None:
        assert repr(JasparDatabase("2024", "all")) == (
            f"JasparDatabase(release='2024', tax_group='all', motifs={FIXTURE_COUNT})"
        )

    def test_an_unknown_release_raises_before_anything_is_fetched(self, served: FakeFetch) -> None:
        with pytest.raises(ValueError, match="no JASPAR release"):
            JasparDatabase("2019")
        assert served.calls == []


class TestJasparDatabaseQueries:
    def test_it_is_a_motif_set(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        assert isinstance(database, MotifSet)
        assert len(database) == FIXTURE_COUNT

    def test_indexing_by_id_base_id_and_unique_name(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        assert database["MA0139.2"].motif_name == "CTCF"
        assert database["MA0139"].motif_id == "MA0139.2"
        assert database["lin-14"].motif_id == "MA0261.1"

    def test_an_ambiguous_name_names_all_three_ctcfs(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        with pytest.raises(AmbiguousMotifNameError) as raised:
            database["CTCF"]
        assert raised.value.motif_ids == ("MA0139.2", "MA1929.2", "MA1930.2")

    def test_name_lookup_returns_a_tuple_either_way(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        assert len(database.by_name("CTCF")) == 3
        assert len(database.by_name("lin-14")) == 1

    def test_filtering_hands_back_a_plain_motif_set_and_not_a_release(
        self, served: FakeFetch
    ) -> None:
        # A filtered release is no longer that release, so it must not claim to be one.
        database = JasparDatabase("2024", "all")
        filtered = database.filter(tax_group="vertebrates")
        assert type(filtered) is MotifSet
        assert not isinstance(filtered, JasparDatabase)
        assert not hasattr(filtered, "release")

    def test_a_filtered_release_still_does_everything_a_set_does(self, served: FakeFetch) -> None:
        filtered = JasparDatabase("2024", "all").filter(tf_class="zinc finger")
        assert filtered["MA0079.5"].motif_name == "SP1"
        assert len(filtered.by_name("CTCF")) == 3

    def test_a_de_novo_set_answers_every_question_a_release_does(self, served: FakeFetch) -> None:
        # What the container abstraction is for: matrices JASPAR never published get the
        # whole API, and neither call below knows which kind of set it is holding.
        database = JasparDatabase("2024", "all")
        de_novo = MotifSet(
            [Motif(f"pattern_{index}", "CTCF-like", np.ones((4, 8))) for index in range(3)]
        )
        for motifs in (database, de_novo):
            assert isinstance(motifs.filter(lambda motif: len(motif) == 8), MotifSet)
            assert isinstance(motifs.motif_ids, tuple)
        assert de_novo["pattern_1"].motif_name == "CTCF-like"
        assert len(de_novo.by_name("CTCF-like")) == 3


class TestJasparDatabaseIntegrity:
    def _cache(self, text: str, release: str = "2024", tax_group: str = "all") -> Path:
        """Put ``text`` where a download would have left it, so the cache path is read."""
        path = jaspar_data_dir() / jaspar_filename(release, tax_group)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def test_a_truncated_file_raises_rather_than_yielding_a_partial_release(
        self, jaspar_counts: None, fixture_text: str
    ) -> None:
        short = "".join(part + "//\n" for part in fixture_text.split("//\n")[:3])
        path = self._cache(short)
        with pytest.raises(JasparReleaseError) as raised:
            JasparDatabase("2024", "all")
        message = str(raised.value)
        assert "holds 3 motifs" in message
        assert f"has {FIXTURE_COUNT}" in message
        assert str(path) in message

    def test_the_count_is_checked_on_a_cache_read_and_not_only_on_a_download(
        self, jaspar_counts: None, fixture_text: str
    ) -> None:
        # Nothing was fetched here at all: the file was already on disk.
        self._cache(fixture_text + fixture_text.split("//\n")[0] + "//\n")
        with pytest.raises(JasparReleaseError, match="holds 11 motifs"):
            JasparDatabase("2024", "all")

    def test_the_real_count_is_what_an_unpatched_read_is_held_to(self, fixture_text: str) -> None:
        # Ten records where the 2024 union file has 2346: the constant is not decoration.
        self._cache(fixture_text)
        with pytest.raises(JasparReleaseError, match="has 2346"):
            JasparDatabase("2024", "all")

    def test_two_versions_of_one_matrix_are_refused(
        self, jaspar_counts: None, fixture_text: str
    ) -> None:
        # A non-redundant release ships one version of each, which is what makes a bare
        # base id address one motif — asserted rather than assumed.
        self._cache(fixture_text.replace("AC MA1929.2", "AC MA0139.1"))
        with pytest.raises(JasparReleaseError, match="two versions of the matrix MA0139"):
            JasparDatabase("2024", "all")

    def test_the_base_ids_of_a_good_file_are_all_distinct(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        bases = [motif_id.split(".")[0] for motif_id in database.motif_ids]
        assert len(set(bases)) == len(bases)

    def test_a_bad_record_raises_a_parse_error_and_not_a_release_error(
        self, jaspar_counts: None, fixture_text: str
    ) -> None:
        self._cache(fixture_text.replace("PO\tA\tC\tG\tT", "PO\tA\tG\tC\tT", 1))
        with pytest.raises(TransfacError):
            JasparDatabase("2024", "all")
