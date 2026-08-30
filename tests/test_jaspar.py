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

A motif download writes a **Completion marker** like every other **Prepared set**, and
what it buys over the atomic rename it replaced gets tests of its own: a release whose
record went missing reads as unfinished rather than as finished, and one rewritten
afterwards is refused by the digest the record holds it to. The count check is tested
beside them, because it says something no digest of ours can — that this is the right
release, whole.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from types import MappingProxyType
from typing import Any, Protocol

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.io import fetch as fetch_mod
from genome.io.completion import RECORD_NAME, UnfinishedRegistrationError, read_record
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
    MotifSetNotDownloadedError,
    TransfacError,
    jaspar_data_dir,
    jaspar_filename,
    jaspar_prepare_command,
    jaspar_set_dir,
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


class PublishText(Protocol):
    """Makes JASPAR publish some text, whatever is in it, for one test."""

    def __call__(self, text: str) -> None:
        """Serve ``text`` as the release's published file from now on."""
        ...


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
    def test_the_fixture_matches_what_the_readme_says_about_it(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        read = tuple((m.motif_id, m.motif_name, len(m), m.tax_group) for m in fixture_motifs)
        assert read == FIXTURE_MOTIFS
        for record in fixture_motifs:
            assert record.tax_group  # the one annotation no record leaves empty
            assert record.pubmed_ids
            assert record.data_type
        names = [m.motif_name for m in fixture_motifs]
        assert names.count("CTCF") == 3
        assert {m.tax_group for m in fixture_motifs} == {
            "vertebrates",
            "plants",
            "nematodes",
            "fungi",
            "diatoms",
        }
        # Two records are below the minimum scannable length, and only one is a diatom —
        # JASPAR really does publish just the one diatom matrix.
        short = [m.motif_id for m in fixture_motifs if len(m) < MIN_MOTIF_LENGTH]
        assert short == ["MA2355.1", "MA0261.1"]
        diatoms = MotifSet(fixture_motifs).filter(tax_group="diatoms")
        assert diatoms.motif_ids == ("MA1407.2",)

    def test_semicolon_splits_a_value_comma_does_not_and_a_blank_is_empty(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        by_id = MotifSet(fixture_motifs)
        # MA0119.1 NFIC::TLX1 — the semicolon-separated case, and why tf_class is plural.
        dimer = by_id["MA0119.1"]
        assert dimer.tf_class == (
            "SMAD/NF-1 DNA-binding domain factors",
            "Homeo domain factors",
        )
        assert dimer.tf_family == ("Nuclear factor 1", "NK")
        assert dimer.uniprot_ids == ("P08651", "P31314")
        assert by_id["MA0789.1"].pubmed_ids == ("8876240", "2350782")
        # Splitting on the comma would make two classes out of one, silently.
        assert by_id["MA2355.1"].tf_class == ("C3H(C),C2HC zinc-fingers like factors",)
        assert by_id["MA0283.1"].data_type == "PBM, CSA and/or DIP-chip"
        # And an annotation the source left blank is an empty tuple, not a stray value.
        assert by_id["MA2355.1"].uniprot_ids == ()
        assert (by_id["MA0261.1"].tf_class, by_id["MA0261.1"].tf_family) == ((), ())
        assert by_id["MA0283.1"].tf_family == ()

    def test_counts_and_information_content_are_kept_as_given(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        # What the transfac serialization keeps and `.jaspar` rounds away.
        by_id = MotifSet(fixture_motifs)
        sp1 = by_id["MA0079.5"]
        assert sp1.counts[0, 0] == 1.05485
        fractional = [
            m.motif_id for m in fixture_motifs if not (m.counts == m.counts.astype(int)).all()
        ]
        assert fractional == ["MA0079.5"]

        # The longest record has the least informative flanks, and trimming acts only on
        # the ends: this matrix carries twelve interior positions under 0.25 bits, and
        # every one of them survives.
        long = by_id["MA1930.2"]
        assert len(long) == 33 == max(len(m) for m in fixture_motifs)
        bits = long.information_content
        assert bits[0] == pytest.approx(0.36, abs=0.005)
        assert bits[-1] == pytest.approx(0.31, abs=0.005)
        trimmed = long.trim(0.4)
        assert (len(trimmed), trimmed.offset) == (30, 1)
        assert (trimmed.information_content < 0.25).sum() == 12


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
    def test_one_record_becomes_a_motif_with_a_4_by_l_matrix(self) -> None:
        (record,) = parse_transfac(ONE_RECORD)
        assert (record.motif_id, record.motif_name) == ("MA0260.1", "che-1")
        assert record.counts.shape == (4, 2)
        # Column zero observed G 37 times; positions are columns, bases are rows.
        assert record.counts[:, 0].tolist() == [0.0, 0.0, 37.0, 0.0]
        assert record.consensus == "GA"

    def test_all_six_annotations_are_read(self) -> None:
        (record,) = parse_transfac(ONE_RECORD)
        assert record.tax_group == "nematodes"
        assert record.tf_class == ("C2H2 zinc finger factors",)
        assert record.tf_family == ("More than 3 adjacent zinc fingers",)
        assert record.uniprot_ids == ("Q966L8",)
        assert record.pubmed_ids == ("17606643",)
        assert record.data_type == "COMPILED"

    def test_values_split_on_semicolon_but_never_on_a_comma(self) -> None:
        split = ONE_RECORD.replace("CC tf_class:C2H2", "CC tf_class:Alpha; Beta ;C2H2")
        (record,) = parse_transfac(split)
        assert record.tf_class == ("Alpha", "Beta", "C2H2 zinc finger factors")

        kept_whole = ONE_RECORD.replace(
            "CC tf_class:C2H2", "CC tf_class:Zinc finger, BED-type; C2H2"
        )
        (record,) = parse_transfac(kept_whole)
        assert record.tf_class == ("Zinc finger, BED-type", "C2H2 zinc finger factors")

    def test_empty_and_trailing_separators_never_yield_a_stray_value(self) -> None:
        blanked = ONE_RECORD.replace(
            "CC tf_family:More than 3 adjacent zinc fingers", "CC tf_family:"
        )
        (record,) = parse_transfac(blanked)
        assert record.tf_family == ()

        trailing = ONE_RECORD.replace("CC uniprot_ids:Q966L8", "CC uniprot_ids:Q966L8; ")
        (record,) = parse_transfac(trailing)
        assert record.uniprot_ids == ("Q966L8",)

    def test_empty_or_terminator_less_text_is_read_rather_than_lost(
        self, fixture_text: str
    ) -> None:
        assert parse_transfac("") == ()
        assert parse_transfac("\n\n") == ()
        # A pure function from text to motifs: no path, no release, no network.
        assert len(parse_transfac(fixture_text)) == FIXTURE_COUNT
        # Losing a final record silently is the truncation the count check exists to
        # catch, so a record missing only its terminator must still be read.
        assert len(parse_transfac(ONE_RECORD.replace("//\n", ""))) == 1

    def test_malformed_records_are_refused_and_named(self) -> None:
        without_accession = ONE_RECORD.replace("AC MA0260.1\n", "")
        with pytest.raises(TransfacError, match="record 1"):
            parse_transfac(without_accession)

        without_matrix = "AC MA0001.1\nXX\nID x\nXX\nCC tax_group:fungi\nXX\n//\n"
        with pytest.raises(TransfacError, match=r"MA0001\.1"):
            parse_transfac(without_matrix)

        transposed_header = ONE_RECORD.replace("PO\tA\tC\tG\tT", "PO\tA\tG\tC\tT")
        with pytest.raises(TransfacError, match="another order"):
            parse_transfac(transposed_header)

        non_numeric_row = ONE_RECORD.replace("01\t0.0\t0.0\t37.0\t0.0", "01\t0.0\tx\t37.0\t0.0")
        with pytest.raises(TransfacError, match="four numbers"):
            parse_transfac(non_numeric_row)


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

    def test_every_fixture_record_has_a_well_formed_matrix_and_preserves_proportions(
        self, fixture_motifs: tuple[Motif, ...]
    ) -> None:
        assert len(fixture_motifs) == FIXTURE_COUNT
        for record in fixture_motifs:
            assert record.counts.ndim == 2
            assert record.counts.shape[0] == 4
            assert len({len(row) for row in record.counts}) == 1
            assert (record.counts.sum(axis=0) > 0).all()
            restored = record.probabilities * record.counts.sum(axis=0)
            assert restored == pytest.approx(record.counts, rel=1e-9)


# ---------------------------------------------------------------------------
# URLs, file names and the arguments they take
# ---------------------------------------------------------------------------


class TestReleaseVocabulary:
    def test_the_vocabulary_shape_tax_groups_releases_and_counts(self) -> None:
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
        assert JASPAR_RELEASES == ("2024", "2026")
        assert set(MOTIF_COUNTS) == set(itertools.product(JASPAR_RELEASES, JASPAR_TAX_GROUPS))
        for release in JASPAR_RELEASES:
            named = [group for group in JASPAR_TAX_GROUPS if group != "all"]
            assert MOTIF_COUNTS[release, "all"] == sum(MOTIF_COUNTS[release, g] for g in named)
        # The degenerate tax group: JASPAR really does publish one diatom matrix, in
        # both releases.
        assert MOTIF_COUNTS["2024", "diatoms"] == MOTIF_COUNTS["2026", "diatoms"] == 1
        assert DEFAULT_RELEASE == "2026"
        assert DEFAULT_TAX_GROUP == "vertebrates"

    def test_an_unknown_release_or_tax_group_names_the_ones_there_are(self) -> None:
        with pytest.raises(ValueError, match="no JASPAR release '2020'"):
            jaspar_url("2020", "vertebrates")
        with pytest.raises(ValueError, match="no JASPAR tax group 'mammals'"):
            jaspar_url("2024", "mammals")


class TestUrlAndFilename:
    def test_the_url_and_cached_name_carry_release_and_tax_group_and_every_name_is_distinct(
        self,
    ) -> None:
        assert jaspar_url("2024", "vertebrates") == (
            "https://jaspar.elixir.no/download/data/2024/CORE/"
            "JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt"
        )
        # "all" is the one tax group whose URL and filename drop the taxon segment.
        assert jaspar_url("2026", "all") == (
            "https://jaspar.elixir.no/download/data/2026/CORE/"
            "JASPAR2026_CORE_non-redundant_pfms_transfac.txt"
        )
        assert jaspar_filename("2024", "vertebrates") == (
            "JASPAR2024_CORE_vertebrates_non-redundant_pfms_transfac.txt"
        )
        # Unlike the published union file, whose own name says neither.
        assert jaspar_filename("2026", "all") == (
            "JASPAR2026_CORE_all_non-redundant_pfms_transfac.txt"
        )
        names = {
            jaspar_filename(release, tax_group)
            for release, tax_group in itertools.product(JASPAR_RELEASES, JASPAR_TAX_GROUPS)
        }
        assert len(names) == len(JASPAR_RELEASES) * len(JASPAR_TAX_GROUPS)


# ---------------------------------------------------------------------------
# The database: downloading, caching, and what it knows about itself
# ---------------------------------------------------------------------------


class TestJasparDatabaseDownload:
    @pytest.mark.parametrize(
        ("release", "tax_group"), [("2024", "vertebrates"), ("2026", "all"), ("2024", "diatoms")]
    )
    def test_the_url_it_asked_for(self, served: FakeFetch, release: str, tax_group: str) -> None:
        # Read off the recorded fetch call, which is the only thing the package built.
        # "all" is the one tax group whose URL drops the taxon segment.
        JasparDatabase(release, tax_group)
        segment = "" if tax_group == "all" else f"{tax_group}_"
        assert served.last.url == (
            f"https://jaspar.elixir.no/download/data/{release}/CORE/"
            f"JASPAR{release}_CORE_{segment}non-redundant_pfms_transfac.txt"
        )

    def test_the_layout_is_one_set_per_directory_beside_the_assembly_tree(
        self, served: FakeFetch, liulab_data: Path
    ) -> None:
        # One directory per release and tax group, because each carries a marker of its
        # own — the shape an Xref set and a Homology set are prepared in.
        database = JasparDatabase("2024", "nematodes")
        expected = (
            liulab_data
            / "motif"
            / "jaspar"
            / "2024"
            / "nematodes"
            / "JASPAR2024_CORE_nematodes_non-redundant_pfms_transfac.txt"
        )
        assert database.path == expected
        assert expected.is_file()
        assert jaspar_set_dir("2024", "nematodes") == expected.parent
        assert jaspar_data_dir() == liulab_data / "motif" / "jaspar"
        assert not (liulab_data / "genome").exists()
        # The marker is written last, beside the file it claims, and it is what says the
        # release is finished — not the file's mere existence.
        record = read_record(expected.parent)
        assert record is not None
        assert record.kind == "motif"
        assert record.source_url == jaspar_url("2024", "nematodes")
        assert record.details["release"] == "2024"
        assert record.details["tax_group"] == "nematodes"
        assert record.files == {expected.name: expected.stat().st_size}
        # Nothing JASPAR published to pin, so what is recorded is the digest of what was
        # stored — which is what a re-read is held to.
        assert record.sha256 == hashlib.sha256(expected.read_bytes()).hexdigest()
        assert not (expected.parent / ".work").exists()

        JasparDatabase("2026", "all")
        prepared = sorted(
            str(path.relative_to(jaspar_data_dir())) for path in jaspar_data_dir().glob("*/*/*.txt")
        )
        assert prepared == [
            "2024/nematodes/JASPAR2024_CORE_nematodes_non-redundant_pfms_transfac.txt",
            "2026/all/JASPAR2026_CORE_all_non-redundant_pfms_transfac.txt",
        ]

    def test_fetching_happens_once_per_release_and_tax_group(self, served: FakeFetch) -> None:
        first = JasparDatabase("2024", "vertebrates")
        second = JasparDatabase("2024", "vertebrates")
        assert len(served.calls) == 1
        assert first.path == second.path
        assert first.motif_ids == second.motif_ids

        JasparDatabase("2024", "plants")
        assert len(served.calls) == 2

    def test_an_explicit_cache_dir_overrides_the_layout_and_the_progress_bar_is_the_callers_choice(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        elsewhere = tmp_path / "somewhere-else"
        database = JasparDatabase("2024", "fungi", cache_dir=elsewhere, progressbar=False)
        assert database.path.parent == elsewhere
        assert served.last.progressbar is False

    def test_an_interrupted_download_leaves_no_release_and_names_the_login_node(
        self, monkeypatch: pytest.MonkeyPatch, jaspar_counts: None, data_dir: Path
    ) -> None:
        # Half a file arrives and the fetch then dies. Nothing is placed and no marker is
        # written, so the release reads as absent rather than as a short one — and what a
        # compute node with no internet gets is the call to make on a login node, where it
        # used to get pooch's own transport error.
        half = (data_dir / FIXTURE).read_text(encoding="utf-8")[:2000]

        def die(url: str, dest_dir: Path, **kwargs: Any) -> Path:
            dest_dir.mkdir(parents=True, exist_ok=True)
            (dest_dir / str(kwargs["fname"])).write_text(half)
            raise ConnectionError("the network went away")

        monkeypatch.setattr(fetch_mod, "fetch_url", die)
        with pytest.raises(MotifSetNotDownloadedError) as raised:
            JasparDatabase("2024", "vertebrates")
        assert "login node" in str(raised.value)
        assert jaspar_prepare_command("2024", "vertebrates") in str(raised.value)

        directory = jaspar_set_dir("2024", "vertebrates")
        assert not (directory / jaspar_filename("2024", "vertebrates")).exists()
        assert read_record(directory) is None

    def test_what_an_interrupted_download_left_behind_is_never_adopted(
        self, monkeypatch: pytest.MonkeyPatch, jaspar_counts: None, data_dir: Path
    ) -> None:
        # The half file is still in the working area when the next construction starts.
        # JASPAR pins no checksum, so nothing could vouch for it — and pooch serves a file
        # already sitting at the destination — which is why the working area is swept
        # before the fetch rather than picked up as though it had finished.
        half = (data_dir / FIXTURE).read_text(encoding="utf-8")[:2000]
        directory = jaspar_set_dir("2024", "vertebrates")
        work = directory / ".work"
        work.mkdir(parents=True)
        stale = work / jaspar_filename("2024", "vertebrates")
        stale.write_text(half)

        fake = FakeFetch(FIXTURE)
        monkeypatch.setattr(fetch_mod, "fetch_url", fake)
        database = JasparDatabase("2024", "vertebrates")
        assert len(database) == FIXTURE_COUNT
        assert not stale.exists()
        # And the working area itself is gone once a normal download is in place too.
        assert not work.exists()

    def test_a_release_left_unfinished_reads_as_unfinished_and_the_repair_names_the_rebuild(
        self, served: FakeFetch
    ) -> None:
        # What the marker buys over an atomic rename: a file whose record went missing is
        # not silently re-read as a finished release, and the message names both halves of
        # the repair rather than leaving the caller with nothing and no way back.
        database = JasparDatabase("2024", "nematodes")
        (database.path.parent / RECORD_NAME).unlink()

        with pytest.raises(UnfinishedRegistrationError, match="rm -rf") as raised:
            JasparDatabase("2024", "nematodes")
        assert jaspar_prepare_command("2024", "nematodes") in str(raised.value)

    def test_a_release_rewritten_after_it_was_prepared_is_refused_by_its_own_digest(
        self, served: FakeFetch, fixture_text: str
    ) -> None:
        # The same ten records with one byte changed: the motif count still accepts this
        # file and the recorded digest does not, which is the whole of what the marker
        # adds. JASPAR publishes no checksum, so the digest is of what was stored.
        database = JasparDatabase("2024", "nematodes")
        edited = fixture_text.replace("ID lin-14", "ID lin-15")
        assert edited != fixture_text
        database.path.write_text(edited, encoding="utf-8")
        marker = database.path.parent / RECORD_NAME
        payload = json.loads(marker.read_text(encoding="utf-8"))
        payload["files"][database.path.name] = database.path.stat().st_size
        marker.write_text(json.dumps(payload), encoding="utf-8")

        with pytest.raises(JasparReleaseError, match="hashes to"):
            JasparDatabase("2024", "nematodes")


class TestJasparDatabaseIdentity:
    def test_identity_reflects_the_release_and_tax_group_explicit_or_default(
        self, served: FakeFetch
    ) -> None:
        database = JasparDatabase("2024", "nematodes")
        assert (database.release, database.tax_group) == ("2024", "nematodes")
        assert database.source_url == jaspar_url("2024", "nematodes")

        defaulted = JasparDatabase()
        assert (defaulted.release, defaulted.tax_group) == (DEFAULT_RELEASE, DEFAULT_TAX_GROUP)

        assert repr(JasparDatabase("2024", "all")) == (
            f"JasparDatabase(release='2024', tax_group='all', motifs={FIXTURE_COUNT})"
        )

    def test_an_unknown_release_raises_before_anything_is_fetched(self, served: FakeFetch) -> None:
        with pytest.raises(ValueError, match="no JASPAR release"):
            JasparDatabase("2019")
        assert served.calls == []


class TestJasparDatabaseQueries:
    def test_it_behaves_as_a_motif_set_for_indexing_and_name_lookup(
        self, served: FakeFetch
    ) -> None:
        database = JasparDatabase("2024", "all")
        assert isinstance(database, MotifSet)
        assert len(database) == FIXTURE_COUNT
        assert database["MA0139.2"].motif_name == "CTCF"
        assert database["MA0139"].motif_id == "MA0139.2"
        assert database["lin-14"].motif_id == "MA0261.1"
        assert len(database.by_name("CTCF")) == 3
        assert len(database.by_name("lin-14")) == 1

    def test_an_ambiguous_name_names_all_three_ctcfs(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        with pytest.raises(AmbiguousMotifNameError) as raised:
            database["CTCF"]
        assert raised.value.motif_ids == ("MA0139.2", "MA1929.2", "MA1930.2")

    def test_filtering_drops_the_release_identity_but_keeps_the_set_behavior(
        self, served: FakeFetch
    ) -> None:
        # A filtered release is no longer that release, so it must not claim to be one.
        database = JasparDatabase("2024", "all")
        filtered = database.filter(tax_group="vertebrates")
        assert type(filtered) is MotifSet
        assert not isinstance(filtered, JasparDatabase)
        assert not hasattr(filtered, "release")
        # But it still does everything a set does.
        zinc_fingers = database.filter(tf_class="zinc finger")
        assert zinc_fingers["MA0079.5"].motif_name == "SP1"
        assert len(zinc_fingers.by_name("CTCF")) == 3

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
    @pytest.fixture
    def publish(self, fake_fetch: FakeFetch, tmp_path: Path) -> PublishText:
        """Return a helper that makes JASPAR publish ``text``, whatever is in it.

        Through the fetch step rather than written into the set's directory by hand: the
        release is then prepared exactly as a real one is, marker and all, so what these
        tests exercise is a *prepared* release being read — which is where a truncated
        publisher file would actually be caught.
        """

        def publish(text: str) -> None:
            path = tmp_path / "published_transfac.txt"
            path.write_text(text, encoding="utf-8")
            fake_fetch.serve(path)

        return publish

    def test_a_wrong_motif_count_raises_rather_than_yielding_a_partial_release(
        self, publish: PublishText, jaspar_counts: None, fixture_text: str
    ) -> None:
        publish("".join(part + "//\n" for part in fixture_text.split("//\n")[:3]))
        with pytest.raises(JasparReleaseError) as raised:
            JasparDatabase("2024", "all")
        message = str(raised.value)
        assert "holds 3 motifs" in message
        assert f"has {FIXTURE_COUNT}" in message
        assert str(jaspar_set_dir("2024", "all")) in message

        # Checked on every read and not only on the download: the release is prepared and
        # recorded, and constructing it again reads what is there and still refuses it.
        with pytest.raises(JasparReleaseError, match="holds 3 motifs"):
            JasparDatabase("2024", "all")

    def test_the_real_count_is_what_an_unpatched_read_is_held_to(
        self, publish: PublishText, fixture_text: str
    ) -> None:
        # Ten records where the 2024 union file has 2346: the constant is not decoration.
        publish(fixture_text)
        with pytest.raises(JasparReleaseError, match="has 2346"):
            JasparDatabase("2024", "all")

    def test_two_versions_of_one_matrix_are_refused(
        self, publish: PublishText, jaspar_counts: None, fixture_text: str
    ) -> None:
        # A non-redundant release ships one version of each, which is what makes a bare
        # base id address one motif — asserted rather than assumed.
        publish(fixture_text.replace("AC MA1929.2", "AC MA0139.1"))
        with pytest.raises(JasparReleaseError, match="two versions of the matrix MA0139"):
            JasparDatabase("2024", "all")

    def test_the_base_ids_of_a_good_file_are_all_distinct(self, served: FakeFetch) -> None:
        database = JasparDatabase("2024", "all")
        bases = [motif_id.split(".")[0] for motif_id in database.motif_ids]
        assert len(set(bases)) == len(bases)

    def test_a_bad_record_raises_a_parse_error_and_not_a_release_error(
        self, publish: PublishText, jaspar_counts: None, fixture_text: str
    ) -> None:
        publish(fixture_text.replace("PO\tA\tC\tG\tT", "PO\tA\tG\tC\tT", 1))
        with pytest.raises(TransfacError):
            JasparDatabase("2024", "all")
