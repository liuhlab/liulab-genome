"""Tests for genome.homology — the Compara slice, the partition guard and the answer.

Offline throughout: the ``fake_fetch`` fixture stands in for the package's one fetch step
and serves the committed Compara subsamples under ``tests/data/homology``, and the URL and
the checksum asserted here are the ones the package *built and demanded*, read back off
the recorded call rather than off the network. The autouse data-root fixture puts the set
under the test's own directory, so the layout is exercised for real.

The fixtures are whole published rows of three real release-116 dumps, and everything
``tests/data/README.md`` says about them is asserted here rather than trusted — the traps
especially: the human dump really does hold **zero** human↔mouse rows, so the partition
guard is exercised against the publisher's own partition and not a staged one; the
paralogy rows really are same-species, so the boundary a **Homology link** draws is held
against real bytes; and both quality scores really are ``NULL`` on every row of either
worm pairing.

**Two things this suite deliberately does not assert.** That a cross-species **Paralogy
link** exists — release 116 publishes none, counted over the whole human dump, so a test
claiming one would claim something about the publisher that is not true. And any number
this package computed: every value checked is a cell of the published file.

Nothing reaches into a set's index: what a caller can see is the answer, the path, the
record and the file, and the stems a test asks about are read out of the committed
fixture rather than out of the object under test.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import ast
import gzip
import json
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from genome.homology import (
    COMPARA_COLUMNS,
    DEFAULT_RELEASE,
    METADATA_COLUMNS,
    QUALITY_SCORE_COLUMNS,
    ComparaFileError,
    ComparaPartitionError,
    HomologyMetadata,
    HomologyMetadataError,
    HomologySet,
    NoHomologyPairError,
    UnknownHomologySpeciesError,
    VersionedGeneIdError,
    compara_url,
    homology_data_dir,
    homology_metadata,
    homology_releases,
    homology_species,
    homology_table,
    read_metadata,
    resolve_homologs,
)
from genome.io.completion import RECORD_NAME, UnfinishedRegistrationError, read_record
from genome.io.gtf import AnnotationRegistry
from genome.io.registration import AssemblyDir

from .conftest import FakeFetch

#: Which committed subsample stands in for which species' published dump.
FIXTURES = {
    "Homo sapiens": "homology/compara116.homo_sapiens.homologies.tsv.gz",
    "Mus musculus": "homology/compara116.mus_musculus.homologies.tsv.gz",
    "Caenorhabditis elegans": "homology/compara116.caenorhabditis_elegans.homologies.tsv.gz",
}

#: Compara's own name for each species, which for these three is the species slug.
SLUGS = {
    "Homo sapiens": "homo_sapiens",
    "Mus musculus": "mus_musculus",
    "Caenorhabditis elegans": "caenorhabditis_elegans",
}

#: What each fixture holds, as ``tests/data/README.md`` states it: total rows, the rows of
#: each cross-species pair in it, and how many same-species paralogy rows ride along.
FIXTURE_SHAPE: dict[str, tuple[int, dict[str, int], int]] = {
    "Homo sapiens": (17, {"caenorhabditis_elegans": 11, "mus_musculus": 0}, 4),
    "Mus musculus": (27, {"homo_sapiens": 10, "caenorhabditis_elegans": 11}, 4),
    "Caenorhabditis elegans": (4, {"homo_sapiens": 0, "mus_musculus": 0}, 4),
}

#: The three pairings among the lab's species, and how many links each set holds in the
#: fixtures. Every one of them must answer — an acceptance criterion of its own.
PAIRS: tuple[tuple[str, str, int], ...] = (
    ("Homo sapiens", "Mus musculus", 10),
    ("Homo sapiens", "Caenorhabditis elegans", 11),
    ("Mus musculus", "Caenorhabditis elegans", 11),
)

#: A human gene the fixture gives three worm orthologs, all ``ortholog_one2many``.
ONE2MANY_HUMAN = "ENSG00000152670"
ONE2MANY_WORMS = ("WBGene00001598", "WBGene00001599", "WBGene00001600")

#: Two human genes the fixture gives the *same* two worm orthologs, ``many2many`` both.
MANY2MANY_HUMAN = ("ENSG00000163331", "ENSG00000112977")
MANY2MANY_WORMS = ("WBGene00010019", "WBGene00012140")

#: A human gene with exactly one worm ortholog, ``ortholog_one2one`` and high confidence.
ONE2ONE_HUMAN, ONE2ONE_WORM = "ENSG00000177479", "WBGene00020462"

#: A mouse gene the fixture gives two human orthologs, ``ortholog_one2many``, with real
#: quality scores on both — the pair Compara does score.
ONE2MANY_MOUSE = "ENSMUSG00000074698"
ONE2MANY_HUMANS = ("ENSG00000101266", "ENSG00000254598")

#: A mouse gene with two human orthologs the publisher calls ``ortholog_many2many``, which
#: has exactly as many partners as the one2many gene above and a different label.
MANY2MANY_MOUSE = "ENSMUSG00000070378"

#: A human gene with two mouse orthologs, ``ortholog_many2many``.
MANY2MANY_HUMAN_OF_MOUSE = "ENSG00000172150"
MANY2MANY_MICE = ("ENSMUSG00000070377", "ENSMUSG00000070378")

#: A stem no fixture mentions, for the unresolved bucket.
ABSENT = "ENSG00000000000"


def _row(species: str, other: str) -> HomologyMetadata:
    """The shipped provenance row for one pair, which every test builds through."""
    found = homology_metadata(species, other, DEFAULT_RELEASE)
    assert found is not None
    return found


def _fixture_rows(data_dir: Path, species: str) -> list[list[str]]:
    """Every row of one committed subsample, split into cells, header dropped."""
    with gzip.open(data_dir / FIXTURES[species], "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t") for line in handle.read().splitlines()[1:]]


def _published(data_dir: Path, species: str, other: str) -> list[list[str]]:
    """The pair's own rows out of whichever fixture the shipped table names for it."""
    wanted = {SLUGS[species], SLUGS[other]}
    return [
        row
        for row in _fixture_rows(data_dir, _row(species, other).holding_species)
        if {row[2], row[7]} == wanted and row[2] != row[7]
    ]


def _stems(data_dir: Path, species: str, other: str) -> list[str]:
    """Every ``species`` gene id stem the fixtures give a homolog in ``other``, sorted."""
    at = 0 if _row(species, other).holding_species == species else 5
    return sorted({row[at] for row in _published(data_dir, species, other)})


def _built(
    fake_fetch: FakeFetch,
    species: str,
    other: str,
    *,
    cache_dir: str | Path | None = None,
) -> HomologySet:
    """Build one set, serving the subsample of whichever dump the shipped table names.

    Which file holds a pair is exactly what the shipped table records, so a test serves
    what the package asked for rather than deciding for itself — the tests that serve the
    *wrong* file do so on purpose and say why.
    """
    fake_fetch.serve(FIXTURES[_row(species, other).holding_species])
    return HomologySet(species, other, progressbar=False, cache_dir=cache_dir)


def _registry_declaring(tmp_path: Path, *gene_ids: str) -> AnnotationRegistry:
    """Register a GTF declaring ``gene_ids`` under ``mine`` and open its registry.

    The annotation hop is :meth:`AnnotationRegistry.resolve_gene_ids` used unchanged, so a
    test of the crossing needs a real registered annotation and nothing else.
    """
    directory = tmp_path / "annotation"
    directory.mkdir(parents=True, exist_ok=True)
    source = directory / "mine.gtf"
    source.write_text(
        "".join(
            f'chrI\ttest\t{feature}\t1\t100\t.\t+\t.\tgene_id "{gene_id}"; '
            f'transcript_id "{gene_id}_t";\n'
            for gene_id in gene_ids
            for feature in ("gene", "transcript", "exon")
        )
    )
    AnnotationRegistry(AssemblyDir(assembly="tiny", path=directory)).register_path(
        source, "mine", disable_infer_genes=True, disable_infer_transcripts=True
    )
    return AnnotationRegistry.locate("tiny", directory)


# ---------------------------------------------------------------------------
# The committed bytes: everything the README says about the fixtures
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    @pytest.mark.parametrize("species", sorted(FIXTURES))
    def test_each_fixture_leads_with_comparas_own_header(
        self, data_dir: Path, species: str
    ) -> None:
        with gzip.open(data_dir / FIXTURES[species], "rt", encoding="utf-8") as handle:
            assert tuple(handle.readline().rstrip("\n").split("\t")) == COMPARA_COLUMNS

    @pytest.mark.parametrize("species", sorted(FIXTURES))
    def test_each_fixture_holds_the_rows_the_readme_counts(
        self, data_dir: Path, species: str
    ) -> None:
        total, pairs, paralogy = FIXTURE_SHAPE[species]
        rows = _fixture_rows(data_dir, species)
        assert len(rows) == total
        for partner, count in pairs.items():
            here = [row for row in rows if row[7] == partner and row[2] != row[7]]
            assert len(here) == count, partner
        assert len([row for row in rows if row[2] == row[7]]) == paralogy

    def test_the_human_dump_holds_no_mouse_row_which_is_the_published_partition(
        self, data_dir: Path
    ) -> None:
        # Not a staged absence: the real release-116 human file holds 0 human/mouse rows
        # and 23,982 human/worm rows, and the mouse file holds the human pair. The guard
        # below is therefore exercised against the publisher's own partition.
        rows = _fixture_rows(data_dir, "Homo sapiens")
        assert [row for row in rows if row[7] == "mus_musculus"] == []
        assert len([row for row in rows if row[7] == "caenorhabditis_elegans"]) == 11

    @pytest.mark.parametrize("species", sorted(FIXTURES))
    def test_every_paralogy_row_relates_two_genes_of_one_species(
        self, data_dir: Path, species: str
    ) -> None:
        # Release 116 publishes no `between_species_paralog` anywhere in the human dump —
        # 4.0M rows over some 200 partner species — so a duplication label always means
        # one species, and these rows are here to hold that boundary.
        for row in _fixture_rows(data_dir, species):
            if not row[4].startswith("ortholog_"):
                assert row[2] == row[7], row

    def test_the_worm_dump_holds_no_row_of_either_pair(self, data_dir: Path) -> None:
        for row in _fixture_rows(data_dir, "Caenorhabditis elegans"):
            assert {row[2], row[7]} == {"caenorhabditis_elegans"}

    def test_both_quality_scores_are_null_on_every_worm_row_and_set_on_the_mouse_pair(
        self, data_dir: Path
    ) -> None:
        worm = [
            row
            for species in ("Homo sapiens", "Mus musculus")
            for row in _fixture_rows(data_dir, species)
            if row[7] == "caenorhabditis_elegans"
        ]
        assert len(worm) == 22
        assert {(row[11], row[12]) for row in worm} == {("NULL", "NULL")}
        scored = [
            row for row in _fixture_rows(data_dir, "Mus musculus") if row[7] == "homo_sapiens"
        ]
        assert all(row[11] != "NULL" and row[12] != "NULL" for row in scored)

    @pytest.mark.parametrize(("species", "other", "_links"), PAIRS)
    def test_every_pair_carries_all_three_cardinalities(
        self, data_dir: Path, species: str, other: str, _links: int
    ) -> None:
        assert {row[4] for row in _published(data_dir, species, other)} == {
            "ortholog_one2one",
            "ortholog_one2many",
            "ortholog_many2many",
        }

    def test_the_high_confidence_flag_keeps_all_three_of_its_states(self, data_dir: Path) -> None:
        rows = _fixture_rows(data_dir, "Homo sapiens")
        assert {row[13] for row in rows if row[2] == row[7]} == {"NULL"}
        assert {row[13] for row in rows if row[7] == "caenorhabditis_elegans"} == {"0", "1"}


# ---------------------------------------------------------------------------
# The shipped provenance table
# ---------------------------------------------------------------------------


class TestShippedTable:
    def test_the_table_pins_one_release_and_the_labs_three_species(self) -> None:
        assert homology_releases() == (DEFAULT_RELEASE,)
        assert homology_species() == ("Caenorhabditis elegans", "Homo sapiens", "Mus musculus")

    def test_every_row_pins_a_publisher_a_url_and_a_real_checksum(self) -> None:
        for row in homology_table():
            assert row.publisher == "Ensembl Compara"
            assert row.source_url.startswith("https://ftp.ensembl.org/pub/")
            assert len(row.md5) == 32
            assert set(row.md5) <= set("0123456789abcdef")
            assert row.pubmed_id > 0

    def test_every_rows_url_is_its_holding_species_own_dump_for_that_release(self) -> None:
        # The URL is written down rather than formatted — release 113 ships these dumps
        # uncompressed, so no template covers every release — and this holds the written
        # row and the message builder to one file for the release actually pinned.
        for row in homology_table():
            assert row.source_url == compara_url(row.holding_species, row.release)

    def test_the_three_pairings_are_pinned_and_two_of_them_share_one_file(self) -> None:
        held = {row.pair: row.holding_species for row in homology_table()}
        assert held == {
            ("Caenorhabditis elegans", "Homo sapiens"): "Homo sapiens",
            ("Caenorhabditis elegans", "Mus musculus"): "Mus musculus",
            ("Homo sapiens", "Mus musculus"): "Mus musculus",
        }

    def test_the_pair_is_unordered_so_either_spelling_finds_one_row(self) -> None:
        assert _row("Homo sapiens", "Mus musculus") is _row("Mus musculus", "Homo sapiens")

    def test_a_pair_the_table_does_not_pin_answers_none_rather_than_guessing(self) -> None:
        assert homology_metadata("Homo sapiens", "Danio rerio", DEFAULT_RELEASE) is None

    def test_a_blank_cell_raises_naming_the_column(self) -> None:
        header = "\t".join(METADATA_COLUMNS)
        cells = ["116", "Homo sapiens", "Mus musculus", "Mus musculus", "Ensembl", "1", "u", ""]

        with pytest.raises(HomologyMetadataError, match="'md5'"):
            read_metadata(f"{header}\n" + "\t".join(cells) + "\n", origin="shipped.tsv")

    def test_a_header_that_is_not_the_tables_raises_naming_the_columns_it_should_be(self) -> None:
        with pytest.raises(HomologyMetadataError, match="holding_species"):
            read_metadata("release\tspecies\n", origin="shipped.tsv")

    def test_an_attribution_line_names_the_publisher_the_release_and_the_file(self) -> None:
        line = _row("Homo sapiens", "Mus musculus").attribution()

        assert line.startswith("Ensembl Compara release 116 (PMID 26896847)")
        assert line.endswith("Compara.116.protein_default.homologies.tsv.gz")


# ---------------------------------------------------------------------------
# Preparing a set: construction fetches once, verifies, records and re-reads
# ---------------------------------------------------------------------------


class TestPreparation:
    def test_constructing_a_set_fetches_the_dump_the_shipped_row_pins(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")

        assert fake_fetch.last.url == homologs.source_url
        assert "mus_musculus/Compara.116.protein_default" in fake_fetch.last.url

    def test_the_publishers_own_checksum_is_demanded_of_the_bytes_as_they_are_fetched(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Load-bearing, not a formality: a resumed download of one of these gzips has been
        # seen to pass `gzip -t` with the wrong md5, so opening cleanly proves nothing.
        _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        assert fake_fetch.last.known_hash == "md5:59857f48bbbdf6812999d58d7a24ccc4"

    def test_a_set_is_filed_beside_the_assembly_tree_under_its_release_and_pair(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")

        assert homology_data_dir() == liulab_data / "homology"
        assert homologs.path.parent == (
            liulab_data / "homology" / "ensembl_compara" / "116" / "homo_sapiens-mus_musculus"
        )

    def test_a_second_construction_re_reads_what_is_there_and_fetches_nothing(
        self, fake_fetch: FakeFetch
    ) -> None:
        first = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        fetched = len(fake_fetch.calls)

        again = HomologySet("Homo sapiens", "Caenorhabditis elegans", progressbar=False)

        assert len(fake_fetch.calls) == fetched
        assert again.path == first.path
        assert len(again) == len(first)

    def test_the_reverse_orientation_reads_the_same_prepared_file_and_fetches_nothing(
        self, fake_fetch: FakeFetch
    ) -> None:
        # A pair is one download and one file; which species a caller asks *about* is a
        # property of the question, so the reverse orientation is a re-read.
        forward = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        fetched = len(fake_fetch.calls)

        backward = HomologySet("Caenorhabditis elegans", "Homo sapiens", progressbar=False)

        assert len(fake_fetch.calls) == fetched
        assert backward.path == forward.path
        assert len(backward) == len(forward)

    def test_cache_dir_overrides_the_homology_root_and_keeps_the_layout_beneath_it(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        elsewhere = tmp_path / "elsewhere"

        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus", cache_dir=elsewhere)

        assert homologs.path.parent == (
            elsewhere / "ensembl_compara" / "116" / "homo_sapiens-mus_musculus"
        )

    def test_the_stored_slice_is_a_plain_gzipped_tsv_of_the_publishers_own_rows(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # A collaborator must be able to read it in R or a shell without this package, and
        # every cell must still be the publisher's — nothing here re-derives one.
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")

        with gzip.open(homologs.path, "rt", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert tuple(lines[0].split("\t")) == COMPARA_COLUMNS
        assert [line.split("\t") for line in lines[1:]] == _published(
            data_dir, "Mus musculus", "Homo sapiens"
        )

    def test_a_completion_marker_carries_the_publishers_checksum_and_the_slices_own(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        record = read_record(homologs.path.parent)
        assert record is not None
        assert record.kind == "homology"
        assert record.source_url == homologs.source_url
        assert record.details["source_md5"] == "59857f48bbbdf6812999d58d7a24ccc4"
        assert record.sha256 is not None
        assert record.sha256 != record.details["source_md5"]
        assert record.details["links"] == len(homologs)
        assert record.files == {homologs.path.name: homologs.path.stat().st_size}

    def test_an_interrupted_preparation_reads_as_unfinished_rather_than_present(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        (homologs.path.parent / RECORD_NAME).unlink()

        with pytest.raises(UnfinishedRegistrationError, match="rm -rf"):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

    def test_a_slice_changed_after_it_was_prepared_raises_rather_than_answering_short(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Same rows, one label edited: only the recorded digest catches this, which is why
        # the slice's own sha256 is re-checked on every read and not only when it is made.
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        with gzip.open(homologs.path, "rt", encoding="utf-8") as handle:
            text = handle.read()
        with gzip.open(homologs.path, "wt", encoding="utf-8") as handle:
            handle.write(text.replace("ortholog_one2one", "ortholog_one2two"))
        marker = homologs.path.parent / RECORD_NAME
        payload = json.loads(marker.read_text())
        payload["files"][homologs.path.name] = homologs.path.stat().st_size
        marker.write_text(json.dumps(payload))

        with pytest.raises(ComparaFileError, match="does not hash"):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

    def test_a_fetched_file_that_is_not_comparas_raises_naming_the_columns_it_should_have(
        self, fake_fetch: FakeFetch
    ) -> None:
        fake_fetch.serve("tiny.gtf.gz")

        with pytest.raises(ComparaFileError, match="gene_stable_id"):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)


# ---------------------------------------------------------------------------
# The partition guard
# ---------------------------------------------------------------------------


class TestPartitionGuard:
    def test_a_pair_taken_from_the_file_that_does_not_hold_it_raises_naming_the_other(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The published partition, not a staged one: the human dump really does hold zero
        # human/mouse rows. Serving it for that pair is what a release that had
        # re-partitioned would look like from here.
        fake_fetch.serve(FIXTURES["Homo sapiens"])

        with pytest.raises(ComparaPartitionError) as raised:
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

        assert "homo_sapiens/Compara.116.protein_default.homologies.tsv.gz" in str(raised.value)
        assert "homology_metadata.tsv" in str(raised.value)

    def test_the_guard_names_the_publishers_own_statement_of_the_partition(
        self, fake_fetch: FakeFetch
    ) -> None:
        fake_fetch.serve(FIXTURES["Homo sapiens"])

        with pytest.raises(ComparaPartitionError, match="arbitrary subset"):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

    def test_a_guarded_pair_writes_no_slice_so_an_empty_set_is_never_re_read(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        fake_fetch.serve(FIXTURES["Homo sapiens"])
        directory = (
            liulab_data / "homology" / "ensembl_compara" / "116" / "homo_sapiens-mus_musculus"
        )

        with pytest.raises(ComparaPartitionError):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

        assert list(directory.glob("*.tsv.gz")) == []
        assert read_record(directory) is None

    def test_a_file_holding_neither_pair_raises_rather_than_answering_empty(
        self, fake_fetch: FakeFetch
    ) -> None:
        # The worm dump holds no cross-species row at all, which is the real file's shape.
        fake_fetch.serve(FIXTURES["Caenorhabditis elegans"])

        with pytest.raises(ComparaPartitionError):
            HomologySet("Caenorhabditis elegans", "Mus musculus", progressbar=False)


# ---------------------------------------------------------------------------
# Answering: all three pairings, in the shape ResolvedGeneIds established
# ---------------------------------------------------------------------------


class TestAnswers:
    @pytest.mark.parametrize(("species", "other", "links"), PAIRS)
    def test_all_three_pairings_among_human_mouse_and_worm_answer(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, links: int
    ) -> None:
        homologs = _built(fake_fetch, species, other)

        answer = homologs.homologs(_stems(data_dir, species, other))

        assert len(homologs) == links
        assert (answer.species, answer.other_species) == (species, other)
        assert answer.release == DEFAULT_RELEASE
        assert len(answer.links) == links
        assert answer.unresolved == ()

    @pytest.mark.parametrize(("species", "other", "links"), PAIRS)
    def test_each_pairing_answers_the_same_links_asked_the_other_way_round(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, links: int
    ) -> None:
        forward = _built(fake_fetch, species, other)
        backward = _built(fake_fetch, other, species)

        ahead = {
            (link.gene_id_stem, link.homolog_gene_id_stem, link.homology_type)
            for link in forward.homologs(_stems(data_dir, species, other)).links
        }
        behind = {
            (link.homolog_gene_id_stem, link.gene_id_stem, link.homology_type)
            for link in backward.homologs(_stems(data_dir, other, species)).links
        }
        assert ahead == behind
        assert len(ahead) == links

    def test_a_stem_answers_with_every_homolog_and_never_a_chosen_one(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs([ONE2MANY_HUMAN])

        partners = tuple(link.homolog_gene_id_stem for link in answer.resolved[ONE2MANY_HUMAN])
        assert partners == ONE2MANY_WORMS

    def test_what_named_nothing_rides_back_and_no_resolved_value_is_ever_empty(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs([ONE2ONE_HUMAN, ABSENT])

        assert answer.unresolved == (ABSENT,)
        assert ABSENT not in answer.resolved
        assert all(found for found in answer.resolved.values())

    def test_the_answer_keeps_ask_order_and_asks_a_repeat_once(self, fake_fetch: FakeFetch) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs([MANY2MANY_HUMAN[1], ABSENT, ONE2ONE_HUMAN, MANY2MANY_HUMAN[1]])

        assert list(answer.resolved) == [MANY2MANY_HUMAN[1], ONE2ONE_HUMAN]
        assert answer.unresolved == (ABSENT,)

    def test_flattening_keeps_a_partner_two_asked_stems_both_name(
        self, fake_fetch: FakeFetch
    ) -> None:
        # What flattening loses is which stem reached which partner, so the two worm genes
        # both human genes are many2many to appear twice rather than once.
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs(list(MANY2MANY_HUMAN))

        assert answer.homolog_gene_id_stems == list(MANY2MANY_WORMS) * 2

    def test_asking_about_nothing_answers_emptily_rather_than_about_every_gene(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")

        answer = homologs.homologs([])

        assert answer.resolved == {}
        assert answer.unresolved == ()

    def test_as_json_writes_the_publishers_type_the_nulls_and_the_flattened_partners(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        payload = homologs.homologs([ONE2ONE_HUMAN, ABSENT]).as_json()

        assert json.loads(json.dumps(payload)) == payload
        assert payload["release"] == DEFAULT_RELEASE
        assert payload["unresolved"] == [ABSENT]
        assert payload["null_quality_scores"] == list(QUALITY_SCORE_COLUMNS)
        assert payload["homolog_gene_id_stems"] == [ONE2ONE_WORM]
        (link,) = payload["resolved"][ONE2ONE_HUMAN]
        assert link["homology_type"] == "ortholog_one2one"
        assert link["is_ortholog"] is True
        assert link["goc_score"] is None
        assert link["is_high_confidence"] is True

    @settings(suppress_health_check=[HealthCheck.function_scoped_fixture], max_examples=40)
    @given(
        asked=st.lists(
            st.sampled_from(
                [ONE2MANY_HUMAN, ONE2ONE_HUMAN, *MANY2MANY_HUMAN, ABSENT, "ENSG00000000001"]
            ),
            max_size=10,
        )
    )
    def test_resolved_and_unresolved_partition_the_asked_stems_exactly_once(
        self, fake_fetch: FakeFetch, asked: list[str]
    ) -> None:
        # Held over generated inputs because it is the invariant every answer shape here
        # promises: nothing is dropped, nothing is counted twice, a repeat is asked once,
        # and the order the caller passed is the order they read back. The set is prepared
        # once for the whole run and read-only, which is why the fixture is shared across
        # examples rather than reset between them.
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs(asked)

        once = list(dict.fromkeys(asked))
        assert sorted([*answer.resolved, *answer.unresolved]) == sorted(once)
        assert list(answer.resolved) == [stem for stem in once if stem in answer.resolved]
        assert list(answer.unresolved) == [stem for stem in once if stem not in answer.resolved]


# ---------------------------------------------------------------------------
# The Homology type, and what a filter may not do to it
# ---------------------------------------------------------------------------


class TestHomologyType:
    @pytest.mark.parametrize(("species", "other", "_links"), PAIRS)
    def test_every_link_carries_the_publishers_own_type_verbatim(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, _links: int
    ) -> None:
        homologs = _built(fake_fetch, species, other)
        at, then = (0, 5) if _row(species, other).holding_species == species else (5, 0)

        answer = homologs.homologs(_stems(data_dir, species, other), paralogs=True)

        assert {
            (link.gene_id_stem, link.homolog_gene_id_stem, link.homology_type)
            for link in answer.links
        } == {(row[at], row[then], row[4]) for row in _published(data_dir, species, other)}

    def test_a_type_is_read_and_never_counted_so_two_partners_may_read_one_to_many(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Counting would call both of these `one2many`: each has exactly two partners here.
        # The publisher's tree says otherwise about one of them, and the label is its say.
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")

        answer = homologs.homologs([ONE2MANY_MOUSE, MANY2MANY_MOUSE])

        assert len(answer.resolved[ONE2MANY_MOUSE]) == len(answer.resolved[MANY2MANY_MOUSE]) == 2
        assert {link.homology_type for link in answer.resolved[ONE2MANY_MOUSE]} == {
            "ortholog_one2many"
        }
        assert {link.homology_type for link in answer.resolved[MANY2MANY_MOUSE]} == {
            "ortholog_many2many"
        }


class TestParalogy:
    @pytest.mark.parametrize(("species", "other", "_links"), PAIRS)
    def test_orthologs_are_the_default_and_every_type_returned_is_a_speciation_label(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, _links: int
    ) -> None:
        homologs = _built(fake_fetch, species, other)

        answer = homologs.homologs(_stems(data_dir, species, other))

        assert answer.links
        assert all(link.is_ortholog for link in answer.links)

    def test_a_same_species_paralogy_row_is_not_in_a_pairs_set_because_a_link_needs_two(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # The fixture's four human `other_paralog` rows are real published rows sitting in
        # the very file this pair was cut from. They relate two *human* genes, so they are
        # not this pair's, and `paralogs=True` does not reach them either.
        paralogy = sorted(
            {row[0] for row in _fixture_rows(data_dir, "Homo sapiens") if row[2] == row[7]}
        )
        assert paralogy

        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        answer = homologs.homologs(paralogy, paralogs=True)
        assert answer.resolved == {}
        assert sorted(answer.unresolved) == paralogy

    def test_release_116_publishes_no_cross_species_paralogy_so_the_switch_changes_nothing(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # The measurement, asserted rather than assumed: counted over the whole human dump
        # there is not one `between_species_paralog` row, and every duplication label is
        # same-species. The switch is kept because it is where such a row would land the
        # release Compara publishes one, and because *not an ortholog* must have somewhere
        # to be distinguishable from *absent* (ADR-0013). On this release it agrees.
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        asked = _stems(data_dir, "Homo sapiens", "Caenorhabditis elegans")

        orthologs = homologs.homologs(asked)
        everything = homologs.homologs(asked, paralogs=True)

        assert orthologs.as_json() == everything.as_json()
        assert orthologs.dropped_partners == ()

    def test_the_set_keeps_every_row_the_publisher_wrote_for_the_pair_unfiltered(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Kept and marked, never excluded: the stored slice is what the publisher wrote for
        # the pair, of every type, and the filtering happens on the way out.
        homologs = _built(fake_fetch, "Mus musculus", "Caenorhabditis elegans")

        with gzip.open(homologs.path, "rt", encoding="utf-8") as handle:
            stored = [line.split("\t") for line in handle.read().splitlines()[1:]]
        assert stored == _published(data_dir, "Mus musculus", "Caenorhabditis elegans")


# ---------------------------------------------------------------------------
# Quality scores, and being told when they are null
# ---------------------------------------------------------------------------


class TestQualityScores:
    @pytest.mark.parametrize("other", ["Homo sapiens", "Mus musculus"])
    def test_a_worm_pairing_says_both_quality_scores_are_null_rather_than_emptying_a_filter(
        self, fake_fetch: FakeFetch, data_dir: Path, other: str
    ) -> None:
        homologs = _built(fake_fetch, other, "Caenorhabditis elegans")

        answer = homologs.homologs(_stems(data_dir, other, "Caenorhabditis elegans"))

        assert homologs.null_quality_scores == QUALITY_SCORE_COLUMNS
        assert answer.null_quality_scores == QUALITY_SCORE_COLUMNS
        assert all(link.goc_score is None for link in answer.links)
        assert all(link.wga_coverage is None for link in answer.links)
        # And the filter a caller would have written really does empty, which is the whole
        # reason the answer says so before they ever run it.
        assert [link for link in answer.links if (link.goc_score or 0) > 50] == []

    def test_the_mouse_human_pairing_carries_comparas_quality_scores_through(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")

        answer = homologs.homologs([ONE2MANY_MOUSE])

        assert homologs.null_quality_scores == ()
        assert answer.null_quality_scores == ()
        assert {
            link.homolog_gene_id_stem: (link.goc_score, link.wga_coverage, link.is_high_confidence)
            for link in answer.resolved[ONE2MANY_MOUSE]
        } == {
            "ENSG00000101266": (100, 100.0, True),
            "ENSG00000254598": (0, 0.0, False),
        }

    def test_a_null_score_is_a_fact_about_the_pair_and_not_about_the_file_it_came_from(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Both sets below are cut from the *same* mouse dump, and one of them is scored.
        scored = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        unscored = _built(fake_fetch, "Mus musculus", "Caenorhabditis elegans")

        assert scored.null_quality_scores == ()
        assert unscored.null_quality_scores == QUALITY_SCORE_COLUMNS


# ---------------------------------------------------------------------------
# Crossing into a registered annotation
# ---------------------------------------------------------------------------


class TestCrossingIntoAnAnnotation:
    def test_resolving_into_an_annotation_leaves_the_type_unchanged_and_counts_the_dropped(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # The one that matters. Two human partners, an annotation that spells one of them:
        # the view now looks one-to-one, the label must still read one2many, and the
        # partner that fell away is counted rather than quietly gone (ADR-0020).
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        answer = homologs.homologs([ONE2MANY_MOUSE])
        registry = _registry_declaring(tmp_path, f"{ONE2MANY_HUMANS[0]}.7")

        crossed = resolve_homologs(answer, registry, "mine")

        (link,) = crossed.resolved[ONE2MANY_MOUSE]
        assert link.homology_type == "ortholog_one2many"
        assert link.homolog_gene_id_stem == ONE2MANY_HUMANS[0]
        assert crossed.dropped_partners == (ONE2MANY_HUMANS[1],)
        assert crossed.gene_ids == {ONE2MANY_HUMANS[0]: (f"{ONE2MANY_HUMANS[0]}.7",)}
        assert crossed.homolog_gene_ids == [f"{ONE2MANY_HUMANS[0]}.7"]

    def test_a_stem_whose_every_partner_the_annotation_is_missing_becomes_unresolved(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE, ABSENT])
        registry = _registry_declaring(tmp_path, "ENSMUSG99999999.1")

        crossed = resolve_homologs(answer, registry, "mine")

        assert crossed.resolved == {}
        assert crossed.unresolved == (MANY2MANY_HUMAN_OF_MOUSE, ABSENT)
        assert crossed.dropped_partners == MANY2MANY_MICE

    def test_a_partner_the_annotation_spells_twice_answers_with_both_gene_ids(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # The pseudoautosomal shape `resolve_gene_ids` already answers with both of, and
        # nothing here picks one either.
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE])
        registry = _registry_declaring(
            tmp_path, f"{MANY2MANY_MICE[0]}.3", f"{MANY2MANY_MICE[0]}.3_PAR_Y"
        )

        crossed = resolve_homologs(answer, registry, "mine")

        assert crossed.gene_ids[MANY2MANY_MICE[0]] == (
            f"{MANY2MANY_MICE[0]}.3",
            f"{MANY2MANY_MICE[0]}.3_PAR_Y",
        )
        assert crossed.homolog_gene_ids == [
            f"{MANY2MANY_MICE[0]}.3",
            f"{MANY2MANY_MICE[0]}.3_PAR_Y",
        ]

    def test_the_crossing_names_the_assembly_and_annotation_and_keeps_the_pair_and_release(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE])
        registry = _registry_declaring(tmp_path, f"{MANY2MANY_MICE[0]}.3")

        crossed = resolve_homologs(answer, registry, "mine")

        assert (crossed.assembly, crossed.annotation) == ("tiny", "mine")
        assert (crossed.species, crossed.other_species) == ("Homo sapiens", "Mus musculus")
        assert crossed.release == DEFAULT_RELEASE
        payload = crossed.as_json()
        assert json.loads(json.dumps(payload)) == payload
        assert payload["dropped_partners"] == [MANY2MANY_MICE[1]]
        assert payload["resolved"][MANY2MANY_HUMAN_OF_MOUSE][0]["homology_type"] == (
            "ortholog_many2many"
        )

    def test_a_crossing_adds_no_link_back_that_the_answer_had_already_filtered_out(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # Whatever the answer was filtered to is what is crossed: the hop never re-reads
        # the set, so it can neither add a partner nor recompute a label.
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE])
        registry = _registry_declaring(
            tmp_path, f"{MANY2MANY_MICE[0]}.3", f"{MANY2MANY_MICE[1]}.4", "ENSMUSG00000051306.1"
        )

        crossed = resolve_homologs(answer, registry, "mine")

        assert list(crossed.gene_ids) == list(MANY2MANY_MICE)
        assert crossed.dropped_partners == ()


# ---------------------------------------------------------------------------
# What is refused, and what the refusal names
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_a_species_with_no_set_raises_naming_the_species_that_have_one(self) -> None:
        with pytest.raises(UnknownHomologySpeciesError) as raised:
            HomologySet("Danio rerio", "Homo sapiens", progressbar=False)

        assert "Homo sapiens" in str(raised.value)
        assert "Mus musculus" in str(raised.value)

    def test_a_release_that_is_not_pinned_raises_naming_the_ones_that_are(self) -> None:
        with pytest.raises(ValueError, match="116"):
            HomologySet("Homo sapiens", "Mus musculus", "115", progressbar=False)

    def test_a_pair_of_one_species_raises_because_a_set_relates_two(self) -> None:
        with pytest.raises(ValueError, match="two different species"):
            HomologySet("Homo sapiens", "homo_sapiens", progressbar=False)

    def test_a_pair_the_table_does_not_pin_raises_naming_the_pairs_it_does(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Both species are prepared and the release is pinned; the *pair* is not, which is
        # a third thing and says so rather than answering empty.
        from genome.homology import metadata as metadata_mod

        kept = tuple(row for row in homology_table() if "Mus musculus" not in row.pair)
        monkeypatch.setattr(metadata_mod, "homology_metadata", lambda *_a, **_k: None)
        monkeypatch.setattr(metadata_mod, "homology_table", lambda: kept)

        with pytest.raises(NoHomologyPairError, match="Caenorhabditis elegans/Homo sapiens"):
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

    def test_a_versioned_gene_id_is_refused_naming_the_stem_to_pass(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Compara writes its ids bare, so the versioned spelling would match nothing and
        # come back unresolved looking exactly like a gene it never placed in a tree.
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        with pytest.raises(VersionedGeneIdError, match=ONE2ONE_HUMAN):
            homologs.homologs([f"{ONE2ONE_HUMAN}.18"])

    def test_either_spelling_of_a_species_is_accepted_and_the_answer_names_one(
        self, fake_fetch: FakeFetch
    ) -> None:
        fake_fetch.serve(FIXTURES["Homo sapiens"])

        homologs = HomologySet("homo_sapiens", "caenorhabditis_elegans", progressbar=False)

        assert (homologs.species, homologs.other_species) == (
            "Homo sapiens",
            "Caenorhabditis elegans",
        )
        assert repr(homologs).startswith("HomologySet(species='Homo sapiens'")


# ---------------------------------------------------------------------------
# Orthology is served and never consumed (ADR-0019)
# ---------------------------------------------------------------------------


def test_nothing_this_package_publishes_is_derived_through_homology() -> None:
    """No module outside ``genome.homology`` may import it, so no shipped table can use it.

    The structural half of ADR-0019: a **TF gene table**, a **Cofactor table** or any list
    this package ships cannot be derived through homology if nothing that builds one can
    reach it — not even ``genome/__init__.py``, which is why the subpackage is addressed
    as ``genome.homology`` and never re-exported at the top level. Read off the source
    rather than observed, since importing a module loads its package and ``sys.modules``
    would then say more than any module asked for.
    """
    package = Path(__file__).resolve().parents[1] / "src" / "genome"
    offenders: list[str] = []
    for path in sorted(package.rglob("*.py")):
        if path.parent.name == "homology":
            continue
        for node in ast.parse(path.read_text(encoding="utf-8")).body:
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            else:
                continue
            if any(name.startswith("genome.homology") for name in names):
                offenders.append(str(path.relative_to(package)))
    assert offenders == []


def test_the_shipped_homology_directory_holds_provenance_and_no_homology_data() -> None:
    """Data is downloaded, never redistributed, which is what keeps its checksum honest.

    It is also what makes *orthology is served and never consumed* a property of the wheel
    rather than a promise: there is no homology table in here for anything to be built
    from.
    """
    directory = Path(__file__).resolve().parents[1] / "src" / "genome" / "data" / "homology"

    assert sorted(entry.name for entry in directory.iterdir()) == [
        "ATTRIBUTION.md",
        "homology_metadata.tsv",
    ]
