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
from dataclasses import replace
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from genome.annotation import AnnotationRegistry
from genome.assembly.registration import AssemblyDir
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
    HomologySetNotDownloadedError,
    NoHomologyPairError,
    UnknownHomologySpeciesError,
    VersionedGeneIdError,
    compara_url,
    homology_data_dir,
    homology_metadata,
    homology_prepare_command,
    homology_releases,
    homology_species,
    homology_table,
    read_metadata,
    resolve_homologs,
)
from genome.store import fetch as fetch_mod
from genome.store.completion import RECORD_NAME, UnfinishedRegistrationError, read_record

from ..conftest import FakeFetch

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

#: A mouse stem no fixture mentions, standing in for a partner a filter removed before an
#: answer ever reached an annotation. Sorts before every mouse stem the fixtures carry.
ABSENT_PARTNER = "ENSMUSG00000000001"

#: The package's own source tree, which the structural bans at the foot of this file read.
PACKAGE = Path(__file__).resolve().parents[2] / "src" / "genome"


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


def _rewrite_cell(path: Path, *, column: str, value: str) -> None:
    """Put ``value`` in one column of a stored slice's first row, record still agreeing.

    The **Completion marker**'s file sizes are re-stated afterwards, so what the re-read
    trips over is the edited cell rather than a directory that no longer looks finished.
    """
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    cells = lines[1].split("\t")
    cells[COMPARA_COLUMNS.index(column)] = value
    lines[1] = "\t".join(cells)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    marker = path.parent / RECORD_NAME
    payload = json.loads(marker.read_text())
    payload["files"][path.name] = path.stat().st_size
    marker.write_text(json.dumps(payload))


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
    def test_the_readmes_claims_about_the_fixtures_hold_on_the_committed_bytes(
        self, data_dir: Path
    ) -> None:
        for species in sorted(FIXTURES):
            with gzip.open(data_dir / FIXTURES[species], "rt", encoding="utf-8") as handle:
                assert tuple(handle.readline().rstrip("\n").split("\t")) == COMPARA_COLUMNS
            total, pairs, paralogy = FIXTURE_SHAPE[species]
            rows = _fixture_rows(data_dir, species)
            assert len(rows) == total, species
            for partner, count in pairs.items():
                here = [row for row in rows if row[7] == partner and row[2] != row[7]]
                assert len(here) == count, (species, partner)
            assert len([row for row in rows if row[2] == row[7]]) == paralogy, species
            # Release 116 publishes no `between_species_paralog` anywhere in the human
            # dump — 4.0M rows over some 200 partner species — so a duplication label
            # always means one species, for every fixture.
            for row in rows:
                if not row[4].startswith("ortholog_"):
                    assert row[2] == row[7], row

        # Not a staged absence: the real release-116 human file holds 0 human/mouse rows
        # and 23,982 human/worm rows, and the mouse file holds the human pair — so the
        # partition guard exercised in TestPartitionGuard is against the publisher's own
        # partition.
        human_rows = _fixture_rows(data_dir, "Homo sapiens")
        assert [row for row in human_rows if row[7] == "mus_musculus"] == []
        assert len([row for row in human_rows if row[7] == "caenorhabditis_elegans"]) == 11

        for row in _fixture_rows(data_dir, "Caenorhabditis elegans"):
            assert {row[2], row[7]} == {"caenorhabditis_elegans"}

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

        assert {row[13] for row in human_rows if row[2] == row[7]} == {"NULL"}
        assert {row[13] for row in human_rows if row[7] == "caenorhabditis_elegans"} == {"0", "1"}

        for species, other, _links in PAIRS:
            assert {row[4] for row in _published(data_dir, species, other)} == {
                "ortholog_one2one",
                "ortholog_one2many",
                "ortholog_many2many",
            }


# ---------------------------------------------------------------------------
# The shipped provenance table
# ---------------------------------------------------------------------------


class TestShippedTable:
    def test_the_shipped_table_pins_releases_species_pairs_publisher_metadata_and_attribution(
        self,
    ) -> None:
        assert homology_releases() == (DEFAULT_RELEASE,)
        assert homology_species() == ("Caenorhabditis elegans", "Homo sapiens", "Mus musculus")

        for row in homology_table():
            assert row.publisher == "Ensembl Compara"
            assert row.source_url.startswith("https://ftp.ensembl.org/pub/")
            assert len(row.md5) == 32
            assert set(row.md5) <= set("0123456789abcdef")
            assert row.pubmed_id > 0
            # The URL is written down rather than formatted — release 113 ships these
            # dumps uncompressed, so no template covers every release.
            assert row.source_url == compara_url(row.holding_species, row.release)

        held = {row.pair: row.holding_species for row in homology_table()}
        assert held == {
            ("Caenorhabditis elegans", "Homo sapiens"): "Homo sapiens",
            ("Caenorhabditis elegans", "Mus musculus"): "Mus musculus",
            ("Homo sapiens", "Mus musculus"): "Mus musculus",
        }
        assert _row("Homo sapiens", "Mus musculus") is _row("Mus musculus", "Homo sapiens")
        assert homology_metadata("Homo sapiens", "Danio rerio", DEFAULT_RELEASE) is None

        line = _row("Homo sapiens", "Mus musculus").attribution()
        assert line.startswith("Ensembl Compara release 116 (PMID 26896847)")
        assert line.endswith("Compara.116.protein_default.homologies.tsv.gz")

    def test_a_malformed_metadata_row_or_header_raises_naming_what_is_wrong(self) -> None:
        header = "\t".join(METADATA_COLUMNS)
        cells = ["116", "Homo sapiens", "Mus musculus", "Mus musculus", "Ensembl", "1", "u", ""]
        with pytest.raises(HomologyMetadataError, match="'md5'"):
            read_metadata(f"{header}\n" + "\t".join(cells) + "\n", origin="shipped.tsv")

        with pytest.raises(HomologyMetadataError, match="holding_species"):
            read_metadata("release\tspecies\n", origin="shipped.tsv")


# ---------------------------------------------------------------------------
# Preparing a set: construction fetches once, verifies, records and re-reads
# ---------------------------------------------------------------------------


class TestPreparation:
    def test_construction_fetches_verifies_and_records_and_every_later_ask_reuses_the_one_file(
        self, fake_fetch: FakeFetch, liulab_data: Path, data_dir: Path, tmp_path: Path
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        assert fake_fetch.last.url == homologs.source_url
        assert "mus_musculus/Compara.116.protein_default" in fake_fetch.last.url
        assert homology_data_dir() == liulab_data / "homology"
        assert homologs.path.parent == (
            liulab_data / "homology" / "ensembl_compara" / "116" / "homo_sapiens-mus_musculus"
        )

        # Load-bearing, not a formality: a resumed download of one of these gzips has been
        # seen to pass `gzip -t` with the wrong md5, so opening cleanly proves nothing.
        worm_pair = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        assert fake_fetch.last.known_hash == "md5:59857f48bbbdf6812999d58d7a24ccc4"

        fetched = len(fake_fetch.calls)
        again = HomologySet("Homo sapiens", "Caenorhabditis elegans", progressbar=False)
        assert len(fake_fetch.calls) == fetched
        assert again.path == worm_pair.path
        assert len(again) == len(worm_pair)

        # A pair is one download and one file; which species a caller asks *about* is a
        # property of the question, so the reverse orientation is a re-read.
        backward = HomologySet("Caenorhabditis elegans", "Homo sapiens", progressbar=False)
        assert len(fake_fetch.calls) == fetched
        assert backward.path == worm_pair.path
        assert len(backward) == len(worm_pair)

        # `cache_dir` means the same thing across the package: an `XrefSet` and a
        # `JasparDatabase` both prepare *in* the directory they are handed, and a set that
        # re-applied a layout beneath it would be the one exception nobody expects.
        elsewhere = tmp_path / "elsewhere"
        elsewhere_homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus", cache_dir=elsewhere)
        assert elsewhere_homologs.path.parent == elsewhere

        # A collaborator must be able to read the slice in R or a shell without this
        # package, and every cell must still be the publisher's — nothing here re-derives
        # one.
        mouse_first = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        with gzip.open(mouse_first.path, "rt", encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        assert tuple(lines[0].split("\t")) == COMPARA_COLUMNS
        assert [line.split("\t") for line in lines[1:]] == _published(
            data_dir, "Mus musculus", "Homo sapiens"
        )

        record = read_record(worm_pair.path.parent)
        assert record is not None
        assert record.kind == "homology"
        assert record.source_url == worm_pair.source_url
        assert record.details["source_md5"] == "59857f48bbbdf6812999d58d7a24ccc4"
        assert record.sha256 is not None
        assert record.sha256 != record.details["source_md5"]
        assert record.details["links"] == len(worm_pair)
        assert record.files == {worm_pair.path.name: worm_pair.path.stat().st_size}

        # Either spelling of a species is accepted, and the answer names one.
        spelled = HomologySet("homo_sapiens", "caenorhabditis_elegans", progressbar=False)
        assert (spelled.species, spelled.other_species) == (
            "Homo sapiens",
            "Caenorhabditis elegans",
        )
        assert repr(spelled).startswith("HomologySet(species='Homo sapiens'")

    def test_an_interrupted_preparation_reads_as_unfinished_and_the_repair_names_the_rebuild(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        (homologs.path.parent / RECORD_NAME).unlink()

        with pytest.raises(UnfinishedRegistrationError, match="rm -rf") as raised:
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

        # Deleting the directory is half a repair: it leaves the caller with nothing and no
        # way back, so the message must also name the call that rebuilds it.
        assert homology_prepare_command("Homo sapiens", "Mus musculus", DEFAULT_RELEASE) in str(
            raised.value
        )

    def test_a_quality_cell_that_cannot_be_read_raises_rather_than_reading_as_a_null(
        self, fake_fetch: FakeFetch
    ) -> None:
        # Two different facts, and the same `None` would have said both: *Compara recorded
        # no score* and *this package could not read the one it recorded*. The second is a
        # hole in a column a caller filters on, so it is named — file, column and cell.
        homologs = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        _rewrite_cell(homologs.path, column="goc_score", value="not-a-score")

        with pytest.raises(ComparaFileError) as raised:
            HomologySet("Mus musculus", "Homo sapiens", progressbar=False)

        message = str(raised.value)
        assert str(homologs.path) in message
        assert "goc_score" in message
        assert "not-a-score" in message

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

    def test_a_set_that_cannot_be_fetched_names_the_call_to_make_on_a_login_node(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # What a compute node looks like from here: fetching is the one step in this
        # package that needs the network, and the lab's CPU cluster nodes have none. A
        # transport error on its own would leave the reader to work out that the repair is
        # to run this somewhere else first.
        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)

        with pytest.raises(HomologySetNotDownloadedError) as raised:
            HomologySet("Mus musculus", "Homo sapiens", progressbar=False)

        assert "login node" in str(raised.value)
        assert homology_prepare_command("Homo sapiens", "Mus musculus", DEFAULT_RELEASE) in str(
            raised.value
        )


# ---------------------------------------------------------------------------
# The partition guard
# ---------------------------------------------------------------------------


class TestPartitionGuard:
    def test_a_pair_taken_from_the_wrong_file_is_refused_and_writes_no_slice(
        self, fake_fetch: FakeFetch, liulab_data: Path
    ) -> None:
        # The published partition, not a staged one: the human dump really does hold zero
        # human/mouse rows. Serving it for that pair is what a release that had
        # re-partitioned would look like from here.
        fake_fetch.serve(FIXTURES["Homo sapiens"])
        directory = (
            liulab_data / "homology" / "ensembl_compara" / "116" / "homo_sapiens-mus_musculus"
        )

        with pytest.raises(ComparaPartitionError) as raised:
            HomologySet("Homo sapiens", "Mus musculus", progressbar=False)

        message = str(raised.value)
        assert "homo_sapiens/Compara.116.protein_default.homologies.tsv.gz" in message
        assert "homology_metadata.tsv" in message
        assert "arbitrary subset" in message
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
    def test_all_three_pairings_answer_and_agree_when_asked_the_other_way_round(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, links: int
    ) -> None:
        forward = _built(fake_fetch, species, other)
        answer = forward.homologs(_stems(data_dir, species, other))
        assert len(forward) == links
        assert (answer.species, answer.other_species) == (species, other)
        assert answer.release == DEFAULT_RELEASE
        assert len(answer.links) == links
        assert answer.unresolved == ()

        backward = _built(fake_fetch, other, species)
        ahead = {
            (link.gene_id_stem, link.homolog_gene_id_stem, link.homology_type)
            for link in answer.links
        }
        behind = {
            (link.homolog_gene_id_stem, link.gene_id_stem, link.homology_type)
            for link in backward.homologs(_stems(data_dir, other, species)).links
        }
        assert ahead == behind

    def test_the_answer_holds_every_homolog_in_order_deduped_and_flattened_without_loss(
        self, fake_fetch: FakeFetch
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")

        every = homologs.homologs([ONE2MANY_HUMAN])
        partners = tuple(link.homolog_gene_id_stem for link in every.resolved[ONE2MANY_HUMAN])
        assert partners == ONE2MANY_WORMS

        missing = homologs.homologs([ONE2ONE_HUMAN, ABSENT])
        assert missing.unresolved == (ABSENT,)
        assert ABSENT not in missing.resolved
        assert all(found for found in missing.resolved.values())

        ordered = homologs.homologs([MANY2MANY_HUMAN[1], ABSENT, ONE2ONE_HUMAN, MANY2MANY_HUMAN[1]])
        assert list(ordered.resolved) == [MANY2MANY_HUMAN[1], ONE2ONE_HUMAN]
        assert ordered.unresolved == (ABSENT,)

        # What flattening loses is which stem reached which partner, so the two worm genes
        # both human genes are many2many to appear twice rather than once.
        flattened = homologs.homologs(list(MANY2MANY_HUMAN))
        assert flattened.homolog_gene_id_stems == list(MANY2MANY_WORMS) * 2

    def test_an_empty_ask_answers_emptily_and_json_carries_the_type_nulls_and_partners(
        self, fake_fetch: FakeFetch
    ) -> None:
        empty_set = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        empty = empty_set.homologs([])
        assert empty.resolved == {}
        assert empty.unresolved == ()

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
    def test_the_publishers_type_is_read_verbatim_and_never_counted(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Every pairing, not just one: the type-vs-published-row check must hold for all
        # three so a link built for the mouse/worm pair is never taken on faith from the
        # human/mouse one.
        for species, other, _links in PAIRS:
            homologs = _built(fake_fetch, species, other)
            at, then = (0, 5) if _row(species, other).holding_species == species else (5, 0)
            answer = homologs.homologs(_stems(data_dir, species, other), paralogs=True)
            assert {
                (link.gene_id_stem, link.homolog_gene_id_stem, link.homology_type)
                for link in answer.links
            } == {(row[at], row[then], row[4]) for row in _published(data_dir, species, other)}, (
                species,
                other,
            )

        # Counting would call both of these `one2many`: each has exactly two partners here.
        # The publisher's tree says otherwise about one of them, and the label is its say.
        mouse_human = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        typed = mouse_human.homologs([ONE2MANY_MOUSE, MANY2MANY_MOUSE])
        assert len(typed.resolved[ONE2MANY_MOUSE]) == len(typed.resolved[MANY2MANY_MOUSE]) == 2
        assert {link.homology_type for link in typed.resolved[ONE2MANY_MOUSE]} == {
            "ortholog_one2many"
        }
        assert {link.homology_type for link in typed.resolved[MANY2MANY_MOUSE]} == {
            "ortholog_many2many"
        }


class TestParalogy:
    def test_orthologs_default_paralogs_true_changes_nothing_here_and_the_slice_is_unfiltered(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Every pairing, not just one: orthologs are the default and every type any of the
        # three pairs returns is a speciation label, unasked.
        for species, other, _links in PAIRS:
            paired = _built(fake_fetch, species, other)
            answer = paired.homologs(_stems(data_dir, species, other))
            assert answer.links, (species, other)
            assert all(link.is_ortholog for link in answer.links), (species, other)

        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        asked = _stems(data_dir, "Homo sapiens", "Caenorhabditis elegans")
        orthologs = homologs.homologs(asked)

        # The fixture's four human `other_paralog` rows are real published rows sitting in
        # the very file this pair was cut from. They relate two *human* genes, so they are
        # not this pair's, and `paralogs=True` does not reach them either.
        paralogy = sorted(
            {row[0] for row in _fixture_rows(data_dir, "Homo sapiens") if row[2] == row[7]}
        )
        assert paralogy
        assert homologs.homologs(paralogy, paralogs=True).resolved == {}

        # The measurement, asserted rather than assumed: counted over the whole human dump
        # there is not one `between_species_paralog` row, so `paralogs=True` changes
        # nothing on this release, which is where such a row would land were there one
        # (ADR-0013).
        everything = homologs.homologs(asked, paralogs=True)
        assert orthologs.as_json() == everything.as_json()
        assert orthologs.dropped_partners == ()

        # Kept and marked, never excluded: the stored slice is what the publisher wrote for
        # a pair, of every type, and the filtering happens on the way out.
        other_pair = _built(fake_fetch, "Mus musculus", "Caenorhabditis elegans")
        with gzip.open(other_pair.path, "rt", encoding="utf-8") as handle:
            stored = [line.split("\t") for line in handle.read().splitlines()[1:]]
        assert stored == _published(data_dir, "Mus musculus", "Caenorhabditis elegans")


# ---------------------------------------------------------------------------
# Quality scores, and being told when they are null
# ---------------------------------------------------------------------------


class TestQualityScores:
    def test_a_null_score_is_a_fact_about_the_pair_and_real_scores_carry_through_unchanged(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        # Both sets below are cut from the *same* mouse dump, and one of them is scored.
        scored = _built(fake_fetch, "Mus musculus", "Homo sapiens")
        unscored = _built(fake_fetch, "Mus musculus", "Caenorhabditis elegans")
        assert scored.null_quality_scores == ()
        assert unscored.null_quality_scores == QUALITY_SCORE_COLUMNS

        answer = unscored.homologs(_stems(data_dir, "Mus musculus", "Caenorhabditis elegans"))
        assert answer.null_quality_scores == QUALITY_SCORE_COLUMNS
        assert all(link.goc_score is None for link in answer.links)
        assert all(link.wga_coverage is None for link in answer.links)
        # The filter a caller would have written really does empty, which is the whole
        # reason the answer says so before they ever run it.
        assert [link for link in answer.links if (link.goc_score or 0) > 50] == []

        scored_answer = scored.homologs([ONE2MANY_MOUSE])
        assert scored_answer.null_quality_scores == ()
        assert {
            link.homolog_gene_id_stem: (link.goc_score, link.wga_coverage, link.is_high_confidence)
            for link in scored_answer.resolved[ONE2MANY_MOUSE]
        } == {
            "ENSG00000101266": (100, 100.0, True),
            "ENSG00000254598": (0, 0.0, False),
        }


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

    def test_a_partner_spelled_twice_answers_with_both_ids_and_crossing_names_assembly_and_release(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE])

        # The pseudoautosomal shape `resolve_gene_ids` already answers with both of, and
        # nothing here picks one either.
        par_registry = _registry_declaring(
            tmp_path / "par", f"{MANY2MANY_MICE[0]}.3", f"{MANY2MANY_MICE[0]}.3_PAR_Y"
        )
        par_crossed = resolve_homologs(answer, par_registry, "mine")
        assert par_crossed.gene_ids[MANY2MANY_MICE[0]] == (
            f"{MANY2MANY_MICE[0]}.3",
            f"{MANY2MANY_MICE[0]}.3_PAR_Y",
        )
        assert par_crossed.homolog_gene_ids == [
            f"{MANY2MANY_MICE[0]}.3",
            f"{MANY2MANY_MICE[0]}.3_PAR_Y",
        ]

        registry = _registry_declaring(tmp_path / "single", f"{MANY2MANY_MICE[0]}.3")
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

    def test_a_crossing_never_adds_back_a_filtered_link_and_counts_what_was_dropped_before_it(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        homologs = _built(fake_fetch, "Homo sapiens", "Mus musculus")
        answer = homologs.homologs([MANY2MANY_HUMAN_OF_MOUSE])

        # Whatever the answer was filtered to is what is crossed: the hop never re-reads
        # the set, so it can neither add a partner nor recompute a label.
        wide_registry = _registry_declaring(
            tmp_path / "wide",
            f"{MANY2MANY_MICE[0]}.3",
            f"{MANY2MANY_MICE[1]}.4",
            "ENSMUSG00000051306.1",
        )
        wide_crossed = resolve_homologs(answer, wide_registry, "mine")
        assert list(wide_crossed.gene_ids) == list(MANY2MANY_MICE)
        assert wide_crossed.dropped_partners == ()

        # A **Dropped partner** is one an answer no longer names, whichever step removed
        # it — a Homology type filter or the annotation — so the crossing adds to the
        # count rather than replacing it. Release 116 publishes no cross-species paralogy
        # for these pairs, so the already-filtered answer is written out here instead of
        # asked for from a set that could not produce one.
        filtered = replace(answer, dropped_partners=(ABSENT_PARTNER,))
        narrow_registry = _registry_declaring(tmp_path / "narrow", f"{MANY2MANY_MICE[0]}.3")
        narrow_crossed = resolve_homologs(filtered, narrow_registry, "mine")
        assert narrow_crossed.dropped_partners == (ABSENT_PARTNER, MANY2MANY_MICE[1])

    def test_the_crossing_carries_the_null_quality_scores_it_was_handed(
        self, fake_fetch: FakeFetch, tmp_path: Path
    ) -> None:
        # The measurement rides on every answer, and a crossing is still an answer: a
        # caller who filters on `goc_score` after resolving must be told it is empty for
        # this pair in the same breath, not left to discover it.
        homologs = _built(fake_fetch, "Homo sapiens", "Caenorhabditis elegans")
        answer = homologs.homologs([ONE2MANY_HUMAN])
        registry = _registry_declaring(tmp_path, f"{ONE2MANY_WORMS[0]}.1")

        crossed = resolve_homologs(answer, registry, "mine")

        assert crossed.null_quality_scores == QUALITY_SCORE_COLUMNS
        assert crossed.as_json()["null_quality_scores"] == list(QUALITY_SCORE_COLUMNS)


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


# ---------------------------------------------------------------------------
# Orthology is served and never consumed (ADR-0019)
# ---------------------------------------------------------------------------


def _imported_modules(path: Path) -> set[str]:
    """Return every module ``path`` names in an import, read off its source.

    The whole tree of the file is walked rather than its top level, so an import deferred
    into a function body or hidden under ``TYPE_CHECKING`` counts too — a ban that only
    the laziest evasion defeats is not a structural guarantee. Read rather than observed,
    since importing a module loads its package and ``sys.modules`` would then say more
    than any one module asked for.
    """
    found: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            found.add(node.module or "")
    return found


def test_nothing_this_package_publishes_is_derived_through_homology() -> None:
    """No module outside ``genome.homology`` may import it, so no shipped table can use it.

    The structural half of ADR-0019: a **TF gene table**, a **Cofactor table** or any list
    this package ships cannot be derived through homology if nothing that builds one can
    reach it — not even ``genome/__init__.py``, which is why the subpackage is addressed
    as ``genome.homology`` and never re-exported at the top level.

    ``cli.py`` is the one exemption, and it is sound because *consume* here means
    *derive a claim of this package's own*. The CLI publishes no table: it parses
    arguments, makes the one call the caller asked for and renders the answer, so
    ``genome homologs`` is the served half of ADR-0019 reaching a shell rather than
    Python. Nothing it prints outlives the process, and no list that ships in the wheel
    passes through it. Every other module keeps the ban whole — a reader of this list
    should be able to say, of any name added to it, which table it could not build.
    """
    served = {PACKAGE / "cli.py"}
    offenders = sorted(
        str(path.relative_to(PACKAGE))
        for path in PACKAGE.rglob("*.py")
        if path.parent.name != "homology"
        and path not in served
        and any(name.startswith("genome.homology") for name in _imported_modules(path))
    )
    assert offenders == []


@pytest.mark.parametrize("context", ["homology", "xref"])
def test_neither_cross_species_context_imports_the_tf_half(context: str) -> None:
    """*Orthology → TF gene* and *Xref → TF gene* both run one way, held structurally.

    The context map says of the first that it is "a prohibition rather than a call …
    neither half imports the other", and of the second that "the xref half reads no
    census". One guard for both edges, and it lives beside its mirror image above — that
    one lets nothing outside ``genome.homology`` import *it*, this one lets neither the
    Orthology half nor the Xref half import ``genome.tf``.

    What made this fail before it was written: ``species_slug``, a file-naming helper that
    happened to be defined in ``genome.tf.gene.census`` and is now in ``genome.assembly.metadata``,
    the Assembly context's module that every context may read. A helper is not a concept,
    and an import edge does not care about the difference.
    """
    offenders = sorted(
        str(path.relative_to(PACKAGE))
        for path in (PACKAGE / context).rglob("*.py")
        if any(name.startswith("genome.tf") for name in _imported_modules(path))
    )
    assert offenders == []


def test_the_shipped_homology_directory_holds_provenance_and_no_homology_data() -> None:
    """Data is downloaded, never redistributed, which is what keeps its checksum honest.

    It is also what makes *orthology is served and never consumed* a property of the wheel
    rather than a promise: there is no homology table in here for anything to be built
    from.
    """
    directory = Path(__file__).resolve().parents[2] / "src" / "genome" / "data" / "homology"

    assert sorted(entry.name for entry in directory.iterdir()) == [
        "ATTRIBUTION.md",
        "homology_metadata.tsv",
    ]
