"""Tests for genome.xref — id normalisation, the Alliance reader, and the Xref set.

Offline throughout: the ``fake_fetch`` fixture stands in for the package's one fetch step
and serves ``tests/data/xref/alliance_genecrossreference_tiny.tsv.gz``, and every URL
asserted here is the one the package *built*, read back off the recorded call rather than
off the network. The autouse data-root fixture puts the set under the test's own directory,
so the layout is exercised for real.

Everything is driven through the two seams the design names: the :class:`XrefSet`
constructor and its two verbs. Nothing here asserts which private helper ran, and nothing
asserts more about the on-disk layout than a caller can observe — the path the set says it
read, and that the stored form really is a plain gzipped TSV, which is a promise to a
collaborator who does not use Python rather than an implementation detail.

The fixture is 14 whole genes of the real Alliance release 9.0.0 file, copied verbatim with
their duplicate rows intact. The curated row's pinned checksum is over the *whole*
publisher file, so a fixture that is 14 genes of it cannot match: the ``pinned_to_fixture``
fixture re-stamps the shipped rows with the fixture's own digest, exactly as the JASPAR
suite re-points its motif counts, and one test deliberately leaves it off to prove that a
file that does not match its pin is refused.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Protocol

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from genome.io import fetch as fetch_mod
from genome.io.completion import (
    RECORD_NAME,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
)
from genome.io.gtf import _gene_id_stem as annotation_gene_id_stem
from genome.io.results import ResolvedStems, ResolvedXrefIds
from genome.xref import (
    ALLIANCE,
    ENSEMBL,
    ENTREZ,
    HGNC,
    MGI,
    NAMESPACES,
    UNIPROT,
    WORMBASE,
    AllianceFileError,
    NamespaceNotCarriedError,
    NoXrefSetError,
    XrefSet,
    XrefSetNotDownloadedError,
    XrefTableError,
    gene_id_stem,
    lookup_xref,
    normalise_id,
    xref_data_dir,
    xref_releases,
    xref_set_dir,
    xref_slice_name,
    xref_sources,
    xref_species,
    xref_table,
)
from genome.xref import metadata as metadata_mod
from genome.xref.alliance import ALLIANCE_COLUMNS, read_alliance
from genome.xref.xref import SLICE_COLUMNS, parse_slice, xref_prepare_command

from .conftest import FakeFetch

#: The committed fixture, cut whole-gene from Alliance release 9.0.0. See tests/data/README.md.
FIXTURE = "xref/alliance_genecrossreference_tiny.tsv.gz"

#: The pinned release every test names, which is the only one the shipped table lists.
RELEASE = "9.0.0"

#: What the fixture holds per species, once sliced: the namespaces it carries, how many
#: ``(namespace, id, stem)`` triples survive deduplication, and how many stems that is.
#: Every number is counted off the committed bytes, so a fixture edited without the table
#: fails loudly rather than quietly changing what the verb tests mean.
FIXTURE_SLICES: dict[str, tuple[tuple[str, ...], int, int]] = {
    "Homo sapiens": ((ENSEMBL, ENTREZ, UNIPROT, HGNC), 24, 6),
    "Mus musculus": ((ENSEMBL, ENTREZ, UNIPROT, MGI), 29, 5),
    "Caenorhabditis elegans": ((ENSEMBL, ENTREZ, UNIPROT, WORMBASE), 28, 6),
}

#: Human genes the fixture carries, by HGNC id, and the **Gene id stem**s each names. Two
#: of the five name two stems, which is the collision the answer shape exists to keep.
HUMAN_GENES: dict[str, tuple[str, ...]] = {
    "HGNC:1100": ("ENSG00000012048",),
    "HGNC:11998": ("ENSG00000141510",),
    "HGNC:13666": ("ENSG00000094914", "ENSG00000291836"),
    "HGNC:7622": ("ENSG00000196132", "ENSG00000276876"),
}

#: One real human gene the Alliance lists with no ``ENSEMBL:`` cross-reference at all, so
#: it has no hub and appears nowhere in the slice.
HUMAN_GENE_WITHOUT_A_HUB = "HGNC:10041"


def _fixture_md5(path: Path) -> str:
    """The md5 of the fixture's **unpacked** bytes, the way Alliance publishes its own."""
    with gzip.open(path, "rb") as handle:
        # md5 because that is the algorithm Alliance publishes, not a choice made here.
        return hashlib.md5(handle.read()).hexdigest()


@pytest.fixture
def fixture_path(data_dir: Path) -> Path:
    """The committed Alliance fixture."""
    return data_dir / FIXTURE


@pytest.fixture
def fixture_lines(fixture_path: Path) -> list[str]:
    """The committed Alliance fixture, unpacked into lines."""
    with gzip.open(fixture_path, "rt", encoding="utf-8") as handle:
        return handle.read().splitlines()


@pytest.fixture
def pinned_to_fixture(monkeypatch: pytest.MonkeyPatch, fixture_path: Path) -> None:
    """Pin the curated rows to the fixture's digest, so any species may be prepared.

    The checksum check is what holds a truncated download to be an error rather than a
    quietly short answer, so it is never switched off — only pointed at what the fake fetch
    actually serves. Every other cell of every shipped row survives, including the real URL.
    """
    digest = f"md5:{_fixture_md5(fixture_path)}"
    rows = tuple(replace(row, source_checksum=digest) for row in xref_table())
    monkeypatch.setattr(metadata_mod, "xref_table", lambda: rows)


@pytest.fixture
def served(fake_fetch: FakeFetch, pinned_to_fixture: None) -> FakeFetch:
    """A fetch step serving the Alliance fixture, with the checksum pinned to it."""
    fake_fetch.serve(FIXTURE)
    return fake_fetch


class ServeSource(Protocol):
    """Serves Alliance-shaped text as the publisher's file, checksum and all."""

    def __call__(self, rows: list[str]) -> None:
        """Gzip ``rows`` under the header, pin the curated rows to it, and serve it."""
        ...


@pytest.fixture
def serve_source(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> ServeSource:
    """Return a helper that serves Alliance-shaped rows as if the publisher had them.

    For the handful of cases the committed fixture cannot show, because the publisher's
    own file does not contain them — a species it dropped, or an id spelled with a version.
    Everything else is driven off the real bytes.
    """

    def serve(rows: list[str]) -> None:
        payload = "".join(f"{row}\n" for row in ["\t".join(ALLIANCE_COLUMNS), *rows]).encode()
        path = tmp_path / "served.tsv.gz"
        with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
            out.write(payload)
        digest = f"md5:{hashlib.md5(payload).hexdigest()}"
        pinned = tuple(replace(row, source_checksum=digest) for row in xref_table())
        monkeypatch.setattr(metadata_mod, "xref_table", lambda: pinned)
        fake_fetch.serve(path)

    return serve


# ---------------------------------------------------------------------------
# The committed bytes: everything tests/data/README.md says about the fixture
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    def test_it_is_the_publishers_own_columns(self, fixture_lines: list[str]) -> None:
        header = next(line for line in fixture_lines if not line.startswith("#"))
        assert tuple(header.split("\t")) == ALLIANCE_COLUMNS

    def test_it_carries_the_three_species_and_no_others(self, fixture_lines: list[str]) -> None:
        taxa = {line.split("\t")[4] for line in fixture_lines if not line.startswith("#")}
        assert taxa == {"NCBITaxon:9606", "NCBITaxon:10090", "NCBITaxon:6239", "TaxonID"}

    def test_the_duplication_is_on_the_key_and_never_on_the_whole_row(
        self, fixture_lines: list[str]
    ) -> None:
        # The same pair recurs once per page the Alliance links it from, and the page is a
        # column, so no two *rows* are identical while a third of the *keys* repeat.
        # Counted on the whole release 9.0.0 file: 2,659,704 rows reduce to 1,811,267
        # distinct (GeneID, GlobalCrossReferenceID, TaxonID) — 31.9% redundant, and a
        # whole-row `uniq` removes none of it. The fixture is real rows, so it says the same.
        rows = [
            line
            for line in fixture_lines
            if not line.startswith("#") and not line.startswith("GeneID\t")
        ]
        keys = [tuple(line.split("\t")[i] for i in (0, 1, 4)) for line in rows]
        assert len(rows) == len(set(rows))
        assert len(keys) > len(set(keys))
        assert ("HGNC:1100", "NCBI_Gene:672", "NCBITaxon:9606") in keys

    def test_it_carries_a_human_gene_with_no_ensembl_cross_reference(
        self, fixture_lines: list[str]
    ) -> None:
        rows = [line for line in fixture_lines if line.startswith(f"{HUMAN_GENE_WITHOUT_A_HUB}\t")]
        assert rows
        assert not any("ENSEMBL:" in row for row in rows)

    def test_it_carries_worms_symbol_rows_under_the_authority_prefix(
        self, fixture_lines: list[str]
    ) -> None:
        # WB:WBGene00000001 -> WB:aap-1 is a *symbol* under the authority's own prefix, and
        # the trap the reader exists to sidestep.
        assert any(line.startswith("WB:WBGene00000001\tWB:aap-1\t") for line in fixture_lines)


# ---------------------------------------------------------------------------
# Normalisation, which is the first thing tested
# ---------------------------------------------------------------------------


class TestGeneIdStem:
    @pytest.mark.parametrize(
        ("gene_id", "stem"),
        [
            ("ENSG00000141510.18", "ENSG00000141510"),
            ("ENSG00000141510", "ENSG00000141510"),
            ("ENSMUSG00000059552.16", "ENSMUSG00000059552"),
            ("WBGene00000001", "WBGene00000001"),
            ("ENSG00000182378.14_PAR_Y", "ENSG00000182378"),
        ],
    )
    def test_the_version_goes_and_nothing_else_does(self, gene_id: str, stem: str) -> None:
        assert gene_id_stem(gene_id) == stem

    def test_a_worm_gene_id_is_its_own_stem(self) -> None:
        # Unversioned as published, so the hop is the identity function.
        assert gene_id_stem("WBGene00000912") == "WBGene00000912"

    @given(st.text())
    def test_it_is_idempotent(self, gene_id: str) -> None:
        assert gene_id_stem(gene_id_stem(gene_id)) == gene_id_stem(gene_id)

    @given(st.text(alphabet=st.characters(blacklist_characters="."), min_size=1))
    def test_an_unversioned_id_is_its_own_stem(self, gene_id: str) -> None:
        assert gene_id_stem(gene_id) == gene_id

    @given(st.text())
    def test_it_agrees_with_the_annotation_half(self, gene_id: str) -> None:
        # The two sides of the crossing must reduce a gene id the same way, or a stem
        # answered here joins to nothing over there and says nothing about it.
        assert gene_id_stem(gene_id) == annotation_gene_id_stem(gene_id)


class TestNormaliseId:
    @pytest.mark.parametrize(
        ("spelled", "namespace", "canonical"),
        [
            ("HGNC:1100", HGNC, "HGNC:1100"),
            ("1100", HGNC, "HGNC:1100"),
            ("hgnc:1100", HGNC, "HGNC:1100"),
            ("MGI:MGI:88276", MGI, "MGI:88276"),
            ("UniProtKB:P38398", UNIPROT, "P38398"),
            ("P38398", UNIPROT, "P38398"),
            ("NCBI_Gene:672", ENTREZ, "672"),
            ("672", ENTREZ, "672"),
            ("  ENSEMBL:ENSG00000141510.18  ", ENSEMBL, "ENSG00000141510"),
            ("WB:WBGene00000001", WORMBASE, "WBGene00000001"),
        ],
    )
    def test_every_published_spelling_lands_on_one(
        self, spelled: str, namespace: str, canonical: str
    ) -> None:
        assert normalise_id(spelled, namespace) == canonical

    def test_a_namespace_it_does_not_know_still_drops_the_version(self) -> None:
        assert normalise_id("SOMETHING.3", "not-a-namespace") == "SOMETHING"

    @given(st.text(), st.sampled_from(NAMESPACES))
    def test_it_is_idempotent(self, identifier: str, namespace: str) -> None:
        once = normalise_id(identifier, namespace)
        assert normalise_id(once, namespace) == once


# ---------------------------------------------------------------------------
# The curated table
# ---------------------------------------------------------------------------


class TestXrefTable:
    def test_the_three_species_ship_with_alliance_as_the_default(self) -> None:
        assert xref_species() == ("Homo sapiens", "Mus musculus", "Caenorhabditis elegans")
        for species in xref_species():
            assert xref_sources(species) == (ALLIANCE,)
            assert xref_releases(species, ALLIANCE) == (RELEASE,)
            assert lookup_xref(species).source == ALLIANCE
            assert lookup_xref(species).default is True

    def test_every_row_pins_a_publisher_a_version_a_url_and_a_checksum(self) -> None:
        for row in xref_table():
            assert row.publisher == "Alliance of Genome Resources"
            assert row.version == RELEASE
            assert row.url.startswith("https://download.alliancegenome.org/")
            algorithm, _, digest = row.source_checksum.partition(":")
            assert algorithm == "md5"
            assert len(digest) == 32
            assert set(digest) <= set("0123456789abcdef")
            assert row.pubmed_id == 38552170

    def test_the_taxids_are_the_ones_the_publishers_file_uses(self) -> None:
        assert {row.species: row.ncbi_taxid for row in xref_table()} == {
            "Homo sapiens": 9606,
            "Mus musculus": 10090,
            "Caenorhabditis elegans": 6239,
        }

    def test_a_species_slug_names_the_same_row_as_the_species(self) -> None:
        assert lookup_xref("homo_sapiens") == lookup_xref("Homo sapiens")

    def test_an_unsupported_species_names_the_species_that_have_a_set(self) -> None:
        with pytest.raises(NoXrefSetError) as raised:
            lookup_xref("Danio rerio")
        for species in xref_species():
            assert species in str(raised.value)

    def test_an_unknown_source_names_the_sources_there_are(self) -> None:
        with pytest.raises(NoXrefSetError, match="alliance"):
            lookup_xref("Homo sapiens", "ncbi")

    def test_an_unknown_release_names_the_releases_there_are(self) -> None:
        with pytest.raises(NoXrefSetError, match=RELEASE):
            lookup_xref("Homo sapiens", ALLIANCE, "1.0")

    def test_no_release_named_answers_with_the_newest(self) -> None:
        rows = (
            replace(lookup_xref("Homo sapiens"), release="8.2.0"),
            replace(lookup_xref("Homo sapiens"), release="9.0.0"),
        )
        assert lookup_xref("Homo sapiens", ALLIANCE, table=rows).release == "9.0.0"

    def test_a_species_with_two_sources_and_no_default_names_them_both(self) -> None:
        rows = (
            replace(lookup_xref("Homo sapiens"), source="alliance", default=False),
            replace(lookup_xref("Homo sapiens"), source="ensembl", default=False),
        )
        with pytest.raises(NoXrefSetError, match="no default xref source"):
            lookup_xref("Homo sapiens", table=rows)

    def test_the_attribution_names_the_publisher_the_release_and_the_paper(self) -> None:
        line = lookup_xref("Mus musculus").attribution()
        assert "Alliance of Genome Resources 9.0.0" in line
        assert "PMID 38552170" in line


# ---------------------------------------------------------------------------
# The Alliance reader
# ---------------------------------------------------------------------------


class TestAllianceReader:
    @pytest.mark.parametrize("species", sorted(FIXTURE_SLICES))
    def test_it_reads_the_fixture_into_the_slice_the_table_records(
        self, fixture_lines: list[str], species: str
    ) -> None:
        namespaces, triples, stems = FIXTURE_SLICES[species]
        read = read_alliance(
            fixture_lines, ncbi_taxid=lookup_xref(species).ncbi_taxid, origin="fixture"
        )
        assert len(read) == triples
        assert {namespace for namespace, _id, _stem in read} == set(namespaces)
        assert len({stem for _ns, _id, stem in read}) == stems

    def test_it_takes_the_authority_id_from_the_gene_column_and_never_from_a_symbol(
        self, fixture_lines: list[str]
    ) -> None:
        # WB:WBGene00000001 -> WB:aap-1 is a symbol. Reading the authority off the
        # cross-reference column would key the gene by 'aap-1' and call it a WormBase id.
        read = read_alliance(fixture_lines, ncbi_taxid=6239, origin="fixture")
        worm_ids = {i for namespace, i, _stem in read if namespace == WORMBASE}
        assert "WBGene00000001" in worm_ids
        assert not any("aap-1" in identifier for identifier in worm_ids)

    def test_the_worm_hop_is_the_identity(self, fixture_lines: list[str]) -> None:
        # Counted on the whole release 9.0.0 file, not assumed: all 46,926 worm genes carry
        # an ENSEMBL cross-reference whose id is the WBGene id, with zero differing. So a
        # **Gene id stem** and a WormBase gene id are the same string, and the hop worm data
        # makes into this package's answers is the identity function.
        read = read_alliance(fixture_lines, ncbi_taxid=6239, origin="fixture")
        worm = [(identifier, stem) for ns, identifier, stem in read if ns == WORMBASE]
        assert worm
        assert all(identifier == stem for identifier, stem in worm)
        assert {identifier for identifier, _stem in worm} == {
            identifier for ns, identifier, _stem in read if ns == ENSEMBL
        }

    def test_a_gene_with_no_ensembl_cross_reference_contributes_nothing(
        self, fixture_lines: list[str]
    ) -> None:
        read = read_alliance(fixture_lines, ncbi_taxid=9606, origin="fixture")
        assert HUMAN_GENE_WITHOUT_A_HUB not in {i for _ns, i, _stem in read}

    def test_it_deduplicates_the_publishers_repeated_rows(self, fixture_lines: list[str]) -> None:
        read = read_alliance(fixture_lines, ncbi_taxid=9606, origin="fixture")
        assert len(read) == len(set(read))

    def test_a_taxon_the_file_does_not_carry_reads_as_nothing(
        self, fixture_lines: list[str]
    ) -> None:
        assert read_alliance(fixture_lines, ncbi_taxid=7955, origin="fixture") == ()

    def test_a_missing_header_names_the_columns_it_wanted(self) -> None:
        with pytest.raises(AllianceFileError, match="GeneID"):
            read_alliance(["# only a comment"], ncbi_taxid=9606, origin="somewhere.tsv")

    def test_a_respelled_header_refuses_rather_than_reading_by_position(self) -> None:
        rows = ["\t".join(("Gene", *ALLIANCE_COLUMNS[1:])), "a\tb\tc\td\te"]
        with pytest.raises(AllianceFileError, match=r"somewhere\.tsv"):
            read_alliance(rows, ncbi_taxid=9606, origin="somewhere.tsv")

    def test_an_unknown_species_authority_names_the_ones_it_knows(self) -> None:
        rows = [
            "\t".join(ALLIANCE_COLUMNS),
            "ZFIN:ZDB-1\tENSEMBL:ENSDARG1\tu\tgeneric\tNCBITaxon:9606",
        ]
        with pytest.raises(AllianceFileError, match="HGNC"):
            read_alliance(rows, ncbi_taxid=9606, origin="somewhere.tsv")


# ---------------------------------------------------------------------------
# Preparing the set: construction fetches once and re-reads thereafter
# ---------------------------------------------------------------------------


class TestXrefSetPreparation:
    def test_the_url_it_asked_for_is_the_curated_rows(self, served: FakeFetch) -> None:
        XrefSet("Homo sapiens")
        assert served.last.url == lookup_xref("Homo sapiens").url

    def test_the_slice_lands_where_the_layout_puts_it(
        self, served: FakeFetch, liulab_data: Path
    ) -> None:
        human = XrefSet("Homo sapiens")
        expected = (
            liulab_data
            / "xref"
            / ALLIANCE
            / RELEASE
            / "homo_sapiens"
            / "homo_sapiens.xref_table.tsv.gz"
        )
        assert human.path == expected
        assert expected.is_file()
        assert xref_set_dir("Homo sapiens", ALLIANCE, RELEASE) == expected.parent
        assert xref_slice_name("Homo sapiens") == expected.name

    def test_xref_data_is_a_sibling_of_the_assembly_tree(
        self, served: FakeFetch, liulab_data: Path
    ) -> None:
        XrefSet("Homo sapiens")
        assert xref_data_dir() == liulab_data / "xref"
        assert (liulab_data / "xref").is_dir()
        assert not (liulab_data / "genome").exists()

    def test_a_second_construction_fetches_nothing(self, served: FakeFetch) -> None:
        first = XrefSet("Homo sapiens")
        second = XrefSet("Homo sapiens")
        assert len(served.calls) == 1
        assert first.path == second.path
        assert first.to_stems(["7157"], ENTREZ) == second.to_stems(["7157"], ENTREZ)

    def test_another_species_is_another_set(self, served: FakeFetch) -> None:
        human = XrefSet("Homo sapiens")
        mouse = XrefSet("Mus musculus")
        assert len(served.calls) == 2
        assert human.path != mouse.path

    def test_an_explicit_cache_dir_names_the_directory_itself(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        elsewhere = tmp_path / "somewhere-else"
        human = XrefSet("Homo sapiens", cache_dir=elsewhere)
        assert human.path.parent == elsewhere
        assert not (xref_data_dir() / ALLIANCE).exists()

    def test_the_progress_bar_is_the_callers_choice(self, served: FakeFetch) -> None:
        XrefSet("Homo sapiens", progressbar=False)
        assert served.last.progressbar is False

    def test_it_knows_which_set_it_is(self, served: FakeFetch) -> None:
        human = XrefSet("Homo sapiens")
        assert (human.species, human.source, human.release) == ("Homo sapiens", ALLIANCE, RELEASE)
        assert human.source_url == lookup_xref("Homo sapiens").url

    def test_repr_names_the_set_and_how_many_stems_it_holds(self, served: FakeFetch) -> None:
        assert repr(XrefSet("Homo sapiens")) == (
            f"XrefSet(species='Homo sapiens', source='alliance', release='{RELEASE}', stems=6)"
        )

    @pytest.mark.parametrize("species", sorted(FIXTURE_SLICES))
    def test_every_species_carries_the_namespaces_the_source_has_for_it(
        self, served: FakeFetch, species: str
    ) -> None:
        namespaces, _triples, stems = FIXTURE_SLICES[species]
        prepared = XrefSet(species)
        assert prepared.namespaces == namespaces
        assert len(prepared) == stems

    def test_an_unsupported_species_raises_before_anything_is_fetched(
        self, served: FakeFetch
    ) -> None:
        with pytest.raises(NoXrefSetError, match="Homo sapiens"):
            XrefSet("Danio rerio")
        assert served.calls == []

    def test_a_file_that_does_not_match_its_pin_is_refused(self, fake_fetch: FakeFetch) -> None:
        # No `pinned_to_fixture` here: the shipped row pins the whole 25 MB publisher file,
        # and 14 genes of it is not that file. A truncated download is not a smaller
        # release, and slicing one would answer with silently fewer genes.
        fake_fetch.serve(FIXTURE)
        with pytest.raises(XrefTableError, match="hashes to"):
            XrefSet("Homo sapiens")

    def test_a_file_carrying_no_row_for_the_species_names_the_url(
        self, serve_source: ServeSource
    ) -> None:
        serve_source(["MGI:98834\tENSEMBL:ENSMUSG00000059552\tu\tgeneric\tNCBITaxon:10090"])
        with pytest.raises(XrefTableError, match="NCBITaxon:9606"):
            XrefSet("Homo sapiens")


class TestVersionsOnTheSourceSide:
    """A publisher that spells a version is joined to one that does not, and neither loses.

    **Alliance publishes none**: 0 of the 43,867 human ids in release 9.0.0 carry a dot, so
    the committed fixture cannot show this and these two tests write the rows themselves.
    They are here because the trap is real elsewhere — NCBI's ``gene2ensembl`` writes the
    gene id bare and the transcript id versioned in the same row — and because the reader is
    where a second source (#150) will meet it.
    """

    def test_the_reader_stems_a_versioned_id_on_ingest(self) -> None:
        rows = [
            "\t".join(ALLIANCE_COLUMNS),
            "HGNC:11998\tENSEMBL:ENSG00000141510.18\tu\tgeneric\tNCBITaxon:9606",
            "HGNC:11998\tUniProtKB:P04637.2\tu\tgeneric\tNCBITaxon:9606",
        ]
        assert read_alliance(rows, ncbi_taxid=9606, origin="somewhere.tsv") == (
            (ENSEMBL, "ENSG00000141510", "ENSG00000141510"),
            (HGNC, "HGNC:11998", "ENSG00000141510"),
            (UNIPROT, "P04637", "ENSG00000141510"),
        )

    def test_a_versioned_source_and_a_bare_ask_still_meet(self, serve_source: ServeSource) -> None:
        serve_source(["HGNC:11998\tENSEMBL:ENSG00000141510.18\tu\tgeneric\tNCBITaxon:9606"])
        human = XrefSet("Homo sapiens")
        assert human.to_stems(["HGNC:11998"], HGNC).resolved == {"HGNC:11998": ("ENSG00000141510",)}
        assert human.from_stems(["ENSG00000141510"], HGNC).resolved == {
            "ENSG00000141510": ("HGNC:11998",)
        }


class TestXrefSetIsNotDownloaded:
    def test_a_set_that_cannot_be_fetched_names_the_next_action(
        self, monkeypatch: pytest.MonkeyPatch, pinned_to_fixture: None
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: Any) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        with pytest.raises(XrefSetNotDownloadedError) as raised:
            XrefSet("Homo sapiens")
        message = str(raised.value)
        assert xref_prepare_command("Homo sapiens", ALLIANCE, RELEASE) in message
        assert "login node" in message

    def test_an_interrupted_download_leaves_a_directory_that_reads_as_unfinished(
        self, monkeypatch: pytest.MonkeyPatch, pinned_to_fixture: None, served: FakeFetch
    ) -> None:
        directory = xref_set_dir("Homo sapiens", ALLIANCE, RELEASE)
        directory.mkdir(parents=True)
        (directory / xref_slice_name("Homo sapiens")).write_bytes(b"half a file")
        with pytest.raises(UnfinishedRegistrationError, match=RECORD_NAME):
            XrefSet("Homo sapiens")

    def test_what_an_interrupted_download_left_in_the_working_area_is_never_adopted(
        self, served: FakeFetch
    ) -> None:
        directory = xref_set_dir("Homo sapiens", ALLIANCE, RELEASE)
        work = directory / ".work"
        work.mkdir(parents=True)
        (work / "GENECROSSREFERENCE_COMBINED_11.tsv.gz").write_bytes(b"not a gzip")
        human = XrefSet("Homo sapiens")
        assert len(human) == FIXTURE_SLICES["Homo sapiens"][2]
        assert not work.exists()


class TestCompletionMarker:
    def test_it_records_where_the_bytes_came_from_and_both_checksums(
        self, served: FakeFetch, fixture_path: Path
    ) -> None:
        human = XrefSet("Homo sapiens")
        record = read_record(human.path.parent)
        assert record is not None
        assert record.source_url == lookup_xref("Homo sapiens").url
        # The publisher's own checksum is provenance: what is stored is a derived slice.
        assert record.details["source_checksum"] == f"md5:{_fixture_md5(fixture_path)}"
        # The derived slice's own checksum is the integrity check.
        with gzip.open(human.path, "rb") as handle:
            assert record.sha256 == hashlib.sha256(handle.read()).hexdigest()
        assert record.details["release"] == RELEASE
        assert record.details["source"] == ALLIANCE
        assert sorted(record.details["namespaces"]) == sorted(FIXTURE_SLICES["Homo sapiens"][0])

    def test_a_marker_that_disagrees_about_the_checksum_means_unfinished(
        self, served: FakeFetch
    ) -> None:
        human = XrefSet("Homo sapiens")
        marker = human.path.parent / RECORD_NAME
        payload = json.loads(marker.read_text(encoding="utf-8"))
        payload["sha256"] = "0" * 64
        marker.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(RegistrationMismatchError, match="hashes to"):
            XrefSet("Homo sapiens")

    def test_a_marker_that_disagrees_about_the_size_means_unfinished(
        self, served: FakeFetch
    ) -> None:
        human = XrefSet("Homo sapiens")
        human.path.write_bytes(b"")
        with pytest.raises(RegistrationMismatchError, match=RECORD_NAME):
            XrefSet("Homo sapiens")


class TestStoredForm:
    def test_the_stored_slice_is_a_plain_gzipped_tsv(self, served: FakeFetch) -> None:
        human = XrefSet("Homo sapiens")
        with gzip.open(human.path, "rt", encoding="utf-8") as handle:
            assert handle.readline() == "\t".join(SLICE_COLUMNS) + "\n"

    def test_a_collaborator_who_does_not_use_python_can_read_it(self, served: FakeFetch) -> None:
        # No library of this package's is involved: three tab-separated columns and a
        # header, which is all R's read.delim or a shell's cut needs.
        human = XrefSet("Homo sapiens")
        frame = pd.read_csv(human.path, sep="\t", dtype=str)
        assert list(frame.columns) == list(SLICE_COLUMNS)
        assert len(frame) == FIXTURE_SLICES["Homo sapiens"][1]

    def test_two_machines_slicing_one_release_write_the_same_bytes(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        first = XrefSet("Homo sapiens", cache_dir=tmp_path / "one")
        second = XrefSet("Homo sapiens", cache_dir=tmp_path / "two")
        assert first.path.read_bytes() == second.path.read_bytes()

    def test_a_slice_this_package_did_not_write_is_refused(self) -> None:
        with pytest.raises(XrefTableError, match="namespace"):
            parse_slice("gene\tid\n", origin="somewhere.tsv")

    def test_a_slice_naming_an_unknown_namespace_is_refused(self) -> None:
        text = "\t".join(SLICE_COLUMNS) + "\nrefseq\tNM_1\tENSG1\n"
        with pytest.raises(XrefTableError, match="refseq"):
            parse_slice(text, origin="somewhere.tsv")

    def test_a_short_row_is_refused(self) -> None:
        text = "\t".join(SLICE_COLUMNS) + "\nentrez\t7157\n"
        with pytest.raises(XrefTableError, match="line 2"):
            parse_slice(text, origin="somewhere.tsv")


# ---------------------------------------------------------------------------
# The two verbs
# ---------------------------------------------------------------------------


@pytest.fixture
def human(served: FakeFetch) -> XrefSet:
    """The human set, prepared from the fixture."""
    return XrefSet("Homo sapiens")


@pytest.fixture
def worm(served: FakeFetch) -> XrefSet:
    """The worm set, prepared from the fixture."""
    return XrefSet("Caenorhabditis elegans")


@pytest.fixture
def mouse(served: FakeFetch) -> XrefSet:
    """The mouse set, prepared from the fixture."""
    return XrefSet("Mus musculus")


class TestToStems:
    @pytest.mark.parametrize(
        ("namespace", "identifier"),
        [(ENTREZ, "7157"), (UNIPROT, "P04637"), (HGNC, "HGNC:11998"), (ENSEMBL, "ENSG00000141510")],
    )
    def test_every_namespace_reaches_the_hub(
        self, human: XrefSet, namespace: str, identifier: str
    ) -> None:
        answer = human.to_stems([identifier], namespace)
        assert answer.resolved == {identifier: ("ENSG00000141510",)}
        assert answer.unresolved == ()
        assert (answer.species, answer.source, answer.release, answer.namespace) == (
            "Homo sapiens",
            ALLIANCE,
            RELEASE,
            namespace,
        )

    def test_mouse_reaches_the_hub_through_its_own_authority(self, mouse: XrefSet) -> None:
        answer = mouse.to_stems(["MGI:98834"], MGI)
        assert answer.resolved == {"MGI:98834": ("ENSMUSG00000059552",)}

    def test_worm_reaches_the_hub_through_its_own_authority(self, worm: XrefSet) -> None:
        answer = worm.to_stems(["WBGene00000001"], WORMBASE)
        # The WormBase gene id *is* the Ensembl stable gene id, so the hop is the identity.
        assert answer.resolved == {"WBGene00000001": ("WBGene00000001",)}

    def test_a_foreign_id_naming_two_stems_answers_with_both(self, human: XrefSet) -> None:
        answer = human.to_stems(["HGNC:13666"], HGNC)
        assert answer.resolved == {"HGNC:13666": HUMAN_GENES["HGNC:13666"]}

    def test_a_uniprot_accession_naming_four_worm_genes_answers_with_all_four(
        self, worm: XrefSet
    ) -> None:
        answer = worm.to_stems(["P05634"], UNIPROT)
        assert answer.resolved["P05634"] == (
            "WBGene00003425",
            "WBGene00003432",
            "WBGene00003449",
            "WBGene00003463",
        )

    def test_ids_that_named_nothing_ride_back_in_ask_order(self, human: XrefSet) -> None:
        answer = human.to_stems(["999999998", "7157", "999999999"], ENTREZ)
        assert answer.unresolved == ("999999998", "999999999")
        assert list(answer.resolved) == ["7157"]

    def test_no_resolved_value_is_ever_empty(self, human: XrefSet) -> None:
        answer = human.to_stems(["7157", "999999999", ""], ENTREZ)
        assert all(answer.resolved.values())

    def test_ask_order_is_preserved(self, human: XrefSet) -> None:
        asked = ["HGNC:13666", "HGNC:11998", "HGNC:1100", "HGNC:7622"]
        assert list(human.to_stems(asked, HGNC).resolved) == asked

    def test_a_versioned_and_an_unversioned_spelling_resolve_identically(
        self, human: XrefSet
    ) -> None:
        answer = human.to_stems(["ENSG00000141510.18", "ENSG00000141510"], ENSEMBL)
        assert answer.resolved == {
            "ENSG00000141510.18": ("ENSG00000141510",),
            "ENSG00000141510": ("ENSG00000141510",),
        }

    def test_the_authoritys_two_spellings_resolve_identically(self, human: XrefSet) -> None:
        assert human.to_stems(["11998"], HGNC).gene_id_stems == (
            human.to_stems(["HGNC:11998"], HGNC).gene_id_stems
        )

    def test_duplicate_source_rows_do_not_produce_duplicate_answers(self, human: XrefSet) -> None:
        # The publisher lists HGNC:1100 -> NCBI_Gene:672 twice, under two pages.
        assert human.to_stems(["672"], ENTREZ).resolved == {"672": ("ENSG00000012048",)}

    def test_an_id_asked_for_twice_is_answered_once(self, human: XrefSet) -> None:
        answer = human.to_stems(["7157", "7157"], ENTREZ)
        assert list(answer.resolved) == ["7157"]

    def test_a_gene_with_no_hub_reaches_nothing(self, human: XrefSet) -> None:
        answer = human.to_stems([HUMAN_GENE_WITHOUT_A_HUB], HGNC)
        assert answer.unresolved == (HUMAN_GENE_WITHOUT_A_HUB,)


class TestFromStems:
    @pytest.mark.parametrize(
        ("namespace", "expected"),
        [
            (ENTREZ, ("7157",)),
            (UNIPROT, ("P04637",)),
            (HGNC, ("HGNC:11998",)),
            (ENSEMBL, ("ENSG00000141510",)),
        ],
    )
    def test_every_namespace_answers_from_the_hub(
        self, human: XrefSet, namespace: str, expected: tuple[str, ...]
    ) -> None:
        answer = human.from_stems(["ENSG00000141510"], namespace)
        assert answer.resolved == {"ENSG00000141510": expected}
        assert answer.namespace == namespace

    def test_a_stem_naming_two_ids_answers_with_both(self, mouse: XrefSet) -> None:
        answer = mouse.from_stems(["ENSMUSG00000059552"], UNIPROT)
        assert len(answer.resolved["ENSMUSG00000059552"]) == 14
        assert "P02340" in answer.resolved["ENSMUSG00000059552"]

    def test_a_versioned_gene_id_is_accepted_as_its_stem(self, human: XrefSet) -> None:
        assert human.from_stems(["ENSG00000141510.18"], HGNC).resolved == {
            "ENSG00000141510.18": ("HGNC:11998",)
        }

    def test_stems_that_named_nothing_ride_back_in_ask_order(self, human: XrefSet) -> None:
        answer = human.from_stems(["ENSG00000000001", "ENSG00000141510", "ENSG00000000002"], HGNC)
        assert answer.unresolved == ("ENSG00000000001", "ENSG00000000002")

    def test_ask_order_is_preserved(self, worm: XrefSet) -> None:
        asked = ["WBGene00003449", "WBGene00000001", "WBGene00003425"]
        assert list(worm.from_stems(asked, ENTREZ).resolved) == asked

    def test_no_resolved_value_is_ever_empty(self, human: XrefSet) -> None:
        answer = human.from_stems(["ENSG00000141510", "ENSG00000000001"], UNIPROT)
        assert all(answer.resolved.values())

    def test_there_is_no_foreign_to_foreign_verb(self, human: XrefSet) -> None:
        # Entrez to HGNC is two calls and the caller owns the join (ADR-0017).
        assert not [name for name in dir(human) if name.startswith("convert")]
        stems = human.to_stems(["7157"], ENTREZ).gene_id_stems
        assert human.from_stems(stems, HGNC).xref_ids == ["HGNC:11998"]


#: One **Gene id stem** per species that the fixture carries in every namespace, so that a
#: round trip through both verbs can be run for all twelve species-and-namespace pairs.
ROUND_TRIP_STEMS: dict[str, str] = {
    "Homo sapiens": "ENSG00000141510",
    "Mus musculus": "ENSMUSG00000059552",
    "Caenorhabditis elegans": "WBGene00000001",
}


class TestBothVerbsForEveryNamespaceAndSpecies:
    """Twelve combinations: three species, four namespaces each, out and back again.

    The acceptance criterion in the flesh — a stem the fixture carries goes out through
    :meth:`XrefSet.from_stems` and comes back through :meth:`XrefSet.to_stems`, so neither
    verb can pass by answering nothing.
    """

    @pytest.mark.parametrize("species", sorted(ROUND_TRIP_STEMS))
    def test_a_stem_goes_out_and_comes_back_in_every_namespace(
        self, served: FakeFetch, species: str
    ) -> None:
        prepared = XrefSet(species)
        stem = ROUND_TRIP_STEMS[species]
        assert prepared.namespaces == FIXTURE_SLICES[species][0]
        for namespace in prepared.namespaces:
            out = prepared.from_stems([stem], namespace)
            assert out.resolved[stem], namespace
            back = prepared.to_stems(list(out.resolved[stem]), namespace)
            assert set(back.resolved) == set(out.resolved[stem]), namespace
            assert all(stem in stems for stems in back.resolved.values()), namespace


class TestNamespaces:
    def test_a_namespace_the_source_does_not_carry_names_the_ones_it_does(
        self, mouse: XrefSet
    ) -> None:
        with pytest.raises(NamespaceNotCarriedError) as raised:
            mouse.to_stems(["HGNC:11998"], HGNC)
        message = str(raised.value)
        for namespace in mouse.namespaces:
            assert namespace in message

    def test_the_reverse_verb_refuses_the_same_namespace(self, worm: XrefSet) -> None:
        with pytest.raises(NamespaceNotCarriedError, match="wormbase"):
            worm.from_stems(["WBGene00000001"], MGI)

    def test_a_namespace_nobody_has_heard_of_refuses(self, human: XrefSet) -> None:
        with pytest.raises(NamespaceNotCarriedError, match="refseq"):
            human.to_stems(["NM_000546"], "refseq")

    def test_case_is_not_significant(self, human: XrefSet) -> None:
        assert human.to_stems(["7157"], "Entrez").namespace == ENTREZ


class TestAnswerShape:
    def test_to_stems_flattens_every_stem_and_not_one_per_id(self, human: XrefSet) -> None:
        answer = human.to_stems(["HGNC:13666", "HGNC:11998"], HGNC)
        assert answer.gene_id_stems == [
            "ENSG00000094914",
            "ENSG00000291836",
            "ENSG00000141510",
        ]

    def test_from_stems_flattens_every_id_and_not_one_per_stem(self, worm: XrefSet) -> None:
        answer = worm.from_stems(["WBGene00000001", "WBGene00000912"], ENTREZ)
        assert answer.xref_ids == ["172141", "172981"]

    def test_to_stems_json_names_what_produced_it(self, human: XrefSet) -> None:
        rendered = human.to_stems(["7157", "999999999"], ENTREZ).as_json()
        assert rendered == {
            "species": "Homo sapiens",
            "source": ALLIANCE,
            "release": RELEASE,
            "namespace": ENTREZ,
            "resolved": {"7157": ["ENSG00000141510"]},
            "unresolved": ["999999999"],
            "gene_id_stems": ["ENSG00000141510"],
        }
        assert json.loads(json.dumps(rendered)) == rendered

    def test_from_stems_json_names_what_produced_it(self, human: XrefSet) -> None:
        rendered = human.from_stems(["ENSG00000141510"], HGNC).as_json()
        assert rendered["namespace"] == HGNC
        assert rendered["xref_ids"] == ["HGNC:11998"]
        assert json.loads(json.dumps(rendered)) == rendered

    def test_the_answers_are_the_shapes_results_declares(self, human: XrefSet) -> None:
        assert isinstance(human.to_stems(["7157"], ENTREZ), ResolvedStems)
        assert isinstance(human.from_stems(["ENSG00000141510"], HGNC), ResolvedXrefIds)


#: The prepared set is read-only and every generated ask is independent of the last, so a
#: set built once and reused across inputs cannot carry state between them.
_REUSES_THE_SET = settings(suppress_health_check=[HealthCheck.function_scoped_fixture])


class TestAnswerProperties:
    """Invariants over generated asks — the partition, and what it never loses."""

    @_REUSES_THE_SET
    @given(
        asked=st.lists(
            st.sampled_from([*HUMAN_GENES, "HGNC:999999", "999999", HUMAN_GENE_WITHOUT_A_HUB]),
            max_size=12,
        )
    )
    def test_resolved_and_unresolved_partition_the_ask_exactly_once(
        self, human: XrefSet, asked: list[str]
    ) -> None:
        answer = human.to_stems(asked, HGNC)
        once = list(dict.fromkeys(asked))
        assert sorted([*answer.resolved, *answer.unresolved]) == sorted(once)
        assert len({*answer.resolved} & {*answer.unresolved}) == 0
        assert len(answer.unresolved) == len(set(answer.unresolved))

    @_REUSES_THE_SET
    @given(asked=st.lists(st.sampled_from([*HUMAN_GENES, "HGNC:999999"]), max_size=12))
    def test_ask_order_survives_both_buckets(self, human: XrefSet, asked: list[str]) -> None:
        answer = human.to_stems(asked, HGNC)
        once = list(dict.fromkeys(asked))
        assert list(answer.resolved) == [key for key in once if key in answer.resolved]
        assert list(answer.unresolved) == [key for key in once if key not in answer.resolved]

    @_REUSES_THE_SET
    @given(
        stems=st.lists(
            st.sampled_from(
                [stem for stems in HUMAN_GENES.values() for stem in stems] + ["ENSG00000000001"]
            ),
            max_size=12,
        )
    )
    def test_the_reverse_verb_partitions_the_same_way(
        self, human: XrefSet, stems: list[str]
    ) -> None:
        answer = human.from_stems(stems, HGNC)
        once = list(dict.fromkeys(stems))
        assert sorted([*answer.resolved, *answer.unresolved]) == sorted(once)
        assert all(answer.resolved.values())
