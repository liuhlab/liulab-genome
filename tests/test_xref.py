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
from genome.io.prepared import PreparedChecksumError
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
    def test_the_committed_bytes_match_everything_the_readme_says(
        self, fixture_lines: list[str]
    ) -> None:
        header = next(line for line in fixture_lines if not line.startswith("#"))
        assert tuple(header.split("\t")) == ALLIANCE_COLUMNS

        taxa = {line.split("\t")[4] for line in fixture_lines if not line.startswith("#")}
        assert taxa == {"NCBITaxon:9606", "NCBITaxon:10090", "NCBITaxon:6239", "TaxonID"}

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

        hub_rows = [
            line for line in fixture_lines if line.startswith(f"{HUMAN_GENE_WITHOUT_A_HUB}\t")
        ]
        assert hub_rows
        assert not any("ENSEMBL:" in row for row in hub_rows)

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
            ("ENSG00000182378.14_PAR_Y", "ENSG00000182378"),
        ],
    )
    def test_the_version_goes_and_nothing_else_does(self, gene_id: str, stem: str) -> None:
        assert gene_id_stem(gene_id) == stem

    def test_a_worm_gene_id_is_its_own_stem(self) -> None:
        # Unversioned as published, so the hop is the identity function.
        assert gene_id_stem("WBGene00000912") == "WBGene00000912"

    @given(st.text())
    def test_it_is_idempotent_and_agrees_with_the_annotation_half(self, gene_id: str) -> None:
        assert gene_id_stem(gene_id_stem(gene_id)) == gene_id_stem(gene_id)
        # The two sides of the crossing must reduce a gene id the same way, or a stem
        # answered here joins to nothing over there and says nothing about it.
        assert gene_id_stem(gene_id) == annotation_gene_id_stem(gene_id)

    @given(st.text(alphabet=st.characters(blacklist_characters="."), min_size=1))
    def test_an_unversioned_id_is_its_own_stem(self, gene_id: str) -> None:
        assert gene_id_stem(gene_id) == gene_id


class TestNormaliseId:
    @pytest.mark.parametrize(
        ("spelled", "namespace", "canonical"),
        [
            ("1100", HGNC, "HGNC:1100"),
            ("hgnc:1100", HGNC, "HGNC:1100"),
            ("MGI:MGI:88276", MGI, "MGI:88276"),
            ("NCBI_Gene:672", ENTREZ, "672"),
            ("UniProtKB:P38398", UNIPROT, "P38398"),
            ("P38398", UNIPROT, "P38398"),
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

    @pytest.mark.parametrize("spelled", ["7157\r.", " 7157 . "])
    def test_whitespace_hidden_behind_the_version_separator_goes_on_the_first_pass(
        self, spelled: str
    ) -> None:
        # Stripping before stemming leaves it behind, so the id only settled on its second
        # pass — and two spellings of one id that settle after a different number of passes
        # join to nothing and say nothing about it.
        assert normalise_id(spelled, ENTREZ) == "7157"

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
            assert ALLIANCE in xref_sources(species)
            assert xref_releases(species, ALLIANCE) == (RELEASE,)
            assert lookup_xref(species).source == ALLIANCE
            assert lookup_xref(species).default is True

    def test_every_row_pins_full_provenance(self) -> None:
        alliance_rows = [row for row in xref_table() if row.source == ALLIANCE]
        assert len(alliance_rows) == len(xref_species())
        for row in alliance_rows:
            assert row.publisher == "Alliance of Genome Resources"
            assert row.version == RELEASE
            assert row.url.startswith("https://download.alliancegenome.org/")
            assert row.pubmed_id == 38552170

        for row in xref_table():
            algorithm, _, digest = row.source_checksum.partition(":")
            assert algorithm == "md5"
            assert len(digest) == 32
            assert set(digest) <= set("0123456789abcdef")

        assert {row.species: row.ncbi_taxid for row in xref_table() if row.source == ALLIANCE} == {
            "Homo sapiens": 9606,
            "Mus musculus": 10090,
            "Caenorhabditis elegans": 6239,
        }

    def test_lookup_finds_by_slug_and_attributes_the_release(self) -> None:
        assert lookup_xref("homo_sapiens") == lookup_xref("Homo sapiens")
        line = lookup_xref("Mus musculus").attribution()
        assert "Alliance of Genome Resources 9.0.0" in line
        assert "PMID 38552170" in line

    def test_lookup_xref_refuses_and_names_the_next_action(self) -> None:
        with pytest.raises(NoXrefSetError) as unsupported:
            lookup_xref("Danio rerio")
        for species in xref_species():
            assert species in str(unsupported.value)

        with pytest.raises(NoXrefSetError, match="alliance"):
            lookup_xref("Homo sapiens", "ncbi")

        with pytest.raises(NoXrefSetError, match=RELEASE):
            lookup_xref("Homo sapiens", ALLIANCE, "1.0")

        # `116` belongs to the other source. Answering it with 9.0.0 would hand back bytes
        # nobody asked for, under a release string that says they did — which is the whole
        # of what pinning is for.
        two_releases = (
            replace(lookup_xref("Homo sapiens"), source=ALLIANCE, release="9.0.0", default=True),
            replace(lookup_xref("Homo sapiens"), source="ensembl", release="116", default=False),
        )
        with pytest.raises(NoXrefSetError) as no_such_release:
            lookup_xref("Homo sapiens", release="116", table=two_releases)
        assert "9.0.0" in str(no_such_release.value)

        two_sources_no_default = (
            replace(lookup_xref("Homo sapiens"), source="alliance", default=False),
            replace(lookup_xref("Homo sapiens"), source="ensembl", default=False),
        )
        with pytest.raises(NoXrefSetError, match="no default xref source"):
            lookup_xref("Homo sapiens", table=two_sources_no_default)

    def test_a_release_alone_finds_the_source_that_has_it(self) -> None:
        newest = (
            replace(lookup_xref("Homo sapiens"), release="8.2.0"),
            replace(lookup_xref("Homo sapiens"), release="9.0.0"),
        )
        assert lookup_xref("Homo sapiens", ALLIANCE, table=newest).release == "9.0.0"

        # Built here rather than read off the shipped table, so this says what `lookup_xref`
        # does rather than what today's rows happen to make it do.
        two_releases = (
            replace(lookup_xref("Homo sapiens"), source=ALLIANCE, release="9.0.0", default=True),
            replace(lookup_xref("Homo sapiens"), source="ensembl", release="116", default=False),
        )
        answered = lookup_xref("Homo sapiens", release="9.0.0", table=two_releases)
        assert (answered.source, answered.release) == (ALLIANCE, "9.0.0")


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

    def test_worm_authority_and_identity_hop_missing_hub_dedup_and_unknown_taxon(
        self, fixture_lines: list[str]
    ) -> None:
        # WB:WBGene00000001 -> WB:aap-1 is a symbol. Reading the authority off the
        # cross-reference column would key the gene by 'aap-1' and call it a WormBase id.
        worm = read_alliance(fixture_lines, ncbi_taxid=6239, origin="fixture")
        worm_ids = {i for namespace, i, _stem in worm if namespace == WORMBASE}
        assert "WBGene00000001" in worm_ids
        assert not any("aap-1" in identifier for identifier in worm_ids)

        # Counted on the whole release 9.0.0 file, not assumed: all 46,926 worm genes carry
        # an ENSEMBL cross-reference whose id is the WBGene id, with zero differing. So a
        # **Gene id stem** and a WormBase gene id are the same string, and the hop worm data
        # makes into this package's answers is the identity function.
        worm_pairs = [(identifier, stem) for ns, identifier, stem in worm if ns == WORMBASE]
        assert worm_pairs
        assert all(identifier == stem for identifier, stem in worm_pairs)
        assert {identifier for identifier, _stem in worm_pairs} == {
            identifier for ns, identifier, _stem in worm if ns == ENSEMBL
        }

        human = read_alliance(fixture_lines, ncbi_taxid=9606, origin="fixture")
        assert HUMAN_GENE_WITHOUT_A_HUB not in {i for _ns, i, _stem in human}
        assert len(human) == len(set(human))

        assert read_alliance(fixture_lines, ncbi_taxid=7955, origin="fixture") == ()

    def test_a_bad_file_is_refused_and_names_what_it_expected(self) -> None:
        with pytest.raises(AllianceFileError, match="GeneID"):
            read_alliance(["# only a comment"], ncbi_taxid=9606, origin="somewhere.tsv")

        respelled = ["\t".join(("Gene", *ALLIANCE_COLUMNS[1:])), "a\tb\tc\td\te"]
        with pytest.raises(AllianceFileError, match=r"somewhere\.tsv"):
            read_alliance(respelled, ncbi_taxid=9606, origin="somewhere.tsv")

        unknown_authority = [
            "\t".join(ALLIANCE_COLUMNS),
            "ZFIN:ZDB-1\tENSEMBL:ENSDARG1\tu\tgeneric\tNCBITaxon:9606",
        ]
        with pytest.raises(AllianceFileError, match="HGNC"):
            read_alliance(unknown_authority, ncbi_taxid=9606, origin="somewhere.tsv")


# ---------------------------------------------------------------------------
# Preparing the set: construction fetches once and re-reads thereafter
# ---------------------------------------------------------------------------


class TestXrefSetPreparation:
    def test_a_prepared_sets_identity_layout_and_construction_options(
        self, served: FakeFetch, liulab_data: Path, tmp_path: Path
    ) -> None:
        human = XrefSet("Homo sapiens")
        assert served.last.url == lookup_xref("Homo sapiens").url

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
        assert xref_data_dir() == liulab_data / "xref"
        assert not (liulab_data / "genome").exists()

        assert (human.species, human.source, human.release) == ("Homo sapiens", ALLIANCE, RELEASE)
        assert human.source_url == lookup_xref("Homo sapiens").url
        assert repr(human) == (
            f"XrefSet(species='Homo sapiens', source='alliance', release='{RELEASE}', stems=6)"
        )

        mouse = XrefSet("Mus musculus")
        assert len(served.calls) == 2
        assert human.path != mouse.path

        again = XrefSet("Homo sapiens")
        assert len(served.calls) == 2
        assert human.path == again.path
        assert human.to_stems(["7157"], ENTREZ) == again.to_stems(["7157"], ENTREZ)

        XrefSet("Homo sapiens", cache_dir=tmp_path / "no-progressbar", progressbar=False)
        assert served.last.progressbar is False

    def test_an_explicit_cache_dir_writes_nothing_into_the_shared_default_tree(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        # A fresh data root: the mega test above already builds a default-location set,
        # so running this check there would find the default tree already populated and
        # prove nothing. Here nothing has touched it yet.
        elsewhere = tmp_path / "somewhere-else"
        cached = XrefSet("Homo sapiens", cache_dir=elsewhere)
        assert cached.path.parent == elsewhere
        assert not (xref_data_dir() / ALLIANCE).exists()

    @pytest.mark.parametrize("species", sorted(FIXTURE_SLICES))
    def test_every_species_carries_the_namespaces_the_source_has_for_it(
        self, served: FakeFetch, species: str
    ) -> None:
        namespaces, _triples, stems = FIXTURE_SLICES[species]
        prepared = XrefSet(species)
        assert prepared.namespaces == namespaces
        assert len(prepared) == stems

    def test_preparation_refuses_before_and_after_the_fetch(
        self, fake_fetch: FakeFetch, serve_source: ServeSource
    ) -> None:
        with pytest.raises(NoXrefSetError, match="Homo sapiens"):
            XrefSet("Danio rerio")
        assert fake_fetch.calls == []

        # No `pinned_to_fixture` here: the shipped row pins the whole 25 MB publisher file,
        # and 14 genes of it is not that file. A truncated download is not a smaller
        # release, and slicing one would answer with silently fewer genes.
        fake_fetch.serve(FIXTURE)
        with pytest.raises(PreparedChecksumError, match="hashes to") as pinned:
            XrefSet("Homo sapiens")
        # The pin covers the publisher's unpacked bytes, and the repair is both halves of
        # preparing the set again — the shared pipeline's message, not one of this module's.
        assert "unpacked bytes" in str(pinned.value)
        assert xref_prepare_command("Homo sapiens", ALLIANCE, RELEASE) in str(pinned.value)

        serve_source(["MGI:98834\tENSEMBL:ENSMUSG00000059552\tu\tgeneric\tNCBITaxon:10090"])
        with pytest.raises(XrefTableError, match="NCBITaxon:9606"):
            XrefSet("Homo sapiens")


class TestVersionsOnTheSourceSide:
    """A publisher that spells a version is joined to one that does not, and neither loses.

    **Alliance publishes none**: 0 of the 43,867 human ids in release 9.0.0 carry a dot, so
    the committed fixture cannot show this and this test writes the rows itself. It is here
    because the trap is real elsewhere — NCBI's ``gene2ensembl`` writes the gene id bare and
    the transcript id versioned in the same row — and because the reader is where a second
    source (#150) will meet it.
    """

    def test_the_reader_stems_on_ingest_and_a_versioned_ask_still_meets_a_bare_source(
        self, serve_source: ServeSource
    ) -> None:
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

    def test_an_interrupted_download_reads_as_unfinished_and_stale_work_is_never_adopted(
        self, monkeypatch: pytest.MonkeyPatch, pinned_to_fixture: None, served: FakeFetch
    ) -> None:
        directory = xref_set_dir("Homo sapiens", ALLIANCE, RELEASE)
        directory.mkdir(parents=True)
        (directory / xref_slice_name("Homo sapiens")).write_bytes(b"half a file")
        with pytest.raises(UnfinishedRegistrationError, match=RECORD_NAME):
            XrefSet("Homo sapiens")
        (directory / xref_slice_name("Homo sapiens")).unlink()

        work = directory / ".work"
        work.mkdir(parents=True)
        (work / "GENECROSSREFERENCE_COMBINED_11.tsv.gz").write_bytes(b"not a gzip")
        human = XrefSet("Homo sapiens")
        assert len(human) == FIXTURE_SLICES["Homo sapiens"][2]
        assert not work.exists()


class TestCompletionMarker:
    def test_the_marker_records_provenance_and_is_verified_on_every_read(
        self, served: FakeFetch, fixture_path: Path, tmp_path: Path
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

        marker = human.path.parent / RECORD_NAME
        payload = json.loads(marker.read_text(encoding="utf-8"))
        payload["sha256"] = "0" * 64
        marker.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(RegistrationMismatchError, match="hashes to"):
            XrefSet("Homo sapiens")

        size_broken = XrefSet("Homo sapiens", cache_dir=tmp_path / "size-broken")
        size_broken.path.write_bytes(b"")
        with pytest.raises(RegistrationMismatchError, match=RECORD_NAME):
            XrefSet("Homo sapiens", cache_dir=tmp_path / "size-broken")


class TestStoredForm:
    def test_the_stored_slice_is_a_plain_gzipped_tsv_a_collaborator_can_read(
        self, served: FakeFetch, tmp_path: Path
    ) -> None:
        human = XrefSet("Homo sapiens")
        with gzip.open(human.path, "rt", encoding="utf-8") as handle:
            assert handle.readline() == "\t".join(SLICE_COLUMNS) + "\n"

        # No library of this package's is involved: three tab-separated columns and a
        # header, which is all R's read.delim or a shell's cut needs.
        frame = pd.read_csv(human.path, sep="\t", dtype=str)
        assert list(frame.columns) == list(SLICE_COLUMNS)
        assert len(frame) == FIXTURE_SLICES["Homo sapiens"][1]

        second = XrefSet("Homo sapiens", cache_dir=tmp_path / "two")
        assert human.path.read_bytes() == second.path.read_bytes()

    def test_parse_slice_refuses_a_bad_file_and_names_the_line(self) -> None:
        with pytest.raises(XrefTableError, match="namespace"):
            parse_slice("gene\tid\n", origin="somewhere.tsv")

        unknown_namespace = "\t".join(SLICE_COLUMNS) + "\nrefseq\tNM_1\tENSG1\n"
        with pytest.raises(XrefTableError, match="refseq"):
            parse_slice(unknown_namespace, origin="somewhere.tsv")

        short_row = "\t".join(SLICE_COLUMNS) + "\nentrez\t7157\n"
        with pytest.raises(XrefTableError, match="line 2"):
            parse_slice(short_row, origin="somewhere.tsv")


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
        [(ENTREZ, "7157"), (ENSEMBL, "ENSG00000141510")],
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

    def test_mouse_and_worm_reach_the_hub_through_their_own_authority(
        self, mouse: XrefSet, worm: XrefSet
    ) -> None:
        answer = mouse.to_stems(["MGI:98834"], MGI)
        assert answer.resolved == {"MGI:98834": ("ENSMUSG00000059552",)}

        # The WormBase gene id *is* the Ensembl stable gene id, so the hop is the identity.
        worm_answer = worm.to_stems(["WBGene00000001"], WORMBASE)
        assert worm_answer.resolved == {"WBGene00000001": ("WBGene00000001",)}

    def test_an_id_naming_several_stems_answers_with_all_of_them(
        self, human: XrefSet, worm: XrefSet
    ) -> None:
        human_answer = human.to_stems(["HGNC:13666"], HGNC)
        assert human_answer.resolved == {"HGNC:13666": HUMAN_GENES["HGNC:13666"]}

        worm_answer = worm.to_stems(["P05634"], UNIPROT)
        assert worm_answer.resolved["P05634"] == (
            "WBGene00003425",
            "WBGene00003432",
            "WBGene00003449",
            "WBGene00003463",
        )

    def test_the_answer_never_empties_a_value_preserves_order_and_treats_equivalent_ids_alike(
        self, human: XrefSet
    ) -> None:
        # Ids that named nothing ride back in ask order, no resolved value is ever empty, a
        # repeated ask is answered once, a versioned and an unversioned spelling and the
        # authority's two spellings resolve identically, duplicate source rows don't
        # duplicate the answer, and a gene the fixture gives no hub reaches nothing.
        named_nothing = human.to_stems(["999999998", "7157", "999999999"], ENTREZ)
        assert named_nothing.unresolved == ("999999998", "999999999")
        assert list(named_nothing.resolved) == ["7157"]

        never_empty = human.to_stems(["7157", "999999999", ""], ENTREZ)
        assert all(never_empty.resolved.values())

        asked = ["HGNC:13666", "HGNC:11998", "HGNC:1100", "HGNC:7622"]
        assert list(human.to_stems(asked, HGNC).resolved) == asked

        versioned = human.to_stems(["ENSG00000141510.18", "ENSG00000141510"], ENSEMBL)
        assert versioned.resolved == {
            "ENSG00000141510.18": ("ENSG00000141510",),
            "ENSG00000141510": ("ENSG00000141510",),
        }
        assert human.to_stems(["11998"], HGNC).gene_id_stems == (
            human.to_stems(["HGNC:11998"], HGNC).gene_id_stems
        )

        # The publisher lists HGNC:1100 -> NCBI_Gene:672 twice, under two pages.
        assert human.to_stems(["672"], ENTREZ).resolved == {"672": ("ENSG00000012048",)}

        twice = human.to_stems(["7157", "7157"], ENTREZ)
        assert list(twice.resolved) == ["7157"]

        assert human.to_stems([HUMAN_GENE_WITHOUT_A_HUB], HGNC).unresolved == (
            HUMAN_GENE_WITHOUT_A_HUB,
        )


class TestFromStems:
    @pytest.mark.parametrize(
        ("namespace", "expected"),
        [(ENTREZ, ("7157",)), (ENSEMBL, ("ENSG00000141510",))],
    )
    def test_every_namespace_answers_from_the_hub(
        self, human: XrefSet, namespace: str, expected: tuple[str, ...]
    ) -> None:
        answer = human.from_stems(["ENSG00000141510"], namespace)
        assert answer.resolved == {"ENSG00000141510": expected}
        assert answer.namespace == namespace

    def test_the_answer_never_empties_a_value_preserves_order_and_meets_a_versioned_ask(
        self, human: XrefSet, mouse: XrefSet, worm: XrefSet
    ) -> None:
        multi = mouse.from_stems(["ENSMUSG00000059552"], UNIPROT)
        assert len(multi.resolved["ENSMUSG00000059552"]) == 14
        assert "P02340" in multi.resolved["ENSMUSG00000059552"]

        versioned = human.from_stems(["ENSG00000141510.18"], HGNC)
        assert versioned.resolved == {"ENSG00000141510.18": ("HGNC:11998",)}

        named_nothing = human.from_stems(
            ["ENSG00000000001", "ENSG00000141510", "ENSG00000000002"], HGNC
        )
        assert named_nothing.unresolved == ("ENSG00000000001", "ENSG00000000002")

        never_empty = human.from_stems(["ENSG00000141510", "ENSG00000000001"], UNIPROT)
        assert all(never_empty.resolved.values())

        asked = ["WBGene00003449", "WBGene00000001", "WBGene00003425"]
        assert list(worm.from_stems(asked, ENTREZ).resolved) == asked

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
    def test_an_unnamed_or_unknown_namespace_refuses_and_names_what_there_is(
        self, mouse: XrefSet, worm: XrefSet, human: XrefSet
    ) -> None:
        with pytest.raises(NamespaceNotCarriedError) as raised:
            mouse.to_stems(["HGNC:11998"], HGNC)
        message = str(raised.value)
        for namespace in mouse.namespaces:
            assert namespace in message

        with pytest.raises(NamespaceNotCarriedError, match="wormbase"):
            worm.from_stems(["WBGene00000001"], MGI)

        with pytest.raises(NamespaceNotCarriedError, match="refseq"):
            human.to_stems(["NM_000546"], "refseq")

    def test_case_is_not_significant(self, human: XrefSet) -> None:
        assert human.to_stems(["7157"], "Entrez").namespace == ENTREZ


class TestAnswerShape:
    def test_both_verbs_flatten_and_serialize_their_answer(
        self, human: XrefSet, worm: XrefSet
    ) -> None:
        to_stems_answer = human.to_stems(["HGNC:13666", "HGNC:11998"], HGNC)
        assert to_stems_answer.gene_id_stems == [
            "ENSG00000094914",
            "ENSG00000291836",
            "ENSG00000141510",
        ]

        from_stems_answer = worm.from_stems(["WBGene00000001", "WBGene00000912"], ENTREZ)
        assert from_stems_answer.xref_ids == ["172141", "172981"]

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

        from_stems_rendered = human.from_stems(["ENSG00000141510"], HGNC).as_json()
        assert from_stems_rendered["namespace"] == HGNC
        assert from_stems_rendered["xref_ids"] == ["HGNC:11998"]
        assert json.loads(json.dumps(from_stems_rendered)) == from_stems_rendered

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
    def test_resolved_and_unresolved_partition_the_ask_exactly_once_in_order(
        self, human: XrefSet, asked: list[str]
    ) -> None:
        answer = human.to_stems(asked, HGNC)
        once = list(dict.fromkeys(asked))
        assert sorted([*answer.resolved, *answer.unresolved]) == sorted(once)
        assert len({*answer.resolved} & {*answer.unresolved}) == 0
        assert len(answer.unresolved) == len(set(answer.unresolved))
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
