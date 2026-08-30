"""Tests for Ensembl as a second **Xref source**, and for the empty-filter trap.

Offline throughout, on committed real bytes. Three fixtures are served here, routed by the
URL the package built — ``ensembl_entrez_human_tiny.tsv.gz`` and
``ensembl_entrez_mouse_tiny.tsv.gz`` cut from Ensembl release 116's per-species dumps, and
``alliance_genecrossreference_disagreeing.tsv.gz`` cut from Alliance 9.0.0 — so that one
test can construct two **Xref set**s over one species and ask both about one id.

Everything is driven through the two seams the design names: the :class:`XrefSet`
constructor and its two verbs. What is asserted is what a caller can observe — the answers,
the source each names, the error and the next action it points at — and never which private
helper ran.

The one deliberate exception is :class:`TestTheFanOutIsStatedWhereTheSourceIsChosen`, which
reads docstrings and shipped prose. Ensembl's fan-out is a fact about the data that a
caller must meet *before* choosing the source rather than after being surprised by it, so
"it is written where the choice is made" is an acceptance criterion and gets a test like
any other.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import gzip
import hashlib
from collections.abc import Callable
from dataclasses import replace
from importlib.resources import files
from pathlib import Path
from typing import Any

import pytest

import genome.xref.ensembl as ensembl_mod
from genome.store import fetch as fetch_mod
from genome.store.completion import read_record
from genome.store.prepared import PreparedChecksumError
from genome.xref import (
    ALLIANCE,
    ENSEMBL,
    ENSEMBL_TSV,
    ENTREZ,
    HGNC,
    UNIPROT,
    EmptyEvidenceFilterError,
    EnsemblTsvFileError,
    EvidenceNotRecordedError,
    NamespaceNotCarriedError,
    NoXrefSetError,
    XrefSet,
    XrefSetNotDownloadedError,
    XrefTableError,
    lookup_xref,
    xref_prepare_command,
    xref_releases,
    xref_set_dir,
    xref_slice_name,
    xref_sources,
)
from genome.xref import metadata as metadata_mod
from genome.xref.ensembl import (
    DEPENDENT,
    DIRECT,
    ENSEMBL_COLUMNS,
    ENTREZ_DB_NAME,
    TRANSCRIPT_NAME_DB_NAME,
    read_ensembl,
)

from ..conftest import FakeFetch

#: The committed fixtures. See tests/data/README.md for what each is cut for.
ENSEMBL_HUMAN = "xref/ensembl_entrez_human_tiny.tsv.gz"
ENSEMBL_MOUSE = "xref/ensembl_entrez_mouse_tiny.tsv.gz"
ALLIANCE_DISAGREEING = "xref/alliance_genecrossreference_disagreeing.tsv.gz"

#: The Ensembl **Release** the shipped rows pin. A number, and Alliance's ``9.0.0`` is not.
RELEASE = "116"

#: The Alliance release the other rows pin, so the two are visibly independent.
ALLIANCE_RELEASE = "9.0.0"

#: The id both sources carry and disagree about: Entrez GeneID 79166. Alliance 9.0.0 lists
#: its gene, ``HGNC:15497``, with **two** ``ENSEMBL:`` cross-references; Ensembl release
#: 116 asserts **72** stems for the same GeneID. Both counted on the publishers' own files.
DISAGREEING_GENEID = "79166"

#: What Alliance 9.0.0 says GeneID 79166 names — the whole of it.
ALLIANCE_STEMS: tuple[str, ...] = ("ENSG00000170858", "ENSG00000293273")

#: How many stems Ensembl release 116 says the same GeneID names.
ENSEMBL_FAN_OUT = 72

#: One human stem the fixture carries that four different GeneIDs name — the fan-out
#: running the other way, which is what makes ``from_stems`` many-valued here.
FAN_IN_STEM = "ENSG00000173213"
FAN_IN_GENEIDS: tuple[str, ...] = ("124908013", "124908015", "124908041", "260334")

#: What each Ensembl fixture comes to once sliced: triples, then stems. Counted off the
#: committed bytes, so a fixture edited without this table fails loudly.
HUMAN_TRIPLES, HUMAN_STEMS = 159, 78
MOUSE_TRIPLES, MOUSE_STEMS = 13, 6

#: The distinct (GeneID, stem) pairs the human fixture's 149 rows collapse to — the same
#: collapse release 116 makes at full size, where 552,633 rows become 36,824 pairs.
HUMAN_PAIRS = 81

#: Every row of the human fixture is ``DEPENDENT`` and not one is ``DIRECT`` — 149 of them,
#: which is release 116's whole shape in miniature: 552,633 human ``EntrezGene`` rows,
#: **zero** direct. Mouse is the same at 358,853.
HUMAN_DEPENDENT_ROWS = 149


def _md5_unpacked(path: Path) -> str:
    """The md5 of a publisher file's **unpacked** bytes, which is what a row pins."""
    with gzip.open(path, "rb") as handle:
        return hashlib.md5(handle.read()).hexdigest()


@pytest.fixture
def sources(fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> FakeFetch:
    """Serve every shipped row's publisher file from the fixture cut out of it.

    Routed by URL rather than switched by hand, so a single test may construct the Alliance
    set and the Ensembl set and ask both about one id. Each row is re-pinned to the digest
    of the fixture its own URL serves: the checksum check is never switched off, only
    pointed at the bytes that actually arrive.
    """
    by_url = {
        lookup_xref("Homo sapiens", ALLIANCE).url: data_dir / ALLIANCE_DISAGREEING,
        lookup_xref("Homo sapiens", ENSEMBL_TSV).url: data_dir / ENSEMBL_HUMAN,
        lookup_xref("Mus musculus", ENSEMBL_TSV).url: data_dir / ENSEMBL_MOUSE,
    }
    rows = tuple(
        replace(row, source_checksum=f"md5:{_md5_unpacked(by_url[row.url])}")
        if row.url in by_url
        else row
        for row in metadata_mod.xref_table()
    )
    monkeypatch.setattr(metadata_mod, "xref_table", lambda: rows)

    def route(url: str, dest_dir: Path, **kwargs: Any) -> Path:
        fake_fetch.serve(by_url[url])
        return fake_fetch(url, dest_dir, **kwargs)

    monkeypatch.setattr(fetch_mod, "fetch_url", route)
    return fake_fetch


@pytest.fixture
def human(sources: FakeFetch) -> XrefSet:
    """The human Ensembl set, prepared from the committed release-116 fixture."""
    return XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)


@pytest.fixture
def mouse(sources: FakeFetch) -> XrefSet:
    """The mouse Ensembl set, prepared from the committed release-116 fixture."""
    return XrefSet("Mus musculus", ENSEMBL_TSV, RELEASE)


@pytest.fixture
def human_lines(data_dir: Path) -> list[str]:
    """The human fixture, unpacked into lines."""
    with gzip.open(data_dir / ENSEMBL_HUMAN, "rt", encoding="utf-8") as handle:
        return handle.read().splitlines()


@pytest.fixture
def attribution_text() -> str:
    """The attribution that ships beside the curated table, as installed."""
    resource = files("genome").joinpath("data/xref/ATTRIBUTION.md")
    return resource.read_text(encoding="utf-8")


@pytest.fixture
def serve_ensembl_rows(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Callable[[list[str]], None]:
    """Serve Ensembl-shaped rows as if the publisher's dump held only those.

    For the handful of cases the committed fixtures cannot show, because the publisher's
    own file does not contain them. Everything else is driven off the real bytes.
    """

    def serve(rows: list[str]) -> None:
        payload = "".join(f"{row}\n" for row in ["\t".join(ENSEMBL_COLUMNS), *rows]).encode()
        path = tmp_path / "served.tsv.gz"
        with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as out:
            out.write(payload)
        digest = f"md5:{hashlib.md5(payload).hexdigest()}"
        pinned = tuple(replace(row, source_checksum=digest) for row in metadata_mod.xref_table())
        monkeypatch.setattr(metadata_mod, "xref_table", lambda: pinned)
        fake_fetch.serve(path)

    return serve


# ---------------------------------------------------------------------------
# The committed bytes
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    def test_the_fixture_carries_the_traps_the_reader_exists_to_catch(
        self, human_lines: list[str]
    ) -> None:
        assert tuple(human_lines[0].split("\t")) == ENSEMBL_COLUMNS

        # The trap, on real bytes: release 116's human dump carries 552,633 EntrezGene rows,
        # every one DEPENDENT and not one DIRECT, and the fixture is cut from it unedited.
        rows = [line.split("\t") for line in human_lines[1:]]
        entrez = [row for row in rows if row[4] == ENTREZ_DB_NAME]
        assert len(entrez) == HUMAN_DEPENDENT_ROWS
        assert {row[5] for row in entrez} == {DEPENDENT}
        assert DIRECT not in {row[5] for row in entrez}

        # `EntrezGene_trans_name` rows put a *transcript name* in the xref column —
        # `KU-MEL-3-201`, not a GeneID — so a reader that did not split on db_name would
        # key the Entrez namespace by a transcript label. They also carry MISC rather than
        # DEPENDENT, which is why the empty-filter message must not offer MISC as a type
        # the caller could have asked for.
        named = [row for row in rows if row[4] == TRANSCRIPT_NAME_DB_NAME]
        assert named
        assert {row[5] for row in named} == {"MISC"}
        assert not any(row[3].isdigit() for row in named)

        entrez_rows = [row for row in rows if row[4] == "EntrezGene"]
        stems = {row[0].split(".")[0] for row in entrez_rows if row[3] == DISAGREEING_GENEID}
        assert len(stems) == ENSEMBL_FAN_OUT
        assert set(ALLIANCE_STEMS) <= stems
        assert {row[3] for row in entrez_rows if row[0].split(".")[0] == FAN_IN_STEM} == set(
            FAN_IN_GENEIDS
        )


# ---------------------------------------------------------------------------
# Ensembl is selectable, and pins to a numbered release of its own
# ---------------------------------------------------------------------------


class TestEnsemblIsSelectable:
    def test_it_is_a_second_source_for_human_and_mouse_pinning_its_own_release_worm_has_none(
        self,
    ) -> None:
        # Which sources a species has grows as rows are added, so what is asserted is that
        # this one is among them and that the first one still leads — not the whole list.
        for species in ("Homo sapiens", "Mus musculus"):
            assert xref_sources(species)[0] == ALLIANCE
            assert ENSEMBL_TSV in xref_sources(species)
            assert lookup_xref(species).source == ALLIANCE
            assert lookup_xref(species, ENSEMBL_TSV).default is False

            row = lookup_xref(species, ENSEMBL_TSV)
            assert row.publisher == "Ensembl"
            assert row.version == RELEASE
            assert row.pubmed_id == 39656687
            algorithm, _, digest = row.source_checksum.partition(":")
            assert algorithm == "md5"
            assert len(digest) == 32
            assert f"/release-{RELEASE}/" in row.url
            assert row.url.endswith(f".{RELEASE}.entrez.tsv.gz")
            # `current_tsv` is a live 404: there is no current-release shortcut for these
            # dumps, so the number is mandatory and travels in the URL.
            assert "current_tsv" not in row.url

        assert xref_releases("Homo sapiens", ENSEMBL_TSV) == (RELEASE,)
        assert xref_releases("Homo sapiens", ALLIANCE) == (ALLIANCE_RELEASE,)

        # Ensembl files worm under Ensembl Genomes' own numbering — release-116's worm
        # directory holds `Caenorhabditis_elegans.WBcel235.63.entrez.tsv.gz` — so "116"
        # would name a file that does not exist. Worm is answered by Alliance, where the
        # hop is the identity, and no worm row is invented here.
        assert ENSEMBL_TSV not in xref_sources("Caenorhabditis elegans")
        with pytest.raises(NoXrefSetError, match=ALLIANCE):
            lookup_xref("Caenorhabditis elegans", ENSEMBL_TSV)

    def test_the_two_sources_are_prepared_and_cached_in_two_directories(
        self, sources: FakeFetch
    ) -> None:
        alliance = XrefSet("Homo sapiens", ALLIANCE)
        ensembl = XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)
        assert alliance.path != ensembl.path
        assert alliance.path.parent == xref_set_dir("Homo sapiens", ALLIANCE, ALLIANCE_RELEASE)
        assert ensembl.path.parent == xref_set_dir("Homo sapiens", ENSEMBL_TSV, RELEASE)
        assert (ensembl.species, ensembl.source, ensembl.release) == (
            "Homo sapiens",
            ENSEMBL_TSV,
            RELEASE,
        )
        assert ensembl.source_url == lookup_xref("Homo sapiens", ENSEMBL_TSV).url

        again = XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)
        assert len(sources.calls) == 2
        assert again.path == ensembl.path


class TestNamingAReleaseWithoutASource:
    """Two sources numbering their releases differently, on the rows that actually ship.

    What `lookup_xref` does with a release and no source is pinned on synthetic rows in
    ``tests/xref/test_xref.py``; what is here is that the *shipped* table now makes the case
    reachable, which it could not while every species listed one release.
    """

    def test_a_release_named_alone_finds_the_source_that_has_it_and_no_other(self) -> None:
        answered = lookup_xref("Homo sapiens", release=ALLIANCE_RELEASE)
        assert (answered.source, answered.release) == (ALLIANCE, ALLIANCE_RELEASE)

        # `116` is Ensembl's release string, not Alliance's. Asking the default source for
        # it must miss and say so rather than quietly answering with 9.0.0.
        with pytest.raises(NoXrefSetError, match=ALLIANCE_RELEASE):
            lookup_xref("Homo sapiens", release=RELEASE)


# ---------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------


class TestEnsemblReader:
    def test_it_reads_the_fixture_gene_level_and_stems_a_versioned_id_on_ingest(
        self, human_lines: list[str]
    ) -> None:
        read = read_ensembl(human_lines, ncbi_taxid=9606, origin="fixture")
        assert len(read) == HUMAN_TRIPLES
        assert {namespace for namespace, _id, _stem in read} == {ENSEMBL, ENTREZ}
        assert len({stem for _ns, _id, stem in read}) == HUMAN_STEMS

        entrez_ids = {identifier for ns, identifier, _stem in read if ns == ENTREZ}
        assert "KU-MEL-3-201" not in entrez_ids
        assert all(identifier.isdigit() for identifier in entrez_ids)

        hub = [(identifier, stem) for ns, identifier, stem in read if ns == ENSEMBL]
        assert len(hub) == HUMAN_STEMS
        assert all(identifier == stem for identifier, stem in hub)

        # The publisher's file is transcript- and protein-grained: 552,633 human rows
        # collapse to 36,824 gene-level pairs. Deduplication is on the pair and happens
        # here, so an answer counts genes rather than transcripts.
        assert len(read) == len(set(read))
        pairs = [triple for triple in read if triple[0] == ENTREZ]
        assert len(pairs) == HUMAN_PAIRS < HUMAN_DEPENDENT_ROWS

        versioned = [
            "\t".join(ENSEMBL_COLUMNS),
            "ENSG00000141510.18\tENST1\t-\t7157\tEntrezGene\tDEPENDENT\t-\t-\t-",
        ]
        assert read_ensembl(versioned, ncbi_taxid=9606, origin="x") == (
            (ENSEMBL, "ENSG00000141510", "ENSG00000141510"),
            (ENTREZ, "7157", "ENSG00000141510"),
        )

    def test_a_broken_file_is_refused_and_names_what_it_expected(self) -> None:
        with pytest.raises(EnsemblTsvFileError, match="gene_stable_id"):
            read_ensembl([], ncbi_taxid=9606, origin="somewhere.tsv")

        respelled = ["\t".join(("gene_id", *ENSEMBL_COLUMNS[1:])), "a\tb\tc\td\te\tf\tg\th\ti"]
        with pytest.raises(EnsemblTsvFileError, match=r"somewhere\.tsv"):
            read_ensembl(respelled, ncbi_taxid=9606, origin="somewhere.tsv")

        too_narrow = ["\t".join(ENSEMBL_COLUMNS), "too\tfew\tfields"]
        with pytest.raises(EnsemblTsvFileError, match="fields"):
            read_ensembl(too_narrow, ncbi_taxid=9606, origin="somewhere.tsv")


# ---------------------------------------------------------------------------
# Two sources, two answers, nothing merged
# ---------------------------------------------------------------------------


class TestTwoSourcesDisagree:
    """One id, two **Xref set**s, two answers — each naming its own source (ADR-0017)."""

    def test_one_geneid_gets_two_unmerged_answers_each_naming_its_own_source(
        self, sources: FakeFetch
    ) -> None:
        alliance = XrefSet("Homo sapiens", ALLIANCE).to_stems([DISAGREEING_GENEID], ENTREZ)
        ensembl = XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE).to_stems(
            [DISAGREEING_GENEID], ENTREZ
        )
        assert alliance.resolved[DISAGREEING_GENEID] == ALLIANCE_STEMS
        assert len(ensembl.resolved[DISAGREEING_GENEID]) == ENSEMBL_FAN_OUT
        assert alliance.resolved != ensembl.resolved
        assert (alliance.source, alliance.release) == (ALLIANCE, ALLIANCE_RELEASE)
        assert (ensembl.source, ensembl.release) == (ENSEMBL_TSV, RELEASE)
        assert alliance.as_json()["source"] == ALLIANCE
        assert ensembl.as_json()["source"] == ENSEMBL_TSV
        assert alliance.as_json()["release"] != ensembl.as_json()["release"]

        # The narrower answer is a strict subset of the wider one and stays that way: no
        # union, no vote, no "best" mapping. A caller who wants both makes both calls.
        assert set(ALLIANCE_STEMS) < set(ensembl.resolved[DISAGREEING_GENEID])

    def test_neither_set_answers_in_a_namespace_only_the_other_carries(
        self, sources: FakeFetch, human: XrefSet
    ) -> None:
        # Alliance reaches human through HGNC; Ensembl's entrez dump carries no authority
        # id at all. A namespace one source has is not one the other answers in.
        alliance = XrefSet("Homo sapiens", ALLIANCE)
        assert HGNC in alliance.namespaces
        assert human.namespaces == (ENSEMBL, ENTREZ)
        with pytest.raises(NamespaceNotCarriedError, match=ENTREZ):
            human.to_stems(["HGNC:15497"], HGNC)

        with pytest.raises(NamespaceNotCarriedError) as raised:
            human.from_stems(["ENSG00000170858"], UNIPROT)
        assert ENSEMBL in str(raised.value)
        assert ENTREZ in str(raised.value)


class TestTheTwoVerbs:
    def test_the_fan_out_and_fan_in_answer_with_everything_named_and_a_mouse_id_is_its_own(
        self, human: XrefSet, mouse: XrefSet
    ) -> None:
        answer = human.to_stems([DISAGREEING_GENEID], ENTREZ)
        assert len(answer.resolved[DISAGREEING_GENEID]) == ENSEMBL_FAN_OUT
        assert answer.unresolved == ()

        assert human.from_stems([FAN_IN_STEM], ENTREZ).resolved[FAN_IN_STEM] == FAN_IN_GENEIDS

        assert human.from_stems(["ENSG00000170858.14"], ENTREZ).resolved == {
            "ENSG00000170858.14": (DISAGREEING_GENEID,)
        }

        rides_back = human.to_stems(["999999998", DISAGREEING_GENEID, "999999999"], ENTREZ)
        assert rides_back.unresolved == ("999999998", "999999999")
        assert list(rides_back.resolved) == [DISAGREEING_GENEID]

        assert len(mouse) == MOUSE_STEMS
        assert mouse.namespaces == (ENSEMBL, ENTREZ)
        assert mouse.to_stems(["22059"], ENTREZ).resolved == {"22059": ("ENSMUSG00000059552",)}
        # A human GeneID is not a mouse one, and the mouse set says so rather than guessing.
        assert mouse.to_stems([DISAGREEING_GENEID], ENTREZ).unresolved == (DISAGREEING_GENEID,)

        record = read_record(human.path.parent)
        assert record is not None
        assert record.details["rows"] == HUMAN_TRIPLES
        assert record.details["source"] == ENSEMBL_TSV
        assert record.details["release"] == RELEASE


# ---------------------------------------------------------------------------
# The empty-filter trap
# ---------------------------------------------------------------------------


class TestTheEmptyFilterTrap:
    """The intuitive quality filter empties the set, so it raises and says why.

    Every human ``EntrezGene`` row Ensembl release 116 publishes carries
    ``info_type=DEPENDENT`` and not one carries ``DIRECT`` — 552,633 rows, zero direct;
    mouse is the same at 358,853. So ``evidence="DIRECT"`` is not a stricter set, it is no
    set at all, and answering nothing would look like a gene list with no matches.
    """

    def test_the_trap_fires_for_human_and_mouse_and_never_offers_a_type_it_dropped(
        self, sources: FakeFetch
    ) -> None:
        with pytest.raises(EmptyEvidenceFilterError) as raised:
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=DIRECT)
        message = str(raised.value)
        assert DIRECT in message
        assert DEPENDENT in message
        assert str(HUMAN_DEPENDENT_ROWS) in message
        # `MISC` belongs to the transcript-name rows, which are a different assertion and
        # are never kept. Offering it would send the caller to a filter that also empties.
        assert "MISC" not in message
        assert "without an evidence filter" in message

        with pytest.raises(EmptyEvidenceFilterError, match=DEPENDENT):
            XrefSet("Mus musculus", ENSEMBL_TSV, RELEASE, evidence=DIRECT)

        with pytest.raises(EmptyEvidenceFilterError, match=DIRECT):
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=[DIRECT, "MISC"])

    def test_an_emptied_filter_leaves_nothing_behind_and_the_set_is_still_preparable(
        self, sources: FakeFetch
    ) -> None:
        directory = xref_set_dir("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=(DIRECT,))
        for _ in range(2):
            with pytest.raises(EmptyEvidenceFilterError, match=DEPENDENT):
                XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=DIRECT)
        assert read_record(directory) is None
        assert not (directory / xref_slice_name("Homo sapiens")).exists()
        assert len(XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)) == HUMAN_STEMS


class TestTheFilterIsARealCapability:
    """A guard on a filter nobody can call is not a guard, so the filter really filters."""

    def test_filtering_to_the_carried_type_answers_in_full_spelled_one_way_beside_the_plain_set(
        self, sources: FakeFetch
    ) -> None:
        plain = XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)
        assert plain.evidence == ()

        filtered = XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=DEPENDENT)
        assert len(filtered) == HUMAN_STEMS
        assert filtered.evidence == (DEPENDENT,)
        assert (
            filtered.to_stems([DISAGREEING_GENEID], ENTREZ).resolved[DISAGREEING_GENEID]
            == plain.to_stems([DISAGREEING_GENEID], ENTREZ).resolved[DISAGREEING_GENEID]
        )
        assert plain.path != filtered.path
        assert plain.path.is_file()
        assert filtered.path.is_file()

        for asked in ("dependent", " DEPENDENT ", ["dependent", "DEPENDENT"]):
            assert XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=asked).evidence == (
                DEPENDENT,
            )

    def test_a_source_recording_no_evidence_type_says_so_and_the_repair_carries_the_filter(
        self, sources: FakeFetch
    ) -> None:
        with pytest.raises(EvidenceNotRecordedError) as raised:
            XrefSet("Homo sapiens", ALLIANCE, evidence=DIRECT)
        message = str(raised.value)
        assert ALLIANCE in message
        assert ENSEMBL_TSV in message

        assert "evidence" in xref_prepare_command(
            "Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=(DEPENDENT,)
        )


# ---------------------------------------------------------------------------
# The fan-out, stated where the source is chosen
# ---------------------------------------------------------------------------


class TestTheFanOutIsStatedWhereTheSourceIsChosen:
    """Ensembl is not the equal of the first source, and must not read as one.

    Measured on human: the two agree on 57.6% of gene-level (GeneID, ENSG) pairs. The cause
    is method — NCBI's mapping is a sequence match at a published overlap threshold and is
    near-one-to-one, while Ensembl's fans out to 72 stems for one GeneID and 208 GeneIDs for
    one stem. A caller meets that where they choose the source, not in a note afterwards.
    """

    def test_the_fan_out_figures_are_written_where_a_caller_chooses_the_source(
        self, attribution_text: str
    ) -> None:
        assert XrefSet.__doc__ is not None
        assert ensembl_mod.__doc__ is not None
        for figure in ("72", "208", "57.6%"):
            assert figure in XrefSet.__doc__
            assert figure in ensembl_mod.__doc__
            assert figure in attribution_text

    def test_the_shipped_attribution_also_states_the_empty_filter_trap_and_checksum_convention(
        self, attribution_text: str
    ) -> None:
        assert "552,633" in attribution_text
        assert DIRECT in attribution_text
        # Alliance publishes an md5 of the *unpacked* TSV; Ensembl publishes a BSD `sum` of
        # the *served* `.gz`. Assuming either convention holds for the other is a live trap.
        assert "CHECKSUMS" in attribution_text
        assert "28782 5952" in attribution_text


# ---------------------------------------------------------------------------
# Who to cite for an answer
# ---------------------------------------------------------------------------


class TestProvenance:
    def test_a_set_cites_the_row_it_actually_resolved_to_and_not_what_was_asked(
        self, human: XrefSet, sources: FakeFetch
    ) -> None:
        assert human.provenance.source == ENSEMBL_TSV
        assert human.provenance.release == RELEASE
        assert human.provenance.url == human.source_url
        line = human.provenance.attribution()
        assert line.startswith(f"Ensembl {RELEASE}")
        assert "PMID 39656687" in line

        # Named no source and no release, so the row is the one the defaults resolved to —
        # which is the point of reading it off the set rather than looking it up again.
        prepared = XrefSet("Homo sapiens")
        assert prepared.provenance.source == ALLIANCE == prepared.source
        assert prepared.provenance.release == prepared.release


# ---------------------------------------------------------------------------
# Preparing, and what the guards still do for a second source
# ---------------------------------------------------------------------------


class TestPreparation:
    def test_a_file_that_does_not_match_its_pin_is_refused(self, fake_fetch: FakeFetch) -> None:
        # No re-pinning: the shipped row pins Ensembl's whole 6 MB human dump, and 151 rows
        # of it is not that file. A truncated download is not a smaller release.
        fake_fetch.serve(ENSEMBL_HUMAN)
        with pytest.raises(PreparedChecksumError, match="hashes to"):
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)

    def test_a_set_that_cannot_be_fetched_names_the_next_action_filter_and_all(
        self, monkeypatch: pytest.MonkeyPatch, sources: FakeFetch
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: Any) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        with pytest.raises(XrefSetNotDownloadedError) as raised:
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=DEPENDENT)
        message = str(raised.value)
        assert "login node" in message
        # The repair prepares the set that was asked for, filter and all.
        assert "evidence" in message

    def test_an_empty_file_names_the_url_rather_than_answering_nothing_filtered_or_not(
        self, serve_ensembl_rows: Callable[[list[str]], None]
    ) -> None:
        serve_ensembl_rows([])
        with pytest.raises(XrefTableError, match="carries no"):
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE)

        # Nothing was graded, so nothing can be said about which gradings exist. The file
        # is the problem, and the message that names the URL is the one that helps.
        serve_ensembl_rows([])
        with pytest.raises(XrefTableError, match="carries no"):
            XrefSet("Homo sapiens", ENSEMBL_TSV, RELEASE, evidence=DIRECT)
