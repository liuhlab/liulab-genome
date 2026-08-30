"""Tests for symbols — the asymmetric hop, and the HGNC and Alliance-BGI sources.

Offline throughout, on committed real bytes. Three fixtures are served here, routed by the
URL the package built: ``hgnc_complete_set_tiny.txt`` cut from HGNC's quarterly archive
file of 2026-07-07, and ``bgi_mgi_tiny.json.gz`` and ``bgi_wb_tiny.json.gz`` cut from the
Alliance's per-species gene submissions — MGI's and WormBase's WS298. Every record is the
publisher's own bytes, whitespace and all.

Everything is driven through the seam the design names: the :class:`XrefSet` constructor
and its verbs. What is asserted is what a caller can observe — the matches, the kind on
each, the source and release that answered, and what a set says it could not have matched —
never which private helper ran.

**The two directions are deliberately not mirror images**, and most of what follows is
about that: away from the hub one approved symbol per stem through ``from_stems``, toward
it every stem any spelling names through ``match_symbols``, and ``to_stems`` refusing the
symbol namespace rather than answering it on approved spellings alone.

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
from hypothesis import given
from hypothesis import strategies as st

from genome.io import fetch as fetch_mod
from genome.io.completion import read_record
from genome.io.results import ResolvedSymbols, SymbolMatch
from genome.xref import (
    ALIAS,
    ALLIANCE,
    ALLIANCE_BGI,
    APPROVED,
    ENSEMBL,
    ENTREZ,
    HGNC,
    HGNC_ARCHIVE,
    MGI,
    PREVIOUS,
    SYMBOL,
    SYMBOL_KINDS,
    UNIPROT,
    WORMBASE,
    BgiFileError,
    EvidenceNotRecordedError,
    HgncFileError,
    NamespaceNotCarriedError,
    NoXrefSetError,
    SymbolDirectionError,
    XrefSet,
    fold_symbol,
    lookup_xref,
    normalise_symbol,
    xref_sources,
)
from genome.xref import metadata as metadata_mod
from genome.xref.bgi import BGI_SYMBOL_LIMIT, read_bgi
from genome.xref.hgnc import HGNC_COLUMNS, read_hgnc
from genome.xref.symbols import KIND_NAMESPACES

from .conftest import FakeFetch

#: The committed fixtures. See tests/data/README.md for what each is cut for.
HGNC_FIXTURE = "xref/hgnc_complete_set_tiny.txt"
BGI_MOUSE = "xref/bgi_mgi_tiny.json.gz"
BGI_WORM = "xref/bgi_wb_tiny.json.gz"

#: The **Release** each new source pins. HGNC's is a date because the archive's files are
#: dated and irregular; the Alliance's is its own release number, as it already was.
HGNC_RELEASE = "2026-07-07"
BGI_RELEASE = "9.0.0"

#: What each fixture reduces to once sliced: triples, then **Gene id stem**s. Counted off
#: the committed bytes, so a fixture edited without this table fails loudly rather than
#: quietly changing what every test below means.
HGNC_TRIPLES, HGNC_STEMS = 73, 9
MOUSE_TRIPLES, MOUSE_STEMS = 68, 7
WORM_TRIPLES, WORM_STEMS = 33, 5

#: One gene of the fixture, spelled every way HGNC spells it. The point of the source.
BMAL1 = "ENSG00000133794"

#: A spelling HGNC retired, and one of the 31 EpiFactors rows the cofactor work measured —
#: `ARNTL` for `BMAL1`, named in that table's own attribution beside `C11orf30` for `EMSY`.
RETIRED_EPIFACTORS_SYMBOL = "ARNTL"

#: The other one, kept so the measured failure is shown twice rather than once.
RETIRED_EPIFACTORS_SYMBOL_TOO = "C11orf30"

#: One spelling that names two different genes in two different ways: HGNC's approved
#: symbol for `ADCY3`, and a symbol it retired from `ADCY8`. Both are real and both come
#: back, which is what makes ambiguity the return type here.
AMBIGUOUS_SYMBOL = "ADCY3"
ADCY3_STEM, ADCY8_STEM = "ENSG00000138031", "ENSG00000155897"

#: The mouse spelling of a human gene, which is what an analyst pastes by accident. Exact
#: matching must miss it: the species is fixed by the set, so this is the wrong authority's
#: spelling rather than a typo to absorb.
MOUSE_CASED_SYMBOL = "Brca1"
BRCA1_STEM = "ENSG00000012048"

#: An alias HGNC writes in lower case, so that "exact" and "insensitive" differ on a real
#: spelling rather than on one invented for the test.
LOWER_CASE_ALIAS = "p53"
TP53_STEM = "ENSG00000141510"

#: The human gene the fixture carries with no Ensembl gene id at all, so it has no hub and
#: appears nowhere in the slice — HGNC's own row for the AAVS1 integration site.
HGNC_ROW_WITHOUT_A_HUB = "AAVS1"

#: The shipped cofactor table the 801-row EpiFactors measurement was made on.
COFACTOR_TABLE = "data/tf_cofactor/homo_sapiens.cofactor_table.tsv.gz"


def _digest(path: Path) -> str:
    """The md5 of a publisher file's **unpacked** bytes, which is what a curated row pins.

    HGNC's archive file is plain text and the Alliance's submissions are gzipped, so the
    two conventions this directory already holds side by side meet here in one helper.
    """
    if path.suffix != ".gz":
        return hashlib.md5(path.read_bytes()).hexdigest()
    with gzip.open(path, "rb") as handle:
        return hashlib.md5(handle.read()).hexdigest()


@pytest.fixture
def sources(fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> FakeFetch:
    """Serve each new source's publisher file from the fixture cut out of it.

    Routed by URL rather than switched by hand, so one test may hold the human, mouse and
    worm symbol sets at once. Each row is re-pinned to the digest of the fixture its own
    URL serves: the checksum check is never switched off, only pointed at the bytes that
    actually arrive.
    """
    by_url = {
        lookup_xref("Homo sapiens", HGNC_ARCHIVE).url: data_dir / HGNC_FIXTURE,
        lookup_xref("Mus musculus", ALLIANCE_BGI).url: data_dir / BGI_MOUSE,
        lookup_xref("Caenorhabditis elegans", ALLIANCE_BGI).url: data_dir / BGI_WORM,
    }
    rows = tuple(
        replace(row, source_checksum=f"md5:{_digest(by_url[row.url])}")
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
    """The human HGNC set, prepared from the committed quarterly-archive fixture."""
    return XrefSet("Homo sapiens", HGNC_ARCHIVE)


@pytest.fixture
def mouse(sources: FakeFetch) -> XrefSet:
    """The mouse set, prepared from the committed MGI submission fixture."""
    return XrefSet("Mus musculus", ALLIANCE_BGI)


@pytest.fixture
def worm(sources: FakeFetch) -> XrefSet:
    """The worm set, prepared from the committed WormBase WS298 submission fixture."""
    return XrefSet("Caenorhabditis elegans", ALLIANCE_BGI)


@pytest.fixture
def alliance_human(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, data_dir: Path
) -> XrefSet:
    """The human **Default xref source** set, which carries ids and no symbol at all."""
    fixture = data_dir / "xref/alliance_genecrossreference_tiny.tsv.gz"
    digest = f"md5:{_digest(fixture)}"
    rows = tuple(replace(row, source_checksum=digest) for row in metadata_mod.xref_table())
    monkeypatch.setattr(metadata_mod, "xref_table", lambda: rows)
    fake_fetch.serve(fixture)
    return XrefSet("Homo sapiens", ALLIANCE)


@pytest.fixture
def hgnc_lines(data_dir: Path) -> list[str]:
    """The HGNC fixture, as lines."""
    return (data_dir / HGNC_FIXTURE).read_text(encoding="utf-8").splitlines()


@pytest.fixture
def worm_lines(data_dir: Path) -> list[str]:
    """The WormBase submission fixture, unpacked into lines."""
    with gzip.open(data_dir / BGI_WORM, "rt", encoding="utf-8") as handle:
        return handle.read().splitlines(keepends=True)


@pytest.fixture
def serve_hgnc_rows(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> Callable[[list[str], tuple[str, ...]], None]:
    """Serve HGNC-shaped rows under a header of the caller's own column order.

    The column order is the caller's because that is the thing under test: a reader that
    reads by name answers the same whichever order the header is in, and HGNC has already
    moved from 52 columns to 54.
    """

    def serve(rows: list[str], columns: tuple[str, ...]) -> None:
        payload = "".join(f"{row}\n" for row in ["\t".join(columns), *rows]).encode()
        path = tmp_path / "hgnc_complete_set_served.txt"
        path.write_bytes(payload)
        digest = f"md5:{hashlib.md5(payload).hexdigest()}"
        pinned = tuple(replace(row, source_checksum=digest) for row in metadata_mod.xref_table())
        monkeypatch.setattr(metadata_mod, "xref_table", lambda: pinned)
        fake_fetch.serve(path)

    return serve


def _shipped_cofactor_stems() -> dict[str, str]:
    """Return the shipped human **Cofactor table**'s symbol-to-stem column, as it ships.

    Read off the wheel's own bytes rather than restated here: the 801-row EpiFactors
    measurement was made on that table, so a test about it must be about that table.
    """
    with files("genome").joinpath(COFACTOR_TABLE).open("rb") as raw, gzip.open(raw, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        stem, symbol = header.index("gene_id_stem"), header.index("symbol")
        return {line.split("\t")[symbol]: line.split("\t")[stem] for line in handle}


# ---------------------------------------------------------------------------
# The committed bytes
# ---------------------------------------------------------------------------


class TestFixtureBytes:
    def test_the_hgnc_fixture_is_the_publishers_own_header(self, hgnc_lines: list[str]) -> None:
        header = hgnc_lines[0].split("\t")
        # 54 columns, which is the current schema; the reader needs seven of them and finds
        # every one by name.
        assert len(header) == 54
        assert set(HGNC_COLUMNS) <= set(header)

    def test_a_multi_valued_hgnc_cell_is_quoted_and_pipe_separated(
        self, hgnc_lines: list[str]
    ) -> None:
        # The trap, on real bytes: `"MOP3|JAP3|PASD3|bHLHe5|ARNTL1"` with the quotes in the
        # file, and `AC3` without them. Splitting before stripping keys the namespace by
        # `"MOP3` and `ARNTL1"`, neither of which anybody types.
        alias = hgnc_lines[0].split("\t").index("alias_symbol")
        cells = [line.split("\t")[alias] for line in hgnc_lines[1:]]
        assert any(cell.startswith('"') and "|" in cell for cell in cells)
        assert any(cell and not cell.startswith('"') for cell in cells)

    def test_the_worm_submission_is_wormbases_final_release(self, worm_lines: list[str]) -> None:
        # WS298 is WormBase's last release, and its own download host answers 403 to an
        # automated client — so this copy, served by the Alliance, is the worm authority's
        # own bytes reachable at a pinned URL.
        assert any('"release" : "WS298"' in line for line in worm_lines)

    def test_a_row_with_no_ensembl_gene_id_contributes_nothing(self, hgnc_lines: list[str]) -> None:
        symbol = hgnc_lines[0].split("\t").index("symbol")
        assert any(line.split("\t")[symbol] == HGNC_ROW_WITHOUT_A_HUB for line in hgnc_lines[1:])
        read = read_hgnc(hgnc_lines, ncbi_taxid=9606, origin="fixture")
        assert not [triple for triple in read if HGNC_ROW_WITHOUT_A_HUB in triple]


# ---------------------------------------------------------------------------
# The readers
# ---------------------------------------------------------------------------


class TestHgncReader:
    def test_it_reads_the_fixture_into_the_slice_the_table_records(
        self, hgnc_lines: list[str]
    ) -> None:
        read = read_hgnc(hgnc_lines, ncbi_taxid=9606, origin="fixture")
        assert len(read) == HGNC_TRIPLES
        assert len({stem for _namespace, _id, stem in read}) == HGNC_STEMS

    def test_reordering_a_column_does_not_change_the_answer(self, hgnc_lines: list[str]) -> None:
        # The named acceptance test. HGNC's header is 52 columns wide in the 2020-07-01
        # snapshot and 54 in this one, so a reader indexing by position reads the wrong
        # column on the older file and says nothing about it.
        header = hgnc_lines[0].split("\t")
        order = list(reversed(range(len(header))))
        shuffled = [
            "\t".join(_padded(line.split("\t"), len(header))[index] for index in order)
            for line in hgnc_lines
        ]
        assert shuffled[0] != hgnc_lines[0]
        assert read_hgnc(shuffled, ncbi_taxid=9606, origin="reordered") == read_hgnc(
            hgnc_lines, ncbi_taxid=9606, origin="fixture"
        )

    def test_a_header_naming_none_of_the_columns_it_needs_raises(self) -> None:
        with pytest.raises(HgncFileError, match="prev_symbol"):
            read_hgnc(["not\ta\theader"], ncbi_taxid=9606, origin="x")

    def test_an_empty_file_raises_rather_than_reading_nothing(self) -> None:
        with pytest.raises(HgncFileError, match="no header line"):
            read_hgnc([], ncbi_taxid=9606, origin="x")

    def test_it_refuses_an_evidence_filter_rather_than_ignoring_one(
        self, hgnc_lines: list[str]
    ) -> None:
        with pytest.raises(EvidenceNotRecordedError, match="grades"):
            read_hgnc(hgnc_lines, ncbi_taxid=9606, origin="x", evidence=("DIRECT",))


class TestBgiReader:
    def test_it_reads_each_submission_into_the_slice_the_table_records(
        self, data_dir: Path
    ) -> None:
        for name, taxid, expected in (
            (BGI_MOUSE, 10090, (MOUSE_TRIPLES, MOUSE_STEMS)),
            (BGI_WORM, 6239, (WORM_TRIPLES, WORM_STEMS)),
        ):
            with gzip.open(data_dir / name, "rt", encoding="utf-8") as handle:
                read = read_bgi(handle, ncbi_taxid=taxid, origin=name)
            assert (len(read), len({stem for _n, _i, stem in read})) == expected

    def test_a_taxon_the_submission_does_not_carry_reads_as_the_wrong_file(
        self, worm_lines: list[str]
    ) -> None:
        assert read_bgi(worm_lines, ncbi_taxid=9606, origin="fixture") == ()

    def test_a_file_with_no_data_array_raises(self) -> None:
        with pytest.raises(BgiFileError, match="data"):
            read_bgi(['{"metaData": {}}'], ncbi_taxid=10090, origin="x")

    def test_it_reads_the_symbol_and_never_the_untyped_synonyms(
        self, worm_lines: list[str]
    ) -> None:
        # `daf-16`'s record lists `daf-17`, `R13H8.1` and `CELE_R13H8.1` in one synonyms
        # list, mixing a genuine former name with two sequence names and saying which is
        # which nowhere. Reading them would put a kind on a claim no publisher made.
        assert any("daf-17" in line for line in worm_lines)
        read = read_bgi(worm_lines, ncbi_taxid=6239, origin="fixture")
        assert ("symbol", "daf-16", "WBGene00000912") in read
        assert not [triple for triple in read if "daf-17" in triple]

    def test_it_refuses_an_evidence_filter_rather_than_ignoring_one(
        self, worm_lines: list[str]
    ) -> None:
        with pytest.raises(EvidenceNotRecordedError, match="grades"):
            read_bgi(worm_lines, ncbi_taxid=6239, origin="x", evidence=("DIRECT",))


def _padded(fields: list[str], width: int) -> list[str]:
    """Return ``fields`` at exactly ``width``, so a reordering may index any column."""
    return (fields + [""] * width)[:width]


# ---------------------------------------------------------------------------
# The curated rows, and what a prepared set holds
# ---------------------------------------------------------------------------


class TestTheNewSourcesAreRows:
    def test_human_gains_hgnc_and_mouse_and_worm_gain_the_submission_source(self) -> None:
        assert HGNC_ARCHIVE in xref_sources("Homo sapiens")
        assert ALLIANCE_BGI in xref_sources("Mus musculus")
        assert ALLIANCE_BGI in xref_sources("Caenorhabditis elegans")

    def test_neither_displaces_the_default_source(self) -> None:
        for species in ("Homo sapiens", "Mus musculus", "Caenorhabditis elegans"):
            assert lookup_xref(species).source == ALLIANCE

    def test_hgnc_pins_a_dated_file_from_the_archive_listing(self) -> None:
        # The archive's dates are irregular — 2024-07-02, 2025-01-06, 2025-10-07, and both
        # 2026-07-03 and 2026-07-07 — so a URL built from "the first of the quarter" 404s
        # about half the time. The release *is* the file's date, and the URL ends in it.
        row = lookup_xref("Homo sapiens", HGNC_ARCHIVE)
        assert row.release == HGNC_RELEASE
        assert row.url.endswith(f"hgnc_complete_set_{HGNC_RELEASE}.txt")
        # The doubled path segment, which the remembered EBI FTP path is a 404 without.
        assert "/hgnc/archive/archive/quarterly/tsv/" in row.url

    def test_each_submission_row_names_the_authority_that_curated_it(self) -> None:
        assert "Mouse Genome Informatics" in lookup_xref("Mus musculus", ALLIANCE_BGI).publisher
        worm = lookup_xref("Caenorhabditis elegans", ALLIANCE_BGI)
        assert "WormBase" in worm.publisher
        # The publisher's own release identifier, which for worm is WormBase's own last one.
        assert worm.version == "WS298"


# ---------------------------------------------------------------------------
# Which set answers when no source is named
# ---------------------------------------------------------------------------


class TestTheDefaultWhenTheQuestionIsSymbols:
    """The **Default xref source** is per species *and per question* (ADR-0021).

    The papercut this closes: human's id default is the Alliance, whose cross-reference file
    carries no human symbol at all, so the epic's most common question — a gene list copied
    out of a paper — failed on the first try until a source was named. The flag that answers
    it is a second column of the curated table, and the path that reads it is
    :meth:`XrefSet.for_symbols`, which is the constructor with the question named.

    Both halves are asserted here, because the two behaviours sitting next to each other is
    the confusing part: the question-named path reaches the symbol-carrying source, and the
    plain constructor still does not — it carries no symbol and says exactly what to pass.

    :meth:`XrefSet.for_namespace` is the same fill-in for a caller holding a **Namespace**
    rather than a verb, which is what every caller reading one off a flag holds. It exists
    so the shell has nothing left to decide.
    """

    def test_every_species_flags_exactly_one_source_for_symbols(self) -> None:
        # Read off the shipped rows, so a fourth species added without the flag fails here
        # rather than at the first `match-symbols` somebody runs.
        flagged = {
            species: [
                row.source
                for row in metadata_mod.xref_table()
                if row.species == species and row.symbol_default
            ]
            for species in ("Homo sapiens", "Mus musculus", "Caenorhabditis elegans")
        }
        assert flagged == {
            "Homo sapiens": [HGNC_ARCHIVE],
            "Mus musculus": [ALLIANCE_BGI],
            "Caenorhabditis elegans": [ALLIANCE_BGI],
        }

    def test_the_two_defaults_are_different_rows_and_the_id_one_did_not_move(self) -> None:
        for species in ("Homo sapiens", "Mus musculus", "Caenorhabditis elegans"):
            assert lookup_xref(species).source == ALLIANCE
            assert lookup_xref(species, for_symbols=True).source != ALLIANCE

    def test_a_human_symbol_needs_no_source_named(self, sources: FakeFetch) -> None:
        # Story 10 of the epic, on the retired spelling the whole feature exists for.
        answer = XrefSet.for_symbols("Homo sapiens").match_symbols([RETIRED_EPIFACTORS_SYMBOL])

        assert answer.source == HGNC_ARCHIVE
        assert answer.resolved[RETIRED_EPIFACTORS_SYMBOL] == (
            SymbolMatch(symbol=RETIRED_EPIFACTORS_SYMBOL, gene_id_stem=BMAL1, kind=PREVIOUS),
        )

    def test_mouse_and_worm_reach_their_third_source_the_same_way(self, sources: FakeFetch) -> None:
        # Symbols come from a third source again for these two, so the flag has to hold for
        # them or it is a human special case wearing a general name.
        mouse, worm = (
            XrefSet.for_symbols("Mus musculus"),
            XrefSet.for_symbols("Caenorhabditis elegans"),
        )

        assert (mouse.source, worm.source) == (ALLIANCE_BGI, ALLIANCE_BGI)
        assert mouse.match_symbols(["Trp53"]).gene_id_stems == ["ENSMUSG00000059552"]
        assert worm.match_symbols(["daf-16"]).gene_id_stems == ["WBGene00000912"]

    def test_a_named_source_is_never_swapped_for_the_flagged_one(self) -> None:
        # Naming a source is how the scientific choice gets made deliberately, so the
        # question changes what an unnamed source resolves to and nothing else.
        assert lookup_xref("Homo sapiens", ALLIANCE, for_symbols=True).source == ALLIANCE

    def test_a_release_named_alone_is_still_honoured_against_the_flagged_source(self) -> None:
        # The release check runs after the source is filled in, exactly as it does for the
        # id default — a release asked for is answered or raises, never quietly swapped.
        assert lookup_xref("Homo sapiens", release=HGNC_RELEASE, for_symbols=True).source == (
            HGNC_ARCHIVE
        )

    def test_naming_the_namespace_reaches_the_same_set_as_naming_the_verb(
        self, sources: FakeFetch
    ) -> None:
        # The fill-in a caller holding a namespace asks for. It is `for_symbols` under
        # another name when that namespace is the symbol one, and must not be a second
        # opinion about which source answers.
        by_namespace = XrefSet.for_namespace("Homo sapiens", SYMBOL)
        by_verb = XrefSet.for_symbols("Homo sapiens")

        assert (by_namespace.source, by_namespace.release) == (by_verb.source, by_verb.release)
        assert by_namespace.source == HGNC_ARCHIVE

    def test_the_labelling_direction_answers_off_it_with_no_source_named(
        self, sources: FakeFetch
    ) -> None:
        # `from_stems(stems, "symbol")` is the other symbol question, so it reaches the same
        # source — which is what makes the rule a rule rather than one verb's special case.
        assert XrefSet.for_namespace("Homo sapiens", SYMBOL).from_stems(
            [BMAL1], SYMBOL
        ).resolved == {BMAL1: ("BMAL1",)}

    def test_every_other_namespace_reaches_the_identifier_default(
        self, alliance_human: XrefSet
    ) -> None:
        # Only the symbol question moves the default; for every other namespace this is the
        # ordinary constructor under a name that says which question is being asked.
        assert XrefSet.for_namespace("Homo sapiens", ENTREZ).source == ALLIANCE

    def test_a_source_named_beside_a_namespace_is_still_never_swapped(
        self, alliance_human: XrefSet
    ) -> None:
        # Naming a source is the deliberate scientific choice, and the question fills a
        # default in rather than overriding one — so this is the Alliance, symbols and all.
        chosen = XrefSet.for_namespace("Homo sapiens", SYMBOL, ALLIANCE)

        assert chosen.source == ALLIANCE
        with pytest.raises(NamespaceNotCarriedError):
            chosen.from_stems([BMAL1], SYMBOL)

    def test_the_plain_constructor_still_answers_from_the_id_default(
        self, alliance_human: XrefSet
    ) -> None:
        assert alliance_human.source == ALLIANCE
        assert alliance_human.symbol_kinds == ()

    def test_a_set_carrying_no_symbols_names_exactly_what_to_pass(
        self, alliance_human: XrefSet
    ) -> None:
        # The set genuinely carries no symbol, and answering from another publisher's bytes
        # is what one query reading exactly one set forbids (ADR-0017) — so it raises, and
        # the message names this species' symbol source and the call that fills it in.
        with pytest.raises(NamespaceNotCarriedError) as raised:
            alliance_human.match_symbols([RETIRED_EPIFACTORS_SYMBOL])

        message = str(raised.value)
        assert HGNC_ARCHIVE in message
        assert "for_symbols" in message
        # Mouse's and worm's source is not what to pass for a human set.
        assert ALLIANCE_BGI not in message

    def test_the_labelling_direction_misses_with_the_same_route(
        self, alliance_human: XrefSet
    ) -> None:
        # Both symbol questions miss on this set, so both must route the same way: naming
        # the namespaces it carries and stopping there sent nobody anywhere.
        with pytest.raises(NamespaceNotCarriedError) as labelling:
            alliance_human.from_stems([BMAL1], SYMBOL)
        with pytest.raises(NamespaceNotCarriedError) as matching:
            alliance_human.match_symbols([RETIRED_EPIFACTORS_SYMBOL])

        assert str(labelling.value) == str(matching.value)
        assert "for_symbols" in str(labelling.value)

    def test_no_other_namespace_is_sent_to_a_symbol_source(self, alliance_human: XrefSet) -> None:
        # The hint belongs to the symbol question alone: a human set asked in mouse's
        # namespace is not helped by being told where human symbols come from.
        with pytest.raises(NamespaceNotCarriedError) as raised:
            alliance_human.from_stems([BMAL1], MGI)

        message = str(raised.value)
        assert "for_symbols" not in message
        assert MGI in message

    def test_a_species_no_row_flags_names_the_sources_it_has(self) -> None:
        rows = tuple(
            replace(row, symbol_default=False)
            for row in metadata_mod.xref_table()
            if row.species == "Homo sapiens"
        )

        with pytest.raises(NoXrefSetError) as raised:
            lookup_xref("Homo sapiens", for_symbols=True, table=rows)

        for source in xref_sources("Homo sapiens"):
            assert source in str(raised.value)

    def test_the_miss_names_a_route_the_caller_can_actually_take(self) -> None:
        # The curated table ships inside the wheel, so "flag it in the table" is not
        # something the person reading this can do. Handing rows of their own to the lookup
        # is, and naming a source is the ordinary answer.
        rows = tuple(
            replace(row, symbol_default=False)
            for row in metadata_mod.xref_table()
            if row.species == "Homo sapiens"
        )

        with pytest.raises(NoXrefSetError) as raised:
            lookup_xref("Homo sapiens", for_symbols=True, table=rows)

        assert "table=" in str(raised.value)

    def test_one_source_is_not_a_symbol_default_just_by_being_the_only_one(self) -> None:
        # The id default answers a species' single source with no flag; this one does not. A
        # source that carries no symbol does not start carrying them by being alone.
        rows = (replace(lookup_xref("Homo sapiens"), symbol_default=False),)

        assert lookup_xref("Homo sapiens", table=rows).source == ALLIANCE
        with pytest.raises(NoXrefSetError):
            lookup_xref("Homo sapiens", for_symbols=True, table=rows)


class TestWhatAPreparedSetHolds:
    def test_each_set_carries_the_symbol_namespace_beside_its_ids(
        self, human: XrefSet, mouse: XrefSet, worm: XrefSet
    ) -> None:
        assert human.namespaces == (ENSEMBL, ENTREZ, UNIPROT, HGNC, SYMBOL)
        assert mouse.namespaces == (ENSEMBL, ENTREZ, UNIPROT, MGI, SYMBOL)
        assert worm.namespaces == (ENSEMBL, ENTREZ, UNIPROT, WORMBASE, SYMBOL)

    def test_the_kinds_a_set_can_match_are_read_off_the_slice(
        self, human: XrefSet, mouse: XrefSet, worm: XrefSet
    ) -> None:
        assert human.symbol_kinds == SYMBOL_KINDS
        assert mouse.symbol_kinds == (APPROVED,)
        assert worm.symbol_kinds == (APPROVED,)

    def test_the_stored_slice_says_which_kind_each_spelling_is(self, human: XrefSet) -> None:
        # A collaborator who does not use Python reads the set in a shell, so the kind is a
        # value in the namespace column and reads for itself there.
        with gzip.open(human.path, "rt", encoding="utf-8") as handle:
            namespaces = {line.split("\t")[0] for line in handle.read().splitlines()[1:]}
        assert set(KIND_NAMESPACES.values()) <= namespaces

    def test_the_completion_marker_records_the_kinds_the_set_carries(
        self, human: XrefSet, mouse: XrefSet
    ) -> None:
        typed, current = read_record(human.path.parent), read_record(mouse.path.parent)
        assert typed is not None
        assert current is not None
        assert typed.details["symbol_kinds"] == list(SYMBOL_KINDS)
        assert current.details["symbol_kinds"] == [APPROVED]

    def test_a_second_construction_re_reads_and_fetches_nothing(self, sources: FakeFetch) -> None:
        first = XrefSet("Homo sapiens", HGNC_ARCHIVE)
        fetched = len(sources.calls)
        again = XrefSet("Homo sapiens", HGNC_ARCHIVE)
        assert len(sources.calls) == fetched
        assert again.path == first.path
        assert again.match_symbols(["TP53"]).resolved == first.match_symbols(["TP53"]).resolved


# ---------------------------------------------------------------------------
# Away from the hub: labelling, and one symbol per stem
# ---------------------------------------------------------------------------


class TestFromAStem:
    def test_a_stem_answers_with_exactly_one_current_approved_symbol(self, human: XrefSet) -> None:
        answer = human.from_stems([BMAL1, TP53_STEM], SYMBOL)
        assert answer.resolved == {BMAL1: ("BMAL1",), TP53_STEM: ("TP53",)}
        assert all(len(symbols) == 1 for symbols in answer.resolved.values())

    def test_it_is_the_current_spelling_and_never_the_retired_one(self, human: XrefSet) -> None:
        # This direction is labelling a plot, so the retired spelling would be wrong on the
        # axis even though it is right to *match*.
        assert RETIRED_EPIFACTORS_SYMBOL not in human.from_stems([BMAL1], SYMBOL).xref_ids

    def test_every_species_labels_a_stem_with_its_own_authoritys_spelling(
        self, human: XrefSet, mouse: XrefSet, worm: XrefSet
    ) -> None:
        assert human.from_stems([BRCA1_STEM], SYMBOL).xref_ids == ["BRCA1"]
        assert mouse.from_stems(["ENSMUSG00000017146"], SYMBOL).xref_ids == ["Brca1"]
        assert worm.from_stems(["WBGene00000912"], SYMBOL).xref_ids == ["daf-16"]

    def test_a_versioned_gene_id_is_accepted_and_reduced_to_its_stem(self, human: XrefSet) -> None:
        assert human.from_stems([f"{BMAL1}.17"], SYMBOL).resolved == {f"{BMAL1}.17": ("BMAL1",)}


# ---------------------------------------------------------------------------
# Toward the hub: matching, and the kind on every match
# ---------------------------------------------------------------------------


class TestTowardAStem:
    def test_a_symbol_matches_approved_previous_and_alias_spellings(self, human: XrefSet) -> None:
        answer = human.match_symbols(["BMAL1", RETIRED_EPIFACTORS_SYMBOL, "MOP3"])
        assert answer.unresolved == ()
        assert {
            asked: {match.kind for match in matches} for asked, matches in answer.resolved.items()
        } == {"BMAL1": {APPROVED}, RETIRED_EPIFACTORS_SYMBOL: {PREVIOUS}, "MOP3": {ALIAS}}
        # All three spellings are one gene, which is why matching them all is worth doing.
        assert set(answer.gene_id_stems) == {BMAL1}

    def test_each_match_records_which_kind_it_was(self, human: XrefSet) -> None:
        matches = human.match_symbols([AMBIGUOUS_SYMBOL]).resolved[AMBIGUOUS_SYMBOL]
        assert matches == (
            SymbolMatch(symbol=AMBIGUOUS_SYMBOL, gene_id_stem=ADCY3_STEM, kind=APPROVED),
            SymbolMatch(symbol=AMBIGUOUS_SYMBOL, gene_id_stem=ADCY8_STEM, kind=PREVIOUS),
        )

    def test_a_symbol_naming_several_genes_answers_with_all_of_them(self, human: XrefSet) -> None:
        # `ADCY3` is HGNC's approved symbol for one gene and a symbol it retired from
        # another. Nothing here picks one, and the caller judges the ambiguity themselves.
        stems = {
            match.gene_id_stem
            for match in human.match_symbols([AMBIGUOUS_SYMBOL]).resolved[AMBIGUOUS_SYMBOL]
        }
        assert stems == {ADCY3_STEM, ADCY8_STEM}

    def test_a_symbol_the_release_knows_nothing_of_rides_back_in_ask_order(
        self, human: XrefSet
    ) -> None:
        answer = human.match_symbols(["nosuchgene", "TP53", "alsonot"])
        assert answer.unresolved == ("nosuchgene", "alsonot")
        assert list(answer.resolved) == ["TP53"]
        assert all(matches for matches in answer.resolved.values())

    def test_the_answer_names_the_source_and_the_release_that_produced_it(
        self, human: XrefSet, worm: XrefSet
    ) -> None:
        assert (human.source, human.release) == (HGNC_ARCHIVE, HGNC_RELEASE)
        answer = human.match_symbols(["TP53"])
        assert (answer.source, answer.release, answer.species) == (
            HGNC_ARCHIVE,
            HGNC_RELEASE,
            "Homo sapiens",
        )
        assert worm.match_symbols(["daf-16"]).source == ALLIANCE_BGI


class TestExactAndInsensitiveMatching:
    def test_matching_is_exact_by_default_and_a_mouse_cased_symbol_does_not_match(
        self, human: XrefSet
    ) -> None:
        # The species is fixed by the set, so `Brca1` asked of a human set is the wrong
        # authority's spelling rather than a typo — and saying so is better than half
        # working.
        answer = human.match_symbols([MOUSE_CASED_SYMBOL])
        assert answer.case_insensitive is False
        assert answer.resolved == {}
        assert answer.unresolved == (MOUSE_CASED_SYMBOL,)

    def test_case_insensitive_matching_is_opt_in(self, human: XrefSet) -> None:
        answer = human.match_symbols([MOUSE_CASED_SYMBOL], case_insensitive=True)
        assert answer.case_insensitive is True
        assert answer.resolved[MOUSE_CASED_SYMBOL] == (
            SymbolMatch(symbol="BRCA1", gene_id_stem=BRCA1_STEM, kind=APPROVED),
        )

    def test_the_insensitive_path_returns_every_match_rather_than_picking_one(
        self, human: XrefSet
    ) -> None:
        # Convenience never costs correctness: folding widens what matches and narrows
        # nothing, so the ambiguity that was there exactly is still there folded.
        exact = human.match_symbols([AMBIGUOUS_SYMBOL]).resolved[AMBIGUOUS_SYMBOL]
        folded = human.match_symbols(["adcy3"], case_insensitive=True).resolved["adcy3"]
        assert {(match.gene_id_stem, match.kind) for match in folded} == {
            (match.gene_id_stem, match.kind) for match in exact
        }

    def test_a_match_carries_the_authoritys_own_spelling_not_the_one_asked_about(
        self, human: XrefSet
    ) -> None:
        matched = human.match_symbols(["P53"], case_insensitive=True).resolved["P53"]
        assert [(match.symbol, match.kind) for match in matched] == [(LOWER_CASE_ALIAS, ALIAS)]
        # …and exactly, `P53` is nobody's spelling of it.
        assert human.match_symbols(["P53"]).unresolved == ("P53",)


class TestTheTwoDirectionsAreNotMirrorImages:
    def test_to_stems_refuses_the_symbol_namespace_and_names_the_verb_that_answers(
        self, human: XrefSet
    ) -> None:
        # Answering here would match approved spellings only, which drops exactly the rows
        # this whole thing exists for.
        with pytest.raises(SymbolDirectionError) as raised:
            human.to_stems([RETIRED_EPIFACTORS_SYMBOL], SYMBOL)
        assert "match_symbols" in str(raised.value)

    def test_the_stored_kind_namespaces_are_not_ones_a_caller_may_name(
        self, human: XrefSet
    ) -> None:
        # `previous_symbol` and `alias_symbol` are how the kind is stored, not identifier
        # systems of their own, so neither is a namespace the set offers.
        for stored in ("previous_symbol", "alias_symbol"):
            assert stored not in human.namespaces
            with pytest.raises(NamespaceNotCarriedError):
                human.to_stems([RETIRED_EPIFACTORS_SYMBOL], stored)

    def test_a_set_whose_source_carries_no_symbols_matches_none(
        self, alliance_human: XrefSet
    ) -> None:
        # The Alliance's cross-reference file carries symbols for worm alone, in rows this
        # package reads as gene ids, and none at all for human — measured on the whole
        # 25 MB file. So a human set built on it answers no symbol at all and raises rather
        # than reaching for another publisher's bytes. What the message must name is
        # `TestTheDefaultWhenTheQuestionIsSymbols`'s.
        assert alliance_human.symbol_kinds == ()
        assert SYMBOL not in alliance_human.namespaces
        with pytest.raises(NamespaceNotCarriedError):
            alliance_human.match_symbols(["TP53"])


# ---------------------------------------------------------------------------
# Mouse and worm: current symbols, and why there are no others
# ---------------------------------------------------------------------------


class TestCurrentSymbolsOnly:
    def test_mouse_answers_current_symbols(self, mouse: XrefSet) -> None:
        answer = mouse.match_symbols(["Trp53", "Bmal1"])
        assert answer.unresolved == ()
        assert answer.resolved["Trp53"] == (
            SymbolMatch(symbol="Trp53", gene_id_stem="ENSMUSG00000059552", kind=APPROVED),
        )
        assert {match.kind for matches in answer.resolved.values() for match in matches} == {
            APPROVED
        }

    def test_mouse_says_why_it_has_no_previous_or_alias_matching(self, mouse: XrefSet) -> None:
        # The explanation is part of the behaviour, not a comment: without it, a spelling
        # MGI retired comes back unresolved and looks exactly like a gene that is absent.
        answer = mouse.match_symbols(["Arntl"])
        assert answer.unresolved == ("Arntl",)
        assert answer.kinds == (APPROVED,)
        assert answer.limits is not None
        assert "MGI" in answer.limits
        assert "no dated archive" in answer.limits
        assert PREVIOUS in answer.limits
        assert ALIAS in answer.limits

    def test_worm_answers_current_symbols_and_says_the_same(self, worm: XrefSet) -> None:
        answer = worm.match_symbols(["daf-16", "daf-17"])
        assert answer.resolved["daf-16"] == (
            SymbolMatch(symbol="daf-16", gene_id_stem="WBGene00000912", kind=APPROVED),
        )
        # `daf-17` is in the submission's synonyms list beside two sequence names, untyped.
        assert answer.unresolved == ("daf-17",)
        assert answer.limits is not None
        assert "WormBase" in answer.limits
        assert "403" in answer.limits

    def test_a_set_that_matches_every_kind_carries_no_such_note(self, human: XrefSet) -> None:
        assert human.match_symbols(["TP53"]).limits is None

    def test_the_note_is_the_sources_and_never_the_species(
        self, mouse: XrefSet, worm: XrefSet
    ) -> None:
        assert mouse.match_symbols(["Trp53"]).limits == BGI_SYMBOL_LIMIT
        assert worm.match_symbols(["daf-16"]).limits == BGI_SYMBOL_LIMIT


# ---------------------------------------------------------------------------
# The measured failure this exists to prevent
# ---------------------------------------------------------------------------


class TestTheMeasuredEpiFactorsFailure:
    def test_a_retired_symbol_from_the_measured_epifactors_set_resolves(
        self, human: XrefSet
    ) -> None:
        # Of EpiFactors v2.0's 801 human rows, 31 spell the gene by a symbol HGNC has since
        # retired — `ARNTL` for `BMAL1` among them, named in the shipped cofactor table's
        # own attribution. A join that knows approved spellings only drops exactly those.
        shipped = _shipped_cofactor_stems()
        assert shipped["BMAL1"] == BMAL1
        matched = human.match_symbols([RETIRED_EPIFACTORS_SYMBOL])
        assert matched.resolved[RETIRED_EPIFACTORS_SYMBOL] == (
            SymbolMatch(symbol=RETIRED_EPIFACTORS_SYMBOL, gene_id_stem=BMAL1, kind=PREVIOUS),
        )

    def test_a_second_retired_spelling_reaches_the_same_shipped_row(self, human: XrefSet) -> None:
        shipped = _shipped_cofactor_stems()
        matched = human.match_symbols([RETIRED_EPIFACTORS_SYMBOL_TOO])
        assert [
            match.gene_id_stem for match in matched.resolved[RETIRED_EPIFACTORS_SYMBOL_TOO]
        ] == [shipped["EMSY"]]

    def test_the_approved_only_join_is_what_would_have_dropped_them(self, human: XrefSet) -> None:
        # The counterfactual, stated as a test rather than as prose: reading the symbol
        # namespace as an id — approved spellings only — finds nothing for either.
        retired = [RETIRED_EPIFACTORS_SYMBOL, RETIRED_EPIFACTORS_SYMBOL_TOO]
        approved_only = {
            symbol
            for symbol in retired
            if any(
                match.kind == APPROVED for match in human.match_symbols([symbol]).resolved[symbol]
            )
        }
        assert approved_only == set()
        assert set(human.match_symbols(retired).resolved) == set(retired)


# ---------------------------------------------------------------------------
# The answer's shape
# ---------------------------------------------------------------------------


class TestAnswerShape:
    def test_no_resolved_value_is_ever_empty(self, human: XrefSet) -> None:
        answer = human.match_symbols(["TP53", "nosuchgene", AMBIGUOUS_SYMBOL])
        assert all(matches for matches in answer.resolved.values())

    def test_resolved_and_unresolved_partition_what_was_asked(self, human: XrefSet) -> None:
        asked = ["TP53", "nosuchgene", AMBIGUOUS_SYMBOL, MOUSE_CASED_SYMBOL]
        answer = human.match_symbols(asked)
        assert len(answer.resolved) + len(answer.unresolved) == len(asked)
        assert set(answer.resolved) | set(answer.unresolved) == set(asked)
        assert set(answer.resolved) & set(answer.unresolved) == set()

    def test_the_flattener_loses_the_kind_and_says_so(self, human: XrefSet) -> None:
        answer = human.match_symbols([AMBIGUOUS_SYMBOL])
        assert answer.gene_id_stems == [ADCY3_STEM, ADCY8_STEM]
        assert "kind" in (ResolvedSymbols.gene_id_stems.__doc__ or "")

    def test_json_carries_every_match_with_its_kind_and_what_matched_nothing(
        self, human: XrefSet
    ) -> None:
        payload = human.match_symbols([AMBIGUOUS_SYMBOL, "nosuchgene"]).as_json()
        assert payload["source"] == HGNC_ARCHIVE
        assert payload["release"] == HGNC_RELEASE
        assert payload["kinds"] == list(SYMBOL_KINDS)
        assert payload["limits"] is None
        assert payload["case_insensitive"] is False
        assert payload["unresolved"] == ["nosuchgene"]
        assert [match["kind"] for match in payload["resolved"][AMBIGUOUS_SYMBOL]] == [
            APPROVED,
            PREVIOUS,
        ]
        assert payload["gene_id_stems"] == [ADCY3_STEM, ADCY8_STEM]

    def test_a_repeated_symbol_is_asked_once(self, human: XrefSet) -> None:
        answer = human.match_symbols(["TP53", "TP53"])
        assert list(answer.resolved) == ["TP53"]


# ---------------------------------------------------------------------------
# Normalisation, over generated input
# ---------------------------------------------------------------------------


class TestSymbolNormalisation:
    @given(st.text())
    def test_it_is_idempotent(self, symbol: str) -> None:
        assert normalise_symbol(normalise_symbol(symbol)) == normalise_symbol(symbol)

    @given(st.text())
    def test_folding_is_idempotent_and_agrees_with_normalisation(self, symbol: str) -> None:
        assert fold_symbol(fold_symbol(symbol)) == fold_symbol(symbol)
        assert fold_symbol(symbol) == normalise_symbol(symbol).casefold()

    @given(st.text())
    def test_it_never_cuts_a_symbol_at_a_dot(self, symbol: str) -> None:
        # The reason a symbol does not go through `normalise_id`: WormBase names thousands
        # of genes by their sequence name, and `Y110A7A.10` stemmed is another gene's
        # spelling or none at all.
        assert normalise_symbol(symbol) == symbol.strip()

    @given(st.sampled_from(["TP53", "daf-16", "Y110A7A.10", "bHLHe5"]))
    def test_a_real_spelling_folds_onto_itself_lower_cased(self, symbol: str) -> None:
        assert fold_symbol(symbol) == symbol.lower()
