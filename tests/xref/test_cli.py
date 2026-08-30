"""Tests for the ``genome xref`` sub-app."""

from __future__ import annotations

import gzip
import hashlib
import json as _json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from typer.testing import Result

from genome.cli import app
from genome.store import fetch as fetch_mod
from genome.xref import (
    ALIAS,
    ALLIANCE,
    ALLIANCE_BGI,
    APPROVED,
    ENSEMBL,
    ENSEMBL_TSV,
    ENTREZ,
    HGNC,
    HGNC_ARCHIVE,
    MGI,
    NAMESPACES,
    PREVIOUS,
    SYMBOL,
    UNIPROT,
    XrefSet,
    lookup_xref,
    xref_prepare_command,
    xref_set_dir,
    xref_slice_name,
    xref_species,
    xref_table,
)
from genome.xref import metadata as xref_metadata_mod

from .._cli import help_text, output, runner
from ..conftest import FakeFetch
from .test_xref import FIXTURE as _XREF_FIXTURE
from .test_xref import HUMAN_GENE_WITHOUT_A_HUB as _XREF_NO_HUB
from .test_xref import RELEASE as _XREF_RELEASE
from .test_xref_symbols import ADCY3_STEM as _ADCY3_STEM
from .test_xref_symbols import ADCY8_STEM as _ADCY8_STEM
from .test_xref_symbols import AMBIGUOUS_SYMBOL as _AMBIGUOUS_SYMBOL
from .test_xref_symbols import BGI_MOUSE as _BGI_MOUSE_FIXTURE
from .test_xref_symbols import BMAL1 as _BMAL1_STEM
from .test_xref_symbols import HGNC_FIXTURE as _HGNC_FIXTURE
from .test_xref_symbols import HGNC_RELEASE as _HGNC_RELEASE
from .test_xref_symbols import MOUSE_CASED_SYMBOL as _MOUSE_CASED_SYMBOL
from .test_xref_symbols import RETIRED_EPIFACTORS_SYMBOL as _RETIRED_SYMBOL
from .test_xref_symbols import RETIRED_EPIFACTORS_SYMBOL_TOO as _RETIRED_SYMBOL_TOO
from .test_xref_symbols import _digest as _fixture_digest

#: The species ``genome xref ids`` is exercised against. Human, because the committed
#: Alliance fixture carries the two things the command's shape exists for: a foreign id
#: naming two **Gene id stem**s, and a real gene with no Ensembl cross-reference to reach.
_XREF_SPECIES = "Homo sapiens"

#: The mouse spelling of the species, which the symbol sources are pinned per species.
_MOUSE = "Mus musculus"

#: The one gene the HGNC fixture spells all three ways, so that a single invocation of
#: ``genome xref symbols`` can be asked for an approved, a previous and an alias match at
#: once. The previous spelling is imported rather than written here: it is one of the 31
#: measured EpiFactors rows and ``tests/xref/test_xref_symbols.py`` is where that is recorded.
_APPROVED_SYMBOL, _ALIAS_SYMBOL = "BMAL1", "MOP3"

#: The stem the second measured EpiFactors spelling reaches — EMSY's, as the shipped human
#: **Cofactor table** records it, which is the table the 801-row measurement was made on.
_EMSY_STEM = "ENSG00000158636"

#: A symbol MGI has retired, asked of the mouse set to show what *cannot* be matched there:
#: that source publishes one current approved spelling per gene, so this comes back
#: unresolved and the answer must say why rather than letting it read as an absent gene.
_MOUSE_RETIRED_SYMBOL = "Arntl"

#: One spelling MGI still approves, and the gene it names. Mouse's symbols come from a third
#: source again, so the source a question fills in is asked for here as well as for human.
_MOUSE_APPROVED_SYMBOL, _TRP53_STEM = "Trp53", "ENSMUSG00000059552"


@pytest.fixture
def xref_pinned(monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> None:
    """Pin the curated **Xref source** rows to the committed Alliance fixture's digest.

    The arrangement ``tests/xref/test_xref.py`` uses: the checksum check that holds a truncated
    download to be an error rather than a quietly short answer is never switched off, only
    pointed at what the fake fetch actually serves. Every other cell survives, the real URL
    among them, so what the command prints as provenance is the shipped row's own.
    """
    with gzip.open(data_dir / _XREF_FIXTURE, "rb") as handle:
        # md5 because that is the algorithm Alliance publishes, not a choice made here.
        digest = f"md5:{hashlib.md5(handle.read()).hexdigest()}"
    rows = tuple(replace(row, source_checksum=digest) for row in xref_table())
    monkeypatch.setattr(xref_metadata_mod, "xref_table", lambda: rows)


@pytest.fixture
def xref_release(fake_fetch: FakeFetch, xref_pinned: None) -> FakeFetch:
    """Serve the committed Alliance slice as the publisher's file, pinned to it."""
    fake_fetch.serve(_XREF_FIXTURE)
    return fake_fetch


@pytest.fixture
def symbol_pinned(monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> dict[str, Path]:
    """Pin the symbol-carrying **Xref source** rows to the fixtures cut out of their files.

    ``tests/xref/test_xref_symbols.py``'s arrangement: the checksum check is never switched off,
    only pointed at the bytes that actually arrive. The mapping it returns is what each
    URL should serve, so a fetch can be routed by URL and the human and the mouse set can
    be held at once — which one test below needs, the two sources matching different kinds.
    """
    by_url = {
        lookup_xref(_XREF_SPECIES, HGNC_ARCHIVE).url: data_dir / _HGNC_FIXTURE,
        lookup_xref(_MOUSE, ALLIANCE_BGI).url: data_dir / _BGI_MOUSE_FIXTURE,
    }
    rows = tuple(
        replace(row, source_checksum=f"md5:{_fixture_digest(by_url[row.url])}")
        if row.url in by_url
        else row
        for row in xref_table()
    )
    monkeypatch.setattr(xref_metadata_mod, "xref_table", lambda: rows)
    return by_url


@pytest.fixture
def symbol_sources(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, symbol_pinned: dict[str, Path]
) -> FakeFetch:
    """Serve each symbol source's publisher file, routed by the URL the package built."""

    def route(url: str, dest_dir: Path, **kwargs: Any) -> Path:
        fake_fetch.serve(symbol_pinned[url])
        return fake_fetch(url, dest_dir, **kwargs)

    monkeypatch.setattr(fetch_mod, "fetch_url", route)
    return fake_fetch


def _match_symbols(*arguments: str) -> Result:
    """Ask the human HGNC set, the one shipped source that carries all three kinds."""
    return runner.invoke(
        app, ["xref", "symbols", _XREF_SPECIES, "--source", HGNC_ARCHIVE, *arguments]
    )


class TestXrefCommand:
    """``genome xref ids`` — the shell surface over an **Xref set**, driven off the fixture.

    Offline throughout, the way ``tests/xref/test_xref.py`` is: the fake fetch serves the
    committed Alliance slice and the curated rows are pinned to it, so the command prepares
    and reads a real set under the test's own data root. What is asserted here is the
    command and not the hop — ``tests/xref/test_xref.py`` owns that — so: the direction being
    named rather than sniffed out of the id strings, the stdout/stderr split that makes the
    output pipe, the ids that resolved to nothing staying visible in *both* renderings, and
    a non-zero exit naming the next action for each way it can fail.

    The strongest claim here is that the command holds no logic the API does not, and it is
    checked the only way that can be: ``--json`` is asserted equal to what the same two
    calls answer in Python, whole, in both directions.
    """

    def test_it_converts_ids_never_infers_the_direction_and_json_matches_the_api(
        self, xref_release: FakeFetch
    ) -> None:
        to_stems = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "672"]
        )
        assert to_stems.exit_code == 0
        assert to_stems.stdout.splitlines() == [
            "7157\tENSG00000141510",
            "672\tENSG00000012048",
        ]

        from_stems = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--from-stems", HGNC, "ENSG00000141510"]
        )
        assert from_stems.exit_code == 0
        assert from_stems.stdout.splitlines() == ["ENSG00000141510\tHGNC:11998"]

        # One string, one namespace, two directions, two different answers. `HGNC:11998`
        # is an HGNC id and is not a **Gene id stem**, and nothing here works that out
        # from the characters: the flag says which way the hop goes, and asking the wrong
        # way answers *nothing found* rather than quietly turning around.
        toward = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", HGNC, "HGNC:11998"]
        )
        away = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--from-stems", HGNC, "HGNC:11998"]
        )
        assert toward.exit_code == away.exit_code == 0
        assert toward.stdout.splitlines() == ["HGNC:11998\tENSG00000141510"]
        assert away.stdout.splitlines() == ["HGNC:11998\t"]

        # The strongest claim here is that the command holds no logic the API does not:
        # `--json` is asserted equal to what the same call answers in Python, whole, in
        # both directions.
        to_stems_asked = ["7157", "8086", "999999999"]
        json_toward = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, *to_stems_asked, "--json"]
        )
        assert json_toward.exit_code == 0
        assert (
            _json.loads(json_toward.stdout)
            == XrefSet(_XREF_SPECIES).to_stems(to_stems_asked, ENTREZ).as_json()
        )

        from_stems_asked = ["ENSG00000141510", "ENSG00000012048", "ENSG00000288541"]
        json_away = runner.invoke(
            app,
            ["xref", "ids", _XREF_SPECIES, "--from-stems", UNIPROT, *from_stems_asked, "--json"],
        )
        assert json_away.exit_code == 0
        assert (
            _json.loads(json_away.stdout)
            == XrefSet(_XREF_SPECIES).from_stems(from_stems_asked, UNIPROT).as_json()
        )

        # Naming zero or two directions exits 2 before anything is resolved.
        neither = runner.invoke(app, ["xref", "ids", _XREF_SPECIES, "7157"])
        assert neither.exit_code == 2
        assert neither.stdout == ""
        assert "--to-stems" in output(neither)
        assert "--from-stems" in output(neither)

        both = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "--from-stems", HGNC, "7157"]
        )
        assert both.exit_code == 2
        assert both.stdout == ""
        assert "exactly one" in output(both)

    def test_the_render_keeps_every_id_the_provenance_beside_it_and_ambiguity_intact(
        self, xref_release: FakeFetch
    ) -> None:
        result = runner.invoke(app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])
        assert result.exit_code == 0
        assert result.stdout == "7157\tENSG00000141510\n"
        # Which publisher said so, and when, is what a reader needs and what a pipeline
        # must not be handed: it goes beside the pairs rather than among them.
        assert f"{ALLIANCE} {_XREF_RELEASE}" in result.stderr
        assert _XREF_SPECIES in result.stderr
        assert lookup_xref(_XREF_SPECIES).url in result.stderr
        assert f"{ENTREZ} ids -> gene id stems" in result.stderr

        # 6.2% of human HGNC ids name two stems in this release, so this is the ordinary
        # case rather than an edge one, and nothing here picks between them.
        ambiguous = runner.invoke(app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "8086"])
        assert ambiguous.exit_code == 0
        assert ambiguous.stdout.splitlines() == [
            "8086\tENSG00000094914",
            "8086\tENSG00000291836",
        ]

        # The one thing a hand-rolled join drops. `HGNC:10041` is a real human gene the
        # Alliance lists with no Ensembl cross-reference at all, so it has no hub to
        # reach, and it stays visible in both renderings rather than being dropped.
        unresolved = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", HGNC, "HGNC:11998", _XREF_NO_HUB]
        )
        assert unresolved.exit_code == 0
        assert unresolved.stdout.splitlines() == [
            "HGNC:11998\tENSG00000141510",
            f"{_XREF_NO_HUB}\t",
        ]
        assert "1 this release names none for" in unresolved.stderr

        unresolved_json = runner.invoke(
            app,
            [
                "xref",
                "ids",
                _XREF_SPECIES,
                "--to-stems",
                HGNC,
                "HGNC:11998",
                _XREF_NO_HUB,
                "--json",
            ],
        )
        assert unresolved_json.exit_code == 0
        payload = _json.loads(unresolved_json.stdout)
        assert payload["unresolved"] == [_XREF_NO_HUB]
        assert _XREF_NO_HUB not in payload["resolved"]

    def test_an_unsupported_source_species_or_namespace_exits_one_but_a_named_source_resolves(
        self, xref_release: FakeFetch
    ) -> None:
        # "ncbi" rather than a source that merely happens to be unlisted today: NCBI Gene
        # is excluded on purpose, being rebuilt in place with no retrievable old release
        # (ADR-0018), so this stays a miss however many sources the table grows.
        bad_source = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--source", "ncbi", "--to-stems", ENTREZ, "7157"]
        )
        assert bad_source.exit_code == 1
        assert bad_source.stdout == ""
        assert ALLIANCE in output(bad_source)
        assert ENSEMBL_TSV in output(bad_source)

        bad_species = runner.invoke(
            app, ["xref", "ids", "Danio rerio", "--to-stems", ENTREZ, "7157"]
        )
        assert bad_species.exit_code == 1
        assert bad_species.stdout == ""
        for species in xref_species():
            assert species in output(bad_species)
        # Refused before anything was fetched: a species with no Ensembl presence has no
        # hub to hang a namespace off and is unanswerable by design, not pending a download.
        assert xref_release.calls == []

        # The three species have three different authorities, so a human set asked for
        # MGI ids fails loudly rather than answering nothing — in either direction — the
        # failure that would otherwise look like a gene list with no matches.
        bad_namespace = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", MGI, "MGI:88276"]
        )
        assert bad_namespace.exit_code == 1
        assert bad_namespace.stdout == ""
        for namespace in (ENSEMBL, ENTREZ, UNIPROT, HGNC):
            assert namespace in output(bad_namespace)

        bad_namespace_reverse = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--from-stems", MGI, "ENSG00000141510"]
        )
        assert bad_namespace_reverse.exit_code == 1
        assert bad_namespace_reverse.stdout == ""
        assert HGNC in output(bad_namespace_reverse)

        # Refused before any of these real hops touched the network — a source may be
        # named, and omitting it defaults exactly as the API does.
        named = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--source", ALLIANCE, "--to-stems", ENTREZ, "7157"]
        )
        assert named.exit_code == 0
        assert named.stdout == "7157\tENSG00000141510\n"
        assert f"{ALLIANCE} {_XREF_RELEASE}" in named.stderr

        # Which default an unnamed source is filled in with is the API's decision, named
        # by handing it the namespace the flags carried — so the command and the Python
        # call cannot resolve differently. The symbol namespace is the case that matters
        # and `TestMatchSymbolsCommand` asserts it; this is the same claim for every
        # other one.
        named_json = runner.invoke(
            app,
            [
                "xref",
                "ids",
                _XREF_SPECIES,
                "--source",
                ALLIANCE,
                "--to-stems",
                ENTREZ,
                "7157",
                "--json",
            ],
        )
        defaulted_json = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "--json"]
        )
        assert named_json.exit_code == defaulted_json.exit_code == 0
        defaulted_payload = _json.loads(defaulted_json.stdout)
        assert defaulted_payload == _json.loads(named_json.stdout)
        assert lookup_xref(_XREF_SPECIES).default is True
        assert (defaulted_payload["source"], defaulted_payload["release"]) == (
            ALLIANCE,
            _XREF_RELEASE,
        )
        assert (defaulted_payload["species"], defaulted_payload["namespace"]) == (
            _XREF_SPECIES,
            ENTREZ,
        )
        assert (
            defaulted_payload
            == XrefSet.for_namespace(_XREF_SPECIES, ENTREZ).to_stems(["7157"], ENTREZ).as_json()
        )

    def test_the_symbol_namespace_toward_the_hub_names_match_symbols_and_the_help_agrees(
        self, xref_release: FakeFetch
    ) -> None:
        # The API refuses this call naming `match_symbols(symbols)`, which is a Python call
        # and no next action at all for someone in a shell. The surface owns which of *its*
        # spellings to reach for, so the refusal happens here and names the command.
        result = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", SYMBOL, _RETIRED_SYMBOL]
        )
        assert result.exit_code == 2
        assert result.stdout == ""
        assert "genome xref symbols" in output(result)
        # Refused before anything was fetched: no set has to be prepared to know that this
        # direction is answered by another command.
        assert xref_release.calls == []

        # Help and behaviour agree or the command advertises a conversion it will not make.
        offered = help_text("xref", "ids")
        assert ", ".join(name for name in NAMESPACES if name != SYMBOL) in offered
        assert ", ".join(NAMESPACES) not in offered
        assert "genome xref symbols" in offered

    def test_a_set_not_downloaded_or_left_unfinished_exits_one_naming_the_fix(
        self, monkeypatch: pytest.MonkeyPatch, xref_pinned: None
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        not_downloaded = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"]
        )
        assert not_downloaded.exit_code == 1
        assert not_downloaded.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE) in output(
            not_downloaded
        )
        assert "login node" in output(not_downloaded)

        # A leftover half-written slice is caught before the (still-broken) fetch ever
        # runs: the completion check runs first, so the two never race.
        directory = xref_set_dir(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE)
        directory.mkdir(parents=True)
        (directory / xref_slice_name(_XREF_SPECIES)).write_bytes(b"half a file")
        unfinished = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"]
        )
        assert unfinished.exit_code == 1
        assert unfinished.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, ALLIANCE, _XREF_RELEASE) in output(unfinished)

    def test_the_progress_display_is_suppressed_under_json(self, xref_release: FakeFetch) -> None:
        runner.invoke(app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157", "--json"])
        assert xref_release.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(self, xref_release: FakeFetch) -> None:
        runner.invoke(app, ["xref", "ids", _XREF_SPECIES, "--to-stems", ENTREZ, "7157"])
        assert xref_release.last.progressbar is True


class TestMatchSymbolsCommand:
    """``genome xref symbols`` — the shell surface over ``XrefSet.match_symbols``.

    Offline throughout, the way ``tests/xref/test_xref_symbols.py`` is: the fake fetch serves the
    committed HGNC archive cut and the Alliance's MGI submission, routed by the URL the
    package built, and the curated rows are pinned to them. What is asserted here is the
    command and not the match — ``tests/xref/test_xref_symbols.py`` owns the three kinds, the
    folding and the asymmetry — so: a symbol reaching its genes from a shell at all, the
    kind of every match surviving the render, exact matching being what happens unless case
    is folded on purpose, the symbols that matched nothing staying visible in *both*
    renderings, and a non-zero exit naming the next action for each way it can fail.

    It is a command of its own rather than a third direction of ``genome xref ids`` for the
    reason ``match_symbols`` is a verb of its own: a symbol matches spellings the authority
    has retired, answers with every gene any of them names, and carries a kind on each
    match, so it is neither of that command's two hops and does not render as one.

    The strongest claim here is that the command holds no logic the API does not, and it is
    checked the only way that can be: ``--json`` is asserted equal to what the same call
    answers in Python, whole, both exactly and folded.
    """

    def test_every_match_says_its_kind_and_matching_is_exact_unless_case_is_folded(
        self, symbol_sources: FakeFetch
    ) -> None:
        # One gene, spelled three ways, in one call. Which kind matched is the answer's and
        # not a detail: without it a retired spelling and a current one read alike — and a
        # single approved symbol on its own comes back the same way, unadorned.
        result = _match_symbols(_APPROVED_SYMBOL, _RETIRED_SYMBOL, _ALIAS_SYMBOL)
        assert result.exit_code == 0
        assert [line.split("\t") for line in result.stdout.splitlines()] == [
            [_APPROVED_SYMBOL, _APPROVED_SYMBOL, _BMAL1_STEM, APPROVED],
            [_RETIRED_SYMBOL, _RETIRED_SYMBOL, _BMAL1_STEM, PREVIOUS],
            [_ALIAS_SYMBOL, _ALIAS_SYMBOL, _BMAL1_STEM, ALIAS],
        ]

        # The failure the whole symbol feature exists to prevent, reached from a shell: of
        # EpiFactors v2.0's 801 human rows, 31 spell the gene by a symbol HGNC has since
        # retired, and a join that knows approved spellings only drops exactly those.
        retired = _match_symbols(_RETIRED_SYMBOL, _RETIRED_SYMBOL_TOO)
        assert retired.exit_code == 0
        assert retired.stdout.splitlines() == [
            f"{_RETIRED_SYMBOL}\t{_RETIRED_SYMBOL}\t{_BMAL1_STEM}\t{PREVIOUS}",
            f"{_RETIRED_SYMBOL_TOO}\t{_RETIRED_SYMBOL_TOO}\t{_EMSY_STEM}\t{PREVIOUS}",
        ]

        # `ADCY3` is HGNC's approved symbol for one gene and a symbol it retired from
        # another, so ambiguity is the ordinary case rather than an edge one, and the
        # shell is handed both with the kind that distinguishes them.
        ambiguous = _match_symbols(_AMBIGUOUS_SYMBOL)
        assert ambiguous.exit_code == 0
        assert ambiguous.stdout.splitlines() == [
            f"{_AMBIGUOUS_SYMBOL}\t{_AMBIGUOUS_SYMBOL}\t{_ADCY3_STEM}\t{APPROVED}",
            f"{_AMBIGUOUS_SYMBOL}\t{_AMBIGUOUS_SYMBOL}\t{_ADCY8_STEM}\t{PREVIOUS}",
        ]

        # The species is fixed by the set, so a mouse-cased spelling asked of a human set is
        # the wrong authority's rather than a typo to absorb — and saying so beats half
        # working.
        exact = _match_symbols(_MOUSE_CASED_SYMBOL)
        assert exact.exit_code == 0
        assert exact.stdout.splitlines() == [f"{_MOUSE_CASED_SYMBOL}\t\t\t"]
        assert "exact" in exact.stderr

        # Convenience never costs correctness: folding widens what matches and narrows
        # nothing, so the ambiguity that was there exactly is still there folded — and the
        # authority's own spelling comes back beside the one that was asked about.
        folded = _match_symbols("adcy3", "--case-insensitive")
        assert folded.exit_code == 0
        assert folded.stdout.splitlines() == [
            f"adcy3\t{_AMBIGUOUS_SYMBOL}\t{_ADCY3_STEM}\t{APPROVED}",
            f"adcy3\t{_AMBIGUOUS_SYMBOL}\t{_ADCY8_STEM}\t{PREVIOUS}",
        ]
        assert "case-insensitive" in folded.stderr

    def test_json_carries_unresolved_symbols_matches_the_api_and_agrees_with_the_text(
        self, symbol_sources: FakeFetch
    ) -> None:
        # What your list holds and this release does not is the one thing a hand-rolled
        # join drops silently, so it gets a row here like everything else — in both
        # renderings.
        result = _match_symbols(_APPROVED_SYMBOL, "nosuchgene")
        assert result.exit_code == 0
        assert result.stdout.splitlines() == [
            f"{_APPROVED_SYMBOL}\t{_APPROVED_SYMBOL}\t{_BMAL1_STEM}\t{APPROVED}",
            "nosuchgene\t\t\t",
        ]
        assert "1 this release matched nothing for" in result.stderr

        json_result = _match_symbols(_APPROVED_SYMBOL, "nosuchgene", "--json")
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["unresolved"] == ["nosuchgene"]
        assert "nosuchgene" not in payload["resolved"]

        # The strongest claim here is that the command holds no logic the API does not:
        # `--json` is asserted equal to what the same call answers in Python, whole, both
        # exactly and folded.
        exact_asked = [_AMBIGUOUS_SYMBOL, _RETIRED_SYMBOL, "nosuchgene"]
        exact = _match_symbols(*exact_asked, "--json")
        assert exact.exit_code == 0
        assert (
            _json.loads(exact.stdout)
            == XrefSet(_XREF_SPECIES, HGNC_ARCHIVE).match_symbols(exact_asked).as_json()
        )

        folded_asked = [_MOUSE_CASED_SYMBOL, "adcy3"]
        folded = _match_symbols(*folded_asked, "--case-insensitive", "--json")
        assert folded.exit_code == 0
        assert (
            _json.loads(folded.stdout)
            == XrefSet(_XREF_SPECIES, HGNC_ARCHIVE)
            .match_symbols(folded_asked, case_insensitive=True)
            .as_json()
        )

        # The matches go to stdout and the provenance to stderr so the output pipes.
        stdout_result = _match_symbols(_APPROVED_SYMBOL)
        assert stdout_result.exit_code == 0
        assert (
            stdout_result.stdout
            == f"{_APPROVED_SYMBOL}\t{_APPROVED_SYMBOL}\t{_BMAL1_STEM}\tapproved\n"
        )
        # Which authority said so, and when, is what a reader needs and what a pipeline
        # must not be handed: it goes beside the matches rather than among them.
        assert f"{HGNC_ARCHIVE} {_HGNC_RELEASE}" in stdout_result.stderr
        assert _XREF_SPECIES in stdout_result.stderr
        assert lookup_xref(_XREF_SPECIES, HGNC_ARCHIVE).url in stdout_result.stderr
        assert "gene symbols -> gene id stems" in stdout_result.stderr

        source_json = _match_symbols(_APPROVED_SYMBOL, "--json")
        source_payload = _json.loads(source_json.stdout)
        assert (source_payload["source"], source_payload["release"]) == (
            HGNC_ARCHIVE,
            _HGNC_RELEASE,
        )
        assert (source_payload["species"], source_payload["case_insensitive"]) == (
            _XREF_SPECIES,
            False,
        )

        # The whole claim that the command holds no logic: every text cell printed is a
        # value the JSON answer put there, in the order the JSON answer writes it.
        text = _match_symbols(_AMBIGUOUS_SYMBOL)
        text_payload = _match_symbols(_AMBIGUOUS_SYMBOL, "--json")
        matches = _json.loads(text_payload.stdout)["resolved"][_AMBIGUOUS_SYMBOL]
        assert text.stdout.splitlines() == [
            "\t".join([_AMBIGUOUS_SYMBOL, match["symbol"], match["gene_id_stem"], match["kind"]])
            for match in matches
        ]

    def test_a_symbol_or_mouse_symbol_needs_no_source_named_and_json_matches_the_api(
        self, symbol_sources: FakeFetch
    ) -> None:
        # MGI's submission publishes one current approved symbol per gene, so a spelling it
        # retired matches nothing — and an answer that did not say why would look exactly
        # like a gene that is not in the release. A set that matches every kind carries no
        # such note.
        limited = runner.invoke(
            app, ["xref", "symbols", _MOUSE, "--source", ALLIANCE_BGI, _MOUSE_RETIRED_SYMBOL]
        )
        assert limited.exit_code == 0
        assert limited.stdout.splitlines() == [f"{_MOUSE_RETIRED_SYMBOL}\t\t\t"]
        assert f"on {APPROVED} spellings" in limited.stderr
        for missing in (PREVIOUS, ALIAS):
            assert missing in limited.stderr

        full = _match_symbols(_APPROVED_SYMBOL)
        assert f"on {APPROVED}, {PREVIOUS}, {ALIAS} spellings" in full.stderr
        assert "limits" not in full.stderr

        # The papercut this command was landed with: human's default xref source is the
        # Alliance, which carries no human symbol, so this exited 1 on the first try. The
        # question now picks the default (ADR-0021) and the retired spelling resolves.
        human = runner.invoke(app, ["xref", "symbols", _XREF_SPECIES, _RETIRED_SYMBOL])
        assert human.exit_code == 0
        assert human.stdout.splitlines() == [
            f"{_RETIRED_SYMBOL}\t{_RETIRED_SYMBOL}\t{_BMAL1_STEM}\t{PREVIOUS}"
        ]
        assert f"{HGNC_ARCHIVE} {_HGNC_RELEASE}" in human.stderr

        # Mouse's symbols come from a third source again, so the flag holds for it or it
        # is a human special case wearing a general name.
        mouse = runner.invoke(app, ["xref", "symbols", _MOUSE, _MOUSE_APPROVED_SYMBOL])
        assert mouse.exit_code == 0
        assert mouse.stdout.splitlines() == [
            f"{_MOUSE_APPROVED_SYMBOL}\t{_MOUSE_APPROVED_SYMBOL}\t{_TRP53_STEM}\t{APPROVED}"
        ]
        assert ALLIANCE_BGI in mouse.stderr

        # The strongest form of "the CLI is a thin client" for this path: the command and
        # the Python call fill the default in one place, so they cannot resolve differently.
        json_asked = [_RETIRED_SYMBOL, "nosuchgene"]
        json_result = runner.invoke(app, ["xref", "symbols", _XREF_SPECIES, *json_asked, "--json"])
        assert json_result.exit_code == 0
        assert (
            _json.loads(json_result.stdout)
            == XrefSet.for_symbols(_XREF_SPECIES).match_symbols(json_asked).as_json()
        )

    def test_a_source_carrying_no_symbols_named_on_purpose_exits_one_naming_what_to_pass(
        self, xref_release: FakeFetch
    ) -> None:
        # The Alliance's cross-reference file carries no human symbol at all, measured on
        # the whole 25 MB file — so asking it is a different failure from a gene that is
        # absent. Naming it is deliberate and is never overridden, so this still exits 1,
        # and the message names this species' symbol source rather than every one there is.
        result = runner.invoke(
            app, ["xref", "symbols", _XREF_SPECIES, "--source", ALLIANCE, _APPROVED_SYMBOL]
        )

        assert result.exit_code == 1
        assert result.stdout == ""
        assert HGNC_ARCHIVE in output(result)
        assert ALLIANCE_BGI not in output(result)

    def test_an_unsupported_species_or_an_unlisted_source_exits_one_naming_the_alternatives(
        self, symbol_sources: FakeFetch
    ) -> None:
        bad_species = runner.invoke(app, ["xref", "symbols", "Danio rerio", _APPROVED_SYMBOL])
        assert bad_species.exit_code == 1
        assert bad_species.stdout == ""
        for species in xref_species():
            assert species in output(bad_species)
        assert symbol_sources.calls == []

        bad_source = runner.invoke(
            app, ["xref", "symbols", _XREF_SPECIES, "--source", "ncbi", _APPROVED_SYMBOL]
        )
        assert bad_source.exit_code == 1
        assert bad_source.stdout == ""
        assert HGNC_ARCHIVE in output(bad_source)

    def test_a_set_not_downloaded_or_left_unfinished_exits_one_naming_the_fix(
        self, monkeypatch: pytest.MonkeyPatch, symbol_pinned: dict[str, Path]
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        not_downloaded = _match_symbols(_APPROVED_SYMBOL)
        assert not_downloaded.exit_code == 1
        assert not_downloaded.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, HGNC_ARCHIVE, _HGNC_RELEASE) in output(
            not_downloaded
        )
        assert "login node" in output(not_downloaded)

        # A leftover half-written slice is caught before the (still-broken) fetch ever
        # runs: the completion check runs first, so the two never race.
        directory = xref_set_dir(_XREF_SPECIES, HGNC_ARCHIVE, _HGNC_RELEASE)
        directory.mkdir(parents=True)
        (directory / xref_slice_name(_XREF_SPECIES)).write_bytes(b"half a file")
        unfinished = _match_symbols(_APPROVED_SYMBOL)
        assert unfinished.exit_code == 1
        assert unfinished.stdout == ""
        assert xref_prepare_command(_XREF_SPECIES, HGNC_ARCHIVE, _HGNC_RELEASE) in output(
            unfinished
        )

    def test_the_labelling_direction_needs_no_source_and_is_the_call_the_api_makes(
        self, symbol_sources: FakeFetch
    ) -> None:
        # The pair reads coherently or neither does: away from the hub a stem takes the
        # authority's one current approved spelling, which is a two-column hop like every
        # other and stays on `genome xref ids`.
        named = runner.invoke(
            app,
            [
                "xref",
                "ids",
                _XREF_SPECIES,
                "--source",
                HGNC_ARCHIVE,
                "--from-stems",
                SYMBOL,
                _BMAL1_STEM,
            ],
        )
        assert named.exit_code == 0
        assert named.stdout.splitlines() == [f"{_BMAL1_STEM}\t{_APPROVED_SYMBOL}"]

        # The same papercut in the other direction: labelling a stem is a symbol question
        # too, so the source the question fills in is the same one, and the rule reads as
        # a rule rather than as one command's special case.
        defaulted = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--from-stems", SYMBOL, _BMAL1_STEM]
        )
        assert defaulted.exit_code == 0
        assert defaulted.stdout.splitlines() == [f"{_BMAL1_STEM}\t{_APPROVED_SYMBOL}"]
        assert f"{HGNC_ARCHIVE} {_HGNC_RELEASE}" in defaulted.stderr

        # The command holds no logic the API does not: it hands the namespace it was
        # given to the constructor that fills the question's default in, and renders. It
        # once chose the opener itself, which made this exact Python call raise where the
        # shell answered — two code paths for one question, which is what this asserts is
        # gone. The mouse stem is asked of the human set on purpose: it resolves to
        # nothing, so the two renderings have to agree about what was missed too.
        asked = [_BMAL1_STEM, _TRP53_STEM]
        api_check = runner.invoke(
            app, ["xref", "ids", _XREF_SPECIES, "--from-stems", SYMBOL, *asked, "--json"]
        )
        assert api_check.exit_code == 0
        assert (
            _json.loads(api_check.stdout)
            == XrefSet.for_namespace(_XREF_SPECIES, SYMBOL).from_stems(asked, SYMBOL).as_json()
        )

    def test_the_progress_display_is_suppressed_under_json(self, symbol_sources: FakeFetch) -> None:
        _match_symbols(_APPROVED_SYMBOL, "--json")

        assert symbol_sources.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(self, symbol_sources: FakeFetch) -> None:
        _match_symbols(_APPROVED_SYMBOL)

        assert symbol_sources.last.progressbar is True
