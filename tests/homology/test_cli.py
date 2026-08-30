"""Tests for the ``genome homology`` sub-app."""

from __future__ import annotations

import json as _json
from pathlib import Path

import pytest

from genome.cli import app
from genome.homology import (
    DEFAULT_RELEASE as HOMOLOGY_RELEASE,
)
from genome.homology import (
    QUALITY_SCORE_COLUMNS,
    HomologySet,
    homology_metadata,
    homology_prepare_command,
    homology_species,
)
from genome.store import fetch as fetch_mod
from genome.store.completion import RECORD_NAME

from .._cli import output, runner
from ..conftest import FakeFetch
from .test_homology import ABSENT as _NO_HOMOLOG
from .test_homology import FIXTURES as _COMPARA_FIXTURES
from .test_homology import ONE2MANY_HUMAN as _THREE_WORM_HOMOLOGS
from .test_homology import ONE2MANY_WORMS as _THE_THREE_WORMS
from .test_homology import ONE2ONE_HUMAN as _ONE_WORM_HOMOLOG
from .test_homology import ONE2ONE_WORM as _THE_ONE_WORM
from .test_homology import PAIRS as _HOMOLOGY_PAIRS
from .test_homology import _stems as _homology_stems

#: The three species ``genome homology links`` is exercised against, spelled as the shipped
#: provenance table spells them. All three pairings among them must answer, which is an
#: acceptance criterion of its own rather than a sample.
_HUMAN, _MOUSE, _WORM = "Homo sapiens", "Mus musculus", "Caenorhabditis elegans"


def _serve_compara(fake_fetch: FakeFetch, species: str, other: str) -> FakeFetch:
    """Serve the subsample of whichever published dump the shipped table names for a pair.

    Which of the two per-species files holds a pair is exactly what that table records, so
    a test serves what the package asked for rather than deciding for itself — the test
    that serves the *wrong* file does so on purpose and says why.
    """
    row = homology_metadata(species, other, HOMOLOGY_RELEASE)
    assert row is not None
    fake_fetch.serve(_COMPARA_FIXTURES[row.holding_species])
    return fake_fetch


class TestHomologsCommand:
    """``genome homology links`` — the shell surface over a **Homology set**, on the fixtures.

    Offline throughout, the way ``tests/homology/test_homology.py`` is: the fake fetch serves the
    committed Compara subsamples and the command prepares and reads a real set under the
    test's own data root. What is asserted here is the command and not the set —
    ``tests/homology/test_homology.py`` owns the slice, the partition and the answer — so: all three
    pairings reaching a shell, the publisher's **Homology type** surviving the render, the
    stdout/stderr split that makes the output pipe, the **Dropped partner**s and the null
    quality scores being said out loud, and a non-zero exit naming the next action for each
    way it can fail.

    The strongest claim here is that the command holds no logic the API does not, and it is
    checked the two ways that can be: ``--json`` is asserted equal to what the same call
    answers in Python, whole, and the text rows are asserted to be that JSON's own values
    in its own key order — so nothing the shell sees was assembled here.
    """

    @pytest.mark.parametrize(("species", "other", "links"), _HOMOLOGY_PAIRS)
    def test_it_answers_for_every_pairing_among_human_mouse_and_worm(
        self, fake_fetch: FakeFetch, data_dir: Path, species: str, other: str, links: int
    ) -> None:
        _serve_compara(fake_fetch, species, other)
        asked = _homology_stems(data_dir, species, other)

        result = runner.invoke(app, ["homology", "links", species, other, *asked])

        assert result.exit_code == 0
        assert len(result.stdout.splitlines()) == links

    def test_the_type_column_carries_the_publishers_label_and_defaults_to_orthologs(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        result = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _THREE_WORM_HOMOLOGS])
        assert result.exit_code == 0
        rows = [line.split("\t") for line in result.stdout.splitlines()]
        assert [row[1] for row in rows] == list(_THE_THREE_WORMS)
        # Verbatim, and the label the publisher's tree assigned rather than a count of the
        # rows in front of you: three partners here and the type still reads one2many.
        assert {row[2] for row in rows} == {"ortholog_one2many"}

        _serve_compara(fake_fetch, _MOUSE, _WORM)
        default_asked = _homology_stems(data_dir, _MOUSE, _WORM)
        default = runner.invoke(app, ["homology", "links", _MOUSE, _WORM, *default_asked])
        assert default.exit_code == 0
        types = {line.split("\t")[2] for line in default.stdout.splitlines()}
        assert types
        assert all(kind.startswith("ortholog_") for kind in types)
        assert "orthologs" in default.stderr

        # *Not an ortholog* stays distinguishable from *absent* because the publisher's
        # own label is a column of every row rather than a filter applied before
        # printing: a duplication label would print in the same place `ortholog_one2one`
        # does now. Release 116 publishes no cross-species paralogy for this pair —
        # counted over the whole human dump — so what is asserted is that the flag
        # reaches the API and the render says which question was asked, not that a
        # paralogy row appeared: claiming one would claim something about the publisher
        # that is not true.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        asked = _homology_stems(data_dir, _HUMAN, _WORM)
        paralogs_included = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, *asked, "--paralogs"]
        )
        assert paralogs_included.exit_code == 0
        assert {line.split("\t")[2] for line in paralogs_included.stdout.splitlines()} == {
            "ortholog_one2one",
            "ortholog_one2many",
            "ortholog_many2many",
        }

        result = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, *asked, "--paralogs", "--json"]
        )
        heading = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, *asked, "--paralogs"])
        assert result.exit_code == heading.exit_code == 0
        assert (
            _json.loads(result.stdout)
            == HomologySet(_HUMAN, _WORM, progressbar=False)
            .homologs(asked, paralogs=True)
            .as_json()
        )
        assert "paralogy included" in heading.stderr

    def test_json_matches_the_api_text_matches_the_json_and_nothing_is_ever_dropped(
        self, fake_fetch: FakeFetch, data_dir: Path
    ) -> None:
        _serve_compara(fake_fetch, _MOUSE, _HUMAN)
        asked = [*_homology_stems(data_dir, _MOUSE, _HUMAN), _NO_HOMOLOG]
        result = runner.invoke(app, ["homology", "links", _MOUSE, _HUMAN, *asked, "--json"])
        assert result.exit_code == 0
        assert (
            _json.loads(result.stdout)
            == HomologySet(_MOUSE, _HUMAN, progressbar=False).homologs(asked).as_json()
        )

        # The whole claim that the command holds no logic: every cell printed is a value
        # the API put in the answer, in the order the API writes it, with the
        # publisher's own `NULL` where it recorded nothing.
        text_asked = "ENSMUSG00000074698"
        text = runner.invoke(app, ["homology", "links", _MOUSE, _HUMAN, text_asked])
        rendered = runner.invoke(app, ["homology", "links", _MOUSE, _HUMAN, text_asked, "--json"])
        links = _json.loads(rendered.stdout)["resolved"][text_asked]
        assert text.stdout.splitlines() == [
            "\t".join("NULL" if value is None else str(value) for value in link.values())
            for link in links
        ]

        _serve_compara(fake_fetch, _HUMAN, _WORM)
        row = homology_metadata(_HUMAN, _WORM, HOMOLOGY_RELEASE)
        assert row is not None
        single = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])
        assert single.exit_code == 0
        assert single.stdout.startswith(f"{_ONE_WORM_HOMOLOG}\t{_THE_ONE_WORM}\tortholog_one2one\t")
        assert len(single.stdout.splitlines()) == 1
        # Who asserted it, and from which release and file, is what a reader needs and
        # what a pipeline must not be handed: it goes beside the links rather than among
        # them.
        assert row.attribution() in single.stderr
        assert f"{_HUMAN} -> {_WORM}" in single.stderr

        # The one thing a hand-rolled join drops. A stem with no link gets a row of its
        # own with every other column empty, so nothing leaves shorter than it arrived —
        # and empty is not `NULL`, which would claim a link the publisher scored nothing
        # on — and it stays visible in both renderings rather than being dropped.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        result = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, _NO_HOMOLOG]
        )
        assert result.exit_code == 0
        rows = result.stdout.splitlines()
        assert len(rows) == 2
        assert rows[1].split("\t")[0] == _NO_HOMOLOG
        assert set(rows[1].split("\t")[1:]) == {""}
        assert "1 this release names no homolog for" in result.stderr

        json_result = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, _NO_HOMOLOG, "--json"]
        )
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload["unresolved"] == [_NO_HOMOLOG]
        assert _NO_HOMOLOG not in payload["resolved"]

        # Counted *and* named, both off the answer. None are dropped on release 116 — it
        # publishes no cross-species paralogy for the ortholog filter to remove — and
        # zero is printed as an answer rather than left as a silence.
        asked = _homology_stems(data_dir, _HUMAN, _WORM)
        dropped = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, *asked])
        dropped_payload = _json.loads(
            runner.invoke(app, ["homology", "links", _HUMAN, _WORM, *asked, "--json"]).stdout
        )
        assert dropped.exit_code == 0
        assert "dropped partners" in dropped.stderr
        assert dropped_payload["dropped_partners"] == []

        # Compara records neither score on any link of *either* worm pairing, so a shell
        # user about to write `awk -F'\t' '$6 > 50'` is told the column is null
        # throughout rather than discovering it when the filter comes back empty.
        null_result = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])
        null_payload = _json.loads(
            runner.invoke(
                app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--json"]
            ).stdout
        )
        assert null_result.exit_code == 0
        # The line that says it, not merely the word appearing somewhere: the column
        # list above it names every column, so a warning nobody wrote would pass a
        # looser check.
        (null_quality,) = [
            line for line in null_result.stderr.splitlines() if line.strip().startswith("quality")
        ]
        for column in QUALITY_SCORE_COLUMNS:
            assert column in null_quality
        assert "empties" in null_quality
        assert null_payload["null_quality_scores"] == list(QUALITY_SCORE_COLUMNS)

        # A pair Compara did score says that too: silence would read the same as a
        # warning nobody printed.
        _serve_compara(fake_fetch, _MOUSE, _HUMAN)
        scored = runner.invoke(app, ["homology", "links", _MOUSE, _HUMAN, "ENSMUSG00000074698"])
        assert scored.exit_code == 0
        (scored_quality,) = [
            line for line in scored.stderr.splitlines() if line.strip().startswith("quality")
        ]
        assert "carry values" in scored_quality
        assert "empties" not in scored_quality

    def test_every_way_homologs_can_be_refused_exits_one_naming_the_fix(
        self, monkeypatch: pytest.MonkeyPatch, fake_fetch: FakeFetch
    ) -> None:
        unsupported = runner.invoke(
            app, ["homology", "links", "Danio rerio", _HUMAN, _ONE_WORM_HOMOLOG]
        )
        assert unsupported.exit_code == 1
        assert unsupported.stdout == ""
        for species in homology_species():
            assert species in output(unsupported)
        # Refused before anything was fetched: nobody pinned this species, which must
        # never read as this species having no homologs.
        assert fake_fetch.calls == []

        same_species = runner.invoke(
            app, ["homology", "links", _HUMAN, "homo_sapiens", _ONE_WORM_HOMOLOG]
        )
        assert same_species.exit_code == 1
        assert same_species.stdout == ""
        assert "two different species" in output(same_species)
        assert fake_fetch.calls == []

        bad_release = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--release", "115"]
        )
        assert bad_release.exit_code == 1
        assert bad_release.stdout == ""
        assert HOMOLOGY_RELEASE in output(bad_release)
        assert fake_fetch.calls == []

        # Compara writes its ids bare, so the versioned spelling would match nothing and
        # come back unresolved looking exactly like a gene it never placed in a tree.
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        versioned = runner.invoke(
            app, ["homology", "links", _HUMAN, _WORM, f"{_ONE_WORM_HOMOLOG}.18"]
        )
        assert versioned.exit_code == 1
        assert versioned.stdout == ""
        assert _ONE_WORM_HOMOLOG in output(versioned)

        # The published partition, not a staged one: the real release-116 human dump
        # holds zero human/mouse rows. Serving it for that pair is what a release that
        # had re-partitioned looks like from a shell, and it must be an error naming the
        # other file rather than an empty answer that reads as *these species share no
        # homologs*.
        fake_fetch.serve(_COMPARA_FIXTURES[_HUMAN])
        wrong_file = runner.invoke(app, ["homology", "links", _HUMAN, _MOUSE, "ENSG00000172150"])
        assert wrong_file.exit_code == 1
        assert wrong_file.stdout == ""
        assert "homo_sapiens/Compara.116.protein_default.homologies.tsv.gz" in output(wrong_file)
        assert "homology_metadata.tsv" in output(wrong_file)

        fake_fetch.serve("tiny.gtf.gz")
        not_comparas = runner.invoke(app, ["homology", "links", _HUMAN, _MOUSE, "ENSG00000172150"])
        assert not_comparas.exit_code == 1
        assert not_comparas.stdout == ""
        assert "gene_stable_id" in output(not_comparas)

        # A set that cannot be downloaded, and one left unfinished by a run killed
        # mid-way, are the two remaining ways this command refuses.
        real_fetch = fetch_mod.fetch_url

        def no_internet(url: str, dest_dir: Path, **kwargs: object) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        not_downloaded = runner.invoke(
            app, ["homology", "links", _HUMAN, _MOUSE, "ENSG00000172150"]
        )
        assert not_downloaded.exit_code == 1
        assert not_downloaded.stdout == ""
        assert homology_prepare_command(_HUMAN, _MOUSE, HOMOLOGY_RELEASE) in output(not_downloaded)
        assert "login node" in output(not_downloaded)

        # Restore the real (fake-serving) fetch to prepare a set, then delete only its
        # completion marker, so this is what a run killed mid-way leaves behind.
        monkeypatch.setattr(fetch_mod, "fetch_url", real_fetch)
        _serve_compara(fake_fetch, _HUMAN, _WORM)
        prepared = HomologySet(_HUMAN, _WORM, progressbar=False)
        (prepared.path.parent / RECORD_NAME).unlink()
        unfinished = runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])
        assert unfinished.exit_code == 1
        assert unfinished.stdout == ""
        assert "rm -rf" in output(unfinished)

    def test_the_progress_display_is_suppressed_under_json(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG, "--json"])

        assert fake_fetch.last.progressbar is False

    def test_the_progress_display_is_drawn_without_it(self, fake_fetch: FakeFetch) -> None:
        _serve_compara(fake_fetch, _HUMAN, _WORM)

        runner.invoke(app, ["homology", "links", _HUMAN, _WORM, _ONE_WORM_HOMOLOG])

        assert fake_fetch.last.progressbar is True
