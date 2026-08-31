"""Tests for the root ``genome`` app: its command tree, and the three commands on it.

Every other command moved under a sub-app and its tests moved with it, into
``tests/assembly/``, ``tests/annotation/``, ``tests/tf/``, ``tests/xref/`` and
``tests/homology/``. What is asserted here is the tree itself and the three commands that
belong to no topic.
"""

from __future__ import annotations

import json as _json
import shutil
from pathlib import Path

import pytest
import typer.main
from typer.core import TyperGroup

from genome import __version__ as genome_version
from genome.cli import _DEPRECATED_ALIASES, app
from genome.external import REQUIRED_TOOLS
from genome.external import doctor as doctor_api
from genome.seq import DNA

from ._cli import help_text, output, runner

_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)

#: The whole of what ``genome --help`` offers: three commands belonging to no topic, and
#: one sub-app per topic named for the module it ships from — which is why the Orthology
#: context's is ``homology``. Spelled out rather than derived, because "and nothing else"
#: is the claim: a command left hanging off the root app fails here.
_ROOT_COMMANDS = ("version", "revcomp", "doctor")
_SUB_APPS = ("assembly", "annotation", "tf", "xref", "homology", "motif")


class TestVersion:
    def test_version_reports_the_same_string_as_text_and_json(self) -> None:
        text_result = runner.invoke(app, ["version"])
        json_result = runner.invoke(app, ["version", "--json"])

        assert text_result.exit_code == 0
        assert json_result.exit_code == 0
        payload = _json.loads(json_result.stdout)
        assert payload == {"version": genome_version}
        assert text_result.stdout.strip() == payload["version"]


class TestRevcomp:
    def test_reverses_complements_preserves_case_and_reports_json(self) -> None:
        assert runner.invoke(app, ["revcomp", "ATCG"]).stdout.strip() == "CGAT"

        cased = runner.invoke(app, ["revcomp", "aTcG"])
        assert cased.exit_code == 0
        assert cased.stdout.strip() == "CgAt"

        json_result = runner.invoke(app, ["revcomp", "ATCG", "--json"])
        assert json_result.exit_code == 0
        assert _json.loads(json_result.stdout) == {"input": "ATCG", "reverse_complement": "CGAT"}

    def test_invalid_input_exits_2_naming_the_base_and_the_alphabet(self) -> None:
        result = runner.invoke(app, ["revcomp", "ATCX"])
        assert result.exit_code == 2
        assert "error" in output(result).lower()
        assert "X" in output(result)  # the base it could not complement
        assert "{ACGT}" in output(result)

    def test_alphabet_comes_from_the_type(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # The defect this pins: the edge used to spell A/C/G/T itself, so changing what a
        # DNA may contain reached the CLI not at all. Widen the type's alphabet and the
        # command follows — both the check and the message naming it.
        monkeypatch.setattr(DNA, "ALPHABET", frozenset("ACGTN"))
        accepted = runner.invoke(app, ["revcomp", "ATCN"])
        assert accepted.exit_code == 0
        assert accepted.stdout.strip() == "NGAT"
        refused = runner.invoke(app, ["revcomp", "ATCX"])
        assert refused.exit_code == 2
        assert "{ACGNT}" in output(refused)  # rendered sorted, a set having no order


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="samtools/bedtools not on PATH")
class TestDoctor:
    def test_reports_every_tool_as_text_and_json(self) -> None:
        text = runner.invoke(app, ["doctor"])
        json_result = runner.invoke(app, ["doctor", "--json"])

        assert text.exit_code == 0
        assert json_result.exit_code == 0
        for tool in REQUIRED_TOOLS:
            assert tool in text.stdout
        assert set(_json.loads(json_result.stdout).keys()) == set(REQUIRED_TOOLS)

    def test_doctor_prints_one_line_per_tool_and_nothing_else(self) -> None:
        result = runner.invoke(app, ["doctor"])

        assert result.exit_code == 0
        assert result.stdout == "".join(f"{name}: {ver}\n" for name, ver in doctor_api().items())


class TestTheCommandTree:
    """What ``genome --help`` offers: six sub-apps, three commands, and nothing else."""

    def test_the_root_offers_the_six_sub_apps_the_three_commands_and_nothing_else(self) -> None:
        group = typer.main.get_command(app)
        assert isinstance(group, TyperGroup)

        # Every alias is hidden, so what a reader is shown is the new tree alone — and the
        # sub-apps and the three commands together are all of it.
        listed = {name for name, command in group.commands.items() if not command.hidden}
        assert listed == {*_ROOT_COMMANDS, *_SUB_APPS}

        offered = help_text()
        for name in (*_ROOT_COMMANDS, *_SUB_APPS):
            assert name in offered

    def test_each_sub_app_offers_the_commands_that_moved_under_it(self) -> None:
        assert set(help_text("assembly").split()) >= {"register", "list", "verify", "table-row"}
        assert set(help_text("annotation").split()) >= {
            "register",
            "register-gtf",
            "list",
            "gene-list",
            "gene-categories",
        }
        assert set(help_text("tf").split()) >= {"gene-list", "cofactor-list"}
        assert set(help_text("xref").split()) >= {"ids", "symbols"}
        assert set(help_text("homology").split()) >= {"links"}
        assert set(help_text("motif").split()) >= {"scan"}


class TestTheDeprecatedFlatSpellings:
    """The old spelling of every renamed command still runs, unlisted, and warns on stderr.

    A hard break would reach nobody: this package publishes to PyPI as ``liulab-genome``
    and its callers are not all reachable. So each alias is the very function object its
    sub-app registered — one implementation, two spellings — and the notice naming the
    replacement goes to stderr, which is what leaves ``--json`` on stdout parseable for a
    script that has not moved yet.
    """

    def test_every_alias_is_hidden_and_deprecated_and_shares_its_command_s_callback(
        self,
    ) -> None:
        group = typer.main.get_command(app)
        assert isinstance(group, TyperGroup)

        for old, command in _DEPRECATED_ALIASES.items():
            alias = group.commands[old]
            assert alias.hidden, old
            assert alias.deprecated, old
            # The same function object, not a copy of it: nothing can drift between the
            # two spellings because there is only one implementation to drift from.
            assert alias.callback is not None
            assert alias.callback.__wrapped__ is command  # type: ignore[attr-defined]

    def test_an_old_spelling_runs_with_the_notice_on_stderr_and_the_json_on_stdout(
        self, liulab_data: Path
    ) -> None:
        # `annotations` needs nothing prepared, downloaded or built to answer, so it is
        # the alias that can be run end to end without a fixture of any kind.
        result = runner.invoke(app, ["annotations", "hg38", "--json"])

        assert result.exit_code == 0
        assert "deprecated" in result.stderr
        assert "annotations" in result.stderr
        # The whole of stdout is still the payload — nothing of the notice leaked into it.
        assert _json.loads(result.stdout)["assembly"] == "hg38"
        assert result.stdout == runner.invoke(app, ["annotation", "list", "hg38", "--json"]).stdout

    def test_the_one_old_spelling_a_sub_app_name_shadows_is_answered_by_that_sub_app(
        self, liulab_data: Path
    ) -> None:
        # `xref` names a sub-app and named a command, and a root app holds one of each
        # under one name — so this alias cannot be in the table and the group answers it.
        # Refused at the flags, which is where the old command refused it too, so nothing
        # is fetched to prove the flat spelling arrived.
        flat = runner.invoke(app, ["xref", "Homo sapiens", "ENSG00000141510"])

        assert flat.exit_code == 2
        assert "deprecated" in flat.stderr
        assert "genome xref ids" in flat.stderr
        assert "name exactly one direction" in flat.stderr
        assert flat.stdout == ""
        assert "xref" not in _DEPRECATED_ALIASES
