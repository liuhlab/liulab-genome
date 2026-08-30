"""Tests for genome.io.results — what a registration answers with.

Nothing here registers anything: these are the values registration *returns*, so they are
built by hand from the fields a record carries. The one exception writes a real
annotation, because reading back a record an older version wrote is a claim about a file
on disk and not about a dataclass.

Two things are pinned rather than merely exercised. The ``as_json`` key order, because
``--json`` is a published surface and a reordered key is a break nobody would notice; and
the module's own imports, because these types are reached from both halves of the ``io``
package and an import back into either would close the cycle #100 opened this seam to
break.
"""

from __future__ import annotations

import json
from pathlib import Path

from genome.io import results as results_module
from genome.io.completion import CompletionRecord, read_record, record_path
from genome.io.gtf import AnnotationRegistry, annotation_dir
from genome.io.results import (
    EXPECTED_FROM_RECORD,
    EXPECTED_FROM_TABLE,
    AnnotationStatus,
    AnnotationStatusRow,
    RegisteredAnnotation,
    RegisteredAssembly,
    VerifiedAssembly,
    annotation_register_command,
    chromosome_check_summary,
)

from .test_source import _module_level_imports

_NAME = "ensgene_v101"


def _record(kind: str, name: str, **details: object) -> CompletionRecord:
    """A completion record with everything filled in, so ``as_json`` has every key."""
    return CompletionRecord(
        kind=kind,
        name=name,
        files={f"{name}.db": 34, f"{name}.gtf": 12},
        source_url="https://example.org/x.gz",
        sha256="1a2b3c",
        tool_versions={"samtools": "1.21"},
        package_version="2026.8.0",
        completed_at="2026-08-12T09:00:00+00:00",
        details=dict(details),
    )


def _row(**overrides: object) -> AnnotationStatusRow:
    """An offered-but-not-registered row, with any field overridden by keyword."""
    fields: dict[str, object] = {
        "name": "gencode_v50",
        "offered": True,
        "registered": False,
        "broken": False,
        "default": True,
        "provider": "GENCODE",
        "version": "v50",
        "url": "https://example.org/gencode_v50.gtf.gz",
        "sha256": None,
        "path": None,
        "problem": None,
        "repair": None,
    }
    fields.update(overrides)
    return AnnotationStatusRow(**fields)  # type: ignore[arg-type]


def _status(default: str | None, *rows: AnnotationStatusRow) -> AnnotationStatus:
    return AnnotationStatus(
        assembly="hg38",
        directory=Path("/data/genome/hg38"),
        default_annotation=default,
        annotations=rows,
    )


class TestTheStateARowIsIn:
    """``AnnotationStatusRow.state`` — the four-way state, in the words a surface prints.

    It was the CLI's to derive, which meant the precedence between ``broken`` and
    ``registered`` was written twice: once as the row's invariant, once as an ``if``
    ordering in a renderer. It is the row's, and the renderer chooses column widths.
    """

    def test_the_four_states_and_brokens_precedence_over_both(self) -> None:
        assert _row(registered=True).state == "registered"
        assert _row().state == "offered, not registered"
        assert _row(offered=False, registered=True, default=False).state == (
            "registered, not offered"
        )
        # Broken beats registered rather than reading as neither: a directory nothing
        # vouches for is not registered, so reporting the absence of a registration would
        # be true and useless. What needs acting on is that it is broken.
        assert _row(broken=True).state == "broken"
        assert _row(offered=False, broken=True).state == "broken"  # and beats unlisted too


class TestTheDefaultAnnotationLine:
    """``AnnotationStatus.default_summary`` — the closing line, and what to do about it.

    Four answers, and the two that name a command take it off an interface: the broken
    one from the row's own ``repair``, the absent one from
    :func:`annotation_register_command`. Neither is concatenated here or in a surface.
    """

    def test_the_four_summary_lines(self) -> None:
        assert _status(None).default_summary == "default: (none)"
        assert (
            _status("gencode_v50", _row(registered=True)).default_summary == "default: gencode_v50"
        )

        absent = _status("gencode_v50", _row())
        assert absent.default_summary == (
            "default: gencode_v50 — not registered here; register it with "
            "`genome register-annotation hg38 gencode_v50`"
        )
        # The command it names is the one the package spells once, not a copy of it.
        assert annotation_register_command("hg38", "gencode_v50") in absent.default_summary

        repair = "genome register-annotation hg38 gencode_v50 --force"
        broken = _status("gencode_v50", _row(broken=True, repair=repair))
        assert broken.default_summary == (
            f"default: gencode_v50 — broken here; repair it with `{repair}`"
        )

        # Named by the table, and no row here is about it: there is no row to read
        # `registered` or `broken` off, and the answer is the one that registers it.
        fresh = _status("gencode_v50")
        assert fresh.default_row is None
        assert "not registered here" in fresh.default_summary


class TestTheJsonKeysAndTheirOrder:
    """``as_json`` — every ``--json`` surface, pinned key for key and in order.

    ``--json`` is what a script parses, so a key renamed, dropped or reordered is a break
    whether or not anything in this suite notices. These assert the whole list rather than
    a key inside it, which is the only form that fails on an addition.
    """

    def test_registered_types_and_verified_assembly_pin_key_order_and_expected_from(self) -> None:
        # A registered assembly and annotation are deliberately identical shapes: both are
        # a record plus the two facts a record does not hold about itself.
        record_keys = [
            "kind",
            "name",
            "files",
            "source_url",
            "sha256",
            "tool_versions",
            "package_version",
            "completed_at",
            "details",
            "assembly",
            "directory",
        ]
        registered = RegisteredAssembly(
            assembly="hg38", directory=Path("/data/genome/hg38"), record=_record("genome", "hg38")
        )
        assert list(registered.as_json()) == record_keys
        assert registered.as_json()["directory"] == "/data/genome/hg38"

        annotation = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_record("annotation", "gencode_v50", chromosomes_checked=True),
        )
        assert list(annotation.as_json()) == record_keys

        checked = VerifiedAssembly(
            assembly="sacCer3",
            fasta=Path("/data/genome/sacCer3/sacCer3.fa"),
            sha256="6ff72f07",
            expected="6ff72f07",
            expected_from=EXPECTED_FROM_TABLE,
            components=None,
        )
        assert list(checked.as_json()) == [
            "assembly",
            "fasta",
            "sha256",
            "expected",
            "expected_from",
            "verified",
            "components",
        ]
        assert checked.as_json()["verified"] is True
        assert checked.as_json()["fasta"] == "/data/genome/sacCer3/sacCer3.fa"

        # `expected_from` is serialized as the constant the CLI keys on.
        def _from(expected_from: str | None) -> object:
            return VerifiedAssembly(
                assembly="sacCer3",
                fasta=Path("/tmp/x.fa"),
                sha256="6ff72f07",
                expected=None if expected_from is None else "6ff72f07",
                expected_from=expected_from,
                components=None,
            ).as_json()["expected_from"]

        assert _from(EXPECTED_FROM_TABLE) == "table"
        assert _from(EXPECTED_FROM_RECORD) == "record"
        assert _from(None) is None

    def test_a_status_row_and_the_status_around_it_serialize_with_no_derived_keys(self) -> None:
        assert list(_row().as_json()) == [
            "name",
            "offered",
            "registered",
            "broken",
            "default",
            "provider",
            "version",
            "url",
            "sha256",
            "path",
            "problem",
            "repair",
        ]
        # Writing `state` out would be a second spelling of the precedence for a parser to
        # disagree with, and the three fields it comes from are all here already.
        assert "state" not in _row(broken=True).as_json()

        status = _status("gencode_v50", _row(), _row(name="mine", offered=False, registered=True))
        assert list(status.as_json()) == [
            "assembly",
            "directory",
            "default_annotation",
            "annotations",
        ]
        assert [row["name"] for row in status.as_json()["annotations"]] == ["gencode_v50", "mine"]
        assert status.as_json()["directory"] == "/data/genome/hg38"

    def test_the_whole_report_survives_json_dumps_unchanged(self) -> None:
        # The last step the CLI takes, taken here: nothing in a report is a type json
        # cannot write, which is what makes the `--json` path total.
        status = _status("gencode_v50", _row(broken=True, repair="x", problem="y"))

        assert json.loads(json.dumps(status.as_json())) == status.as_json()


class TestReadingBackWhatWasChecked:
    """``chromosome_check_summary`` — one sentence per state, and never the wrong one.

    The states differ in what a reader should do about them, which is why they are told
    apart at all: an annotation registered before its assembly is waiting for the
    assembly, and one whose check the caller stood down is waiting for nothing.
    """

    _ADVICE = "register the assembly first"

    def test_each_check_state_reads_as_its_own_sentence_and_a_registration_reads_off_its_record(
        self,
    ) -> None:
        # Silence is not how a pass is reported: a surface printing nothing about the
        # check reads exactly like one printing that it passed.
        checked = chromosome_check_summary(
            {"chromosomes_checked": True, "chromosomes_unchecked_because": None}
        )
        assert "chromosomes checked" in checked
        assert self._ADVICE not in checked

        no_sizes = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "no-chrom-sizes"}
        )
        assert "chromosomes not checked" in no_sizes
        assert self._ADVICE in no_sizes

        # An override is never told to register the assembly: it may well be registered,
        # and the caller turned the check off on purpose.
        overridden = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "caller-override"}
        )
        assert "stood down" in overridden
        assert self._ADVICE not in overridden

        # A registration answers off its own record's details — the surface never spells
        # the two `details` keys itself.
        registered = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_record(
                "annotation",
                "gencode_v50",
                chromosomes_checked=False,
                chromosomes_unchecked_because="caller-override",
            ),
        )
        assert registered.chromosome_check == chromosome_check_summary(registered.record.details)
        assert "stood down" in registered.chromosome_check

    def test_every_known_state_is_distinct_and_an_unknown_reason_reads_as_unknown(self) -> None:
        summaries = {
            chromosome_check_summary(details)
            for details in (
                {"chromosomes_checked": True, "chromosomes_unchecked_because": None},
                {"chromosomes_checked": False, "chromosomes_unchecked_because": "no-chrom-sizes"},
                {"chromosomes_checked": False, "chromosomes_unchecked_because": "caller-override"},
                {"chromosomes_checked": False},
            )
        }
        assert len(summaries) == 4

        # Forward as well as backward: a record from a later version claiming some third
        # reason is one this version cannot report, which is the same as not knowing.
        future = chromosome_check_summary(
            {"chromosomes_checked": False, "chromosomes_unchecked_because": "some-later-reason"}
        )
        assert "does not say why" in future

    def test_a_record_written_before_the_reason_existed_reads_as_unknown(
        self, tmp_path: Path, data_dir: Path
    ) -> None:
        # The real back-compatibility case, on a record that is on disk: an older version
        # wrote the bare bool, and which of the two reasons it stood for is not knowable.
        # It must read as neither, and reading it must not raise.
        AnnotationRegistry.locate("tiny", tmp_path).register_path(data_dir / "tiny.gtf", _NAME)
        path = record_path(annotation_dir(tmp_path, _NAME))
        written = json.loads(path.read_text())
        written["details"] = {"chromosomes_checked": False}
        path.write_text(json.dumps(written))

        record = read_record(annotation_dir(tmp_path, _NAME))

        assert record is not None
        assert record.details == {"chromosomes_checked": False}
        summary = chromosome_check_summary(record.details)
        assert "does not say why" in summary
        assert self._ADVICE not in summary  # nor is it claimed to be the override
        assert "stood down" not in summary


class TestARecordIsCarriedWholeRatherThanCopiedOut:
    """The properties that exist so a surface never re-reads a directory."""

    def test_a_record_is_carried_whole_and_not_copied_out(self) -> None:
        registered = RegisteredAssembly(
            assembly="hg38", directory=Path("/data/genome/hg38"), record=_record("genome", "hg38")
        )
        assert registered.file_names == ["hg38.db", "hg38.gtf"]
        first = registered.file_names
        first.append("intruder")
        assert registered.file_names == ["hg38.db", "hg38.gtf"]  # a fresh list each call
        assert registered.chimera is None  # no build merged this one

        annotation = RegisteredAnnotation(
            assembly="hg38",
            directory=Path("/data/genome/hg38/gtf/gencode_v50"),
            record=_record("annotation", "gencode_v50"),
        )
        assert annotation.name == "gencode_v50"
        assert annotation.source_url == "https://example.org/x.gz"
        assert annotation.sha256 == "1a2b3c"


def test_an_answer_imports_nothing_that_produces_one() -> None:
    # The seam, in the one direction it runs. Both halves of the `io` package import
    # these types — `gtf` for the three an annotation answers with, `download` for the
    # two an assembly does — so an import back into either would make this module part of
    # the cycle #100 closed rather than the leaf both can reach.
    forbidden = {"genome.io.gtf", "genome.io.download", "genome.io.chimera", "genome.genome"}

    assert _module_level_imports(results_module) & forbidden == set()
