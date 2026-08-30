"""Tests for genome.io.components — what a chimera's record says it is made of.

Every record here is written by hand, so the whole module is exercised without a native
tool, a network call, or a real concatenation: the questions it answers are answered from
a handful of small JSON files, which is exactly what it reads in production.

``test_chimera_digest`` makes the same record-against-record comparison end to end,
through a real build and behind a tool skip. These place the records directly instead, so
each answer a caller branches on — unchanged, unknown, and the refusal between them — is
asserted on its own rather than inferred from a build that happened to reach it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from genome.io.completion import (
    RegistrationError,
    RegistrationMismatchError,
    build_record,
    read_record,
    record_path,
    write_record,
)
from genome.io.components import (
    COMPONENTS_UNCHANGED,
    COMPONENTS_UNKNOWN,
    ChimeraDetails,
    ComponentDetails,
    components_status,
    merged_annotation_name,
    read_chimera_details,
)
from genome.io.registration import AssemblyDir

_PAIR = "tinyCe_tinySc"


def _record(directory: Path, name: str, **kwargs: Any) -> None:
    """Write a completion record in ``directory`` claiming ``name`` finished there."""
    directory.mkdir(parents=True, exist_ok=True)
    claimed = directory / f"{name}.fa"
    claimed.write_text(">chrI\nACGT\n")
    write_record(
        directory,
        build_record(directory, kind="genome", name=name, files=[claimed], **kwargs),
    )


def _chimera(root: Path, *, details: dict[str, Any] | None) -> AssemblyDir:
    """Return the **Assembly dir** of a chimera recorded under ``root`` with ``details``."""
    here = AssemblyDir.locate(_PAIR, root / "genome" / _PAIR)
    _record(here.path, _PAIR, details=details)
    return here


def _component(chimera: AssemblyDir, name: str, digest: str | None) -> AssemblyDir:
    """Record component ``name`` beside ``chimera``, pinning ``digest``."""
    here = chimera.sibling(name)
    _record(here.path, name, sha256=digest)
    return here


def _annotation(component: AssemblyDir, name: str, digest: str | None) -> Path:
    """Record annotation ``name`` under ``component``, pinning ``digest``, and return its dir."""
    directory = component.annotation_dir(name)
    directory.mkdir(parents=True, exist_ok=True)
    gtf = directory / f"{name}.gtf"
    gtf.write_text("#\n")
    write_record(
        directory,
        build_record(directory, kind="annotation", name=name, files=[gtf], sha256=digest),
    )
    return directory


def _details(*entries: ComponentDetails, merged: bool = False) -> dict[str, Any]:
    """Return the ``details`` a build of ``entries`` writes down."""
    return ChimeraDetails(separator="__", component_details=entries).as_details(merged=merged)


# ---------------------------------------------------------------------------
# What a build writes down, read back
# ---------------------------------------------------------------------------


def test_the_details_a_build_writes_are_the_details_it_reads_back() -> None:
    details = ChimeraDetails(
        separator="___",
        component_details=(
            ComponentDetails("tinyCe", "1a2b", "genes", "3c4d"),
            ComponentDetails("tinySc", "5e6f", None, None),
        ),
    )

    assert ChimeraDetails.from_details(details.as_details(merged=True), assembly=_PAIR) == details


def test_absence_of_a_merge_and_absence_of_a_contributor_are_written_differently() -> None:
    # An absent key is a fact nobody gathered, which is not `null` — the reading a
    # component that contributed nothing to a merge that did happen gets.
    no_merge = ChimeraDetails(
        separator="__",
        component_details=(ComponentDetails("tinyCe", "1a2b", "genes", "3c4d"),),
    )
    written = no_merge.as_details(merged=False)
    assert written == {"separator": "__", "components": [{"name": "tinyCe", "sha256": "1a2b"}]}
    read_back = ChimeraDetails.from_details(written, assembly=_PAIR)
    assert read_back is not None
    assert read_back.component_details[0].annotation is None

    # A merge that happened but nothing contributed to is a third thing again: merged no
    # annotation rather than an empty one.
    empty_merge = ChimeraDetails(
        separator="__",
        component_details=(
            ComponentDetails("tinyCe", None, None, None),
            ComponentDetails("tinySc", None, None, None),
        ),
    )
    assert empty_merge.merged_annotation is None
    assert empty_merge.components == ["tinyCe", "tinySc"]


def test_a_merged_annotation_name_is_the_contributors_joined() -> None:
    # One spelling for the name a merge is written under and the name a record is read
    # back as, so a rebuild can never look for a database under a name it did not write.
    contributors = ["wormbase_ws298", "refseq_rs_2025_06_26"]
    details = ChimeraDetails(
        separator="__",
        component_details=tuple(
            ComponentDetails(name, None, annotation, None)
            for name, annotation in zip(("ce11", "ecHT115"), contributors, strict=True)
        ),
    )

    assert merged_annotation_name(contributors) == "wormbase_ws298+refseq_rs_2025_06_26"
    assert details.merged_annotation == merged_annotation_name(contributors)


# ---------------------------------------------------------------------------
# Not a chimera — the ordinary answer, from every angle
# ---------------------------------------------------------------------------


def test_details_or_a_record_that_say_nothing_about_components_are_not_a_chimeras(
    tmp_path: Path,
) -> None:
    assert ChimeraDetails.from_details({}, assembly="hg38") is None

    # An ordinary registration's details carry keys of their own, and none of them is a
    # half-written chimera: saying nothing about components is saying *not a chimera*.
    ordinary = {"chromosomes_checked": True, "features_inferred": []}
    assert ChimeraDetails.from_details(ordinary, assembly="hg38") is None

    assert ChimeraDetails.from_record(None) is None
    assert read_chimera_details(tmp_path) is None


# ---------------------------------------------------------------------------
# Malformed — a broken registration, which raises and names its repair
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("details", "wrong"),
    [
        pytest.param({"separator": "__"}, "components", id="separator-alone"),
        pytest.param(
            {"components": [{"name": "tinyCe", "sha256": None}]}, "separator", id="components-alone"
        ),
        pytest.param(
            {"separator": "__", "components": {"tinyCe": None}}, "components", id="components-map"
        ),
        pytest.param(
            {"separator": "__", "components": [{"sha256": "1a2b"}]},
            "components",
            id="entry-unnamed",
        ),
    ],
)
def test_a_record_that_claims_to_be_a_chimera_and_is_not_readable_raises(
    details: dict[str, Any], wrong: str
) -> None:
    # The defect this closes: `None` used to mean both *not a chimera* and *this record
    # is broken*, and `verify_assembly` and `components_status` both read the second as
    # the first — a chimera silently demoted to an ordinary assembly nobody checks.
    with pytest.raises(RegistrationError) as excinfo:
        ChimeraDetails.from_details(details, assembly=_PAIR)

    message = str(excinfo.value)
    assert wrong in message  # which half of the shape was wrong
    assert f"genome register {_PAIR} --force" in message  # and the repair


def test_a_broken_record_raises_where_a_chimeras_details_are_asked_for(tmp_path: Path) -> None:
    # Through the record, so the assembly the message names is the one the record does
    # rather than one a caller had to remember to pass.
    chimera = _chimera(tmp_path, details={"separator": "__"})

    with pytest.raises(RegistrationError, match=f"genome register {_PAIR} --force"):
        read_chimera_details(chimera.path)


# ---------------------------------------------------------------------------
# Comparing a chimera's components against their own records
# ---------------------------------------------------------------------------


def test_an_assembly_with_no_components_has_nothing_to_compare(liulab_data: Path) -> None:
    # What makes an ordinary assembly pay nothing: it records no components, so there is
    # nothing to iterate rather than a question to answer about it.
    here = AssemblyDir.locate("hg38", liulab_data / "genome" / "hg38")
    _record(here.path, "hg38")

    assert components_status(here) is None


def test_components_that_are_still_themselves_read_as_unchanged(tmp_path: Path) -> None:
    chimera = _chimera(
        tmp_path,
        details=_details(
            ComponentDetails("tinyCe", "1a2b", None, None),
            ComponentDetails("tinySc", "5e6f", None, None),
        ),
    )
    _component(chimera, "tinyCe", "1a2b")
    _component(chimera, "tinySc", "5e6f")

    assert components_status(chimera) == COMPONENTS_UNCHANGED


def test_a_component_registered_again_underneath_the_chimera_is_refused(tmp_path: Path) -> None:
    # The one failure a digest of the chimera's own bytes cannot show: those bytes are
    # untouched and still agree with their record, and the component they were copied
    # from is gone.
    chimera = _chimera(
        tmp_path,
        details=_details(
            ComponentDetails("tinyCe", "1a2b", None, None),
            ComponentDetails("tinySc", "5e6f", None, None),
        ),
    )
    _component(chimera, "tinyCe", "1a2b")
    _component(chimera, "tinySc", "differentnow")

    with pytest.raises(RegistrationMismatchError) as excinfo:
        components_status(chimera)

    message = str(excinfo.value)
    assert "5e6f" in message  # what the chimera recorded
    assert "differentnow" in message  # and what is pinned now
    assert f"genome register {_PAIR} --force" in message


def test_an_absent_current_or_recorded_digest_reads_as_unknown_rather_than_a_mismatch(
    tmp_path: Path,
) -> None:
    # Nothing to compare against is not a disagreement — the reading a tool that never
    # answered gets in `tool_versions`. Absent on the current side...
    chimera = _chimera(tmp_path, details=_details(ComponentDetails("tinyCe", "1a2b", None, None)))
    component = _component(chimera, "tinyCe", "1a2b")
    record_path(component.path).unlink()
    assert components_status(chimera) == COMPONENTS_UNKNOWN

    # ...and absent on the recorded side: unproven, not refused, either way.
    unguarded = _chimera(
        tmp_path / "other", details=_details(ComponentDetails("tinyCe", None, None, None))
    )
    _component(unguarded, "tinyCe", "1a2b")
    assert components_status(unguarded) == COMPONENTS_UNKNOWN


def test_an_annotation_registered_again_after_the_merge_is_refused(tmp_path: Path) -> None:
    # The failure one level down, which no digest of any FASTA can see: every base of the
    # chimera's sequence is what it was, and the gene models merged into it are not.
    chimera = _chimera(
        tmp_path,
        details=_details(ComponentDetails("tinyCe", "1a2b", "genes", "3c4d"), merged=True),
    )
    component = _component(chimera, "tinyCe", "1a2b")
    _annotation(component, "genes", "rebuiltsince")

    with pytest.raises(RegistrationMismatchError) as excinfo:
        components_status(chimera)

    message = str(excinfo.value)
    assert "genes" in message
    assert "3c4d" in message
    assert f"genome register {_PAIR} --force" in message


def test_an_annotation_with_no_record_leaves_the_component_unproven(tmp_path: Path) -> None:
    chimera = _chimera(
        tmp_path,
        details=_details(ComponentDetails("tinyCe", "1a2b", "genes", "3c4d"), merged=True),
    )
    _component(chimera, "tinyCe", "1a2b")

    assert components_status(chimera) == COMPONENTS_UNKNOWN


def test_the_comparison_reads_records_and_never_the_sequence(tmp_path: Path) -> None:
    # Record against record, with nothing rehashed: the chimera's FASTA is replaced with
    # bytes that hash to something else entirely and the answer does not move, because no
    # base of it is read to reach one.
    chimera = _chimera(tmp_path, details=_details(ComponentDetails("tinyCe", "1a2b", None, None)))
    component = _component(chimera, "tinyCe", "1a2b")
    (component.path / "tinyCe.fa").write_text(">chrI\nTTTTTTTT\n")

    assert components_status(chimera) == COMPONENTS_UNCHANGED
    record = read_record(component.path)
    assert record is not None
    assert record.sha256 == "1a2b"
