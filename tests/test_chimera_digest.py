"""Tests for the digests a chimera is held to: the golden one, and the recorded ones.

Two jobs, both about bytes nobody downloaded. The **golden digest** pins what this
package's concatenation produces, which the metadata table deliberately does not: the
table pins what was downloaded, and a test pins what was derived, so a change in our own
concatenation code fails here rather than on every user's disk. The rest are the
**record-against-record** comparisons that catch a component registered again underneath
a chimera — the one failure no digest of the chimera's own bytes can show, since those
bytes do not change when the component they were copied from does.

These are the tests the shared ``chimera_component`` fixture had to be right for: the
comparison resolves a component by name beside the chimera, so it only has anything to
compare when the components sit where the layout puts them. That is now what the shared
fixture does, and the near-copy of it this module used to carry is gone. Nothing here
reaches the network.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genome import Genome
from genome.io.completion import (
    RegistrationMismatchError,
    read_record,
    record_path,
    write_record,
)
from genome.io.components import COMPONENTS_UNKNOWN, components_status
from genome.io.download import register_assembly, verify_assembly
from genome.io.utils import sha256_file

from .conftest import CHIMERA_COMPONENTS, CHIMERA_EVERYDAY, COMPONENT_ANNOTATION, ComponentFactory

#: The sha256 of the FASTA the everyday three concatenate to. Committed here rather than
#: pinned in the shipped table on purpose (ADR-0008) — see the test that asserts it.
_EVERYDAY_FASTA_SHA256 = "1b04166c34b8cf3e080cfa588fb0c9e7d9fa4c191ec6f8f2b6e155aefe00d964"

#: The pair most of these use — the smallest chimera there is, and the one whose
#: components collide, so the names being right is not incidental to the digests.
_PAIR = "tinyCe_tinySc"


def _other_fasta(destination: Path) -> Path:
    """Write a small valid FASTA carrying ``tinySc``'s chromosome names and other bases.

    What a corrected component looks like: the same sequences by name, not by content, so
    a chimera built from the old ones holds a copy of bytes that are no longer anywhere.
    """
    destination.write_text(">I\nACGTACGTAC\n>II\nACGTACGTAC\n>III\nACGTACGTAC\n")
    return destination


def _other_gtf(source: Path, destination: Path) -> Path:
    """Write ``source`` back one data line short — a corrected annotation, same chromosomes."""
    lines = source.read_text().splitlines()
    destination.write_text("\n".join(lines[:-1]) + "\n")
    return destination


def _blank_component_digest(directory: Path, component: str) -> None:
    """Rewrite the record in ``directory`` with ``component``'s digest unrecorded."""
    record = read_record(directory)
    assert record is not None
    for entry in record.details["components"]:
        if entry["name"] == component:
            entry["sha256"] = None
    write_record(directory, record)


# --------------------------------------------------------------------------------------
# The golden digest: what the table does not pin, a test does
# --------------------------------------------------------------------------------------


def test_the_everyday_chimera_hashes_to_the_digest_this_concatenation_pins(
    chimera_component: ComponentFactory,
) -> None:
    # The everyday three, and not the escalating fourth: these are the components every
    # other build test uses, they derive the `__` separator every shipped assembly
    # derives, and between them they disagree about wrap width (60 against 80) and about
    # soft-masking — the two things "bytes copied verbatim" is a claim about. Adding
    # tinyEcDub would pin the escalated-separator path instead of the ordinary one.
    chimera = Genome.chimera(*(chimera_component(name) for name in CHIMERA_EVERYDAY))
    try:
        digest = sha256_file(chimera.fasta_path)
    finally:
        chimera.close()

    assert digest == _EVERYDAY_FASTA_SHA256, (
        "the everyday chimera's FASTA no longer hashes to the digest committed here. The "
        "components are committed bytes and cannot have changed, so what changed is how "
        "this package concatenates them — a header token, a line ending, the order the "
        "components are written in, or the separator between a chromosome and its "
        "component. That is a change to what every chimera on every disk is, which is "
        "why it is caught here rather than by a pin in the shipped table. If the new "
        "rule is the intended one, re-pin this digest in the commit that changes it."
    )


# --------------------------------------------------------------------------------------
# A component registered again underneath a chimera
# --------------------------------------------------------------------------------------


def test_a_component_registered_again_underneath_a_chimera_refuses_to_reopen(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # The failure the ticket calls the one no digest can see: the chimera's own bytes are
    # untouched and still agree with its record, and its tinySc sequences are a copy of
    # bytes that are no longer anywhere.
    worm, yeast = chimera_component("tinyCe"), chimera_component("tinySc")
    chimera = Genome.chimera(worm, yeast)
    written = chimera.fasta_path.stat().st_mtime_ns
    chimera.close()
    register_assembly(
        "tinySc", source=_other_fasta(tmp_path / "corrected.fa"), force=True, progressbar=False
    )

    with pytest.raises(RegistrationMismatchError, match=f"genome register {_PAIR} --force"):
        Genome.chimera(worm, chimera_component("tinySc"))

    # Nothing was rewritten on the way to refusing, and the stale bytes are still there
    # to be looked at rather than replaced behind the caller's back.
    assert chimera.fasta_path.stat().st_mtime_ns == written


def test_a_component_registered_again_underneath_a_chimera_fails_verification(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    worm, yeast = chimera_component("tinyCe"), chimera_component("tinySc")
    chimera = Genome.chimera(worm, yeast)
    chimera.close()
    register_assembly(
        "tinySc", source=_other_fasta(tmp_path / "corrected.fa"), force=True, progressbar=False
    )

    with pytest.raises(RegistrationMismatchError) as excinfo:
        verify_assembly(_PAIR)

    # The message names both digests and the command that rebuilds the chimera.
    assert f"genome register {_PAIR} --force" in str(excinfo.value)
    assert "tinySc" in str(excinfo.value)


def test_a_component_annotation_registered_again_underneath_a_chimera_refuses(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # One level down, and the reason the FASTA digests are not enough: not a base of the
    # chimera changes when a component's gene models are re-registered, so its own digest
    # agrees with its record while the merged annotation is of models nobody has now.
    worm = chimera_component("tinyCe", with_annotation=True)
    yeast = chimera_component("tinySc", with_annotation=True)
    chimera = Genome.chimera(worm, yeast)
    fasta_digest = sha256_file(chimera.fasta_path)
    chimera.close()
    source = CHIMERA_COMPONENTS["tinyCe"].gtf
    assert source is not None
    corrected = _other_gtf(source, tmp_path / "corrected.gtf")
    worm.annotations.register_path(corrected, COMPONENT_ANNOTATION, force=True)

    with pytest.raises(RegistrationMismatchError) as excinfo:
        Genome.chimera(worm, yeast)

    assert f"genome register {_PAIR} --force" in str(excinfo.value)
    assert COMPONENT_ANNOTATION in str(excinfo.value)
    # The sequence half really is untouched: a failure only the annotation digests see.
    assert sha256_file(chimera.fasta_path) == fasta_digest


# --------------------------------------------------------------------------------------
# An absent digest is unknown, on either side
# --------------------------------------------------------------------------------------


def test_a_component_with_no_record_of_its_own_reads_as_unknown(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # The current side unknown: nothing to compare against is not a disagreement, the
    # same reading a tool that never answered gets in `tool_versions`.
    worm, yeast = chimera_component("tinyCe"), chimera_component("tinySc")
    chimera = Genome.chimera(worm, yeast)
    chimera.close()
    register_assembly(
        "tinySc", source=_other_fasta(tmp_path / "corrected.fa"), force=True, progressbar=False
    )
    record_path(yeast.fasta_path.parent).unlink()

    assert components_status(chimera.assembly_dir) == COMPONENTS_UNKNOWN


def test_a_chimera_that_recorded_no_digest_for_a_component_reads_as_unknown(
    chimera_component: ComponentFactory, tmp_path: Path
) -> None:
    # And the recorded side: a chimera built before its component pinned anything leaves
    # that component unguarded rather than unopenable.
    worm, yeast = chimera_component("tinyCe"), chimera_component("tinySc")
    chimera = Genome.chimera(worm, yeast)
    directory = chimera.fasta_path.parent
    chimera.close()
    _blank_component_digest(directory, "tinySc")
    register_assembly(
        "tinySc", source=_other_fasta(tmp_path / "corrected.fa"), force=True, progressbar=False
    )

    reopened = Genome.chimera(worm, chimera_component("tinySc"))
    reopened.close()


def test_an_assembly_with_no_components_has_nothing_to_compare(
    chimera_component: ComponentFactory,
) -> None:
    # What makes a plain assembly pay nothing: it records no components, so the
    # comparison has nothing to iterate rather than a question to ask about it.
    yeast = chimera_component("tinySc")

    assert components_status(yeast.assembly_dir) is None

    record = read_record(yeast.fasta_path.parent)
    assert record is not None
    assert "components" not in record.details
