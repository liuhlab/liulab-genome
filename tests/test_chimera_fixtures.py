"""The tiny component assemblies a chimera is tested against — and what each exercises.

These are **inputs**, so what is asserted here is the fixture set itself, never the
contract it feeds: the properties are read off the committed bytes with local helpers
rather than through the naming contract, so the set stays pinned independently of the code
that will consume it. A chimera build that stops honouring one of these properties fails
in its own test; a fixture that quietly stops carrying one fails here.

This is the one test module that mirrors a fixture directory rather than a module of
``src/``, because what it guards is ``tests/data/chimera/``.

Provenance is mechanised rather than promised: every component chromosome is asserted to
be the exact slice of ``tests/data/tiny.fa`` that ``tests/data/README.md`` says it is.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from genome.metadata import lookup_assembly

from .conftest import (
    CHIMERA_COMPONENTS,
    CHIMERA_ESCALATION,
    CHIMERA_EVERYDAY,
    COMPONENT_ANNOTATION,
    DATA_DIR,
    ChimeraComponent,
    ComponentFactory,
)

_ALL = list(CHIMERA_COMPONENTS.values())


def _read_fasta(path: Path) -> dict[str, str]:
    """Read a (tiny) FASTA into ``{name: sequence}``, preserving case and file order."""
    sequences: dict[str, list[str]] = {}
    current: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            current = sequences.setdefault(line[1:].split()[0], [])
        else:
            current.append(line)
    return {name: "".join(chunks) for name, chunks in sequences.items()}


def _sequence_line_widths(path: Path) -> list[int]:
    """Return the length of every sequence line, the last line of each record dropped.

    A record's final line is short whenever the sequence does not divide by the wrap
    width, so it says nothing about the wrapping.
    """
    widths: list[int] = []
    lines = path.read_text().splitlines()
    for index, line in enumerate(lines):
        last = index + 1 == len(lines) or lines[index + 1].startswith(">")
        if not line.startswith(">") and not last:
            widths.append(len(line))
    return widths


def _longest_underscore_run(name: str) -> int:
    """Return the length of the longest run of underscores in ``name`` (0 when there is none)."""
    return max((len(run) for run in re.findall(r"_+", name)), default=0)


def _chromosome_names(gtf: Path) -> set[str]:
    """Return the set of chromosome names a GTF's first column carries."""
    return {line.split("\t")[0] for line in gtf.read_text().splitlines() if line.strip()}


def _strict_prefix_pairs(names: set[str]) -> set[tuple[str, str]]:
    """Return every ``(short, long)`` where ``short`` is a strict prefix of ``long``."""
    return {
        (short, longer)
        for short in names
        for longer in names
        if longer != short and longer.startswith(short)
    }


@pytest.fixture(scope="module")
def source() -> dict[str, str]:
    """The committed sacCer3 subsample every component chromosome is cut from."""
    return _read_fasta(DATA_DIR / "tiny.fa")


# --------------------------------------------------------------------------------------
# Provenance and shape: real bytes, the exact ones the README claims, and heterogeneous
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("component", _ALL, ids=lambda c: c.name)
def test_each_component_matches_its_declared_shape(
    component: ChimeraComponent, source: dict[str, str]
) -> None:
    sequences = _read_fasta(component.fasta)

    # Every chromosome is the exact slice of tiny.fa the README claims, in the declared
    # names and order. Case is excluded here on purpose: soft-masking is a separate claim,
    # asserted below.
    for chromosome, (cut_from, start, end) in component.slices.items():
        assert sequences[chromosome].upper() == source[cut_from][start - 1 : end]
    assert list(sequences) == component.chromosomes

    # Wrapped at its declared width...
    assert set(_sequence_line_widths(component.fasta)) == {component.line_width}

    # ...and carrying exactly the masking it declares. ce11 is soft-masked and ecHT115 is
    # not, so a chimera is heterogeneously masked; a build that upper-cased or
    # lower-cased anything loses one component's case and not the others'.
    if component.soft_masked is None:
        assert all(sequence.isupper() for sequence in sequences.values())
        return
    chromosome, bases = component.soft_masked
    assert sequences[chromosome][:bases].islower()
    assert sequences[chromosome][bases:].isupper()
    assert all(seq.isupper() for name, seq in sequences.items() if name != chromosome)


def test_the_fixture_set_is_disjoint_in_bytes_and_heterogeneous_in_shape() -> None:
    # Disjoint slices, so a test can always tell which component a chromosome came from —
    # and so no gene id is carried by two components, which would hand the annotation
    # merge a collision the shipped pair does not have.
    slices = [
        (cut_from, start, end)
        for component in _ALL
        for cut_from, start, end in component.slices.values()
    ]
    for index, (cut_from, start, end) in enumerate(slices):
        for other_cut_from, other_start, other_end in slices[index + 1 :]:
            assert cut_from != other_cut_from or end < other_start or other_end < start

    # Every component chromosome length is distinct too, so a test can name a chimera's
    # chromosome by its length alone.
    lengths = [length for component in _ALL for length in component.lengths.values()]
    assert len(set(lengths)) == len(lengths)

    # And the set disagrees with itself on wrapping and on masking: ce11 wraps at a
    # different column from ecHT115, and only ce11 is soft-masked — so a build that
    # rewrapped or re-cased anything would be visible, which nothing here would notice if
    # every component agreed.
    assert len({component.line_width for component in _ALL}) > 1
    assert len({component.soft_masked is None for component in _ALL}) > 1


# --------------------------------------------------------------------------------------
# What the naming contract needs a fixture for, and no shipped assembly supplies
# --------------------------------------------------------------------------------------


def test_the_collision_and_the_two_prefix_traps_the_contract_needs_a_fixture_for() -> None:
    # The shipped pair does not collide, so no fixture derived from it makes the suffix
    # do real work. tinyCe and tinySc do, twice.
    worm, yeast = CHIMERA_COMPONENTS["tinyCe"], CHIMERA_COMPONENTS["tinySc"]
    assert set(worm.chromosomes) & set(yeast.chromosomes) == {"I", "II"}

    # `I` inside `II` inside `III` is what the attribution regex is anchored at both ends
    # to survive.
    chromosome_names = {name for component in _ALL for name in component.chromosomes}
    assert {("I", "II"), ("I", "III"), ("II", "III")} <= _strict_prefix_pairs(chromosome_names)

    # The same trap one level up: a suffix matched without its right anchor would read
    # `NZ_TINY02000001.1__tinyEcDub` as tinyEc's. Nothing in the package can be fooled —
    # the split is a right split against a known component list — but the names that
    # would expose a hand-rolled regex are here rather than absent.
    component_names = {c.name for c in _ALL}
    assert ("tinyEc", "tinyEcDub") in _strict_prefix_pairs(component_names)


def test_the_separator_and_name_space_properties_the_contract_relies_on() -> None:
    # Both real shapes: ecHT115's `NZ_SMTD01000001.1`, an underscore and a dot in one
    # name, so no split may be a first-occurrence one; and hg38's
    # `chr1_KI270706v1_random`, which under a single-underscore separator would be
    # indistinguishable from a suffixed name. The doubled separator is what makes the
    # difference decidable, and this is the name that would catch it not being.
    ec_names = CHIMERA_COMPONENTS["tinyEc"].chromosomes
    assert all(_longest_underscore_run(name) == 1 for name in ec_names)
    assert "chr1_KI270706v1_random" in ec_names
    assert sum("_" in name for name in ec_names) == len(ec_names)
    assert not any("__" in name for name in ec_names)

    # The separator is the shortest run of underscores, minimum two, strictly longer than
    # the longest run any component's chromosome names carry. The everyday set therefore
    # asks for `__`, as all seven shipped assemblies do, and only the escalation
    # component moves it.
    everyday = [name for key in CHIMERA_EVERYDAY for name in CHIMERA_COMPONENTS[key].chromosomes]
    escalation = CHIMERA_COMPONENTS[CHIMERA_ESCALATION].chromosomes
    assert max(_longest_underscore_run(name) for name in everyday) == 1
    assert max(_longest_underscore_run(name) for name in escalation) == 2

    # Component names are alphanumeric — enforced by the model, not accidental: it is
    # what makes the last-occurrence split sound and the derived chimera name injective —
    # and none collides with a curated row, which would silently pull that row's pinned
    # source and checksum into a fixture that is seeded from a local file.
    assert all(component.name.isalnum() for component in _ALL)
    assert all(lookup_assembly(component.name) is None for component in _ALL)


# --------------------------------------------------------------------------------------
# The annotations the merge is written against
# --------------------------------------------------------------------------------------


def test_the_escalation_ships_no_annotation_and_no_gene_id_is_shared() -> None:
    # tinyEcDub exists for a naming rule; an annotation would only be a second thing to
    # maintain.
    assert CHIMERA_COMPONENTS[CHIMERA_ESCALATION].gtf is None

    # An id collision is a problem the shipped pair does not have, so the fixture set does
    # not invent one for the merge to answer.
    owners: dict[str, str] = {}
    for component in _ALL:
        if component.gtf is None:
            continue
        for gene in set(re.findall(r'gene_id "([^"]+)"', component.gtf.read_text())):
            assert gene not in owners, f"{gene} is in both {owners.get(gene)} and {component.name}"
            owners[gene] = component.name


@pytest.mark.parametrize("component", _ALL, ids=lambda c: c.name)
def test_an_annotation_names_every_chromosome_of_its_component_and_no_other(
    component: ChimeraComponent,
) -> None:
    # Set-equal, as both shipped GTFs are against their assemblies — so a merge that drops
    # or mis-spells one name is caught by the chromosome check rather than passing.
    if component.gtf is None:
        pytest.skip(f"{component.name} ships no annotation")
    assert _chromosome_names(component.gtf) == set(component.chromosomes)


# --------------------------------------------------------------------------------------
# They register — offline, with the real tools, as any other assembly does
# --------------------------------------------------------------------------------------


def test_a_component_registers_as_an_assembly_of_the_expected_sizes(
    chimera_component: ComponentFactory,
) -> None:
    # tinyCe alone: the others are registered as a matter of course by every build test,
    # so this is the one place proving the registration path itself — chromosomes and
    # sizes agreeing with what the table declares — rather than assuming it.
    component = CHIMERA_COMPONENTS["tinyCe"]
    genome = chimera_component(component.name)
    assert genome.chromosomes == component.chromosomes
    assert dict(genome.chrom_sizes) == component.lengths


def test_an_annotation_registers_against_its_own_component_and_reaches_no_network(
    chimera_component: ComponentFactory,
) -> None:
    # On its merits: the chromosome check is left standing, which is the point of the
    # names being set-equal.
    genome = chimera_component("tinyCe", with_annotation=True)
    assert genome.annotations.registered == [COMPONENT_ANNOTATION]
    assert genome.default_gtf == COMPONENT_ANNOTATION

    # The autouse network guard would raise; a component seeded from a local path under a
    # name that is not a UCSC one says so deliberately rather than by luck.
    assert lookup_assembly(genome.assembly) is None
