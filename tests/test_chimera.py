"""Tests for genome.chimera — the naming contract the build and the merge share.

Pure name arithmetic, so nothing here needs a native binary, a temporary directory or a
byte of sequence: what the committed component fixtures supply is their *names*, which
carry the collision, the two prefix traps and the underscore runs no shipped assembly can
demonstrate. What those components are cut from is pinned in ``test_chimera_fixtures``;
this module asserts the contract, never the fixtures.

Hypothesis pins the invariants the fixtures can only sample: suffixing then splitting is
the identity, the derived name does not depend on the order its components arrived in,
the derived separator always beats every run in its inputs, and the published pattern
reads *any* string exactly as the split does — the fixture names all being well formed is
what let those two drift apart once.
"""

from __future__ import annotations

import itertools
import re
import string
from collections.abc import Iterable

import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from genome.chimera import (
    ChimeraNamingError,
    _suffix_pattern,
    check_roundtrip,
    derive_name,
    derive_separator,
    split_name,
    split_suffixed,
    suffixed,
)

from .conftest import CHIMERA_COMPONENTS, CHIMERA_ESCALATION, CHIMERA_EVERYDAY

#: Every component set a chimera could be built from the fixtures — pairs, triples and
#: the whole four. Eleven of them, which is cheap enough to check every name of every one.
_COMBINATIONS = [
    combination
    for size in range(2, len(CHIMERA_COMPONENTS) + 1)
    for combination in itertools.combinations(CHIMERA_COMPONENTS, size)
]


def _chromosomes(*components: str) -> dict[str, list[str]]:
    """Map each named component to its chromosome names, as a build would."""
    return {name: CHIMERA_COMPONENTS[name].chromosomes for name in components}


def _longest_run(names: Iterable[str]) -> int:
    """Length of the longest run of underscores across ``names`` — measured, not imported."""
    return max((len(run) for name in names for run in re.findall(r"_+", name)), default=0)


def _ids(combination: tuple[str, ...]) -> str:
    return "+".join(combination)


# --------------------------------------------------------------------------------------
# The derived name
# --------------------------------------------------------------------------------------


def test_the_name_is_the_component_names_sorted_and_joined() -> None:
    assert derive_name(CHIMERA_EVERYDAY) == "tinyCe_tinyEc_tinySc"


def test_the_shipped_pair_derives_the_name_the_map_names() -> None:
    assert derive_name(["ecHT115", "ce11"]) == "ce11_ecHT115"


@pytest.mark.parametrize("order", list(itertools.permutations(CHIMERA_EVERYDAY)))
def test_the_name_does_not_depend_on_the_order_the_components_arrived_in(
    order: tuple[str, ...],
) -> None:
    # Identity is the component set, so every one of the six orderings is the same chimera.
    assert derive_name(order) == derive_name(CHIMERA_EVERYDAY)


def test_a_chimera_of_one_is_refused() -> None:
    with pytest.raises(ChimeraNamingError, match="at least 2 components"):
        derive_name(["tinyCe"])


def test_a_repeated_component_is_refused() -> None:
    with pytest.raises(ChimeraNamingError, match="must not repeat"):
        derive_name(["tinyCe", "tinySc", "tinyCe"])


@pytest.mark.parametrize(
    "component", ["tiny Ce", "tiny-Ce", "tiny_Ce", "tiny.Ce", "", "tinyCé", "tinyCe\n"]
)
def test_a_non_alphanumeric_component_name_is_refused(component: str) -> None:
    # ASCII alphanumeric, which is narrower than str.isalnum: 'tinyCé' passes that and
    # would then be unreadable by every consumer holding the published regex.
    with pytest.raises(ChimeraNamingError, match="is not alphanumeric"):
        derive_name([component, "tinySc"])


def test_a_chimera_cannot_be_a_component_by_name() -> None:
    # Nesting is forbidden by the model. The half of it a pure module can see is the
    # spelling: a derived name always carries '_', which no component name may.
    with pytest.raises(ChimeraNamingError, match="is not alphanumeric"):
        derive_name([derive_name(CHIMERA_EVERYDAY), "tinyEcDub"])


@pytest.mark.parametrize("components", [["a_b", "c"], ["a", "b_c"]])
def test_the_ambiguity_that_would_break_injectivity_cannot_be_constructed(
    components: list[str],
) -> None:
    # Both would derive 'a_b_c'. The alphanumeric rule is what stops either from existing,
    # which is what makes the derived name injective rather than merely usually unique.
    with pytest.raises(ChimeraNamingError, match="is not alphanumeric"):
        derive_name(components)


def test_one_name_is_not_a_component_list() -> None:
    # Iterating a str would silently derive '1_1_c_e' from 'ce11'.
    with pytest.raises(ChimeraNamingError, match="not one name"):
        derive_name("ce11")


# --------------------------------------------------------------------------------------
# Its inverse
# --------------------------------------------------------------------------------------


def test_the_name_splits_back_into_its_components() -> None:
    assert split_name("tinyCe_tinyEc_tinySc") == ("tinyCe", "tinyEc", "tinySc")


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_every_derived_name_splits_back_to_the_set_it_came_from(
    combination: tuple[str, ...],
) -> None:
    assert split_name(derive_name(combination)) == tuple(sorted(combination))


def test_a_mis_ordered_name_can_be_told_its_canonical_spelling() -> None:
    # Candidates come back in the order the name spells them, so the round trip through
    # derive_name is what detects — and names — a mis-ordering.
    given_name = "ecHT115_ce11"
    assert split_name(given_name) == ("ecHT115", "ce11")
    assert derive_name(split_name(given_name)) == "ce11_ecHT115"


@pytest.mark.parametrize("name", ["hg38", "", "_tinyCe", "tinyCe_", "tinyCe__tinySc", "tiny-Ce_x"])
def test_a_name_not_spelled_like_a_chimeras_is_refused(name: str) -> None:
    with pytest.raises(ChimeraNamingError, match="not spelled like a chimera"):
        split_name(name)


def test_splitting_is_syntactic_and_claims_nothing_about_existence() -> None:
    # 'my_ref' is spelled like a chimera and is not one. Which it is, is settled by asking
    # whether each candidate is prepared or listed — not here.
    assert split_name("my_ref") == ("my", "ref")


# --------------------------------------------------------------------------------------
# The derived separator
# --------------------------------------------------------------------------------------


def test_the_everyday_components_derive_two_underscores() -> None:
    # As all seven shipped assemblies do: the longest run any of them carries is one.
    assert derive_separator(_chromosomes(*CHIMERA_EVERYDAY)) == "__"


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_only_the_escalation_component_pushes_the_separator_past_two(
    combination: tuple[str, ...],
) -> None:
    separator = derive_separator(_chromosomes(*combination))
    if CHIMERA_ESCALATION in combination:
        assert separator == "___"
        assert len(separator) > 2
    else:
        assert separator == "__"


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_the_separator_is_strictly_longer_than_every_run_its_components_carry(
    combination: tuple[str, ...],
) -> None:
    chromosomes = _chromosomes(*combination)
    separator = derive_separator(chromosomes)
    every_name = [name for names in chromosomes.values() for name in names]
    assert len(separator) > _longest_run(every_name)
    assert len(separator) >= 2


def test_the_escalation_component_alone_is_what_forces_the_longer_run() -> None:
    # The doubled underscore is in a chromosome name, not in the component's own name:
    # NZ_TINY02__000002.1 would read as though it were already suffixed under '__'.
    assert _longest_run(CHIMERA_COMPONENTS[CHIMERA_ESCALATION].chromosomes) == 2
    assert derive_separator(_chromosomes("tinyEc", CHIMERA_ESCALATION)) == "___"


def test_a_separator_is_derived_from_a_whole_component_set_not_from_one_component() -> None:
    with pytest.raises(ChimeraNamingError, match="at least 2 components"):
        derive_separator(_chromosomes("tinyCe"))


def test_a_component_name_in_the_mapping_is_checked_too() -> None:
    with pytest.raises(ChimeraNamingError, match="is not alphanumeric"):
        derive_separator({"tiny Ce": ["I"], "tinySc": ["II"]})


# --------------------------------------------------------------------------------------
# Suffixing and splitting
# --------------------------------------------------------------------------------------


def test_a_real_collision_is_resolved_by_the_suffix() -> None:
    # tinyCe and tinySc both carry I and II — the collision the shipped pair never has.
    chromosomes = _chromosomes("tinyCe", "tinySc")
    separator = derive_separator(chromosomes)
    shared = set(chromosomes["tinyCe"]) & set(chromosomes["tinySc"])
    assert shared == {"I", "II"}
    for chromosome in sorted(shared):
        worm = suffixed(chromosome, "tinyCe", separator)
        yeast = suffixed(chromosome, "tinySc", separator)
        assert worm != yeast
        assert split_suffixed(worm, separator) == (chromosome, "tinyCe")
        assert split_suffixed(yeast, separator) == (chromosome, "tinySc")


def test_the_suffix_is_unconditional() -> None:
    # MtDNA collides with nothing and is suffixed anyway, so attribution is one operation
    # for every name in the reference rather than two.
    assert suffixed("MtDNA", "tinyCe", "__") == "MtDNA__tinyCe"


@pytest.mark.parametrize(("short", "longer"), [("I", "II"), ("II", "III")])
def test_the_prefix_trap_between_chromosome_names(short: str, longer: str) -> None:
    # I inside II, and II inside III. An unanchored search for the shorter name's suffix
    # finds the longer name; the split, and the anchored regex, do not.
    name = suffixed(longer, "tinySc", "__")
    assert split_suffixed(name) == (longer, "tinySc")
    assert re.search(f"{short}__tinySc", name) is not None
    assert re.fullmatch(f"^{short}__tinySc$", name) is None


def test_the_prefix_trap_between_component_names() -> None:
    # The same trap one level up: tinyEc inside tinyEcDub. A regex missing its right
    # anchor reads NZ_TINY02000001.1___tinyEcDub as tinyEc's.
    chromosomes = _chromosomes("tinyEc", CHIMERA_ESCALATION)
    separator = derive_separator(chromosomes)
    name = suffixed("NZ_TINY02000001.1", CHIMERA_ESCALATION, separator)
    assert name == "NZ_TINY02000001.1___tinyEcDub"
    assert re.search(f"{separator}tinyEc", name) is not None  # the trap
    assert split_suffixed(name, separator) == ("NZ_TINY02000001.1", "tinyEcDub")
    match = re.match(_suffix_pattern(separator), name)
    assert match is not None
    assert match["component"] == "tinyEcDub"


def test_a_chromosome_name_already_carrying_an_underscore_and_a_dot_survives() -> None:
    # ecHT115's accession shape, so no split may be a first-occurrence one.
    assert split_suffixed("NZ_TINY01000001.1__tinyEc") == ("NZ_TINY01000001.1", "tinyEc")


def test_the_decoy_is_what_a_single_underscore_could_not_tell_from_a_suffix() -> None:
    # chr1_KI270706v1_random is a real hg38 name. Under '_' it is indistinguishable from a
    # suffixed one; under the derived separator it announces that it is not suffixed, and
    # it still reads back whole once it is.
    decoy = "chr1_KI270706v1_random"
    assert "__" not in decoy
    with pytest.raises(ChimeraNamingError, match="carries no component suffix"):
        split_suffixed(decoy)
    assert split_suffixed(suffixed(decoy, "tinyEc", "__")) == (decoy, "tinyEc")


def test_the_split_is_at_the_last_run_not_the_first() -> None:
    # A component whose own chromosome is already named like a suffixed one still splits
    # to the component that actually contributed it.
    assert split_suffixed("bar__tinyCe__tinySc") == ("bar__tinyCe", "tinySc")


@pytest.mark.parametrize("name", ["I", "I_tinyCe", "I__", "I__tiny-Ce", "__", "I__tinyCe\n"])
def test_a_name_with_no_component_suffix_is_refused(name: str) -> None:
    with pytest.raises(ChimeraNamingError, match="carries no component suffix"):
        split_suffixed(name)


def test_reading_defaults_to_two_underscores() -> None:
    assert split_suffixed("I__tinyCe") == ("I", "tinyCe")


def test_writing_takes_no_default_separator() -> None:
    # Deliberate asymmetry: a build that guessed a constant would lose the self-announcing
    # property silently, so suffixed() makes the caller pass the derived separator.
    with pytest.raises(TypeError):
        suffixed("I", "tinyCe")  # type: ignore[call-arg]


@pytest.mark.parametrize("separator", ["_", "", "-", "___x", " __"])
def test_an_illegal_separator_is_refused_at_both_ends(separator: str) -> None:
    with pytest.raises(ChimeraNamingError, match="run of two or more underscores"):
        suffixed("I", "tinyCe", separator)
    with pytest.raises(ChimeraNamingError, match="run of two or more underscores"):
        split_suffixed("I__tinyCe", separator)
    with pytest.raises(ChimeraNamingError, match="run of two or more underscores"):
        _suffix_pattern(separator)


def test_a_non_alphanumeric_component_cannot_be_suffixed_onto_anything() -> None:
    with pytest.raises(ChimeraNamingError, match="is not alphanumeric"):
        suffixed("I", "tiny_Ce", "__")


def test_suffixing_nothing_spells_the_tail_alone_and_reads_back_as_no_chromosome() -> None:
    # The one call that is not a chromosome name: the merged-GTF writer spells the tail
    # every seqname of one component gains with it, once instead of per line. Writing it
    # is legal and reading it back as a name is not, which is the whole asymmetry.
    assert suffixed("", "tinyCe", "__") == "__tinyCe"
    with pytest.raises(ChimeraNamingError, match="suffix and nothing else"):
        split_suffixed("__tinyCe")


@pytest.mark.parametrize(("name", "separator"), [("__tinyCe", "__"), ("___tinyCe", "___")])
def test_a_name_that_is_suffix_and_nothing_else_is_refused(name: str, separator: str) -> None:
    # The one disagreement the published pattern and the split used to have: the regex
    # spells the chromosome `.+` and matched neither of these, while the split returned an
    # empty chromosome for both.
    with pytest.raises(ChimeraNamingError, match="suffix and nothing else"):
        split_suffixed(name, separator)
    assert re.match(_suffix_pattern(separator), name) is None


# --------------------------------------------------------------------------------------
# The published regex
# --------------------------------------------------------------------------------------


def test_the_published_pattern_is_the_documented_one() -> None:
    assert _suffix_pattern() == r"^(?P<chromosome>.+)__(?P<component>[A-Za-z0-9]+)$"


def test_the_published_pattern_is_generated_from_the_separator() -> None:
    assert _suffix_pattern("___") == r"^(?P<chromosome>.+)___(?P<component>[A-Za-z0-9]+)$"


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_the_published_pattern_and_the_split_agree_on_every_name(
    combination: tuple[str, ...],
) -> None:
    # The regex is the contract awk, R and shell consumers attribute names by, so it has
    # to say exactly what the package does.
    chromosomes = _chromosomes(*combination)
    separator = derive_separator(chromosomes)
    pattern = _suffix_pattern(separator)
    for component, names in chromosomes.items():
        for chromosome in names:
            match = re.match(pattern, suffixed(chromosome, component, separator))
            assert match is not None
            assert (match["chromosome"], match["component"]) == (chromosome, component)


# --------------------------------------------------------------------------------------
# The round-trip assertion the build runs before writing
# --------------------------------------------------------------------------------------


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_every_name_of_every_combination_survives_the_round_trip(
    combination: tuple[str, ...],
) -> None:
    chromosomes = _chromosomes(*combination)
    check_roundtrip(chromosomes, derive_separator(chromosomes))


@pytest.mark.parametrize("combination", _COMBINATIONS, ids=_ids)
def test_the_check_is_what_the_build_would_otherwise_hand_roll(
    combination: tuple[str, ...],
) -> None:
    # Same claim, spelled out name by name, so a check that silently checked nothing would
    # be visible here.
    chromosomes = _chromosomes(*combination)
    separator = derive_separator(chromosomes)
    for component, names in chromosomes.items():
        for chromosome in names:
            assert split_suffixed(suffixed(chromosome, component, separator), separator) == (
                chromosome,
                component,
            )


def test_a_constant_separator_is_refused_when_a_component_forces_a_longer_one() -> None:
    # The failure the round trip alone provably cannot catch: splitting at the last run is
    # correct under '__' too, so only comparing against the derived separator catches it.
    chromosomes = _chromosomes("tinyEc", CHIMERA_ESCALATION)
    name = suffixed("NZ_TINY02__000002.1", CHIMERA_ESCALATION, "__")
    assert split_suffixed(name) == ("NZ_TINY02__000002.1", "tinyEcDub")  # round trip is fine
    with pytest.raises(ChimeraNamingError, match="not the one these components derive"):
        check_roundtrip(chromosomes, "__")


def test_a_separator_longer_than_the_derived_one_is_refused_too() -> None:
    # The separator is *the shortest* run that works, so one chimera has exactly one.
    with pytest.raises(ChimeraNamingError, match="not the one these components derive"):
        check_roundtrip(_chromosomes(*CHIMERA_EVERYDAY), "____")


def test_the_check_refuses_an_illegal_component_set() -> None:
    with pytest.raises(ChimeraNamingError, match="at least 2 components"):
        check_roundtrip(_chromosomes("tinyCe"), "__")


def test_the_check_refuses_a_component_chromosome_with_no_name() -> None:
    # A name the round trip used to wave through: it suffixed to '__tinyCe' and split back
    # to ('', 'tinyCe'), while the published pattern matched nothing at all.
    with pytest.raises(ChimeraNamingError, match="suffix and nothing else"):
        check_roundtrip({"tinyCe": [""], "tinySc": ["I"]}, "__")


# --------------------------------------------------------------------------------------
# The invariants, over generated names
# --------------------------------------------------------------------------------------

_ALPHANUMERIC = string.ascii_letters + string.digits

#: A legal component assembly name.
_component = st.text(alphabet=_ALPHANUMERIC, min_size=1, max_size=8)

#: A chromosome name of the shapes real assemblies use — underscores and dots included,
#: since those are what the contract has to survive.
_chromosome = st.text(alphabet=_ALPHANUMERIC + "_.-", min_size=1, max_size=16)

_separator = st.integers(min_value=2, max_value=6).map(lambda length: "_" * length)

_component_set = st.lists(_component, min_size=2, max_size=5, unique=True)

_chromosome_map = st.dictionaries(
    _component, st.lists(_chromosome, min_size=1, max_size=4), min_size=2, max_size=4
)


@st.composite
def _name_under_a_separator(draw: st.DrawFn) -> tuple[str, str]:
    """Draw a ``(name, separator)`` pair — any name at all, legal or not.

    The separator is one of the pieces the name is built from, because a name drawn from
    letters and lone underscores lands on the shapes the two readers could disagree about
    far too rarely — a whole name that is one separator run and a component, say. No
    newline in the alphabet: :func:`_suffix_pattern` anchors with ``$`` and
    :func:`split_suffixed` at the end of the string, so a trailing newline is a documented
    difference between them rather than drift.
    """
    separator = draw(_separator)
    pieces = draw(st.lists(st.sampled_from(["a", "1", ".", "_", separator]), max_size=5))
    return "".join(pieces), separator


@given(chromosome=_chromosome, component=_component, separator=_separator)
def test_suffixing_then_splitting_is_the_identity(
    chromosome: str, component: str, separator: str
) -> None:
    name = suffixed(chromosome, component, separator)
    assert split_suffixed(name, separator) == (chromosome, component)


@given(chromosome=_chromosome, component=_component, separator=_separator)
def test_the_published_pattern_reads_what_suffixing_wrote(
    chromosome: str, component: str, separator: str
) -> None:
    match = re.match(_suffix_pattern(separator), suffixed(chromosome, component, separator))
    assert match is not None
    assert (match["chromosome"], match["component"]) == (chromosome, component)


@given(pair=_name_under_a_separator())
@example(pair=("__a", "__"))
@example(pair=("___a", "___"))
def test_the_published_pattern_and_the_split_agree_on_any_name_at_all(
    pair: tuple[str, str],
) -> None:
    # Not only on names a build wrote: the pattern is the contract awk, R and shell
    # consumers hold, so any string either splits the same way in both or is refused by
    # both. Reading a name one way here and another way there is the drift this module
    # exists to prevent, and the fixtures cannot show it — the two pinned examples are the
    # whole disagreement class, kept explicit so the property is never left to luck.
    name, separator = pair
    match = re.match(_suffix_pattern(separator), name)
    try:
        split = split_suffixed(name, separator)
    except ChimeraNamingError:
        assert match is None
        return
    assert match is not None
    assert (match["chromosome"], match["component"]) == split


@given(components=_component_set, data=st.data())
def test_the_derived_name_does_not_depend_on_order(
    components: list[str], data: st.DataObject
) -> None:
    shuffled = data.draw(st.permutations(components))
    assert derive_name(shuffled) == derive_name(components)


@given(components=_component_set)
def test_the_derived_name_splits_back_to_its_component_set(components: list[str]) -> None:
    # Which is injectivity: two different sets cannot reach the same name, because the
    # name determines the set.
    assert split_name(derive_name(components)) == tuple(sorted(components))


@given(chromosomes=_chromosome_map)
def test_the_derived_separator_beats_every_run_in_its_inputs(
    chromosomes: dict[str, list[str]],
) -> None:
    separator = derive_separator(chromosomes)
    every_name = [name for names in chromosomes.values() for name in names]
    assert set(separator) == {"_"}
    assert len(separator) >= 2
    assert len(separator) > _longest_run(every_name)
    assert len(separator) == max(2, _longest_run(every_name) + 1)


@given(chromosomes=_chromosome_map)
def test_the_check_accepts_whatever_the_derivation_produced(
    chromosomes: dict[str, list[str]],
) -> None:
    check_roundtrip(chromosomes, derive_separator(chromosomes))
