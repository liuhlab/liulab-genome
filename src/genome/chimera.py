"""The naming contract a chimera's build and its annotation merge share.

A chimera is an assembly whose FASTA is concatenated from two or more prepared canonical
assemblies — its components (ADR-0008). Two spellings have to be agreed on by everything
that touches one: the chimera's own assembly name, derived from its component names, and
every chromosome name, which carries the component it came from as a suffix (ADR-0009).
Both are decided here, in one pure module, so that the code writing a chimera's FASTA and
the code merging its annotation cannot drift apart about them.

Pure means names in, names out: nothing here opens a file, and nothing here imports from
``genome.io`` or ``genome.genome``.

Two joins at two levels, spelled with deliberately different characters, so that reading
a name tells you which level you are at:

- **The assembly name** joins component names with a single ``_``, sorted:
  ``ce11_ecHT115``. :func:`derive_name` builds it; :func:`split_name` takes it apart.
- **A chromosome name** joins the chromosome to its component with a *run* of
  underscores, never fewer than two: ``I__ce11``. :func:`suffixed` builds it;
  :func:`split_suffixed` takes it apart; :func:`derive_separator` says how long the run
  must be for a given set of components, and :func:`check_roundtrip` proves the whole set
  survives before a byte is written.

The doubled separator is not what makes the split work — :func:`split_suffixed` splits at
the *last* run and would be correct with any separator at all. It is what makes a suffix
**announce itself**: under a single ``_``, the real hg38 chromosome
``chr1_KI270706v1_random`` is indistinguishable from a suffixed name, so a name arriving
off a BAM header could not be classified without already knowing which reference it came
from.

Component names are alphanumeric, and this module enforces that rather than assuming it.
One rule buys two things: the derived assembly name is injective (otherwise components
``a_b`` and ``c`` would derive the same name as ``a`` and ``b_c``), and a suffixed
chromosome name can be split back to the component that contributed it.

Examples
--------
>>> from genome.chimera import derive_name, derive_separator, split_suffixed, suffixed
>>> derive_name(["ecHT115", "ce11"])          # order in, sorted out
'ce11_ecHT115'
>>> chromosomes = {"ce11": ["I", "MtDNA"], "ecHT115": ["NZ_SMTD01000001.1"]}
>>> separator = derive_separator(chromosomes)
>>> separator
'__'
>>> suffixed("I", "ce11", separator)
'I__ce11'
>>> split_suffixed("NZ_SMTD01000001.1__ecHT115")
('NZ_SMTD01000001.1', 'ecHT115')
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence

#: What a component assembly name may be. ASCII on purpose: this is the character class
#: the published attribution regex carries, and ``str.isalnum`` is wider than it — a name
#: like ``café`` passes that and would then be unreadable by every consumer holding the
#: regex. Matched with ``fullmatch`` rather than anchored with ``$``, which in Python also
#: matches before a trailing newline — and a newline in a name is a broken FASTA header.
_COMPONENT_RE = re.compile(r"[A-Za-z0-9]+")

#: What a separator may be: a run of underscores, never fewer than two. Also ``fullmatch``.
_SEPARATOR_RE = re.compile(r"_{2,}")

#: One run of underscores, however long — used to measure what a component's chromosome
#: names already carry.
_UNDERSCORE_RUN_RE = re.compile(r"_+")

#: The join between component names in a chimera's own assembly name.
_NAME_JOIN = "_"

#: The fewest components a chimera can be built from. One is that assembly under another
#: name and a wasted copy of its bytes.
_MIN_COMPONENTS = 2

#: The shortest separator, whatever the components turn out to carry.
_MIN_SEPARATOR = 2


class ChimeraNamingError(ValueError):
    """A name does not obey the chimera naming contract.

    One type covers every rejection this module makes — a component name that is not
    alphanumeric, a component set that is too small or repeats itself, an assembly name
    that is not spelled like a chimera's, a chromosome name carrying no component suffix,
    a separator that was not derived from the components it is about to be written with.
    The message says which, and what to do about it.

    It is a :class:`ValueError` because every one of those is a bad value, and it is named
    so that a caller deciding *whether* a string is a chimera's name can catch exactly
    this and nothing else.

    Examples
    --------
    >>> from genome.chimera import ChimeraNamingError, split_name
    >>> try:
    ...     split_name("hg38")
    ... except ChimeraNamingError:
    ...     print("not spelled like a chimera")
    not spelled like a chimera
    """


def derive_name(components: Iterable[str]) -> str:
    """Derive a chimera's assembly name from the names of its components.

    The name is the component names sorted lexicographically and joined by ``_``. It is
    derived rather than given and cannot be overridden: one set of components means
    exactly one name, and therefore one directory and one index, whatever order the
    caller happened to list them in.

    A chimera's own name always carries a ``_``, and a component's name never may, so
    this is also where a chimera is refused as a component of another chimera. That is
    only the spelling half of the rule: whether a *prepared* assembly is itself a chimera
    is answered by the record on its disk, which this module cannot see.

    Parameters
    ----------
    components : Iterable[str]
        Two or more component assembly names, in any order. Each must match
        ``[A-Za-z0-9]+`` and appear exactly once.

    Returns
    -------
    str
        The derived assembly name.

    Raises
    ------
    ChimeraNamingError
        If fewer than two components are given, a component repeats, or a component name
        is not alphanumeric.

    Examples
    --------
    >>> derive_name(["ecHT115", "ce11"])
    'ce11_ecHT115'
    >>> derive_name(["tinySc", "tinyCe", "tinyEc"])
    'tinyCe_tinyEc_tinySc'
    """
    if isinstance(components, str):
        raise ChimeraNamingError(
            f"derive_name takes the component names, not one name: {components!r} was "
            f"iterated character by character. Pass a list, as in "
            f"derive_name(['ce11', 'ecHT115'])."
        )
    names = tuple(components)
    _check_component_set(names)
    return _NAME_JOIN.join(sorted(names))


def split_name(name: str) -> tuple[str, ...]:
    """Split a chimera's assembly name back into its candidate component names.

    The inverse of :func:`derive_name`, and syntactic only: it says that ``name`` is
    *spelled* the way a chimera's name is spelled, never that those components exist.
    Deciding that belongs to the caller, which asks whether each candidate is prepared on
    this machine or listed in the shipped table — the step that separates ``ce11_ecHT115``
    from an ordinary assembly someone happened to call ``my_ref``.

    Candidates come back in the order the name spells them, not sorted, so that a caller
    can hand them straight to :func:`derive_name` and compare: a name whose components are
    real but mis-ordered is detectable, and the canonical spelling can be named back.

    Parameters
    ----------
    name : str
        An assembly name to read as a chimera's.

    Returns
    -------
    tuple[str, ...]
        The candidate component names, in the order ``name`` spells them.

    Raises
    ------
    ChimeraNamingError
        If ``name`` does not split into two or more alphanumeric parts.

    Examples
    --------
    >>> split_name("ce11_ecHT115")
    ('ce11', 'ecHT115')
    >>> derive_name(split_name("ecHT115_ce11"))   # the canonical spelling of a mis-ordered name
    'ce11_ecHT115'
    """
    parts = tuple(name.split(_NAME_JOIN))
    if len(parts) < _MIN_COMPONENTS or not all(_COMPONENT_RE.fullmatch(part) for part in parts):
        raise ChimeraNamingError(
            f"{name!r} is not spelled like a chimera's name; a chimera is named by its "
            f"component assembly names — at least {_MIN_COMPONENTS}, each alphanumeric — "
            f"sorted and joined by '_', as in 'ce11_ecHT115'. If this assembly is not a "
            f"chimera, open it under its own name."
        )
    return parts


def derive_separator(chromosomes: Mapping[str, Sequence[str]]) -> str:
    """Derive the separator that spells one chimera's chromosome names.

    The separator is the shortest run of underscores that is at least two long *and*
    strictly longer than the longest run of underscores any component's chromosome names
    already carry. Two is the floor because that is what makes a suffix announce itself;
    strictly longer is what keeps it announcing itself when a component's own chromosome
    is named something like ``NZ_TINY02__000002.1``, which under a two-underscore
    separator would read as though it were already suffixed.

    Every assembly the lab ships derives ``__``: the longest run any of them carries is
    one, in names like ``NZ_SMTD01000001.1`` and ``chr1_KI270706v1_random``.

    The separator belongs to one chimera rather than to this package: it is derived here,
    used to write every name, and recorded, so that the merge and every later reader use
    the separator the chimera was actually built with instead of a constant.

    Parameters
    ----------
    chromosomes : Mapping[str, Sequence[str]]
        Each component assembly name mapped to that component's chromosome names.

    Returns
    -------
    str
        A run of two or more underscores.

    Raises
    ------
    ChimeraNamingError
        If fewer than two components are given, or a component name is not alphanumeric.

    Examples
    --------
    >>> derive_separator({"ce11": ["I", "MtDNA"], "ecHT115": ["NZ_SMTD01000001.1"]})
    '__'
    >>> derive_separator({"tinyEc": ["chr1_KI270706v1_random"], "tinyEcDub": ["NZ_TINY02__000002.1"]})
    '___'
    """
    _check_component_set(tuple(chromosomes))
    longest = _longest_underscore_run(name for names in chromosomes.values() for name in names)
    return "_" * max(_MIN_SEPARATOR, longest + 1)


def suffixed(chromosome: str, component: str, separator: str) -> str:
    """Spell one chimera chromosome name: the chromosome, the separator, the component.

    Unconditional — a chromosome is suffixed whether or not another component carries the
    same name — so attribution is the same operation for every name in the reference and
    no mapping of any kind has to be stored to perform it.

    ``separator`` is required and deliberately has no default, unlike
    :func:`split_suffixed`: it belongs to one chimera, comes from
    :func:`derive_separator`, and a build that wrote a constant instead would quietly lose
    the self-announcing property the derivation exists to preserve, in a way no round trip
    can detect.

    Parameters
    ----------
    chromosome : str
        The chromosome name as the component itself spells it.
    component : str
        The component assembly name, alphanumeric.
    separator : str
        The run of underscores this chimera derived.

    Returns
    -------
    str
        The suffixed chromosome name.

    Raises
    ------
    ChimeraNamingError
        If ``component`` is not alphanumeric, or ``separator`` is not a run of two or more
        underscores.

    Examples
    --------
    >>> suffixed("I", "ce11", "__")
    'I__ce11'
    >>> suffixed("NZ_TINY02__000002.1", "tinyEcDub", "___")
    'NZ_TINY02__000002.1___tinyEcDub'
    """
    _check_component(component)
    _check_separator(separator)
    return f"{chromosome}{separator}{component}"


def split_suffixed(name: str, separator: str = "__") -> tuple[str, str]:
    """Split a suffixed chromosome name back into ``(chromosome, component)``.

    The split is at the **last** run of the separator, which is unconditionally the right
    one: what follows the suffix is an alphanumeric component name, so no later run of the
    separator can exist to be mistaken for it. A component that already spells a
    chromosome ``bar_ce11`` therefore still reads back correctly — ``bar_ce11__ecHT115``
    is ``('bar_ce11', 'ecHT115')``.

    ``separator`` defaults to ``'__'``, which is what a chimera derives whenever no
    component carries a doubled underscore of its own — every assembly the lab ships.
    When one does, the chimera's record is the authority on its separator and the caller
    passes it. Writing has no such default, on purpose: see :func:`suffixed`.

    The same contract in the form to hand a tool that cannot import this module is what
    :func:`suffix_pattern` returns.

    Parameters
    ----------
    name : str
        A suffixed chromosome name.
    separator : str, default ``"__"``
        The run of underscores this chimera recorded.

    Returns
    -------
    tuple[str, str]
        ``(chromosome, component)`` — the chromosome as its component spells it, and the
        component assembly name.

    Raises
    ------
    ChimeraNamingError
        If ``separator`` is not a run of two or more underscores, or ``name`` carries no
        component suffix under it.

    Examples
    --------
    >>> split_suffixed("I__ce11")
    ('I', 'ce11')
    >>> split_suffixed("chr1_KI270706v1_random__tinyEc")
    ('chr1_KI270706v1_random', 'tinyEc')
    >>> split_suffixed("NZ_TINY02__000002.1___tinyEcDub", "___")
    ('NZ_TINY02__000002.1', 'tinyEcDub')
    """
    _check_separator(separator)
    chromosome, found, component = name.rpartition(separator)
    if not found or not _COMPONENT_RE.fullmatch(component):
        raise ChimeraNamingError(
            f"chromosome name {name!r} carries no component suffix under separator "
            f"{separator!r}; a chimera's chromosome is spelled "
            f"<chromosome>{separator}<component> with an alphanumeric component, as in "
            f"'I{separator}ce11'. Check the separator against the one this chimera "
            f"recorded, and check that the name came from a chimera at all."
        )
    return chromosome, component


def suffix_pattern(separator: str = "__") -> str:
    """Return the published regex that reads a suffixed chromosome name.

    The same contract :func:`split_suffixed` performs, in the form to hand something that
    cannot import this module — an awk field split, an R ``sub``, a shell one-liner. It is
    generated from the separator rather than written down, so a chimera whose components
    forced a longer run gets its own pattern instead of the two-underscore one.

    Greedy on the left, so the split falls at the last separator run. **Anchored
    independently at each end**, which is what stops ``I`` from matching inside ``II`` and
    ``tinyEc`` from claiming a name that ends ``tinyEcDub``: an unanchored search for a
    particular component's suffix finds both.

    The right anchor is ``$``, an end-of-line one, because the tools this is for read a
    line at a time and a chromosome name is one field of one line. :func:`split_suffixed`
    anchors at the end of the string instead, since a Python caller can hold a name with a
    trailing newline in it that a line-oriented tool never could.

    Parameters
    ----------
    separator : str, default ``"__"``
        The run of underscores this chimera recorded.

    Returns
    -------
    str
        A regex with named groups ``chromosome`` and ``component``.

    Raises
    ------
    ChimeraNamingError
        If ``separator`` is not a run of two or more underscores.

    Examples
    --------
    >>> suffix_pattern()
    '^(?P<chromosome>.+)__(?P<component>[A-Za-z0-9]+)$'
    >>> import re
    >>> re.match(suffix_pattern("___"), "NZ_TINY02000001.1___tinyEcDub")["component"]
    'tinyEcDub'
    """
    _check_separator(separator)
    return rf"^(?P<chromosome>.+){separator}(?P<component>[A-Za-z0-9]+)$"


def check_roundtrip(chromosomes: Mapping[str, Sequence[str]], separator: str) -> None:
    """Prove every name a chimera would write reads back, before anything is written.

    Two things over the whole set. That ``separator`` is the one these components derive —
    a separator shorter than they need is invisible to the round trip, because splitting
    at the last run is correct with any separator, so nothing else would catch it. And
    that suffixing each component's chromosome and splitting the result returns exactly
    the pair it started from.

    Cheap enough to be unconditional: a chimera of the two assemblies the lab ships is 94
    names.

    Parameters
    ----------
    chromosomes : Mapping[str, Sequence[str]]
        Each component assembly name mapped to that component's chromosome names.
    separator : str
        The separator the caller is about to write with — what :func:`derive_separator`
        returned for this same mapping.

    Raises
    ------
    ChimeraNamingError
        If the component set is illegal, the separator is not the one these components
        derive, or a name does not survive the round trip. The message names the
        chromosome and the component that failed.

    Examples
    --------
    >>> chromosomes = {"ce11": ["I", "MtDNA"], "ecHT115": ["NZ_SMTD01000001.1"]}
    >>> check_roundtrip(chromosomes, derive_separator(chromosomes))
    """
    derived = derive_separator(chromosomes)
    if separator != derived:
        raise ChimeraNamingError(
            f"separator {separator!r} is not the one these components derive ({derived!r}); "
            f"the separator is the shortest run of underscores, at least two, strictly "
            f"longer than the longest run these components' chromosome names carry. Pass "
            f"derive_separator(chromosomes) rather than a constant."
        )
    for component, names in chromosomes.items():
        for chromosome in names:
            spelled = suffixed(chromosome, component, separator)
            recovered = split_suffixed(spelled, separator)
            if recovered != (chromosome, component):
                raise ChimeraNamingError(
                    f"chromosome {chromosome!r} of component {component!r} does not "
                    f"survive suffixing: it is spelled {spelled!r} and reads back as "
                    f"{recovered!r}. Do not write this chimera; report the two names "
                    f"above, since a derived separator is meant to make this impossible."
                )


def _check_component_set(names: Sequence[str]) -> None:
    """Raise unless ``names`` is a legal component set: two or more, distinct, alphanumeric."""
    if len(names) < _MIN_COMPONENTS:
        listed = ", ".join(repr(name) for name in names) or "none"
        raise ChimeraNamingError(
            f"a chimera needs at least {_MIN_COMPONENTS} components, got {len(names)} "
            f"({listed}); a chimera of one is that assembly under another name and a "
            f"second copy of its bytes — open that assembly directly instead."
        )
    repeated = sorted(name for name, count in Counter(names).items() if count > 1)
    if repeated:
        listed = ", ".join(repr(name) for name in repeated)
        raise ChimeraNamingError(
            f"component assembly names must not repeat, got {listed} more than once; the "
            f"same component twice doubles its sequence and makes the suffix ambiguous. "
            f"List each component exactly once."
        )
    for name in names:
        _check_component(name)


def _check_component(name: str) -> None:
    """Raise unless ``name`` is a legal component assembly name."""
    if not _COMPONENT_RE.fullmatch(name):
        raise ChimeraNamingError(
            f"component assembly name {name!r} is not alphanumeric; a component's name "
            f"must match [A-Za-z0-9]+, which is what makes a chimera's derived name "
            f"unambiguous and lets a suffixed chromosome name be split back to it. "
            f"Register these bytes under an alphanumeric assembly name — and note that a "
            f"chimera, whose name always carries '_', can never be a component."
        )


def _check_separator(separator: str) -> None:
    """Raise unless ``separator`` is a run of two or more underscores."""
    if not _SEPARATOR_RE.fullmatch(separator):
        raise ChimeraNamingError(
            f"separator {separator!r} is not a run of two or more underscores; a single "
            f"underscore cannot announce a suffix, since a real name like "
            f"'chr1_KI270706v1_random' already carries one. Derive it with "
            f"derive_separator(chromosomes)."
        )


def _longest_underscore_run(names: Iterable[str]) -> int:
    """Return the length of the longest run of underscores across ``names``, 0 when none."""
    return max((len(run) for name in names for run in _UNDERSCORE_RUN_RE.findall(name)), default=0)
