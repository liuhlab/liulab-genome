"""Which shipped table answers for an **Assembly**'s species, and the three ways none can.

Both gene-keyed halves of this context choose their shipped file the same way: the
**Assembly**'s own curated metadata row names a species, and that species names the file.
The choice is never a caller's — asking for human transcription factors while holding a
mouse assembly is not expressible (ADR-0003) — so what a caller gets when the choice
cannot be made is an error, and the three of them live here because they are one fact told
three ways rather than one half's business.

**Nothing published, and nothing said, are different facts.**
:class:`NoTFCensusError` and :class:`NoCofactorTableError` are about the literature: a
species nobody has assessed is unanswered here rather than answered with none, and the
message names the species somebody *has* assessed, since asking about one of those is what
a caller can do instead. :class:`UnknownSpeciesError` is about the assembly: nothing says
what species it is, so no file could be chosen for it whatever ships. All three are
:class:`LookupError`s, so a caller indifferent to which one it was catches that and a
caller who cares still tells them apart.

**And the two halves do not raise for the same assemblies.** Worm has a **Cofactor table**
and no census, so ``ce11`` is answered by one half and refused by the other — the
publishers' shape, and not a defect here.

Examples
--------
>>> from genome.tf.species import NoCofactorTableError, NoTFCensusError, UnknownSpeciesError
>>> [issubclass(error, LookupError) for error in
...  (NoTFCensusError, NoCofactorTableError, UnknownSpeciesError)]
[True, True, True]
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

#: How many species an error lists before saying how many it left out. Nothing ships for
#: enough species today for it to bind, and a message that grew past it unbounded would be
#: one nobody reads.
_MAX_LISTED_SPECIES = 10


def _elide(names: Sequence[str], limit: int = _MAX_LISTED_SPECIES) -> str:
    """Return ``names`` comma-joined, cut to ``limit`` and counting what was cut."""
    listed = ", ".join(names[:limit])
    hidden = len(names) - limit
    return listed if hidden <= 0 else f"{listed} (and {hidden} more)"


class NoTFCensusError(LookupError):
    """No census has been published for this **Assembly**'s species, so none can answer.

    The first of the two absences a **TF gene list** has, and the one a caller must never
    read as *this species has no transcription factors*: nobody has published a census
    for it, which is a fact about the literature and not about the genome. Worm, yeast and
    *E. coli* land here today. The message names the species that do have one, since asking
    about one of those is the thing a caller can do instead.

    A :class:`LookupError`, so it may be caught together with :class:`UnknownSpeciesError`
    and still told apart — exactly as the **Curated gene list**'s two absences are.

    Parameters
    ----------
    assembly : str
        The **Assembly** asked about.
    species : str
        The species its metadata row names, which is what no census ships for.
    censused : iterable of str
        The species a census does ship for.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    species : str
        Its species.
    censused : tuple of str
        The species that do have a census.

    Examples
    --------
    >>> try:
    ...     raise NoTFCensusError("ce11", "Caenorhabditis elegans", ["Homo sapiens"])
    ... except LookupError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """

    def __init__(self, assembly: str, species: str, censused: Iterable[str]) -> None:
        self.assembly = assembly
        self.species = species
        self.censused: tuple[str, ...] = tuple(censused)
        super().__init__(
            f"no TF census ships for {species!r}, which is the species the assembly table "
            f"names for {assembly!r} — so which of its genes are transcription factors is "
            f"unanswered here rather than answered with none. Censuses ship for: "
            f"{_elide(self.censused) or '(none)'}. Answering for {species!r} means shipping a "
            f"census for it, which scripts/build_tf_census.py writes."
        )


class NoCofactorTableError(LookupError):
    """Nobody has listed this **Assembly**'s species' cofactors, so no table can answer.

    :class:`NoTFCensusError`'s counterpart on the cofactor half, and the same kind of
    fact: nobody has published a cofactor list for this species, which is about the
    literature and not about the genome, and must never be read as *this species has no
    transcription cofactors*. Yeast and *E. coli* land here today.

    **Worm does not**, although :class:`NoTFCensusError` raises for the same assembly. A
    publisher assessed worm cofactors and none has released a worm TF census, so the two
    halves answer differently for one species — the publishers' shape, and not a defect.

    A :class:`LookupError`, so it may be caught together with :class:`UnknownSpeciesError`
    and still told apart, exactly as the census pair is.

    Parameters
    ----------
    assembly : str
        The **Assembly** asked about.
    species : str
        The species its metadata row names, which is what no table ships for.
    shipped_for : iterable of str
        The species a **Cofactor table** does ship for.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    species : str
        Its species.
    shipped_for : tuple of str
        The species that do have a cofactor table.

    Examples
    --------
    >>> try:
    ...     raise NoCofactorTableError("sacCer3", "Saccharomyces cerevisiae", ["Mus musculus"])
    ... except LookupError as error:
    ...     print("Mus musculus" in str(error))
    True
    """

    def __init__(self, assembly: str, species: str, shipped_for: Iterable[str]) -> None:
        self.assembly = assembly
        self.species = species
        self.shipped_for: tuple[str, ...] = tuple(shipped_for)
        super().__init__(
            f"no cofactor table ships for {species!r}, which is the species the assembly "
            f"table names for {assembly!r} — so which of its genes are transcription "
            f"cofactors is unanswered here rather than answered with none. Cofactor tables "
            f"ship for: {_elide(self.shipped_for) or '(none)'}. Answering for {species!r} "
            f"means shipping a table for it, which scripts/build_tf_cofactor.py writes."
        )


class UnknownSpeciesError(LookupError):
    """Nothing says what species this **Assembly** is, so no shipped table can be chosen.

    The second absence both gene-keyed halves have, and a different fact from
    :class:`NoTFCensusError` and :class:`NoCofactorTableError`: the question was not *has
    anyone published for this species* but *which species is this*, and nothing answered
    it. Two ways in, and neither is a mistake — a **Chimera** is more than one species by
    construction and nothing published answers for one, and an assembly the curated table
    does not list carries no species at all, which is the ordinary state of a free-form
    local key.

    The species is read from the assembly's own metadata and never passed in, so this is
    what a caller gets instead of quietly being handed another species' answer.

    One class for both halves because the fact is one fact: only what a caller can ask
    about instead differs, and ``shipped_table`` is what says which of them was asked for.

    Parameters
    ----------
    assembly : str
        The **Assembly** whose species nothing names.
    shipped_for : iterable of str
        The species the table asked for does ship for.
    shipped_table : str
        What could not be chosen, named in the message — ``"TF census"`` or ``"cofactor
        table"``. Keyword-only and required, so a message never says census over an
        answer that was about cofactors.

    Attributes
    ----------
    assembly : str
        The assembly asked about.
    shipped_for : tuple of str
        The species the table asked for does ship for.
    shipped_table : str
        What could not be chosen.

    Examples
    --------
    >>> try:
    ...     raise UnknownSpeciesError(
    ...         "ce11_ecHT115", ["Homo sapiens"], shipped_table="TF census"
    ...     )
    ... except LookupError as error:
    ...     print("Homo sapiens" in str(error))
    True
    """

    def __init__(self, assembly: str, shipped_for: Iterable[str], *, shipped_table: str) -> None:
        self.assembly = assembly
        self.shipped_table = shipped_table
        self.shipped_for: tuple[str, ...] = tuple(shipped_for)
        super().__init__(
            f"nothing says what species {assembly!r} is, so no {shipped_table} can be chosen "
            f"for it. The species comes from the assembly's own metadata row and is never "
            f"passed in, which is what makes asking about one species while holding another's "
            f"assembly impossible. A chimera is more than one species and no {shipped_table} "
            f"answers for one; anything else needs a row naming its species in the assembly "
            f"metadata table. One ships for: {_elide(self.shipped_for) or '(none)'}."
        )
