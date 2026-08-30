"""Tests for genome.tf.link — the join between a TF gene and the motifs that answer for it.

External behaviour only: what the answer contains, what order it arrives in, which error
is raised, and whether its message names the next action. Nothing here reaches for a
private helper or asserts an internal shape — the malformed cases go through
``parse_motif_link_table``, which is public precisely so that every way a shipped file can
be wrong is reachable without writing a broken one into the package.

The shipped tables answer here; no fixture stands in for them, and nothing touches the
network or the **Data dir**. ``tests/test_tf_link.py`` guards the *files* and pins their
counts — this module is about the API over them, so the two overlap nowhere on purpose.

Most tests that are about one table run against a representative pair — one human table
and one mouse table, on different releases — rather than the full cross-product: species
and release are covariates of the same code path here, and the pinned-count tests in
``tests/test_tf_link.py`` are what actually guards every shipped file.
"""

from __future__ import annotations

import pytest

from genome.tf import (
    COMPLEX,
    LINK_COLUMNS,
    LINK_TAX_GROUP,
    MONOMER,
    GeneNotAssessedError,
    MotifLink,
    MotifLinkTable,
    MotifLinkTableError,
    NoMotifLinkTableError,
    TranscriptionCofactorError,
    VersionedGeneIdError,
    motif_link_table,
    motif_links,
    parse_motif_link_table,
    shipped_link_tables,
)
from genome.tf.gene import tf_gene_table
from genome.tf.motif.jaspar import JASPAR_RELEASES
from tests._tables import table_text

#: Every shipped table, as ``(species slug, release)``. Discovered rather than listed, so
#: adding a species or a release is a file drop rather than an edit here.
_TABLES = shipped_link_tables()

#: One human table and one mouse table, on different releases: enough to cross the species
#: boundary and the release boundary at once, without re-running every table-level check
#: against the full cross-product that :mod:`tests.test_tf_link` already pins.
_REPRESENTATIVE_TABLES = (("homo_sapiens", "2026"), ("mus_musculus", "2024"))

#: A gene whose motifs are all monomers, and the ids they answer with on 2026. Named
#: because it is the plainest possible answer: three matrices of one factor, most
#: specifically attributable first.
_CTCF = ("CTCF", "Homo sapiens", ("MA1930.2", "MA1929.2", "MA0139.2"))

#: The AP-1 case, which is what **Attribution specificity** was written for. JASPAR's
#: canonical AP-1 matrix is the *complex* ``MA0099.4 FOS::JUN``, it describes JUN's
#: binding better than either JUN monomer does, and it still ranks below both — because
#: the order says what a matrix is attributable to and explicitly not which is better.
_JUN_MONOMER, _JUN_CROSS_SPECIES, _AP1 = "MA0488.2", "MA0489.3", "MA0099.4"

#: A gene whose only motif is a complex. A join that read ``Ahr::Arnt`` as an AHR monomer
#: would report it motif-less, which is the failure **Role** exists to make impossible.
_COMPLEX_ONLY = ("AHR", "Homo sapiens", "MA0006.2")

#: Mouse CTCF: every profile JASPAR has for it was measured on human, so asking for
#: species-matched profiles empties the answer. The case that proves an empty answer still
#: says which **Release** found nothing.
_ALL_CROSS_SPECIES = ("Ctcf", "Mus musculus")

#: A gene Lambert assessed and turned down. Unlinked because only assessed-positive genes
#: receive links — not because anything is missing.
_ASSESSED_NEGATIVE = ("SMAD2", "Homo sapiens")

#: A gene Lambert judged a transcription factor that JASPAR has no profile for. The same
#: empty answer as the one above, and a different fact — 763 of human's 1,639 positives
#: are this on the 2026 release.
_NO_PROFILE = ("TFAP2D", "Homo sapiens")

#: Lambert's census carries clone-style symbols shaped exactly like a versioned gene id.
#: This one is why a symbol is looked for before a versioned id is diagnosed.
_CLONE_STYLE_SYMBOL = ("AC023509.3", "Homo sapiens", "ENSG00000267281")

#: One gene, in the three spellings the entry point takes.
_JUN_STEM, _JUN_SYMBOL, _JUN_LOWER = "ENSG00000177606", "JUN", "jun"

#: A gene id of one species handed to another species' census. Absent from that census,
#: which is a different answer from judged not to be a TF.
_FOREIGN_STEM = "ENSMUSG00000005698"

#: The dual-classified case: a gene Lambert judged a transcription factor that a publisher
#: also lists as a **Transcription cofactor**. 151 human genes are both, and TBP is the one
#: of them with a profile — so it is what proves a second table never suppresses an answer
#: the census already reached.
_DUAL_CLASSIFIED = ("TBP", "Homo sapiens", ("MA0108.3",))

#: Two more of the 151, judged transcription factors and linked to no profile on this
#: release. They are the pair a lookup that asked the cofactor table first would break:
#: their answer is empty, and an empty answer is not an absence.
_DUAL_UNLINKED = ("KMT2A", "DNMT1")

#: A cofactor Lambert assessed and turned *down*. Assessed is assessed, so it answers with
#: the census's verdict and never with the cofactor error.
_ASSESSED_NEGATIVE_COFACTOR = ("EP300", "Homo sapiens")

#: Cofactors no census assessed, with the **Gene id stem** one of them is also named by and
#: the publishers each is listed by. Human's table is built from two publishers, so a gene
#: both of them list names both and a gene one lists names one — which is what the message
#: has to get right. Mouse's table is one publisher's throughout, and that publisher is also
#: mouse's census: it assessed the gene as a cofactor and never as a transcription factor,
#: which is exactly the sentence the message makes.
_COFACTOR_BOTH = ("WDR5", "Homo sapiens", "AnimalTFDB and EpiFactors")
_COFACTOR_BOTH_STEM = "ENSG00000196363"
_COFACTOR_ONE = ("CCNC", "Homo sapiens", "AnimalTFDB")
_MOUSE_COFACTOR = ("Smarcb1", "Mus musculus", "AnimalTFDB")

#: A name no census assessed and no publisher lists as a cofactor, in any species.
_KNOWN_TO_NOBODY = "not-a-gene-in-any-table"

#: One well-formed row, and the header above it. Built by the tests that are about a
#: malformed one, so that every way a table can be wrong is one changed cell.
_ROW = (
    "2026",
    "Homo sapiens",
    "ENSG00000141510",
    "TP53",
    "MA0106.3",
    "TP53",
    "monomer",
    "",
    "9606",
    "no",
    "20.6607",
    "1",
)


def _table_text(*rows: tuple[str, ...]) -> str:
    """Return a link table's text: the twelve columns, then one line per row."""
    return table_text(LINK_COLUMNS, *rows)


def _row(**changes: str) -> tuple[str, ...]:
    """Return :data:`_ROW` with the named columns replaced."""
    cells = dict(zip(LINK_COLUMNS, _ROW, strict=True))
    return tuple((changes.get(column, cell)) for column, cell in cells.items())


def _specificity_key(link: MotifLink) -> tuple[int, int, float, str]:
    """Return one link's place under **Attribution specificity**, lowest first.

    Derived from the attributes the link carries and nothing else — which is the claim
    being checked, and also what a caller who disagrees with the order would re-sort on.
    """
    return (
        0 if link.role == MONOMER else 1,
        1 if link.is_cross_species else 0,
        -link.total_information_content,
        link.motif_id,
    )


def _shipped(slug: str, release: str) -> MotifLinkTable:
    """Return one shipped table, failing by name when none ships."""
    table = motif_link_table(slug, release)
    assert table is not None, f"no link table ships for {slug} {release}"
    return table


# ---------------------------------------------------------------------------------------
# A gene answers with its links
# ---------------------------------------------------------------------------------------


def test_the_package_ships_link_tables_at_all() -> None:
    # The guard under every parametrized test below: with no tables each would collect
    # zero cases and pass, which is the silent zero this module must not have.
    assert _TABLES


def test_a_gene_answers_with_the_motifs_that_link_to_it() -> None:
    gene, species, expected = _CTCF

    assert motif_links(gene, species).motif_ids == expected


def test_a_link_says_what_the_matrix_is_a_motif_of_and_what_it_was_measured_on() -> None:
    # Everything one link carries, on one row whose values are known: the ids, the role
    # and its partners, the profile's tax ids, the species flag and the information
    # content. A caller who disagrees with the order re-sorts on exactly these.
    (link,) = motif_links("TP53", "Homo sapiens")

    assert (link.motif_id, link.motif_name) == ("MA0106.3", "TP53")
    assert (link.role, link.partners) == (MONOMER, ())
    assert (link.motif_tax_ids, link.is_cross_species) == (("9606",), False)
    assert link.total_information_content == pytest.approx(20.6607)
    assert (link.gene_id_stem, link.symbol) == ("ENSG00000141510", "TP53")
    assert (link.release, link.species) == ("2026", "Homo sapiens")


def test_a_complex_names_its_partners_and_a_monomer_names_none() -> None:
    # The whole reason **Role** exists: a heterodimer matrix must never be read as a
    # monomer's. `FOS::JUN` answers for JUN as a complex naming FOS.
    links = {link.motif_id: link for link in motif_links("JUN", "Homo sapiens")}

    assert (links[_JUN_MONOMER].role, links[_JUN_MONOMER].partners) == (MONOMER, ())
    assert (links[_AP1].role, links[_AP1].partners) == (COMPLEX, ("FOS",))
    assert links[_AP1].is_complex
    assert not links[_JUN_MONOMER].is_complex


def test_a_gene_whose_only_motifs_are_complexes_is_still_linked() -> None:
    # AHR has no monomer matrix at all. A join that read `Ahr::Arnt` as an AHR motif would
    # report it motif-less; a join that dropped complexes would too.
    gene, species, motif_id = _COMPLEX_ONLY
    (link,) = motif_links(gene, species)

    assert (link.motif_id, link.role, link.partners) == (motif_id, COMPLEX, ("ARNT",))


# ---------------------------------------------------------------------------------------
# Attribution specificity: the order the links arrive in, and how filtering treats it
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _REPRESENTATIVE_TABLES)
def test_every_genes_links_are_complete_ordered_stable_and_filterable(
    slug: str, release: str
) -> None:
    # Four claims about one gene's answer, checked together because they all walk the
    # same `table.links_for(stem)` for every gene the table carries: every gene answers
    # with at least one link, the links arrive in **Attribution specificity** order,
    # asking twice gives the same order back with the rank each row was shipped with, and
    # dropping cross-species links removes some without ever adding or reordering one.
    table = _shipped(slug, release)

    for stem in table.gene_id_stems:
        links = table.links_for(stem)
        ordered = list(links)

        assert links.motif_ids
        assert ordered == sorted(ordered, key=_specificity_key), stem
        assert links.links == table.links_for(stem).links
        assert [link.rank for link in ordered] == sorted(link.rank for link in ordered)

        kept = table.links_for(stem, cross_species=False).links
        assert list(kept) == [link for link in ordered if not link.is_cross_species]


def test_the_canonical_ap1_ranks_below_a_jun_monomer_and_a_caller_can_resort_it() -> None:
    # The order states what a matrix is attributable to and explicitly not which motif is
    # better: MA0099.4 FOS::JUN is the canonical AP-1 matrix and describes JUN's binding
    # better than either JUN monomer, and it still comes after both because it is a motif
    # of a complex. The documented escape is that no attribute is hidden behind the order,
    # so a caller who wants the most informative matrix regardless of what it is a motif
    # of can re-sort on what the links already carry.
    links = motif_links("JUN", "Homo sapiens")
    ids = links.motif_ids

    assert ids.index(_AP1) > ids.index(_JUN_MONOMER)
    assert ids.index(_AP1) > ids.index(_JUN_CROSS_SPECIES)

    by_information = sorted(links, key=lambda link: -link.total_information_content)
    reordered = tuple(link.motif_id for link in by_information)
    assert reordered[0] == _JUN_MONOMER
    assert reordered != ids


# ---------------------------------------------------------------------------------------
# Cross-species links are kept and marked, and the caller may drop them (ADR-0013)
# ---------------------------------------------------------------------------------------


def test_cross_species_links_are_kept_and_marked_by_default() -> None:
    gene, species = _ALL_CROSS_SPECIES
    links = motif_links(gene, species)

    assert links.motif_ids
    assert all(link.is_cross_species for link in links)


def test_a_caller_can_ask_for_species_matched_profiles_and_keeps_the_shipped_ranks() -> None:
    kept = motif_links("Jun", "Mus musculus", cross_species=False)

    assert kept.motif_ids == (_JUN_CROSS_SPECIES,)
    assert not any(link.is_cross_species for link in kept)

    # The rank is the row's own place among *all* of this gene's links, not its position
    # in whatever survived a filter — so the gaps are the point. A rank that renumbered
    # itself would stop saying how specific the attribution was.
    human_kept = motif_links("JUN", "Homo sapiens", cross_species=False)
    assert [link.rank for link in human_kept][:3] == [1, 3, 4]


# ---------------------------------------------------------------------------------------
# Provenance is captured before any filtering
# ---------------------------------------------------------------------------------------


def test_an_answer_says_which_release_tax_group_and_file_it_came_from() -> None:
    links = motif_links("CTCF", "Homo sapiens", release="2024")

    assert (links.species, links.release, links.tax_group) == (
        "Homo sapiens",
        "2024",
        LINK_TAX_GROUP,
    )
    assert links.source.endswith("homo_sapiens.jaspar2024.motif_link_table.tsv.gz")


def test_an_answer_a_filter_emptied_still_says_which_release_found_nothing() -> None:
    # The reason provenance is captured before the filter runs: no link survives to carry
    # the release on its row, and an answer that could not say which release it asked
    # would be indistinguishable from one asked of the other.
    gene, species = _ALL_CROSS_SPECIES
    links = motif_links(gene, species, release="2024", cross_species=False)

    assert len(links) == 0
    assert (links.release, links.tax_group, links.species) == (
        "2024",
        LINK_TAX_GROUP,
        "Mus musculus",
    )
    assert (links.gene_id_stem, links.symbol, links.is_tf) == ("ENSMUSG00000005698", "Ctcf", True)


# ---------------------------------------------------------------------------------------
# Pinning a release, and asking for one or a tax group that has no table
# ---------------------------------------------------------------------------------------


def test_every_release_the_tables_cover_can_be_pinned() -> None:
    for release in sorted({release for _, release in _TABLES}):
        links = motif_links("CTCF", "Homo sapiens", release=release)
        assert links.release == release
        assert all(link.release == release for link in links)


def test_the_releases_that_ship_are_the_ones_the_package_prepares() -> None:
    # Covering only the default release would be a hole discovered at runtime, by whoever
    # pinned the other one for reproducibility.
    assert sorted({release for _, release in _TABLES}) == sorted(JASPAR_RELEASES)


def test_pinning_a_release_changes_the_answer_rather_than_being_ignored() -> None:
    # FOXM1 gains its first profile in 2026. If the release argument were cosmetic the two
    # calls would agree, which is exactly the silent hole a per-release table exists to
    # close.
    assert len(motif_links("FOXM1", "Homo sapiens", release="2024")) == 0
    assert motif_links("FOXM1", "Homo sapiens", release="2026").motif_ids == ("MA2469.1",)


def test_no_table_errors_name_what_does_ship_for_release_tax_group_and_species() -> None:
    # Three ways to ask for a table nothing built: a release, a tax group, and a species
    # each raise the same error and each names what does ship instead.
    with pytest.raises(NoMotifLinkTableError) as by_release:
        motif_links("CTCF", "Homo sapiens", release="2020")
    assert all(release in str(by_release.value) for release in JASPAR_RELEASES)

    # Plants is answered, not silently empty: no table was ever built for it, which is not
    # the same fact as this gene having no motifs.
    with pytest.raises(NoMotifLinkTableError) as by_tax_group:
        motif_links("CTCF", "Homo sapiens", tax_group="plants")
    assert LINK_TAX_GROUP in str(by_tax_group.value)

    with pytest.raises(NoMotifLinkTableError) as by_species:
        motif_links("daf-16", "Caenorhabditis elegans")
    assert "homo_sapiens" in str(by_species.value)


@pytest.mark.parametrize(
    ("species", "release", "tax_group"),
    [
        ("Danio rerio", "2026", LINK_TAX_GROUP),
        ("Homo sapiens", "2020", LINK_TAX_GROUP),
        ("Homo sapiens", "2026", "plants"),
    ],
)
def test_the_loader_below_says_that_absence_with_none(
    species: str, release: str, tax_group: str
) -> None:
    # The one place absence is spelled `None`, exactly as the censuses and the curated
    # gene lists spell it: this is the layer below the one a caller touches.
    assert motif_link_table(species, release, tax_group) is None


def test_a_species_name_shaped_like_a_path_reaches_nothing() -> None:
    assert motif_link_table("../tf_gene/homo_sapiens") is None


# ---------------------------------------------------------------------------------------
# The other absence: a gene the census never assessed
# ---------------------------------------------------------------------------------------


def test_a_gene_no_census_assessed_raises_rather_than_answering_emptily() -> None:
    # A mouse gene id handed to the human census. Absent from a census is not the same
    # answer as judged not to be a TF, and neither is an empty tuple.
    with pytest.raises(GeneNotAssessedError) as raised:
        motif_links(_FOREIGN_STEM, "Homo sapiens")

    message = str(raised.value)
    assert "Lambert et al. 2018" in message
    assert "tf_gene_table" in message


def test_a_gene_its_census_turned_down_answers_with_no_links_and_says_so() -> None:
    gene, species = _ASSESSED_NEGATIVE
    links = motif_links(gene, species)

    assert len(links) == 0
    assert links.is_tf is False


def test_a_gene_its_census_judged_a_tf_and_jaspar_has_no_profile_for_answers_the_same_way() -> None:
    # The same empty answer as the one above and a different fact, which is why the
    # verdict is on the answer: without it the two would be one indistinguishable silence.
    gene, species = _NO_PROFILE
    links = motif_links(gene, species)

    assert len(links) == 0
    assert links.is_tf is True


def test_the_two_absences_are_told_apart_by_the_error_they_raise() -> None:
    # Both are LookupErrors, so a caller may catch the pair and still know which happened.
    with pytest.raises(LookupError) as no_table:
        motif_links("CTCF", "Danio rerio")
    with pytest.raises(LookupError) as no_gene:
        motif_links(_FOREIGN_STEM, "Homo sapiens")

    assert isinstance(no_table.value, NoMotifLinkTableError)
    assert isinstance(no_gene.value, GeneNotAssessedError)
    assert not isinstance(no_table.value, GeneNotAssessedError)


# ---------------------------------------------------------------------------------------
# A transcription cofactor, told apart from a gene nothing assessed
# ---------------------------------------------------------------------------------------


def test_a_dual_classified_gene_answers_from_its_census_verdict_whether_linked_or_not() -> None:
    # The property the lookup order exists for: the census is asked first, so the 151
    # human genes that are both a TF and a **Transcription cofactor** keep whatever their
    # census verdict was. TBP is the one of them with a profile; KMT2A and DNMT1 are two
    # more that JASPAR has none for — their answer is empty and their verdict positive, an
    # empty answer being a real answer and never one of the absences.
    gene, species, expected = _DUAL_CLASSIFIED
    links = motif_links(gene, species)
    assert links.motif_ids == expected
    assert links.is_tf is True

    for gene in _DUAL_UNLINKED:
        links = motif_links(gene, "Homo sapiens")
        assert len(links) == 0
        assert links.is_tf is True


def test_a_cofactor_the_census_turned_down_answers_with_that_verdict() -> None:
    # Assessed is assessed, whichever way the verdict went: a gene the census looked at and
    # rejected comes back with its verdict rather than with anything the cofactor table says.
    gene, species = _ASSESSED_NEGATIVE_COFACTOR
    links = motif_links(gene, species)

    assert len(links) == 0
    assert links.is_tf is False


def test_a_cofactor_no_census_assessed_says_so_rather_than_that_nothing_knows_it() -> None:
    # The whole point: *no census assessed this* reads as *nothing here knows this gene*,
    # where the honest answer is that a cofactor recognises no sequence of its own and there
    # is no motif to look for. The message carries all three facts and the next action.
    gene, species, publishers = _COFACTOR_BOTH
    with pytest.raises(TranscriptionCofactorError) as raised:
        motif_links(gene, species)

    message = str(raised.value)
    assert "Lambert et al. 2018" in message
    assert publishers in message
    assert "no motif" in message
    assert f"cofactor_table({species!r})" in message


@pytest.mark.parametrize(
    ("gene", "species", "publishers"), [_COFACTOR_BOTH, _COFACTOR_ONE, _MOUSE_COFACTOR]
)
def test_the_message_names_the_publishers_that_list_the_gene(
    gene: str, species: str, publishers: str
) -> None:
    # Whoever listed it, and only them: both of human's publishers for a gene both list, one
    # for a gene one lists. Mouse's census and its cofactor table are the same publisher,
    # which is a sentence about two of its lists rather than a contradiction.
    with pytest.raises(TranscriptionCofactorError) as raised:
        motif_links(gene, species)

    assert f"transcription cofactor by {publishers}." in str(raised.value)


def test_a_cofactor_error_is_reachable_by_stem_and_still_caught_by_the_older_except_clause() -> (
    None
):
    # The cofactor table is keyed by stem and carries the symbol beside it, as the census
    # does, so both spellings of one gene reach the same answer. And the reason the error
    # subclasses rather than sitting beside the older one: a caller who wrote
    # `except GeneNotAssessedError` before any cofactor table shipped keeps covering every
    # gene they covered, and the narrower type is what tells them which absence this is.
    gene, species, publishers = _COFACTOR_BOTH
    with pytest.raises(TranscriptionCofactorError) as raised:
        motif_links(_COFACTOR_BOTH_STEM, species)
    assert publishers in str(raised.value)

    with pytest.raises(GeneNotAssessedError) as older:
        motif_links(gene, species)
    assert isinstance(older.value, TranscriptionCofactorError)


@pytest.mark.parametrize(("slug", "release"), _REPRESENTATIVE_TABLES)
def test_a_gene_neither_table_knows_raises_the_error_it_always_did(slug: str, release: str) -> None:
    # With nobody listing the gene the lookup falls through to the census's own silence,
    # unnarrowed by either table.
    with pytest.raises(GeneNotAssessedError) as raised:
        motif_links(_KNOWN_TO_NOBODY, slug, release=release)

    assert type(raised.value) is GeneNotAssessedError
    assert "never assessed" in str(raised.value)


# ---------------------------------------------------------------------------------------
# How a gene is named
# ---------------------------------------------------------------------------------------


def test_a_gene_answers_to_its_stem_symbol_and_lowercased_symbol() -> None:
    for gene in (_JUN_STEM, _JUN_SYMBOL, _JUN_LOWER):
        links = motif_links(gene, "Homo sapiens")
        assert links.gene_id_stem == _JUN_STEM
        assert links.symbol == _JUN_SYMBOL


def test_each_census_answers_to_its_own_spelling_of_one_factors_symbol() -> None:
    # The two censuses spell one factor `JUN` and `Jun`, and each spelling is its
    # publisher's. The symbol on the answer is the census's own, never JASPAR's.
    assert motif_links("Jun", "Mus musculus").symbol == "Jun"
    assert motif_links("JUN", "Homo sapiens").symbol == "JUN"


def test_a_versioned_gene_id_is_refused_by_value_error_naming_the_stem_and_is_not_an_absence() -> (
    None
):
    # A stem may name more than one gene id in one annotation, so answering the versioned
    # id would answer for the stem — which names a gene the caller did not. And nothing is
    # missing — the gene is assessed and its links are here — so a caller catching
    # LookupError for a gene that is not there does not swallow this.
    with pytest.raises(VersionedGeneIdError) as raised:
        motif_links(f"{_JUN_STEM}.6", "Homo sapiens")
    assert _JUN_STEM in str(raised.value)

    with pytest.raises(ValueError, match="gene id stem") as value_error:
        motif_links(f"{_JUN_STEM}.6", "Homo sapiens")
    assert not isinstance(value_error.value, LookupError)


def test_a_clone_style_symbol_is_read_as_a_symbol_and_not_as_a_versioned_gene_id() -> None:
    # Lambert's census carries symbols shaped exactly like a versioned gene id, so a
    # symbol is looked for before a versioned id is diagnosed.
    symbol, species, stem = _CLONE_STYLE_SYMBOL
    links = motif_links(symbol, species)

    assert (links.gene_id_stem, links.symbol) == (stem, symbol)


# ---------------------------------------------------------------------------------------
# The whole table, and the frame a collaborator reads it as
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _REPRESENTATIVE_TABLES)
def test_a_table_names_itself_lists_only_the_censuss_positives_and_frames_cleanly(
    slug: str, release: str
) -> None:
    table = _shipped(slug, release)
    census = tf_gene_table(slug)
    assert census is not None

    assert (table.release, table.tax_group) == (release, LINK_TAX_GROUP)
    assert table.species
    assert table.source.endswith(f"{slug}.jaspar{release}.motif_link_table.tsv.gz")
    # Only assessed-positive genes receive links, which is the census speaking and never
    # this package.
    assert set(table.gene_id_stems) <= set(census.assessed_positive)

    frame = table.frame()
    assert list(frame.columns) == list(LINK_COLUMNS)
    assert len(frame) == len(table)


def test_the_frame_is_rebuilt_fresh_and_a_table_is_reachable_however_species_is_spelled() -> None:
    table = _shipped("homo_sapiens", "2026")
    frame = table.frame()
    frame.loc[0, "motif_id"] = "MA9999.9"
    assert table.frame().loc[0, "motif_id"] != "MA9999.9"

    assert motif_link_table("homo_sapiens", "2026") == motif_link_table("Homo sapiens", "2026")


# ---------------------------------------------------------------------------------------
# A shipped file that cannot be trusted never answers
# ---------------------------------------------------------------------------------------


def test_a_well_formed_table_reads_back_as_its_links() -> None:
    table = parse_motif_link_table(_table_text(_ROW), source="one.tsv")

    assert (table.species, table.release, table.tax_group) == (
        "Homo sapiens",
        "2026",
        LINK_TAX_GROUP,
    )
    assert table.links_for("TP53").motif_ids == ("MA0106.3",)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "\t".join(LINK_COLUMNS) + "\n",
    ],
    ids=["empty", "header-only"],
)
def test_a_file_that_is_not_a_link_table_raises_and_names_the_file(text: str) -> None:
    with pytest.raises(MotifLinkTableError, match=r"broken\.tsv"):
        parse_motif_link_table(text, source="broken.tsv")


@pytest.mark.parametrize(
    "changes",
    [
        {"role": "dimer"},
        {"partners": "FOS"},
        {"is_cross_species": "true"},
        {"rank": "first"},
        {"total_information_content": "high"},
        {"gene_id_stem": ""},
    ],
    ids=[
        "a-third-role",
        "a-monomer-naming-one",
        "a-flag-no-census-spells",
        "a-rank-that-is-not-a-number",
        "an-information-content-that-is-not-a-number",
        "no-gene-id-stem",
    ],
)
def test_a_row_a_generator_would_never_write_raises_and_names_its_line(
    changes: dict[str, str],
) -> None:
    with pytest.raises(MotifLinkTableError, match="line 2"):
        parse_motif_link_table(_table_text(_row(**changes)), source="broken.tsv")


def test_a_short_row_raises_a_bad_value_and_not_a_lookup() -> None:
    with pytest.raises(MotifLinkTableError, match="cells"):
        parse_motif_link_table(_table_text(_ROW[:-1]), source="broken.tsv")

    # A packaging defect, not an absence: a caller catching LookupError for a species with
    # no table must not swallow a file that ships broken.
    with pytest.raises(ValueError) as raised:  # noqa: PT011 - the class is the assertion
        parse_motif_link_table("nonsense\n", source="broken.tsv")
    assert isinstance(raised.value, MotifLinkTableError)
    assert not isinstance(raised.value, LookupError)


@pytest.mark.parametrize("changes", [{"release": "2024"}, {"species": "Mus musculus"}])
def test_one_file_naming_two_releases_or_two_species_raises(changes: dict[str, str]) -> None:
    # The release and the species are on every row so that two tables concatenate into one
    # frame that still says where each row came from — which promises nothing unless they
    # are uniform within a file.
    text = _table_text(_ROW, _row(motif_id="MA0106.4", **changes))

    with pytest.raises(MotifLinkTableError, match=r"broken\.tsv"):
        parse_motif_link_table(text, source="broken.tsv")


@pytest.mark.parametrize(
    "text",
    [
        "",
        _table_text(_ROW[:-1]),
        _table_text(_ROW, _row(release="2024", motif_id="MA0106.4")),
    ],
    ids=["empty", "a-short-row", "two-releases"],
)
def test_every_message_a_broken_file_raises_names_the_command_that_rewrites_it(text: str) -> None:
    # A shipped file that cannot be trusted is a packaging defect, and regenerating it is
    # the only repair — so every message names the file and the generator that writes it.
    with pytest.raises(MotifLinkTableError) as raised:
        parse_motif_link_table(text, source="broken.tsv")

    assert "broken.tsv" in str(raised.value)
    assert "scripts/build_tf_links.py" in str(raised.value)
