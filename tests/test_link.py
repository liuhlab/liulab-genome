"""Tests for genome.tf.link — the join between a TF gene and the motifs that answer for it.

External behaviour only: what the answer contains, what order it arrives in, which error
is raised, and whether its message names the next action. Nothing here reaches for a
private helper or asserts an internal shape — the malformed cases go through
``parse_motif_link_table``, which is public precisely so that every way a shipped file can
be wrong is reachable without writing a broken one into the package.

The shipped tables answer here; no fixture stands in for them, and nothing touches the
network or the **Data dir**. ``tests/test_tf_link.py`` guards the *files* and pins their
counts — this module is about the API over them, so the two overlap nowhere on purpose.

Every test that is about one table iterates over *every* shipped table. A second species
or a third release is a file dropped into ``data/tf_link/``, and dropping one in must not
mean rewriting these.
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
    VersionedGeneIdError,
    motif_link_table,
    motif_links,
    parse_motif_link_table,
    shipped_link_tables,
)
from genome.tf.gene import tf_gene_table
from genome.tf.motif.jaspar import JASPAR_RELEASES

#: Every shipped table, as ``(species slug, release)``. Discovered rather than listed, so
#: adding a species or a release is a file drop rather than an edit here.
_TABLES = shipped_link_tables()

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
    lines = ["\t".join(LINK_COLUMNS), *("\t".join(row) for row in rows)]
    return "\n".join(lines) + "\n"


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


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_every_gene_the_table_links_answers_with_at_least_one_link(slug: str, release: str) -> None:
    table = _shipped(slug, release)

    for stem in table.gene_id_stems:
        assert table.links_for(stem).motif_ids


# ---------------------------------------------------------------------------------------
# Attribution specificity: the order the links arrive in
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_links_arrive_in_attribution_specificity_order(slug: str, release: str) -> None:
    # Monomer before complex, species-matched before cross-species, then higher
    # information content, then motif id. Four keys, so the order is total: re-deriving it
    # from the attributes the links carry has to reproduce the order they arrived in.
    table = _shipped(slug, release)

    for stem in table.gene_id_stems:
        links = list(table.links_for(stem))

        assert links == sorted(links, key=_specificity_key), stem


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_the_order_is_stable_and_says_how_specific_each_attribution_is(
    slug: str, release: str
) -> None:
    # Stable: the same question answered twice is the same answer, in the same order, so
    # "the motif for this factor" means one thing on two machines. The rank rises with the
    # position, which is what makes it readable as *how specific* rather than *how good*.
    table = _shipped(slug, release)

    for stem in table.gene_id_stems:
        links = table.links_for(stem)

        assert links.links == table.links_for(stem).links
        assert [link.rank for link in links] == sorted(link.rank for link in links)


def test_the_canonical_ap1_complex_ranks_below_a_jun_monomer() -> None:
    # The order states what a matrix is attributable to and explicitly not which motif is
    # better: MA0099.4 FOS::JUN is the canonical AP-1 matrix and describes JUN's binding
    # better than either JUN monomer, and it still comes after both because it is a motif
    # of a complex.
    links = motif_links("JUN", "Homo sapiens")
    ids = links.motif_ids

    assert ids.index(_AP1) > ids.index(_JUN_MONOMER)
    assert ids.index(_AP1) > ids.index(_JUN_CROSS_SPECIES)


def test_a_caller_who_disagrees_can_re_sort_on_what_the_links_already_carry() -> None:
    # The documented escape: no attribute is hidden behind the order, so a caller who
    # wants the most informative matrix regardless of what it is a motif of has it.
    links = motif_links("JUN", "Homo sapiens")
    by_information = sorted(links, key=lambda link: -link.total_information_content)
    reordered = tuple(link.motif_id for link in by_information)

    assert reordered[0] == _JUN_MONOMER
    assert reordered != links.motif_ids


# ---------------------------------------------------------------------------------------
# Cross-species links are kept and marked, and the caller may drop them (ADR-0013)
# ---------------------------------------------------------------------------------------


def test_cross_species_links_are_kept_and_marked_by_default() -> None:
    gene, species = _ALL_CROSS_SPECIES
    links = motif_links(gene, species)

    assert links.motif_ids
    assert all(link.is_cross_species for link in links)


def test_a_caller_can_ask_for_species_matched_profiles_only() -> None:
    links = motif_links("Jun", "Mus musculus", cross_species=False)

    assert links.motif_ids == (_JUN_CROSS_SPECIES,)
    assert not any(link.is_cross_species for link in links)


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_dropping_cross_species_links_never_adds_or_reorders_one(slug: str, release: str) -> None:
    table = _shipped(slug, release)

    for stem in table.gene_id_stems:
        kept = table.links_for(stem, cross_species=False).links
        all_links = table.links_for(stem).links

        assert list(kept) == [link for link in all_links if not link.is_cross_species]


def test_a_filtered_answer_keeps_the_rank_the_shipped_table_gave_each_link() -> None:
    # The rank is the row's own place among *all* of this gene's links, not its position
    # in whatever survived a filter — so the gaps are the point. A rank that renumbered
    # itself would stop saying how specific the attribution was.
    kept = motif_links("JUN", "Homo sapiens", cross_species=False)

    assert [link.rank for link in kept][:3] == [1, 3, 4]


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


@pytest.mark.parametrize("release", sorted({release for _, release in _TABLES}))
def test_every_release_the_tables_cover_can_be_pinned(release: str) -> None:
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


def test_a_release_no_table_covers_raises_and_names_the_ones_that_ship() -> None:
    with pytest.raises(NoMotifLinkTableError) as raised:
        motif_links("CTCF", "Homo sapiens", release="2020")

    assert all(release in str(raised.value) for release in JASPAR_RELEASES)


def test_a_tax_group_with_no_table_raises_rather_than_answering_emptily() -> None:
    # Plants is answered, not silently empty: no table was ever built for it, which is not
    # the same fact as this gene having no motifs.
    with pytest.raises(NoMotifLinkTableError) as raised:
        motif_links("CTCF", "Homo sapiens", tax_group="plants")

    assert LINK_TAX_GROUP in str(raised.value)


def test_a_species_no_table_ships_for_raises_and_names_those_that_do() -> None:
    with pytest.raises(NoMotifLinkTableError) as raised:
        motif_links("daf-16", "Caenorhabditis elegans")

    assert "homo_sapiens" in str(raised.value)


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
# How a gene is named
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize("gene", [_JUN_STEM, _JUN_SYMBOL, _JUN_LOWER])
def test_a_gene_answers_to_its_stem_and_to_its_censuss_own_symbol(gene: str) -> None:
    links = motif_links(gene, "Homo sapiens")

    assert links.gene_id_stem == _JUN_STEM
    assert links.symbol == _JUN_SYMBOL


def test_each_census_answers_to_its_own_spelling_of_one_factors_symbol() -> None:
    # The two censuses spell one factor `JUN` and `Jun`, and each spelling is its
    # publisher's. The symbol on the answer is the census's own, never JASPAR's.
    assert motif_links("Jun", "Mus musculus").symbol == "Jun"
    assert motif_links("JUN", "Homo sapiens").symbol == "JUN"


def test_a_versioned_gene_id_is_refused_and_the_message_names_the_stem_to_pass() -> None:
    # A stem may name more than one gene id in one annotation, so answering the versioned
    # id would answer for the stem — which names a gene the caller did not.
    with pytest.raises(VersionedGeneIdError) as raised:
        motif_links(f"{_JUN_STEM}.6", "Homo sapiens")

    assert _JUN_STEM in str(raised.value)


def test_refusing_a_versioned_gene_id_is_not_one_of_the_absences() -> None:
    # Nothing is missing — the gene is assessed and its links are here — so a caller
    # catching LookupError for a gene that is not there does not swallow it.
    with pytest.raises(ValueError, match="gene id stem") as raised:
        motif_links(f"{_JUN_STEM}.6", "Homo sapiens")

    assert not isinstance(raised.value, LookupError)


def test_a_clone_style_symbol_is_read_as_a_symbol_and_not_as_a_versioned_gene_id() -> None:
    # Lambert's census carries symbols shaped exactly like a versioned gene id, so a
    # symbol is looked for before a versioned id is diagnosed.
    symbol, species, stem = _CLONE_STYLE_SYMBOL
    links = motif_links(symbol, species)

    assert (links.gene_id_stem, links.symbol) == (stem, symbol)


# ---------------------------------------------------------------------------------------
# The whole table, and the frame a collaborator reads it as
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_a_table_says_which_species_release_and_tax_group_it_is(slug: str, release: str) -> None:
    table = _shipped(slug, release)

    assert (table.release, table.tax_group) == (release, LINK_TAX_GROUP)
    assert table.species
    assert table.source.endswith(f"{slug}.jaspar{release}.motif_link_table.tsv.gz")


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_a_tables_genes_are_the_censuss_and_every_one_was_judged_a_tf(
    slug: str, release: str
) -> None:
    # Only assessed-positive genes receive links, which is the census speaking and never
    # this package.
    table = _shipped(slug, release)
    census = tf_gene_table(slug)
    assert census is not None
    positive = set(census.assessed_positive)

    assert set(table.gene_id_stems) <= positive


@pytest.mark.parametrize(("slug", "release"), _TABLES)
def test_the_frame_carries_the_files_own_columns_and_one_row_per_link(
    slug: str, release: str
) -> None:
    table = _shipped(slug, release)
    frame = table.frame()

    assert list(frame.columns) == list(LINK_COLUMNS)
    assert len(frame) == len(table)


def test_the_frame_is_built_fresh_so_mutating_it_cannot_reach_the_table() -> None:
    table = _shipped("homo_sapiens", "2026")
    frame = table.frame()
    frame.loc[0, "motif_id"] = "MA9999.9"

    assert table.frame().loc[0, "motif_id"] != "MA9999.9"


def test_the_same_table_is_answered_however_the_species_is_spelled() -> None:
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
        "release\tspecies\n",
        "\t".join(LINK_COLUMNS) + "\n",
        "\t".join((*LINK_COLUMNS, "quality_score")) + "\n",
    ],
    ids=["empty", "two-columns", "header-only", "an-extra-column"],
)
def test_a_file_that_is_not_a_link_table_raises_and_names_the_file(text: str) -> None:
    with pytest.raises(MotifLinkTableError, match=r"broken\.tsv"):
        parse_motif_link_table(text, source="broken.tsv")


@pytest.mark.parametrize(
    "changes",
    [
        {"role": "dimer"},
        {"role": "complex"},
        {"partners": "FOS"},
        {"is_cross_species": "true"},
        {"rank": "first"},
        {"total_information_content": "high"},
        {"gene_id_stem": ""},
        {"motif_id": ""},
    ],
    ids=[
        "a-third-role",
        "a-complex-naming-no-partner",
        "a-monomer-naming-one",
        "a-flag-no-census-spells",
        "a-rank-that-is-not-a-number",
        "an-information-content-that-is-not-a-number",
        "no-gene-id-stem",
        "no-motif-id",
    ],
)
def test_a_row_a_generator_would_never_write_raises_and_names_its_line(
    changes: dict[str, str],
) -> None:
    with pytest.raises(MotifLinkTableError, match="line 2"):
        parse_motif_link_table(_table_text(_row(**changes)), source="broken.tsv")


def test_a_row_with_the_wrong_number_of_cells_raises() -> None:
    with pytest.raises(MotifLinkTableError, match="cells"):
        parse_motif_link_table(_table_text(_ROW[:-1]), source="broken.tsv")


@pytest.mark.parametrize("changes", [{"release": "2024"}, {"species": "Mus musculus"}])
def test_one_file_naming_two_releases_or_two_species_raises(changes: dict[str, str]) -> None:
    # The release and the species are on every row so that two tables concatenate into one
    # frame that still says where each row came from — which promises nothing unless they
    # are uniform within a file.
    text = _table_text(_ROW, _row(motif_id="MA0106.4", **changes))

    with pytest.raises(MotifLinkTableError, match=r"broken\.tsv"):
        parse_motif_link_table(text, source="broken.tsv")


def test_a_broken_table_is_a_bad_value_and_not_a_lookup() -> None:
    # A packaging defect, not an absence: a caller catching LookupError for a species with
    # no table must not swallow a file that ships broken.
    with pytest.raises(ValueError) as raised:  # noqa: PT011 - the class is the assertion
        parse_motif_link_table("nonsense\n", source="broken.tsv")

    assert isinstance(raised.value, MotifLinkTableError)
    assert not isinstance(raised.value, LookupError)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "release\tspecies\n",
        "\t".join(LINK_COLUMNS) + "\n",
        _table_text(_ROW[:-1]),
        _table_text(_row(role="dimer")),
        _table_text(_row(rank="first")),
        _table_text(_ROW, _row(release="2024", motif_id="MA0106.4")),
    ],
    ids=[
        "empty",
        "two-columns",
        "header-only",
        "a-short-row",
        "a-third-role",
        "a-rank-that-is-not-a-number",
        "two-releases",
    ],
)
def test_every_message_a_broken_file_raises_names_the_command_that_rewrites_it(text: str) -> None:
    # A shipped file that cannot be trusted is a packaging defect, and regenerating it is
    # the only repair — so every message names the file and the generator that writes it.
    with pytest.raises(MotifLinkTableError) as raised:
        parse_motif_link_table(text, source="broken.tsv")

    assert "broken.tsv" in str(raised.value)
    assert "scripts/build_tf_links.py" in str(raised.value)
