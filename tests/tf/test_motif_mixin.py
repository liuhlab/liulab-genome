"""Tests for genome.tf.motif.mixin — **Region**s scanned against a prepared assembly.

The highest-value tests of the whole motif feature are here, because this is the one
place a region-local position becomes a **Chromosome** coordinate and the one place a
`-` strand region's hits are flipped back into the forward frame. So nothing here settles
for *a hit was found*: the planted sites are asserted at their exact intervals and on
their exact strands, in the assembly's frame, from a plus-strand region and a
minus-strand one alike.

The assembly is `tests/data/planted_motifs.fa` prepared as one — two 600-base
chromosomes, `plantedI` and `plantedII`, with three consensus words written into them at
positions `tests/data/README.md` documents and `tests/tf/test_scan.py` asserts against the
committed bytes. That is what makes an exact expected coordinate available: the site at
`plantedI[100, 115)` is there because it was put there.

The other half of what is asserted here is that a region scan is not a *smaller* scan:
every argument a **Motif set** scan takes is forwarded — the **Background** and the worker
count included, which is what the HPC case is made of — and the one that cannot be,
an output path, is refused by name rather than by an unexpected keyword.

Preparing it needs the native tools, so the fixture skips when they are absent; nothing
here needs STAR or chromap, so these are unit-lane tests and unmarked. Nothing reaches
the network: the assembly is seeded from a committed file behind the same guard the
suite's autouse fixture installs.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st

import genome.assembly.genome as genome_module
from genome import Genome, Region
from genome.assembly.fasta import PREPARATION_TOOLS
from genome.tf.motif import (
    HIT_COLUMNS,
    HIT_DTYPES,
    REGION_HIT_PROVENANCE,
    UNIFORM_BACKGROUND,
    Motif,
    MotifScanMixin,
    MotifSet,
    parse_transfac,
)
from genome.tf.motif import mixin as mixin_module

from ..assembly.test_source import _module_level_imports
from ..conftest import DATA_DIR, install_network_guard

# Every test here spawns its own worker processes. Under `--dist=loadgroup` that pins
# them to ONE xdist worker, so they run one at a time rather than eight of them forking
# pools at once — which oversubscribed the box and made the lane's wall bimodal, 13.5 s
# or 16.1 s depending purely on how they happened to be scheduled.
pytestmark = pytest.mark.xdist_group("spawns_mixin")

#: The committed FASTA prepared as an assembly, and the name it is prepared under. See
#: ``tests/data/README.md`` for what was planted in it.
FIXTURE = "planted_motifs.fa"
ASSEMBLY = "planted"

#: The committed motif records, the same ten every other motif test uses.
MOTIF_FIXTURE = "tiny_jaspar_transfac.txt"

#: Both chromosomes of the prepared assembly are this long.
LENGTH = 600

#: The two CTCF sites planted in ``plantedI``, as ``(start, end, strand)`` in the
#: assembly's own frame — the answer this whole module exists to check a scan against.
FORWARD_SITE = (100, 115, "+")
REVERSE_SITE = (300, 315, "-")

#: The **Motif id** both of them are the consensus of.
CTCF = "MA0139.2"


@pytest.fixture(scope="module")
def assembly(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Genome]:
    """The planted FASTA prepared as an assembly, opened once for the module.

    Module-scoped for two reasons. Hypothesis refuses to drive a test from a
    function-scoped fixture, and the property-based test below needs this one; and
    preparing it is three native-tool runs over 1.2 kb of committed bytes that cannot
    differ between tests. Nothing here writes to it.

    Skips when the preparation tools are not on ``PATH``, so a test using it needs no
    skip marker of its own — the same arrangement ``chimera_component`` uses.
    """
    missing = [tool for tool in PREPARATION_TOOLS if shutil.which(tool) is None]
    if missing:
        pytest.skip(f"not on PATH: {', '.join(missing)}")
    # Built before any test's setup, so the autouse guard is not up yet: stand one up
    # here rather than leave the one registration nobody watches.
    with pytest.MonkeyPatch.context() as guard:
        install_network_guard(guard)
        opened = Genome(
            ASSEMBLY,
            path_or_url=DATA_DIR / FIXTURE,
            cache_dir=tmp_path_factory.mktemp(ASSEMBLY),
            progressbar=False,
        )
    yield opened
    opened.close()


@pytest.fixture(scope="module")
def ctcf() -> MotifSet:
    """The one committed CTCF matrix, which both planted sites are the consensus of."""
    records = parse_transfac((DATA_DIR / MOTIF_FIXTURE).read_text(encoding="utf-8"))
    return MotifSet(records).filter(lambda motif: motif.motif_id == CTCF)


@pytest.fixture(scope="module")
def dense() -> MotifSet:
    """A 7-mer common enough in this yeast that a random window usually holds one.

    The committed matrices are the right thing to assert exact coordinates with and the
    wrong thing to drive a property test with — nine hits in 1 200 bases means most
    generated regions would hold none, and a property that mostly compares two empty
    tables proves nothing.
    """
    return MotifSet([word_motif("AAAAAAA")])


def word_motif(bases: str, motif_id: str = "MA9999.1", name: str = "Testin") -> Motif:
    """A motif fixed on every base of ``bases`` — 2 bits a position."""
    counts = np.zeros((4, len(bases)))
    for column, base in enumerate(bases.upper()):
        counts["ACGT".index(base), column] = 100.0
    return Motif(motif_id, name, counts)


def rows(frame: pd.DataFrame) -> list[tuple[Any, ...]]:
    """Every row as a plain tuple, so two tables compare without their category dtypes."""
    return [tuple(row) for row in frame.itertuples(index=False, name=None)]


def sites(frame: pd.DataFrame) -> set[tuple[Any, ...]]:
    """Every hit as ``(sequence_name, start, end, strand)`` — the identity and score dropped."""
    return {
        (sequence_name, int(start), int(end), strand)
        for _id, _name, sequence_name, start, end, strand, _score in rows(frame)
    }


def into_local_frame(frame: pd.DataFrame, region: Region) -> set[tuple[Any, ...]]:
    """Map every hit back into ``region``'s own frame — the inverse of what the mixin did.

    Written out here rather than imported so the test states the arithmetic independently
    of the code under test: adding the region start is undone by subtracting it, and a
    reverse region's swap of the two ends is undone by swapping them back.
    """
    back: set[tuple[Any, ...]] = set()
    for _id, _name, chrom, start, end, strand, _score in rows(frame):
        assert chrom == region.chrom
        if region.strand == "-":
            local = (region.end - int(end), region.end - int(start))
            local_strand = "-" if strand == "+" else "+"
        else:
            local = (int(start) - region.start, int(end) - region.start)
            local_strand = strand
        back.add((region.chrom, *local, local_strand))
    return back


# ---------------------------------------------------------------------------
# The planted sites, at the exact coordinates they were planted at
# ---------------------------------------------------------------------------


class TestThePlantedSites:
    def test_a_plus_strand_region_lifts_local_coordinates_into_chromosome_ones(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # The everyday case, and the whole point: the scan saw bases starting at zero and
        # the answer is where those bases are — not the region-local number (50, not 100)
        # a plain scan of the fetched bases would have answered with.
        region = Region("plantedI", 50, 200, "+")
        local = ctcf.scan(assembly.fetch_sequence(region), region.chrom)

        assert sites(local) == {("plantedI", 50, 65, "+")}
        assert sites(assembly.scan_regions(ctcf, region)) == {("plantedI", *FORWARD_SITE)}

    def test_a_minus_strand_region_flips_the_interval_and_the_strand(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # The off-by-one this feature exists to centralise. Fetched, the region is the
        # reverse complement of [50, 200), and the scan finds the site at [85, 100) on '-'
        # in that frame; lifted, it is the same [100, 115) on '+' the plus-strand region
        # answered with. 200 - 100 = 100 and 200 - 85 = 115: the ends swap.
        region = Region("plantedI", 50, 200, "-")
        local = ctcf.scan(assembly.fetch_sequence(region), region.chrom)

        assert sites(local) == {("plantedI", 85, 100, "-")}
        assert sites(assembly.scan_regions(ctcf, region)) == {("plantedI", *FORWARD_SITE)}

    def test_the_planted_reverse_site_is_reported_correctly_from_either_region_strand(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # The other planted word is the same consensus reverse-complemented, so from a
        # plus-strand region — where nothing flips — it is a '-' hit over those bases.
        plus = Region("plantedI", 250, 400, "+")
        assert sites(assembly.scan_regions(ctcf, plus)) == {("plantedI", *REVERSE_SITE)}

        # And flipping a region that holds it flips it twice, which is not a flip: a '+'
        # hit in the fetched frame is a '-' hit on the chromosome.
        minus = Region("plantedI", 250, 400, "-")
        local = ctcf.scan(assembly.fetch_sequence(minus), minus.chrom)
        assert sites(local) == {("plantedI", 85, 100, "+")}
        assert sites(assembly.scan_regions(ctcf, minus)) == {("plantedI", *REVERSE_SITE)}

    def test_a_region_holding_both_sites_answers_with_both(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        whole = Region("plantedI", 0, LENGTH, "+")

        assert sites(assembly.scan_regions(ctcf, whole)) == {
            ("plantedI", *FORWARD_SITE),
            ("plantedI", *REVERSE_SITE),
        }


# ---------------------------------------------------------------------------
# The arithmetic itself
# ---------------------------------------------------------------------------


class TestTheCoordinateArithmetic:
    def test_the_answer_does_not_depend_on_which_strand_the_region_was_fetched_as(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # The invariant behind the two assertions above, stated over every hit: which way
        # round a region was fetched is not a property of the sites in it.
        forward = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "+"))
        reverse = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "-"))
        assert sites(forward) == sites(reverse)

        # '.' means nobody knows, and the fetch hands back forward bases for it, so there
        # is nothing to flip — the coordinates match the '+' region's rather than being
        # decided by calling the strand '+'.
        unknown = assembly.scan_regions(ctcf, Region("plantedI", 50, 200, "."))
        plus = assembly.scan_regions(ctcf, Region("plantedI", 50, 200, "+"))
        assert sites(unknown) == sites(plus) == {("plantedI", *FORWARD_SITE)}

    def test_the_lifted_length_is_the_length_the_scan_found(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # Half-open is what makes the reverse case symmetric; a lifted hit that changed
        # length would be the classic 1-based-inclusive mistake.
        hits = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "-"))

        assert set(hits["end"] - hits["start"]) == {15}

    @given(
        chrom=st.sampled_from(["plantedI", "plantedII"]),
        bounds=st.tuples(st.integers(0, LENGTH), st.integers(0, LENGTH)),
        strand=st.sampled_from(["+", "-", "."]),
    )
    def test_a_lifted_hit_maps_back_to_where_the_scan_found_it(
        self, assembly: Genome, dense: MotifSet, chrom: str, bounds: tuple[int, int], strand: str
    ) -> None:
        # The round trip, over any region of either chromosome: lift a hit into the
        # assembly's frame, map it back into the region's, and it is where the scan of
        # that region's own bases put it. Degenerate regions are generated too — an empty
        # one, and one shorter than the motif — and both answer with no hits rather than
        # raising.
        low, high = sorted(bounds)
        region = Region(chrom, low, high, strand)
        local = dense.scan(assembly.fetch_sequence(region), chrom)

        lifted = assembly.scan_regions(dense, region)

        assert into_local_frame(lifted, region) == sites(local)

    @given(
        bounds=st.tuples(st.integers(0, LENGTH), st.integers(0, LENGTH)),
        strand=st.sampled_from(["+", "-", "."]),
    )
    def test_every_lifted_hit_lands_inside_the_region_it_came_from(
        self, assembly: Genome, dense: MotifSet, bounds: tuple[int, int], strand: str
    ) -> None:
        low, high = sorted(bounds)
        region = Region("plantedII", low, high, strand)

        for _chrom, start, end, _strand in sites(assembly.scan_regions(dense, region)):
            assert region.start <= start < end <= region.end


# ---------------------------------------------------------------------------
# The table: one schema, whatever the scan was handed, plus the assembly
# ---------------------------------------------------------------------------


class TestTheTable:
    def test_the_schema_and_chromosome_names_are_correct(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        hits = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "-"))
        assert list(hits.columns) == list(HIT_COLUMNS)
        assert [str(dtype) for dtype in hits.dtypes] == list(HIT_DTYPES.values())

        # Not the region key the scan was driven with, and not a locus string: the name
        # this assembly spells the chromosome with, which is why an unknown one raises
        # rather than being carried through.
        regions = [Region("plantedI", 0, LENGTH, "-"), Region("plantedII", 0, LENGTH, "+")]
        by_two = assembly.scan_regions(ctcf, regions)
        assert set(by_two["sequence_name"]) <= set(assembly.chromosomes)
        assert "plantedI" in set(by_two["sequence_name"])

    def test_the_assembly_and_the_scan_s_own_provenance_both_travel(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        hits = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "+"), threshold=1e-3)

        # Extended, not replaced: the background, threshold, release, tax group and the
        # two motif lists are all still there beside the assembly.
        assert set(hits.attrs) == set(REGION_HIT_PROVENANCE)
        assert hits.attrs["assembly"] == ASSEMBLY
        assert hits.attrs["threshold"] == 1e-3
        assert hits.attrs["motifs_scanned"] == (CTCF,)
        assert hits.attrs["release"] is None

    def test_no_regions_or_no_hits_is_an_empty_table_that_still_says_what_it_was(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        no_regions = assembly.scan_regions(ctcf, [])
        assert no_regions.empty
        assert list(no_regions.columns) == list(HIT_COLUMNS)
        assert no_regions.attrs["assembly"] == ASSEMBLY

        no_hits = assembly.scan_regions(ctcf, Region("plantedI", 400, 500, "-"))
        assert no_hits.empty
        assert no_hits.attrs["assembly"] == ASSEMBLY


# ---------------------------------------------------------------------------
# What a region scan may be asked for, and the one thing it may not
# ---------------------------------------------------------------------------


class TestTheScanArgumentsReachTheScan:
    """Everything a **Motif set** scan takes is forwarded, so a region scan is not a
    smaller scan: the HPC case — a derived **Background** and the whole allocation — is
    reachable from here or it is reachable nowhere.
    """

    def test_background_modes_are_forwarded_to_the_scan(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # 600 bases is far under the derivation floor, so 'auto' stays uniform — which is
        # what makes the derived answer below a different number and not a coincidence.
        default = assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "+"))
        assert default.attrs["background"] == UNIFORM_BACKGROUND

        # The mode reaches the scan, and what it derives from is the very bases the scan
        # saw: the same region scanned region-locally records the same background.
        region = Region("plantedI", 0, LENGTH, "+")
        local = ctcf.scan(assembly.fetch_sequence(region), region.chrom, background="derive")
        derived = assembly.scan_regions(ctcf, region, background="derive")
        assert derived.attrs["background"] == local.attrs["background"] != UNIFORM_BACKGROUND

        given_background = assembly.scan_regions(
            ctcf, Region("plantedI", 0, LENGTH, "+"), background=[0.4, 0.1, 0.1, 0.4]
        )
        assert given_background.attrs["background"] == (0.4, 0.1, 0.1, 0.4)

    def test_a_bad_background_or_worker_count_is_refused(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        with pytest.raises(ValueError, match="auto, uniform, derive"):
            assembly.scan_regions(
                ctcf,
                Region("plantedI", 0, LENGTH, "+"),
                background="gc",  # type: ignore[arg-type]
            )
        with pytest.raises(ValueError, match="at least 1"):
            assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "+"), workers=0)

    def test_two_workers_answer_with_the_identical_table(
        self, assembly: Genome, dense: MotifSet
    ) -> None:
        # The property the whole parallel design is verified by, held after the lift as
        # well as before it: identical means row for row, dtypes and provenance included.
        regions = [Region("plantedI", 0, LENGTH, "-"), Region("plantedII", 0, LENGTH, "+")]

        shared = assembly.scan_regions(dense, regions, workers=2)
        serial = assembly.scan_regions(dense, regions, workers=1)

        assert rows(shared) == rows(serial) != []
        assert [str(dtype) for dtype in shared.dtypes] == list(HIT_DTYPES.values())
        assert shared.attrs == serial.attrs


class TestSeveralRegions:
    def test_two_regions_are_both_scanned_overlapping_ones_answer_separately_and_one_needs_no_list(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        # A mapping keyed by chromosome could hold only one of these, which is why the
        # scan is keyed by position and the name is put back afterwards.
        regions = [Region("plantedI", 0, 200, "+"), Region("plantedI", 250, 400, "+")]
        assert sites(assembly.scan_regions(ctcf, regions)) == {
            ("plantedI", *FORWARD_SITE),
            ("plantedI", *REVERSE_SITE),
        }

        # The same site seen twice is two rows, one per region — deduplicating would be
        # deciding for the caller which of two peaks a hit belongs to.
        overlapping = [Region("plantedI", 0, 200, "+"), Region("plantedI", 50, 250, "-")]
        hits = assembly.scan_regions(ctcf, overlapping)
        assert len(hits) == 2
        assert sites(hits) == {("plantedI", *FORWARD_SITE)}

        # And a single region needs no collection built around it.
        one = Region("plantedI", 0, LENGTH, "+")
        assert rows(assembly.scan_regions(ctcf, one)) == rows(assembly.scan_regions(ctcf, [one]))


# ---------------------------------------------------------------------------
# The raw form, which this method adds to rather than replaces
# ---------------------------------------------------------------------------


class TestTheRawFormIsUntouched:
    def test_scanning_a_mapping_of_sequences_is_still_region_local_and_names_no_assembly(
        self, assembly: Genome, ctcf: MotifSet
    ) -> None:
        region = Region("plantedI", 50, 200, "+")
        local = ctcf.scan_sequences({"peak1": assembly.fetch_sequence(region)})
        assert sites(local) == {("peak1", 50, 65, "+")}

        # A motif belongs to no assembly, and neither does a table scanned from bases
        # nobody said where they came from.
        assert "assembly" not in local.attrs


# ---------------------------------------------------------------------------
# Refusals, each naming what to do instead
# ---------------------------------------------------------------------------


class TestRefusals:
    def test_bad_arguments_are_refused_and_name_the_next_action(
        self, assembly: Genome, ctcf: MotifSet, tmp_path: Path
    ) -> None:
        # A string is iterable, so this would otherwise be read letter by letter; and it
        # carries no strand, which is the one thing this method is about.
        with pytest.raises(TypeError, match=r"Region\.from_string"):
            assembly.scan_regions(ctcf, "plantedI:0-600")  # type: ignore[arg-type]

        with pytest.raises(ValueError, match="unknown chromosome"):
            assembly.scan_regions(ctcf, Region("chrNope", 0, 100, "+"))

        with pytest.raises(ValueError, match="exceeds"):
            assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH + 1, "+"))

        with pytest.raises(ValueError, match="p-value"):
            assembly.scan_regions(ctcf, Region("plantedI", 0, LENGTH, "+"), threshold=5.0)

        # The one scan argument that cannot be forwarded: Parquet hands back a path, and a
        # path holds nothing to lift. Refused by name rather than by an unexpected keyword,
        # and refused before anything is fetched or scanned.
        with pytest.raises(TypeError, match="write_hits"):
            assembly.scan_regions(
                ctcf, Region("plantedI", 0, LENGTH, "+"), output=tmp_path / "hits.parquet"
            )

        # A region this assembly could not fetch would raise its own error first if the
        # refusal came later, which would tell the caller about the wrong problem.
        with pytest.raises(TypeError, match="scan_regions cannot stream"):
            assembly.scan_regions(
                ctcf, Region("chrNope", 0, 100, "+"), output=tmp_path / "hits.parquet"
            )


# ---------------------------------------------------------------------------
# The edge, and which way it runs
# ---------------------------------------------------------------------------


class TestTheDependencyDirection:
    def test_the_dependency_edge_runs_from_genome_to_motif_and_never_back(self) -> None:
        # A motif belongs to no assembly and a motif set is usable with no genome open,
        # so the edge runs Genome to motif and never back. An import at module level here
        # would be the back edge, whatever the type annotations said; the other half is
        # that a base class cannot itself be imported lazily, so that edge is a plain
        # module-level import and is meant to be.
        assert "genome.assembly.genome" not in _module_level_imports(mixin_module)
        assert "genome.tf.motif.mixin" in _module_level_imports(genome_module)

    def test_a_genome_is_one(self, assembly: Genome) -> None:
        assert isinstance(assembly, MotifScanMixin)
