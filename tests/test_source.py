"""Tests for genome.io.source — what one assembly name resolves to, and in what order.

The resolution used to live under the downloader, where the only way to reach it was to
register something. It is a value now, so these are the four ordered checks asserted
directly: a record here, then a source the caller named, then the name, then the fetch
that every other name falls through to (ADR-0003, ADR-0008).

Nothing here needs a native tool or a network: a record is written by hand, and the whole
question is answered from that record, the shipped table and the name.
"""

from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import pytest

from genome.io import components as components_module
from genome.io import download as download_module
from genome.io import source as source_module
from genome.io.completion import build_record, write_record
from genome.io.components import ChimeraDetails, ComponentDetails
from genome.io.registration import AssemblyDir
from genome.io.source import (
    ComponentSource,
    FetchedSource,
    SeededSource,
    fetched_source,
    is_prepared,
    resolve_source,
)
from genome.metadata import AssemblyMetadata

_GOLDEN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

_PINNED_URL = "https://mirror.example.org/references/tiny.fa.gz"


def _row(source_url: str | None) -> AssemblyMetadata:
    """Return a complete curated row for ``tiny``, pinning ``source_url`` or nothing."""
    return AssemblyMetadata(
        assembly_name="tiny",
        species="Testus minimus",
        ucsc_name="tiny",
        ncbi_name="TINY.1",
        ncbi_assembly_id="GCF_000000000.0",
        ncbi_taxid=1,
        source_url=source_url,
    )


def _record_a_genome(directory: Path, assembly: str, details: dict[str, object] | None) -> None:
    """Write a completion record in ``directory`` claiming ``assembly`` finished there."""
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / f"{assembly}.fa"
    fasta.write_text(">chrI\nACGT\n")
    write_record(
        directory,
        build_record(directory, kind="genome", name=assembly, files=[fasta], details=details),
    )


def _chimera_details(*components: str) -> dict[str, object]:
    """Return the ``details`` a chimera build of ``components`` writes."""
    return ChimeraDetails(
        separator="__",
        component_details=tuple(ComponentDetails(name, None, None, None) for name in components),
    ).as_details(merged=False)


# ---------------------------------------------------------------------------
# Check four: the fetch every ordinary name falls through to
# ---------------------------------------------------------------------------


def test_the_fetch_a_name_falls_through_to_is_pinned_derived_or_blank_falls_back(
    tmp_path: Path,
) -> None:
    # A name the shipped table does not list is fetched from the derived golden-path URL...
    source = resolve_source(
        AssemblyDir.locate("hg38", tmp_path), metadata=None, golden_path_url=_GOLDEN
    )
    assert source == FetchedSource(url=_GOLDEN, derived=True)

    # ...a row that pins a source answers outright and is not derived, `derived` being
    # the whole difference that decides whether the name is checked at UCSC first...
    pinned = resolve_source(
        AssemblyDir.locate("tiny", tmp_path), metadata=_row(_PINNED_URL), golden_path_url=_GOLDEN
    )
    assert pinned == FetchedSource(url=_PINNED_URL, derived=False)

    # ...and a blank pin is the same as no pin at all.
    assert fetched_source(_row(""), _GOLDEN).url == _GOLDEN


# ---------------------------------------------------------------------------
# Check three: the name
# ---------------------------------------------------------------------------


def test_a_name_is_read_as_a_component_set_only_when_both_halves_are_assemblies(
    tmp_path: Path,
) -> None:
    # Both parts are in the shipped table, which is what tells a chimera's name from a
    # free-form local key on a machine holding neither.
    chimera = resolve_source(
        AssemblyDir.locate("ce11_ecHT115", tmp_path), metadata=None, golden_path_url=_GOLDEN
    )
    assert chimera == ComponentSource(("ce11", "ecHT115"))

    # Neither half being nobody's assembly reads as one plain name instead...
    plain = resolve_source(
        AssemblyDir.locate("my_ref", tmp_path), metadata=None, golden_path_url=_GOLDEN
    )
    assert plain == FetchedSource(url=_GOLDEN, derived=True)

    # ...and a name with no underscore never reaches the split in the first place.
    assert isinstance(
        resolve_source(
            AssemblyDir.locate("hg38", tmp_path), metadata=None, golden_path_url=_GOLDEN
        ),
        FetchedSource,
    )


def test_a_mis_ordered_chimera_name_is_refused_and_told_its_spelling(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="genome register ce11_ecHT115"):
        resolve_source(
            AssemblyDir.locate("ecHT115_ce11", tmp_path), metadata=None, golden_path_url=_GOLDEN
        )


def test_only_a_prepared_component_counts_when_the_table_lists_neither(
    tmp_path: Path, liulab_data: Path
) -> None:
    # Neither half is in the shipped table, so nothing but a record of its own can make
    # `tinyCe_tinySc` read as two assemblies rather than as one name somebody chose.
    here = AssemblyDir.locate("tinyCe_tinySc", tmp_path / "elsewhere")

    assert isinstance(resolve_source(here, metadata=None, golden_path_url=_GOLDEN), FetchedSource)

    for name in ("tinyCe", "tinySc"):
        _record_a_genome(liulab_data / "genome" / name, name, None)
        assert is_prepared(name)

    assert resolve_source(here, metadata=None, golden_path_url=_GOLDEN) == ComponentSource(
        ("tinyCe", "tinySc")
    )


# ---------------------------------------------------------------------------
# Check one: a record here, believed outright
# ---------------------------------------------------------------------------


def test_an_own_record_is_believed_outright_over_the_name_and_over_no_record_at_all(
    tmp_path: Path,
) -> None:
    # A chimera's own record says what it is made of...
    chimera_dir = tmp_path / "chimera"
    _record_a_genome(chimera_dir, "ce11_ecHT115", _chimera_details("ce11", "ecHT115"))
    assert resolve_source(
        AssemblyDir.locate("ce11_ecHT115", chimera_dir), metadata=None, golden_path_url=_GOLDEN
    ) == ComponentSource(("ce11", "ecHT115"))

    # ...and the reason the record comes first: `ce11_ecHT115` seeded years ago from
    # somebody's own FASTA is not a chimera, and no amount of the name looking like one
    # may change what a finished registration already is.
    plain_dir = tmp_path / "plain"
    _record_a_genome(plain_dir, "ce11_ecHT115", None)
    assert resolve_source(
        AssemblyDir.locate("ce11_ecHT115", plain_dir), metadata=None, golden_path_url=_GOLDEN
    ) == FetchedSource(url=_GOLDEN, derived=True)

    # A record that is lost, though, falls all the way through to the name.
    (plain_dir / ".completion.json").unlink()
    assert resolve_source(
        AssemblyDir.locate("ce11_ecHT115", plain_dir), metadata=None, golden_path_url=_GOLDEN
    ) == ComponentSource(("ce11", "ecHT115"))


# ---------------------------------------------------------------------------
# The kinds themselves
# ---------------------------------------------------------------------------


def test_a_seeded_source_is_declared_and_never_resolved() -> None:
    # The one kind nothing works out: the caller said it, so it is carried as given —
    # a Path stays a Path, since it is quoted back into the repair command verbatim.
    assert SeededSource(Path("/data/my ref.fa")).location == Path("/data/my ref.fa")


# ---------------------------------------------------------------------------
# The layer this module exists to be
# ---------------------------------------------------------------------------


def _module_level_imports(module: ModuleType) -> set[str]:
    """Return the modules ``module`` imports when it is loaded, by reading its source.

    Statements at the top level of the file only — an import inside a function is what
    this is measuring the absence of, so counting one would defeat the purpose. Read
    rather than observed, because importing a submodule loads its package first and
    ``sys.modules`` would then hold half the tree whatever this module itself asked for.

    A ``from x import y`` contributes both ``x`` and ``x.y``, which is what lets a caller
    tell *holds the module* from *holds the name it exports* — the distinction the fetch
    step's one patch point rests on (see test_fetch).

    Shared with the guards in tests/annotation and test_fetch: the assertions differ per
    module, the reading of an import does not.
    """
    assert module.__file__ is not None
    imported: set[str] = set()
    for node in ast.parse(Path(module.__file__).read_text()).body:
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported.add(node.module)
            imported.update(f"{node.module}.{alias.name}" for alias in node.names)
    return imported


def test_the_seam_closes_the_download_annotation_chimera_cycle_in_one_direction_only() -> None:
    # Answering *what is this name* must cost none of what acting on the answer costs,
    # which is what lets the downloader import this at the top of the file. Were `source`
    # to reach the annotation package (which imports the downloader) or `genome` (which
    # imports everything), the resolution would go back behind a deferred import and take
    # the cycle it was split out to close with it. A prefix for the annotation half, since
    # every module of that package is out of bounds and a new one must not arrive
    # unnoticed.
    forbidden = {"genome.io.chimera", "genome.io.download", "genome.genome"}
    reached = _module_level_imports(source_module)
    assert reached & forbidden == set()
    assert {name for name in reached if name.startswith("genome.io.annotation")} == set()

    # The seam the split left, in the one direction it runs. Resolving a name believes an
    # existing record before it consults the name, so it has to tell a chimera's record
    # from any other — one call into `components`. Nothing in `components` needs a name
    # resolved, and an import back would make the two one module again by another route.
    assert "genome.io.components" in _module_level_imports(source_module)
    assert "genome.io.source" not in _module_level_imports(components_module)

    # The other half: four deferred imports paid for the old arrangement, and the only
    # one left is the build, which opens a whole Genome per component.
    assert "genome.io.source" in _module_level_imports(download_module)
    assert "genome.io.chimera" not in _module_level_imports(download_module)
