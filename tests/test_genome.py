"""Tests for genome.Genome — sequence retrieval over a prepared assembly.

A real tiny assembly is built with the native tools, so the whole module skips when
those are absent. Nothing reaches the network by either of two routes: most tests
monkeypatch ``UCSCGenomeDownloader.fetch_genome`` to hand back the prebuilt files,
and the registration tests take the whole path for real behind the shared
``fake_fetch`` fixture, which serves ``tests/data`` in place of a download.
"""

from __future__ import annotations

import shutil
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pytest

import genome.genome as genome_mod
from genome import (
    DNA,
    AnnotationRegistry,
    Genome,
    NoCofactorTableError,
    NoGeneCategoriesError,
    NoTFCensusError,
    Region,
)
from genome.gene_list import curated_gene_list
from genome.io.annotation import AnnotationNotRegisteredError, annotation_dir, register_gtf
from genome.io.completion import (
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
)
from genome.io.download import register_assembly
from genome.io.fasta import prepare_fasta
from genome.metadata import (
    METADATA_FIELDS,
    AnnotationMetadata,
    AssemblyMetadata,
    lookup_assembly,
)

from .conftest import FakeFetch

_REQUIRED = ("samtools", "faToTwoBit", "twoBitInfo")
_TOOLS_PRESENT = all(shutil.which(t) is not None for t in _REQUIRED)
pytestmark = pytest.mark.skipif(
    not _TOOLS_PRESENT, reason="samtools/faToTwoBit/twoBitInfo not on PATH"
)

# chrA: 20 bp with a soft-masked stretch at [8, 12); chrB: 8 bp.
_CHR_A = "ACGTACGTacgtACGTACGT"
_CHR_B = "TTTTGGGG"


def _prepare_cache_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Prepare ``tmp_path/tiny.fa`` and hand its files to every ``Genome`` built here.

    Skips the network/UCSC validation entirely: every ``Genome(...)`` built with this
    cache dir is handed back the prebuilt files, whatever it is named.
    """
    files = prepare_fasta(tmp_path / "tiny.fa")
    monkeypatch.setattr(
        genome_mod.UCSCGenomeDownloader,
        "fetch_genome",
        lambda self, **kwargs: files,
    )
    return tmp_path


@pytest.fixture
def prepared_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Prepare a two-chromosome assembly of known bases and return its cache dir.

    Synthetic on purpose: chrA carries a soft-masked stretch at ``[8, 12)`` and a
    non-palindromic prefix, so the sequence tests can assert exact bases. It is not an
    organism, so nothing is registered against it — the annotation tests use
    :func:`yeast_dir`.
    """
    (tmp_path / "tiny.fa").write_text(f">chrA\n{_CHR_A}\n>chrB\n{_CHR_B}\n")
    return _prepare_cache_dir(tmp_path, monkeypatch)


@pytest.fixture
def yeast_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, data_dir: Path) -> Path:
    """Prepare an assembly and return the cache dir an annotation actually fits.

    The committed ``tiny.fa`` is real subsampled sacCer3 — ``chrI``, ``chrII``,
    ``chrIII`` — which are the very chromosomes the committed ``tiny.gtf`` names, so a
    GTF registered against this assembly passes the chromosome check on its merits
    rather than needing it stood down.
    """
    shutil.copy2(data_dir / "tiny.fa", tmp_path / "tiny.fa")
    return _prepare_cache_dir(tmp_path, monkeypatch)


@pytest.fixture
def genome(prepared_dir: Path) -> Iterator[Genome]:
    g = Genome("tiny", cache_dir=prepared_dir)
    yield g
    g.close()


def test_fetch_sequence_handles_ranges_strings_regions_and_strand(genome: Genome) -> None:
    result = genome.fetch_sequence("chrA:0-8")
    assert result == DNA("ACGTACGT")
    assert isinstance(result, DNA)
    assert genome.fetch_sequence("chrA:8-12") == DNA("acgt")  # soft masking preserved
    assert genome.fetch_sequence("chrB") == DNA(_CHR_B)  # bare chromosome: whole sequence
    assert genome["chrA:0-4"] == DNA("ACGT")  # __getitem__ is sugar for fetch_sequence
    assert genome.fetch_sequence(Region("chrA", 0, 8)) == DNA("ACGTACGT")  # accepts a Region
    assert genome.fetch_sequence("chrA:0-20") == DNA(_CHR_A)  # end == size is allowed

    # chrA[0:6] == "ACGTAC" — not a palindrome, so the reverse complement differs.
    plus = genome.fetch_sequence(Region("chrA", 0, 6, "+"))
    minus = genome.fetch_sequence(Region("chrA", 0, 6, "-"))
    assert plus == DNA("ACGTAC")
    assert minus == plus.reverse_complement() == DNA("GTACGT")


def test_genome_exposes_chrom_sizes_chromosomes_paths_and_repr(genome: Genome) -> None:
    sizes = genome.chrom_sizes
    assert isinstance(sizes, pd.Series)
    assert sizes.to_dict() == {"chrA": 20, "chrB": 8}
    sizes["chrA"] = 0  # mutate the returned copy...
    assert genome.chrom_sizes["chrA"] == 20  # ...the genome is unaffected

    assert genome.chromosomes == ["chrA", "chrB"]

    assert genome.fasta_path == genome.files.fasta
    assert genome.twobit_path == genome.files.twobit
    assert genome.chrom_sizes_path == genome.files.chrom_sizes
    # chrom_sizes_path is the on-disk file, distinct from the in-memory sizes Series.
    assert isinstance(genome.chrom_sizes_path, Path)
    assert genome.chrom_sizes_path.is_file()
    assert isinstance(genome.chrom_sizes, pd.Series)

    assert repr(genome) == "Genome('tiny', 2 sequences)"


@pytest.mark.parametrize(
    ("region", "message"),
    [
        ("chrZ:0-5", "unknown chromosome"),
        ("chrA:0-21", "exceeds chrA length"),
        ("chrA:10-5", "is past end"),
        ("chrA:bad", "malformed region"),
    ],
)
def test_invalid_region_raises(genome: Genome, region: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        genome.fetch_sequence(region)


def test_context_manager_closes_handle(genome: Genome) -> None:
    with genome as g:
        assert g.fetch_sequence("chrA:0-4") == DNA("ACGT")
    with pytest.raises(ValueError, match="closed"):
        g.fetch_sequence("chrA:0-4")


_OVERRIDE = AssemblyMetadata(
    assembly_name="tinyAsm",
    species="Testus minimus",
    ucsc_name="tinyUcsc",
    ncbi_name="TINY.1",
    ncbi_assembly_id="GCF_000000000.0",
    ncbi_taxid=1,
    source_url="https://mirror.example.org/references/tiny.fa.gz",
    sha256="00ff" * 16,
    intron_length_cap=1234,
    intron_length_cap_rationale="a number this test made up",
)


def test_metadata_resolves_override_curated_table_and_unknown_fallback(
    prepared_dir: Path,
) -> None:
    # sacCer3 is listed in the shipped table; a record given here wins over every field of it.
    with Genome("sacCer3", cache_dir=prepared_dir, metadata=_OVERRIDE) as g:
        assert g.metadata == _OVERRIDE
        assert [getattr(g.metadata, field) for field in METADATA_FIELDS] == [
            "tinyAsm",
            "Testus minimus",
            "tinyUcsc",
            "TINY.1",
            "GCF_000000000.0",
            1,
            "https://mirror.example.org/references/tiny.fa.gz",
            "00ff" * 16,
            1234,
            "a number this test made up",
        ]

    # Absent an override, it falls back to the curated table...
    with Genome("sacCer3", cache_dir=prepared_dir) as g:
        assert g.metadata == lookup_assembly("sacCer3")
        assert g.metadata.assembly_name == "sacCer3"
        assert g.metadata.species == "Saccharomyces cerevisiae"
        assert g.metadata.ucsc_name == "sacCer3"
        assert g.metadata.ncbi_name == "R64-1-1"
        assert g.metadata.ncbi_assembly_id == "GCF_000146045.2"
        assert g.metadata.ncbi_taxid == 559292
        assert g.metadata.source_url == (
            "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz"
        )
        assert (
            g.metadata.sha256 == "6ff72f079c3268431fc514a1a88730f8290e717663d343fa8a3590af65c422c3"
        )

    # ...and an assembly the table does not list still constructs: the table is a
    # cross-reference, not an allow-list, so every identifier is simply unknown — and
    # unknown is a record whose fields are unknown, never a missing record, so reading
    # one off the genome needs no guard.
    with Genome("tiny", cache_dir=prepared_dir) as g:
        assert g.metadata == AssemblyMetadata.unknown("tiny")
        assert g.metadata.assembly_name == "tiny"
        assert all(getattr(g.metadata, field) is None for field in METADATA_FIELDS[1:])
        assert g.chromosomes == ["chrA", "chrB"]
        # The genome's metadata is total; *does the curated table list this name?* is
        # not, and it is a different question — the one that separates a chimera's
        # derived name from a free-form local key on a machine holding neither
        # (ADR-0003). An unknown record answering the second would make every name a
        # listed one.
        assert lookup_assembly(g.assembly) is None


#: A record injected instead of the shipped table's row: it pins a URL, so nothing
#: contacts UCSC to validate the name, and no checksum, so the fixture FASTA served in
#: place of a download has nothing to disagree with.
_UNPINNED = AssemblyMetadata(
    assembly_name="hg38",
    species="Homo sapiens",
    ucsc_name="hg38",
    ncbi_name="GRCh38",
    ncbi_assembly_id="GCF_000001405.40",
    ncbi_taxid=9606,
    source_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
)


def test_registering_an_assembly_records_it_and_reopening_costs_no_fetch(
    fake_fetch: FakeFetch, tmp_path: Path
) -> None:
    # The one end-to-end pass: a real download (from the fixture), the real native
    # tools, and the record that says it finished.
    fake_fetch.serve("tiny.fa.gz")
    unpinned = _UNPINNED

    with Genome("hg38", cache_dir=tmp_path, progressbar=False, metadata=unpinned) as first:
        assert first.chromosomes == ["chrI", "chrII", "chrIII"]
        source_url, opening = first.metadata.source_url, first.fetch_sequence("chrI:0-4")

    record = read_record(tmp_path)
    assert record is not None
    assert record.source_url == source_url
    assert sorted(record.files) == ["hg38.2bit", "hg38.chrom.sizes", "hg38.fa", "hg38.fa.fai"]
    assert all((tmp_path / name).stat().st_size == size for name, size in record.files.items())
    assert list(tmp_path.rglob("*.gz")) == []  # the archive went with the working area

    with Genome("hg38", cache_dir=tmp_path, progressbar=False, metadata=unpinned) as again:
        assert again.fetch_sequence("chrI:0-4") == opening

    assert len(fake_fetch.calls) == 1  # the record answered; nothing was fetched twice


#: An annotation row injected instead of the shipped table's: it points at the committed
#: fixture the fake fetch serves and pins that file's digest.
_ANNOTATION = AnnotationMetadata(
    assembly="tiny",
    name="ensgene_v101",
    provider="UCSC",
    version="ensGene.v101",
    url="https://mirror.example.invalid/annotations/tiny.gtf.gz",
    sha256="255f43bd9abef76424d1c2d89a40cccc1a36215409bbc8f32dcead49ca3baf5e",
    default=True,
)


def test_registering_an_annotation_by_name_adopts_it_and_survives_reopening(
    fake_fetch: FakeFetch, yeast_dir: Path
) -> None:
    fake_fetch.serve("tiny.gtf.gz")

    with Genome("tiny", cache_dir=yeast_dir) as g:
        assert g.annotations.registered == []
        annotation = g.annotations.register("ensgene_v101", progressbar=False, metadata=_ANNOTATION)
        assert g.annotations.registered == ["ensgene_v101"]
        assert g.default_gtf == "ensgene_v101"
        assert g.default_gtf_path == annotation.gtf

    # A record was written, so the annotation is registered for anyone who opens the
    # genome next — and opening it fetches nothing.
    with Genome("tiny", cache_dir=yeast_dir) as again:
        assert again.annotations.registered == ["ensgene_v101"]
    assert len(fake_fetch.calls) == 1


def test_annotations_is_the_registry_itself_deliberately_not_a_list(
    prepared_dir: Path,
) -> None:
    # Callers hold one now, so it is a package export rather than an io internal.
    with Genome("tiny", cache_dir=prepared_dir) as g:
        assert isinstance(g.annotations, AnnotationRegistry)
        assert g.annotations.assembly == "tiny"

        # No list protocol, argued rather than overlooked: a registry settles a
        # four-way state — registered, broken, offered, nothing — and a dunder over
        # any one of them would hide which set a reader is walking. So every stale
        # `.annotations` use fails loudly...
        registry = g.annotations
        for stale in (len, iter, sorted, list):
            with pytest.raises(TypeError):
                stale(registry)  # type: ignore[arg-type,call-overload]

        # ...except truthiness, which is the one silent break in the change: an object is
        # always truthy, so `if genome.annotations:` reads as "has annotations" and is
        # not. Pinned here because the CHANGELOG names it and nothing else can catch it.
        assert bool(registry) is True
        assert registry.registered == []

    # A gffutils build killed part-way leaves a database and no record. Opening the
    # genome must not report it as an annotation, whatever files are lying there.
    directory = prepared_dir / "gtf" / "halfway"
    directory.mkdir(parents=True)
    (directory / "halfway.db").write_bytes(b"half a database")

    with Genome("tiny", cache_dir=prepared_dir) as g:
        assert g.annotations.registered == []
        assert g.default_gtf is None


def test_opening_a_genome_over_a_damaged_registration_raises_until_forced(
    fake_fetch: FakeFetch, tmp_path: Path
) -> None:
    # A file truncated after registration surfaces here, naming the file and the
    # command that fixes it — rather than as a confusing failure from deep inside
    # py2bit or an aligner much later.
    fake_fetch.serve("tiny.fa.gz")
    with Genome("hg38", cache_dir=tmp_path, progressbar=False, metadata=_UNPINNED):
        pass
    (tmp_path / "hg38.2bit").write_text("")

    with pytest.raises(RegistrationMismatchError) as excinfo:
        Genome("hg38", cache_dir=tmp_path, progressbar=False, metadata=_UNPINNED)

    message = str(excinfo.value)
    assert "hg38.2bit" in message
    assert "genome register hg38 --force" in message

    register_assembly("hg38", cache_dir=tmp_path, force=True, progressbar=False, metadata=_UNPINNED)

    with Genome("hg38", cache_dir=tmp_path, progressbar=False, metadata=_UNPINNED) as repaired:
        assert repaired.fetch_sequence("chrI:0-4") == DNA("CCAC")


class TestOfferedAgainstRegistered:
    """What the assembly offers, what is registered here, and which one is the default.

    sacCer3 is the assembly for the offered half: the shipped table lists it exactly one
    annotation and flags it as the default, so these assert what a lab member actually
    gets rather than what a table stood up for the test would give. ``tiny`` is in
    neither table, which is how the no-flag fallback is exercised.
    """

    def test_the_tables_flag_is_the_default_and_registering_the_default_starts_no_fetch(
        self, fake_fetch: FakeFetch, prepared_dir: Path
    ) -> None:
        with Genome("sacCer3", cache_dir=prepared_dir) as g:
            offered = g.annotations.offered
            assert [record.name for record in offered] == ["ensgene_v101"]
            assert offered[0].provider == "UCSC"
            # The lab supports it; this machine does not have it.
            assert g.annotations.registered == []
            assert g.default_gtf == "ensgene_v101"

        # A GENCODE registration is a gigabyte download and a database build running
        # many minutes. Naming the default must never start one.
        assert fake_fetch.calls == []
        assert not (prepared_dir / "gtf").exists()

        with (
            Genome("sacCer3", cache_dir=prepared_dir) as g,
            pytest.raises(AnnotationNotRegisteredError) as excinfo,
        ):
            _ = g.default_gtf_path

        assert "genome register-annotation sacCer3 ensgene_v101" in str(excinfo.value)

    def test_an_explicit_default_beats_the_flag_and_need_not_start_registered(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        with Genome("sacCer3", cache_dir=yeast_dir, default_gtf="mine") as g:
            g.annotations.register_path(data_dir / "tiny.gtf", "mine")

            assert [record.name for record in g.annotations.offered] == ["ensgene_v101"]
            assert g.default_gtf == "mine"
            assert g.default_gtf_path == g.annotations.path("mine")

        # One rule for both sources: naming a default is an intention, and the path is
        # where it has to exist. Registering it afterwards is all it takes.
        with Genome("sacCer3", cache_dir=yeast_dir / "elsewhere", default_gtf="mine") as g:
            assert g.default_gtf == "mine"
            with pytest.raises(AnnotationNotRegisteredError, match="register"):
                _ = g.default_gtf_path

            g.annotations.register_path(data_dir / "tiny.gtf", "mine")

            assert g.default_gtf_path == g.annotations.path("mine")

    def test_the_sole_registered_annotation_is_default_only_until_a_second_is_registered(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        # "tiny" is in neither table, so no flag decides anything and the older rule
        # stands — both as it is registered and on reopening.
        with Genome("tiny", cache_dir=yeast_dir) as g:
            assert g.annotations.offered == []
            g.annotations.register_path(data_dir / "tiny.gtf", "mine")
            assert g.default_gtf == "mine"

        with Genome("tiny", cache_dir=yeast_dir) as again:
            assert again.default_gtf == "mine"
            assert again.default_gtf_path == again.annotations.path("mine")

        with Genome("tiny", cache_dir=yeast_dir / "two") as g:
            g.annotations.register_path(data_dir / "tiny.gtf", "one")
            g.annotations.register_path(data_dir / "tiny.gtf", "two")

        with Genome("tiny", cache_dir=yeast_dir / "two") as again_two:
            assert sorted(again_two.annotations.registered) == ["one", "two"]
            assert again_two.default_gtf is None
            assert again_two.default_gtf_path is None

    def test_a_broken_annotation_never_stops_the_genome_opening_but_names_its_repair(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        # The invariant: one annotation nobody can trust must not cost the genome, nor
        # the annotations beside it. It is reported, not raised over.
        register_gtf("sacCer3", data_dir / "tiny.gtf", "healthy", cache_dir=yeast_dir)
        register_gtf("sacCer3", data_dir / "tiny.gtf", "ensgene_v101", cache_dir=yeast_dir)
        record_path(annotation_dir(yeast_dir, "ensgene_v101")).unlink()

        with Genome("sacCer3", cache_dir=yeast_dir) as g:
            assert g.annotations.registered == ["healthy"]
            assert [broken.name for broken in g.annotations.broken] == ["ensgene_v101"]
            assert g.fetch_sequence("chrI:0-4") == DNA("CCAC")

            # Not `genome register-annotation sacCer3 ensgene_v101`, which would itself
            # raise and demand --force: the command named here is the one that works.
            with pytest.raises(AnnotationNotRegisteredError) as excinfo:
                _ = g.default_gtf_path

        assert "genome register-annotation sacCer3 ensgene_v101 --force" in str(excinfo.value)

    def test_repairing_a_broken_annotation_stops_reporting_it(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        register_gtf("tiny", data_dir / "tiny.gtf", "mine", cache_dir=yeast_dir)
        record_path(annotation_dir(yeast_dir, "mine")).unlink()

        with Genome("tiny", cache_dir=yeast_dir) as g:
            assert [broken.name for broken in g.annotations.broken] == ["mine"]

            g.annotations.register_path(data_dir / "tiny.gtf", "mine", force=True)

            assert g.annotations.broken == []
            assert g.annotations.registered == ["mine"]
            assert g.annotations.path("mine").is_file()

    def test_a_broken_directory_and_an_unknown_name_each_name_their_own_repair(
        self, yeast_dir: Path, data_dir: Path, prepared_dir: Path
    ) -> None:
        # A genome knows which assembly it is, so the repair its registry names is a
        # command a shell can run rather than a call with the assembly left blank.
        directory = annotation_dir(yeast_dir, "mine")
        directory.mkdir(parents=True)
        (directory / "mine.db").write_bytes(b"half a database")
        source = data_dir / "tiny.gtf"

        with (
            Genome("tiny", cache_dir=yeast_dir) as g,
            pytest.raises(UnfinishedRegistrationError) as excinfo,
        ):
            g.annotations.register_path(source, "mine")

        assert f"genome register-gtf tiny {source} mine --force" in str(excinfo.value)

        with Genome("sacCer3", cache_dir=prepared_dir) as g, pytest.raises(KeyError) as unknown:
            g.annotations.path("no_such_annotation")

        message = str(unknown.value)
        assert "no_such_annotation" in message
        assert "ensgene_v101" in message  # what the table does offer
        assert "register_path" in message  # and the way in for one it does not


def test_path_or_url_seeds_from_local_fasta(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A local FASTA seed prepares the assembly end-to-end (copy + native tools)
    # without ever contacting UCSC.
    src = tmp_path / "src.fa"
    src.write_text(f">chrA\n{_CHR_A}\n>chrB\n{_CHR_B}\n")

    def boom(self: object, **_kwargs: object) -> None:
        raise AssertionError("UCSC fetch_genome must not run when path_or_url is given")

    monkeypatch.setattr(genome_mod.UCSCGenomeDownloader, "fetch_genome", boom)

    with Genome("tiny", path_or_url=src, cache_dir=tmp_path / "cache") as g:
        assert g.fetch_sequence("chrA:0-8") == DNA("ACGTACGT")
        assert g.chromosomes == ["chrA", "chrB"]


class TestGeneCategories:
    """``Genome.gene_list`` / ``.gene_lists`` — the everyday way to a category's genes.

    Both delegate to the registry, so what is asserted here is that they reach it and
    that the answer arrives whole; what the answer *is* belongs to the annotation tests. The
    shipped
    curated list answers, so the categories are read off it rather than named.
    """

    def _declared(self, annotation: str) -> tuple[str, ...]:
        """The categories the shipped curated list declares for ``annotation``."""
        listed = curated_gene_list(annotation)
        assert listed is not None, f"no curated gene list ships for {annotation}"
        return tuple(listed.categories)

    def test_a_genome_answers_for_its_default_annotation_and_says_what_may_be_asked_for(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        with Genome("sacCer3", cache_dir=yeast_dir) as g:
            g.annotations.register_path(data_dir / "tiny.gtf", "ensgene_v101")
            category = self._declared("ensgene_v101")[0]

            answer = g.gene_list(category)

            assert g.default_gtf == "ensgene_v101"
            assert (answer.assembly, answer.annotation, answer.category) == (
                "sacCer3",
                "ensgene_v101",
                category,
            )
            assert answer.gene_ids == g.annotations.gene_list(category).gene_ids

            answers = g.gene_lists()

            assert [entry.category for entry in answers] == list(self._declared("ensgene_v101"))
            assert all(entry.gene_ids for entry in answers)

    def test_an_annotation_that_cannot_answer_raises_rather_than_answering_emptily(
        self, yeast_dir: Path, data_dir: Path
    ) -> None:
        # The whole point of the surface: a caller must be able to tell *nothing is known
        # about this annotation's categories* from *there are none of these genes*.
        with Genome("tiny", cache_dir=yeast_dir) as g:
            g.annotations.register_path(data_dir / "tiny.gtf", "mine")

            with pytest.raises(NoGeneCategoriesError):
                g.gene_lists()
            with pytest.raises(NoGeneCategoriesError):
                g.gene_list("rRNA")


def test_tf_gene_list_and_tf_cofactor_list_delegate_to_the_registry(yeast_dir: Path) -> None:
    # One line each on the genome object, and both reach the registry bound to this
    # assembly: no census and no cofactor table ships for yeast, which is the
    # registry's answer and not the genome's.
    with Genome("sacCer3", cache_dir=yeast_dir) as g:
        with pytest.raises(NoTFCensusError):
            g.tf_gene_list()
        with pytest.raises(NoCofactorTableError):
            g.tf_cofactor_list()
