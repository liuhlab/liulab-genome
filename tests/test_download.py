"""Tests for genome.io.download.

Every download in the package goes through ``fetch.fetch_url``, so the suite stays
offline by replacing that one function with the shared ``fake_fetch`` fixture (see
tests/conftest.py) and asserting the arguments each caller wires through. The fetch step
itself lives in its own module now and is tested there — see test_fetch.

Metadata is injected as an in-memory :class:`AssemblyMetadata` record, so no test here
reads or fakes the shipped TSV; the table itself is tested in test_metadata.
"""

from __future__ import annotations

import shutil
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

import pooch
import pytest
import requests

from genome import __version__
from genome.io import download as download_mod
from genome.io.completion import (
    CompletionRecord,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
    record_path,
    work_dir,
)
from genome.io.download import (
    EXPECTED_FROM_RECORD,
    EXPECTED_FROM_TABLE,
    RegisteredAssembly,
    UCSCGenomeDownloader,
    VerifiedAssembly,
    assembly_data_dir,
    assembly_table_row,
    liulab_data_dir,
    register_assembly,
    verify_assembly,
)
from genome.io.fasta import PREPARATION_TOOLS, GenomeFiles
from genome.io.source import FetchedSource
from genome.io.utils import ChecksumMismatchError, sha256_file
from genome.metadata import AssemblyMetadata, format_table_row

from .conftest import FakeFetch
from .test_source import _module_level_imports

#: sha256 of the committed ``tiny.fa`` — the *unpacked* bytes ``tiny.fa.gz`` yields.
_TINY_FA_SHA256 = "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"

#: A URL that is nothing like the golden path, so using it can only come from a row.
_PINNED_URL = "https://mirror.example.org/references/tiny.fa.gz"


@dataclass
class _FakeResponse:
    """Minimal stand-in for :class:`requests.Response` for ``head`` stubs."""

    status_code: int = 200

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code} error")


@dataclass
class _HeadRecorder:
    """Records ``requests.head`` calls and returns a configurable response."""

    status_code: int = 200
    calls: list[dict[str, object]] = field(default_factory=list)

    def __call__(self, url: str, **kwargs: object) -> _FakeResponse:
        self.calls.append({"url": url, **kwargs})
        return _FakeResponse(self.status_code)


@pytest.fixture(autouse=True)
def head_recorder(monkeypatch: pytest.MonkeyPatch) -> _HeadRecorder:
    """Patch ``requests.head`` so assembly validation stays offline (200 by default).

    Autouse: every test in this module runs without real network I/O. Tests that
    care about validation request this fixture to inspect calls or set the
    returned status code.
    """
    recorder = _HeadRecorder()
    monkeypatch.setattr(download_mod.requests, "head", recorder)
    return recorder


def _row(*, source_url: str | None = None, sha256: str | None = None) -> AssemblyMetadata:
    """An in-memory metadata record for the ``tiny`` assembly."""
    return AssemblyMetadata(
        assembly_name="tiny",
        species="Testus minimus",
        ucsc_name="tiny",
        ncbi_name="TINY.1",
        ncbi_assembly_id="GCF_000000000.0",
        ncbi_taxid=1,
        source_url=source_url,
        sha256=sha256,
    )


def _derive(fasta: Path) -> GenomeFiles:
    """Write what a real preparation derives from ``fasta`` and return the whole set.

    Stands in for the native tools: the three companion files exist and are named as
    ``prepare_fasta`` names them, which is all a completion record is claiming.
    """
    files = GenomeFiles(
        fasta=fasta,
        fai=fasta.with_name(fasta.name + ".fai"),
        twobit=fasta.with_name(fasta.stem + ".2bit"),
        chrom_sizes=fasta.with_name(fasta.stem + ".chrom.sizes"),
    )
    for derived in (files.fai, files.twobit, files.chrom_sizes):
        derived.write_text(f"derived from {fasta.name}\n")
    return files


@pytest.fixture
def no_native_prepare(monkeypatch: pytest.MonkeyPatch) -> None:
    """Stub the FASTA preparation step so registration runs without native binaries."""

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        return _derive(Path(fasta_path))

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)


def test_liulab_data_dir_cache_dir_and_ucsc_urls_are_derived_correctly(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, liulab_data: Path
) -> None:
    assert assembly_data_dir("hg38") == liulab_data / "genome" / "hg38"
    assert UCSCGenomeDownloader("mm39").cache_dir == liulab_data / "genome" / "mm39"

    assert (
        UCSCGenomeDownloader("hg38").fasta_url
        == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    )
    assert UCSCGenomeDownloader("mm39").fasta_url.endswith("mm39/bigZips/mm39.fa.gz")
    assert (
        UCSCGenomeDownloader("hg38").assembly_url
        == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
    )

    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "lab"))
    assert liulab_data_dir() == tmp_path / "lab"

    monkeypatch.delenv("LIULAB_DATA", raising=False)
    assert liulab_data_dir() == Path.home() / "liulab_data"

    monkeypatch.setenv("LIULAB_DATA", "")  # blank counts as unset
    assert liulab_data_dir() == Path.home() / "liulab_data"


def test_validate_assembly_checks_the_directory_url_and_tells_404_from_other_errors(
    head_recorder: _HeadRecorder,
) -> None:
    dl = UCSCGenomeDownloader("hg38")
    dl.validate_assembly()  # 200 by default -> no raise
    assert head_recorder.calls[0]["url"] == dl.assembly_url

    head_recorder.status_code = 404
    with pytest.raises(ValueError, match="Unknown UCSC assembly 'nope99'"):
        UCSCGenomeDownloader("nope99").validate_assembly()

    head_recorder.status_code = 500
    with pytest.raises(requests.exceptions.HTTPError):
        UCSCGenomeDownloader("hg38").validate_assembly()


def test_the_fetch_and_prepare_pipeline_validates_downloads_decompresses_and_forwards_kwargs(
    fake_fetch: FakeFetch,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    data_dir: Path,
    head_recorder: _HeadRecorder,
) -> None:
    # "tiny" is in no table, so the URL is derived and the name is checked against UCSC —
    # and validation happens before a download, so a bad assembly name aborts rather
    # than fetching a second time.
    fake_fetch.serve("tiny.fa.gz")
    UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "validate").fetch_fasta()
    assert len(head_recorder.calls) == 1

    head_recorder.status_code = 404
    with pytest.raises(ValueError, match="Unknown UCSC assembly"):
        UCSCGenomeDownloader("bad", cache_dir=tmp_path / "bad").fetch_fasta()
    assert len(fake_fetch.calls) == 1
    head_recorder.status_code = 200

    # Decompresses by default — both the archive and its unpacked form land in the
    # working area, never among the assembly's own files — but can keep the archive.
    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path / "decompress")
    result = dl.fetch_fasta()
    assert fake_fetch.last.url == dl.fasta_url
    assert result == work_dir(tmp_path / "decompress") / "hg38.fa"
    assert result.read_text() == (data_dir / "tiny.fa").read_text()
    assert (work_dir(tmp_path / "decompress") / "hg38.fa.gz").is_file()
    assert not (tmp_path / "decompress" / "hg38.fa").exists()
    archived = dl.fetch_fasta(decompress=False)
    assert fake_fetch.last.processor is None
    assert archived.name.endswith(".fa.gz")

    # fetch_genome chains fetch_fasta (network) and prepare_fasta (native tools, stubbed
    # here) and forwards both overwrite and known_hash through to each.
    prepared: dict[str, object] = {}

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        prepared["fasta"] = fasta_path
        prepared["overwrite"] = overwrite
        return _derive(Path(fasta_path))

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    # An unpinned record: this is about the pipeline, not about verification, and every
    # shipped row pins a checksum the fixture deliberately is not.
    dl2 = UCSCGenomeDownloader("hg38", cache_dir=tmp_path / "genome", metadata=_row())
    files = dl2.fetch_genome()

    assert fake_fetch.last.url == dl2.fasta_url
    assert prepared["fasta"] == tmp_path / "genome" / "hg38.fa"
    assert prepared["overwrite"] is False  # default: caches reused
    assert files.fasta == tmp_path / "genome" / "hg38.fa"
    assert files.fai == tmp_path / "genome" / "hg38.fa.fai"
    assert files.twobit == tmp_path / "genome" / "hg38.2bit"
    assert files.chrom_sizes == tmp_path / "genome" / "hg38.chrom.sizes"

    dl2.fetch_genome(overwrite=True)
    assert prepared["overwrite"] is True

    dl3 = UCSCGenomeDownloader("hg38", cache_dir=tmp_path / "again", metadata=_row())
    dl3.fetch_genome(known_hash="md5:abc")
    call = fake_fetch.last
    assert call.url == dl3.fasta_url
    assert call.known_hash == "md5:abc"
    # the pipeline always decompresses, so a Decompress processor is selected.
    assert isinstance(call.processor, pooch.Decompress)


# --- a row's pinned source and checksum -------------------------------------


def test_a_pinned_source_wins_over_no_row_and_no_row_still_uses_the_golden_path(
    fake_fetch: FakeFetch, tmp_path: Path, head_recorder: _HeadRecorder, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    dl.fetch_fasta()

    assert dl.fasta_url == _PINNED_URL
    assert fake_fetch.last.url == _PINNED_URL
    # Validation is a property of the source, and a pinned URL *is* the source, so
    # there is nothing left to guess about the name (ADR-0003).
    assert head_recorder.calls == []

    # hg38 has a row of its own; an explicit record replaces it wholesale.
    dl2 = UCSCGenomeDownloader(
        "hg38", cache_dir=tmp_path / "hg38", metadata=_row(source_url=_PINNED_URL)
    )
    dl2.fetch_fasta()
    assert fake_fetch.last.url == _PINNED_URL

    # The table is a cross-reference, not an allow-list: no row takes nothing away. What
    # the downloader holds is total — an unlisted assembly knows its own name and nothing
    # else — so no step here asks whether there is a record before reading a field off it.
    dl3 = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "no-row")
    assert dl3.metadata == AssemblyMetadata.unknown("tiny")

    files = dl3.fetch_genome()

    assert dl3.fasta_url == "https://hgdownload.soe.ucsc.edu/goldenPath/tiny/bigZips/tiny.fa.gz"
    assert fake_fetch.last.url == dl3.fasta_url
    assert head_recorder.calls[-1]["url"] == dl3.assembly_url  # still validated
    assert files.fasta == tmp_path / "no-row" / "tiny.fa"

    # The two questions stay two, though. What is *known about* "my_ref" is total —
    # unknown, not missing — but whether the curated table *lists* it is still `None`,
    # and that is what tells `my_ref` — one name somebody chose — from a chimera of
    # `my` and `ref` (ADR-0003, ADR-0008).
    other = UCSCGenomeDownloader("my_ref", cache_dir=tmp_path / "other")
    assert other.metadata == AssemblyMetadata.unknown("my_ref")
    assert other._source() == FetchedSource(
        url="https://hgdownload.soe.ucsc.edu/goldenPath/my_ref/bigZips/my_ref.fa.gz",
        derived=True,
    )


def test_registering_checks_the_fasta_against_the_pinned_checksum_and_verify_fasta_agrees(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader(
        "tiny",
        cache_dir=tmp_path,
        metadata=_row(source_url=_PINNED_URL, sha256=_TINY_FA_SHA256),
    )
    files = dl.fetch_genome()
    assert files.fasta.is_file()

    wrong = "0" * 64
    dl_wrong = UCSCGenomeDownloader(
        "tiny", cache_dir=tmp_path / "wrong", metadata=_row(source_url=_PINNED_URL, sha256=wrong)
    )
    with pytest.raises(ChecksumMismatchError) as excinfo:
        dl_wrong.fetch_genome()
    message = str(excinfo.value)
    assert wrong in message  # what the row expected...
    assert _TINY_FA_SHA256 in message  # ...and what actually arrived

    # The whole point of hashing unpacked content: pinning the .fa.gz's digest — which is
    # what pooch's known_hash would check — fails against the FASTA inside it.
    archive_digest = sha256_file(data_dir / "tiny.fa.gz")
    assert archive_digest != _TINY_FA_SHA256
    dl_archive = UCSCGenomeDownloader(
        "tiny",
        cache_dir=tmp_path / "archive",
        metadata=_row(source_url=_PINNED_URL, sha256=archive_digest),
    )
    with pytest.raises(ChecksumMismatchError):
        dl_archive.fetch_genome()

    # And verify_fasta reports that same computed digest — defaulting to the FASTA the
    # assembly registered, in the assembly dir, not whatever a download left behind —
    # exactly as the table row's own checksum column does.
    assert files.fasta.is_file()
    assert dl.verify_fasta(files.fasta) == _TINY_FA_SHA256
    assert dl.verify_fasta() == _TINY_FA_SHA256

    row = assembly_table_row("hg38", cache_dir=tmp_path / "hg38", progressbar=False)
    assert row.sha256 == _TINY_FA_SHA256
    assert row.source_url == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    assert row.ncbi_name == "GRCh38"  # the curated identifiers are carried through
    # It is the table's own row type, so the line to paste is the table's own formatting
    # of it and the columns cannot drift from the ones the table is parsed through.
    assert format_table_row(asdict(row)).split("\t")[0] == "hg38"

    # An assembly the table has never heard of still gets a checksum; only the
    # identifiers a person alone could supply are left blank.
    new_row = assembly_table_row("newAsm", cache_dir=tmp_path / "new", progressbar=False)
    assert new_row.assembly_name == "newAsm"
    assert new_row.species is None
    assert str(new_row.source_url).endswith("/newAsm/bigZips/newAsm.fa.gz")
    assert new_row.sha256 == _TINY_FA_SHA256


# --- the record a finished registration writes -------------------------------


def test_a_finished_registration_writes_a_record_of_what_it_did(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    files = dl.fetch_genome()

    record = read_record(tmp_path)
    assert record is not None
    assert (record.kind, record.name) == ("genome", "tiny")
    assert record.source_url == _PINNED_URL
    assert record.sha256 == _TINY_FA_SHA256
    assert record.package_version == __version__
    assert datetime.fromisoformat(record.completed_at).utcoffset() is not None
    # every file it claims, with the size it is on disk, addressed relative to the
    # assembly dir so the whole directory can be moved
    assert record.files == {
        "tiny.fa": files.fasta.stat().st_size,
        "tiny.fa.fai": files.fai.stat().st_size,
        "tiny.2bit": files.twobit.stat().st_size,
        "tiny.chrom.sizes": files.chrom_sizes.stat().st_size,
    }
    assert set(record.tool_versions) <= set(PREPARATION_TOOLS)


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_the_record_names_the_versions_of_the_tools_that_prepared_it(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")

    UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row()).fetch_genome()

    record = read_record(tmp_path)
    assert record is not None
    assert record.tool_versions["samtools"].startswith("samtools")


def test_the_record_governs_reopening_when_its_written_and_the_archive_it_replaces(
    fake_fetch: FakeFetch,
    tmp_path: Path,
    no_native_prepare: None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))
    first = dl.fetch_genome()

    again = UCSCGenomeDownloader(
        "tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL)
    ).fetch_genome()

    assert again == first
    assert len(fake_fetch.calls) == 1  # the record answered; nothing was fetched

    dl.fetch_genome(overwrite=True)
    assert len(fake_fetch.calls) == 2

    # A record is written only once every derived file it will claim already exists.
    seen: dict[str, bool] = {}

    def watching_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        seen["record_exists_during_preparation"] = record_path(tmp_path / "watched").exists()
        return _derive(Path(fasta_path))

    monkeypatch.setattr(download_mod, "prepare_fasta", watching_prepare)
    UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "watched", metadata=_row()).fetch_genome()

    assert seen["record_exists_during_preparation"] is False
    assert record_path(tmp_path / "watched").is_file()

    # A preempted job must repair from what it already downloaded rather than pull a
    # whole genome down again, so nothing is discarded until the record is written.
    def failing_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        raise RuntimeError("samtools fell over")

    monkeypatch.setattr(download_mod, "prepare_fasta", failing_prepare)
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "unfinished", metadata=_row())

    with pytest.raises(RuntimeError, match="fell over"):
        dl.fetch_genome()

    assert (work_dir(tmp_path / "unfinished") / "tiny.fa.gz").is_file()
    assert read_record(tmp_path / "unfinished") is None  # unfinished, nothing claims otherwise

    # Every prepared assembly used to keep its plain-gzip download forever, which no
    # external tool can even read in place — once the record IS written, it is gone.
    def working_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        return _derive(Path(fasta_path))

    monkeypatch.setattr(download_mod, "prepare_fasta", working_prepare)
    UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "finished", metadata=_row()).fetch_genome()

    assert list((tmp_path / "finished").rglob("*.gz")) == []
    assert not work_dir(tmp_path / "finished").exists()


# --- a broken registration raises; a forced one repairs ----------------------


def test_registration_raises_when_broken_and_succeeds_from_absent_empty_or_interrupted(
    fake_fetch: FakeFetch,
    tmp_path: Path,
    data_dir: Path,
    no_native_prepare: None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A run killed between placing the FASTA and writing the record. A complete build
    # whose record never landed and the wreckage of one killed half-way look identical
    # on disk, so neither is trusted and neither is quietly rebuilt.
    (tmp_path / "no-record" / "tiny.fa").parent.mkdir(parents=True)
    (tmp_path / "no-record" / "tiny.fa").write_text(">chrI\nACGT\n")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "no-record", metadata=_row())

    with pytest.raises(UnfinishedRegistrationError) as excinfo:
        dl.fetch_genome()
    message = str(excinfo.value)
    assert "tiny.fa" in message  # what is there
    assert "genome register tiny --force" in message  # and what to do about it
    assert fake_fetch.calls == []  # nothing was fetched over the top of it

    # The opposite defect: a record exists, but disk disagrees with what it claims.
    fake_fetch.serve("tiny.fa.gz")
    dl2 = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "disagrees", metadata=_row())
    files = dl2.fetch_genome()
    recorded = files.twobit.stat().st_size
    files.twobit.write_text("")  # truncated behind our back

    with pytest.raises(RegistrationMismatchError) as excinfo2:
        dl2.fetch_genome()
    message2 = str(excinfo2.value)
    assert f"tiny.2bit: recorded {recorded} bytes, found 0" in message2
    assert "genome register tiny --force" in message2
    assert len(fake_fetch.calls) == 1  # not silently registered again either

    fresh = tmp_path / "never-registered"
    files = UCSCGenomeDownloader("tiny", cache_dir=fresh, metadata=_row()).fetch_genome()
    assert files.fasta == fresh / "tiny.fa"
    assert read_record(fresh) is not None

    empty = tmp_path / "empty"
    files = UCSCGenomeDownloader("tiny", cache_dir=empty, metadata=_row()).fetch_genome()
    assert files.fasta.is_file()

    # The working area holds working state rather than claimed outputs, so a preempted
    # download is not a broken registration.
    interrupted = tmp_path / "interrupted"
    work_dir(interrupted).mkdir(parents=True)
    shutil.copy2(data_dir / "tiny.fa.gz", work_dir(interrupted) / "tiny.fa.gz")
    files = UCSCGenomeDownloader("tiny", cache_dir=interrupted, metadata=_row()).fetch_genome()
    assert files.fasta.is_file()

    # And a forced re-registration keeps a FASTA that still matches its pin, rebuilding
    # only the derived files a whole genome is not downloaded twice to replace.
    rebuilt: list[bool] = []

    def watching_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        rebuilt.append(overwrite)
        return _derive(Path(fasta_path))

    monkeypatch.setattr(download_mod, "prepare_fasta", watching_prepare)
    pinned_dir = tmp_path / "pinned"
    dl = UCSCGenomeDownloader(
        "tiny",
        cache_dir=pinned_dir,
        metadata=_row(source_url=_PINNED_URL, sha256=_TINY_FA_SHA256),
    )
    calls_before = len(fake_fetch.calls)
    dl.fetch_genome()
    (pinned_dir / "tiny.2bit").unlink()  # a derived file lost; the FASTA is still good

    kept = dl.fetch_genome(overwrite=True)

    assert len(fake_fetch.calls) == calls_before + 1  # a whole genome is not downloaded twice
    assert rebuilt == [False, True]  # only the derived files were rebuilt
    assert kept.twobit.is_file()


@pytest.mark.parametrize(
    ("pinned", "damage"),
    [
        (_TINY_FA_SHA256, "delete"),
        (_TINY_FA_SHA256, "corrupt"),
        (None, "none"),
    ],
    ids=["fasta-missing", "fasta-wrong", "nothing-pinned"],
)
def test_a_forced_re_registration_fetches_again_unless_the_fasta_is_provably_good(
    fake_fetch: FakeFetch,
    tmp_path: Path,
    no_native_prepare: None,
    pinned: str | None,
    damage: str,
) -> None:
    # With nothing pinned there is no way to prove the FASTA on disk is the right one,
    # so it is not trusted — the source is read again.
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader(
        "tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL, sha256=pinned)
    )
    files = dl.fetch_genome()
    if damage == "delete":
        files.fasta.unlink()
    elif damage == "corrupt":
        files.fasta.write_text(">chrI\nNNNN\n")

    dl.fetch_genome(overwrite=True)

    assert len(fake_fetch.calls) == 2
    assert files.fasta.read_text().startswith(">chrI")


def test_a_seeded_assembly_names_its_own_source_in_the_repair_and_records_it_when_it_succeeds(
    tmp_path: Path, data_dir: Path, no_native_prepare: None
) -> None:
    # `genome register tiny --force` would fetch from the golden path, which is not
    # where this assembly came from, so the repair names the source it was seeded from.
    source = data_dir / "tiny.fa.gz"
    broken = tmp_path / "broken"
    broken.mkdir()
    (broken / "tiny.fa").write_text(">chrI\nACGT\n")
    dl = UCSCGenomeDownloader("tiny", cache_dir=broken)

    with pytest.raises(UnfinishedRegistrationError) as excinfo:
        dl.fetch_genome_from(source)
    assert f"genome register tiny --force --source {source}" in str(excinfo.value)

    # And on a clean directory, that same seeding succeeds and records where it came from.
    cache = tmp_path / "cache"
    UCSCGenomeDownloader("tiny", cache_dir=cache).fetch_genome_from(source)
    record = read_record(cache)
    assert record is not None
    assert record.source_url == str(source)
    assert record.sha256 == _TINY_FA_SHA256  # what arrived, recorded rather than compared
    assert sorted(record.files) == ["tiny.2bit", "tiny.chrom.sizes", "tiny.fa", "tiny.fa.fai"]
    assert not work_dir(cache).exists()


# --- registering and verifying an assembly by name ---------------------------


def test_register_assembly_reports_serializes_repairs_and_seeds_from_a_source(
    fake_fetch: FakeFetch,
    tmp_path: Path,
    data_dir: Path,
    no_native_prepare: None,
    head_recorder: _HeadRecorder,
) -> None:
    fake_fetch.serve("tiny.fa.gz")

    registered = register_assembly("tiny", cache_dir=tmp_path / "normal", progressbar=False)

    assert registered.assembly == "tiny"
    assert registered.directory == tmp_path / "normal"
    assert registered.sha256 == _TINY_FA_SHA256
    assert registered.source_url == UCSCGenomeDownloader("tiny").fasta_url
    assert registered.record.files == {
        name: (tmp_path / "normal" / name).stat().st_size
        for name in ("tiny.fa", "tiny.fa.fai", "tiny.2bit", "tiny.chrom.sizes")
    }
    assert registered.file_names == ["tiny.2bit", "tiny.chrom.sizes", "tiny.fa", "tiny.fa.fai"]
    assert registered.chimera is None  # nothing was concatenated from anything

    # The `--json` payload is the completion record under its own on-disk key names,
    # with the two facts a record does not hold about itself. A type wraps those names;
    # it never renames them, because lab directories on shared storage are read by both.
    assert registered.as_json() == {
        **asdict(registered.record),
        "assembly": "tiny",
        "directory": str(tmp_path / "normal"),
    }
    assert list(registered.as_json())[-2:] == ["assembly", "directory"]

    # The command the error message names has to be the command that fixes it.
    (tmp_path / "broken" / "tiny.fa").parent.mkdir(parents=True)
    (tmp_path / "broken" / "tiny.fa").write_text("half a genome\n")
    with pytest.raises(UnfinishedRegistrationError, match="genome register tiny --force"):
        register_assembly("tiny", cache_dir=tmp_path / "broken", progressbar=False)
    repaired = register_assembly(
        "tiny", cache_dir=tmp_path / "broken", force=True, progressbar=False
    )
    assert repaired.sha256 == _TINY_FA_SHA256

    calls_before = len(head_recorder.calls)
    source = data_dir / "tiny.fa.gz"
    seeded = register_assembly(
        "tiny", source=source, cache_dir=tmp_path / "seeded", progressbar=False
    )
    assert seeded.source_url == str(source)
    # UCSC is never consulted about a seeded assembly, unlike the golden-path repair above.
    assert len(head_recorder.calls) == calls_before


def test_verify_assembly_expected_digest_comes_from_table_record_or_nothing_and_catches_bad_ones(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None, data_dir: Path
) -> None:
    # No row lists "tiny", so nothing pins a digest for it — and its own registration
    # wrote one down. A fallback rather than a question about what kind of assembly this
    # is: nothing was downloaded for a chimera either, and it takes the same path.
    fake_fetch.serve("tiny.fa.gz")
    register_assembly("tiny", cache_dir=tmp_path, progressbar=False)

    checked = verify_assembly("tiny", cache_dir=tmp_path)

    assert checked.fasta == tmp_path / "tiny.fa"
    assert checked.sha256 == _TINY_FA_SHA256
    assert checked.expected == _TINY_FA_SHA256
    assert checked.expected_from == "record"
    assert checked.verified is True

    # A row that pins its own digest is a stronger source than the record beside it,
    # even though both would answer here: a pin is what a digest is held to when there
    # is one.
    row = _row(source_url=_PINNED_URL, sha256=_TINY_FA_SHA256)
    register_assembly("tiny", cache_dir=tmp_path / "pinned", progressbar=False, metadata=row)

    checked_pinned = verify_assembly("tiny", cache_dir=tmp_path / "pinned", metadata=row)

    assert checked_pinned.verified is True
    assert checked_pinned.expected == _TINY_FA_SHA256
    assert checked_pinned.expected_from == "table"
    assert checked_pinned.sha256 == _TINY_FA_SHA256
    # And the payload a surface serializes is those same answers, `verified` written out
    # beside the field it is read from.
    assert checked_pinned.as_json() == {
        "assembly": "tiny",
        "fasta": str(tmp_path / "pinned" / "tiny.fa"),
        "sha256": _TINY_FA_SHA256,
        "expected": _TINY_FA_SHA256,
        "expected_from": "table",
        "verified": True,
        "components": None,
    }

    # The third state: no row lists "tiny", and nothing is registered here whose record
    # could answer either. Reported rather than raised — a digest with nothing to compare
    # against is still worth having, which is what `verified` says.
    checked_none = verify_assembly("tiny", fasta=data_dir / "tiny.fa", cache_dir=tmp_path / "none")

    assert checked_none.sha256 == _TINY_FA_SHA256
    assert checked_none.expected is None
    assert checked_none.expected_from is None
    assert checked_none.verified is False

    # The fallback is checked, not merely reported. One base flipped is the same number
    # of bytes, so the registration check — which compares sizes — passes, and what
    # catches it is the digest the record itself pins.
    register_assembly("tiny", cache_dir=tmp_path / "tampered", progressbar=False)
    fasta = tmp_path / "tampered" / "tiny.fa"
    raw = bytearray(fasta.read_bytes())
    raw[-2] = ord("A") if raw[-2] != ord("A") else ord("T")
    fasta.write_bytes(raw)

    with pytest.raises(ChecksumMismatchError, match=_TINY_FA_SHA256):
        verify_assembly("tiny", cache_dir=tmp_path / "tampered")

    # sacCer3's shipped row pins the real genome's digest and the fixture is a subsample
    # of it, so a hand-copied file checked against that row is caught before anything is
    # built on it. Nothing needs to be registered for this.
    with pytest.raises(ChecksumMismatchError) as excinfo:
        verify_assembly("sacCer3", fasta=data_dir / "tiny.fa", cache_dir=tmp_path / "other")
    assert _TINY_FA_SHA256 in str(excinfo.value)

    # And a directory whose files disagree with a claimed registration, or that claims
    # none at all, each raises naming what to run.
    (tmp_path / "broken" / "tiny.fa").parent.mkdir(parents=True)
    (tmp_path / "broken" / "tiny.fa").write_text("half a genome\n")
    with pytest.raises(UnfinishedRegistrationError, match="genome register tiny --force"):
        verify_assembly("tiny", cache_dir=tmp_path / "broken")

    with pytest.raises(FileNotFoundError, match="genome register tiny"):
        verify_assembly("tiny", cache_dir=tmp_path / "elsewhere")


# --- seeding from a user-provided FASTA (path_or_url) -----------------------


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("https://x/y.fa.gz", True),
        ("ftp://x/y.fa", True),
        ("ftps://x/y.fa", False),  # pooch ships no ftps downloader
        ("/data/ref.fa", False),
    ],
)
def test_looks_like_url(source: str, expected: bool) -> None:
    assert download_mod._looks_like_url(source) is expected


def test_materialize_fasta_raises_for_a_missing_local_file_or_an_unfetchable_scheme(
    tmp_path: Path,
) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    with pytest.raises(FileNotFoundError, match="local FASTA source not found"):
        dl._materialize_fasta(tmp_path / "nope.fa")
    with pytest.raises(ValueError, match="no downloader for the 'ftps' scheme"):
        dl._materialize_fasta("ftps://example.org/tiny.fa")


def test_materialize_fasta_copies_plain_decompresses_gz_and_reuses_unless_overwritten(
    tmp_path: Path, data_dir: Path
) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "cache")

    out = dl._materialize_fasta(data_dir / "tiny.fa")
    assert out == tmp_path / "cache" / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()

    dl_gz = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "cache_gz")
    out_gz = dl_gz._materialize_fasta(data_dir / "tiny.fa.gz")
    assert out_gz == tmp_path / "cache_gz" / "tiny.fa"
    assert out_gz.read_text() == (data_dir / "tiny.fa").read_text()
    # the compressed source stays in the working area, never beside the prepared FASTA
    assert (work_dir(tmp_path / "cache_gz") / "tiny.fa.gz").is_file()
    assert not (tmp_path / "cache_gz" / "tiny.fa.gz").exists()

    src = tmp_path / "src.fa"
    src.write_text("ONE")
    dl_reuse = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "reuse")
    dl_reuse._materialize_fasta(src)
    src.write_text("TWO")
    # a fresh <assembly>.fa is reused without re-reading the source...
    assert dl_reuse._materialize_fasta(src).read_text() == "ONE"
    # ...unless overwrite forces a refresh.
    assert dl_reuse._materialize_fasta(src, overwrite=True).read_text() == "TWO"


def test_materialize_fasta_from_a_url_fetches_decompresses_and_overwrite_refetches(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)

    out = dl._materialize_fasta("https://example.org/whatever.fa", progressbar=False)

    assert out == tmp_path / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()
    call = fake_fetch.last
    assert call.url == "https://example.org/whatever.fa"
    assert call.dest_dir == work_dir(tmp_path)  # downloaded into the working area...
    assert call.fname == "tiny.fa"  # ...as <assembly>.fa, then moved into place
    assert call.progressbar is False

    fake_fetch.serve("tiny.fa.gz")
    dl_gz = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "gz")
    out_gz = dl_gz._materialize_fasta("https://example.org/whatever.fa.gz")
    assert fake_fetch.last.fname == "tiny.fa.gz"
    assert out_gz == tmp_path / "gz" / "tiny.fa"
    assert out_gz.read_text() == (data_dir / "tiny.fa").read_text()

    calls_before = len(fake_fetch.calls)
    dl_overwrite = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "overwrite")
    url = "https://example.org/whatever.fa"
    dl_overwrite._materialize_fasta(url)
    dl_overwrite._materialize_fasta(url, overwrite=True)
    # overwrite discards the kept download, so the source is fetched again rather than
    # served from disk.
    assert len(fake_fetch.calls) == calls_before + 2


def test_fetch_genome_from_chains_materialize_and_prepare(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    seen: dict[str, object] = {}

    def fake_materialize(
        self: UCSCGenomeDownloader,
        source: object,
        *,
        progressbar: bool = True,
        overwrite: bool = False,
    ) -> Path:
        seen["source"] = source
        seen["overwrite"] = overwrite
        fasta = tmp_path / "tiny.fa"
        fasta.write_text(">chr1\nACGT\n")
        return fasta

    def fake_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        seen["prepared"] = fasta_path
        return _derive(Path(fasta_path))

    monkeypatch.setattr(UCSCGenomeDownloader, "_materialize_fasta", fake_materialize)
    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare)

    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    files = dl.fetch_genome_from("/some/ref.fa", overwrite=True)

    assert seen["source"] == "/some/ref.fa"
    assert seen["overwrite"] is True
    assert seen["prepared"] == tmp_path / "tiny.fa"
    assert files.fasta == tmp_path / "tiny.fa"
    # UCSC is never contacted when seeding from a user-provided FASTA.
    assert head_recorder.calls == []


# --- what registering and verifying answer with -----------------------------

# Nothing below registers anything: these are the values registration *returns*, so they
# are built by hand from the fields a record carries.


def _completion(kind: str, name: str, **details: object) -> CompletionRecord:
    """A completion record with everything filled in, so ``as_json`` has every key."""
    return CompletionRecord(
        kind=kind,
        name=name,
        files={f"{name}.fa.fai": 21, f"{name}.fa": 12},
        source_url="https://example.org/x.gz",
        sha256="1a2b3c",
        tool_versions={"samtools": "1.21"},
        package_version="2026.8.0",
        completed_at="2026-08-12T09:00:00+00:00",
        details=dict(details),
    )


class TestTheJsonKeysAndTheirOrder:
    """``as_json`` — every ``--json`` surface here, pinned key for key and in order.

    ``--json`` is what a script parses, so a key renamed, dropped or reordered is a break
    whether or not anything in this suite notices. These assert the whole list rather than
    a key inside it, which is the only form that fails on an addition.
    """

    def test_a_registered_assembly_is_a_record_plus_what_a_record_does_not_hold(self) -> None:
        # The same shape a registered annotation serializes in, deliberately: a record
        # plus the two facts a record does not hold about itself. test_gtf pins that half.
        registered = RegisteredAssembly(
            assembly="hg38",
            directory=Path("/data/genome/hg38"),
            record=_completion("genome", "hg38"),
        )

        assert list(registered.as_json()) == [
            "kind",
            "name",
            "files",
            "source_url",
            "sha256",
            "tool_versions",
            "package_version",
            "completed_at",
            "details",
            "assembly",
            "directory",
        ]
        assert registered.as_json()["directory"] == "/data/genome/hg38"

    def test_a_verified_assembly_pins_its_keys_and_serializes_what_supplied_the_digest(
        self,
    ) -> None:
        checked = VerifiedAssembly(
            assembly="sacCer3",
            fasta=Path("/data/genome/sacCer3/sacCer3.fa"),
            sha256="6ff72f07",
            expected="6ff72f07",
            expected_from=EXPECTED_FROM_TABLE,
            components=None,
        )

        assert list(checked.as_json()) == [
            "assembly",
            "fasta",
            "sha256",
            "expected",
            "expected_from",
            "verified",
            "components",
        ]
        assert checked.as_json()["verified"] is True
        assert checked.as_json()["fasta"] == "/data/genome/sacCer3/sacCer3.fa"

        # `expected_from` is serialized as the constant the CLI keys on.
        def _from(expected_from: str | None) -> object:
            return VerifiedAssembly(
                assembly="sacCer3",
                fasta=Path("/tmp/x.fa"),
                sha256="6ff72f07",
                expected=None if expected_from is None else "6ff72f07",
                expected_from=expected_from,
                components=None,
            ).as_json()["expected_from"]

        assert _from(EXPECTED_FROM_TABLE) == "table"
        assert _from(EXPECTED_FROM_RECORD) == "record"
        assert _from(None) is None


def test_a_registered_assembly_is_carried_whole_and_not_copied_out() -> None:
    # The properties that exist so a surface never re-reads a directory.
    registered = RegisteredAssembly(
        assembly="hg38", directory=Path("/data/genome/hg38"), record=_completion("genome", "hg38")
    )

    assert registered.source_url == "https://example.org/x.gz"
    assert registered.sha256 == "1a2b3c"
    assert registered.file_names == ["hg38.fa", "hg38.fa.fai"]
    first = registered.file_names
    first.append("intruder")
    assert registered.file_names == ["hg38.fa", "hg38.fa.fai"]  # a fresh list each call
    assert registered.chimera is None  # no build merged this one


# ---------------------------------------------------------------------------------------
# The edge this module must not grow back
# ---------------------------------------------------------------------------------------


def test_downloading_an_assembly_imports_nothing_that_registers_an_annotation() -> None:
    # The assembly half of the guard test_gtf holds for the annotation half, and the
    # reason the leaf both used to reach through could retire: what registering an
    # assembly answers with lives here and what registering an annotation answers with
    # lives in `io.gtf`, and neither module imports the other. Were this edge to open,
    # the two would be one module by another route and the cycle `io.chimera` closes
    # through `io.gtf` would come back with it.
    forbidden = {"genome.io.gtf", "genome.genome"}

    assert _module_level_imports(download_mod) & forbidden == set()
