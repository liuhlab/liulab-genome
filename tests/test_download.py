"""Tests for genome.io.download.

Every download in the package goes through ``download.fetch_url``, so the suite stays
offline by replacing that one function with the shared ``fake_fetch`` fixture (see
tests/conftest.py) and asserting the arguments each caller wires through. ``fetch_url``
itself is exercised for real against an already-present file, which pooch serves without
touching the network. Nothing here monkeypatches pooch's own retrieve function.

Metadata is injected as an in-memory :class:`AssemblyMetadata` record, so no test here
reads or fakes the shipped TSV; the table itself is tested in test_metadata.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import pooch
import pytest
import requests

from genome.io import download as download_mod
from genome.io.download import (
    Downloader,
    UCSCGenomeDownloader,
    assembly_data_dir,
    assembly_table_row,
    fetch_url,
    liulab_data_dir,
)
from genome.io.fasta import GenomeFiles
from genome.io.utils import ChecksumMismatchError, sha256_file
from genome.metadata import AssemblyMetadata

from .conftest import FakeFetch

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


def _sha256(path: Path) -> str:
    """Return the sha256 of ``path`` in the ``algorithm:hexdigest`` form pooch accepts."""
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


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


@pytest.fixture
def no_native_prepare(monkeypatch: pytest.MonkeyPatch) -> None:
    """Stub the FASTA preparation step so registration runs without native binaries."""

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)


def test_liulab_data_dir_from_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path / "lab"))
    assert liulab_data_dir() == tmp_path / "lab"


def test_liulab_data_dir_defaults_to_home(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("LIULAB_DATA", raising=False)
    assert liulab_data_dir() == Path.home() / "liulab_data"


def test_liulab_data_dir_empty_env_falls_back(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("LIULAB_DATA", "")
    assert liulab_data_dir() == Path.home() / "liulab_data"


def test_assembly_data_dir_layout(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
    assert assembly_data_dir("hg38") == tmp_path / "genome" / "hg38"


def test_ucsc_default_cache_dir_is_assembly_data_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
    dl = UCSCGenomeDownloader("mm39")
    assert dl.cache_dir == tmp_path / "genome" / "mm39"


def test_default_cache_dir_is_under_pooch_os_cache() -> None:
    dl = Downloader()
    assert dl.cache_dir == Path(pooch.os_cache("genome"))


def test_explicit_cache_dir_is_used(tmp_path: Path) -> None:
    dl = Downloader(cache_dir=tmp_path / "cache")
    assert dl.cache_dir == tmp_path / "cache"


# --- the one fetch step ------------------------------------------------------


def test_fetch_url_serves_a_matching_local_file_without_downloading(
    tmp_path: Path, data_dir: Path
) -> None:
    # A file already at the destination whose hash matches is handed back as-is:
    # no downloader is ever constructed, so this exercises fetch_url offline.
    dest = tmp_path / "tiny.fa"
    dest.write_bytes((data_dir / "tiny.fa").read_bytes())

    result = fetch_url(
        "https://example.org/tiny.fa",
        tmp_path,
        known_hash=_sha256(dest),
        fname="tiny.fa",
        progressbar=False,
    )

    assert result.resolve() == dest.resolve()
    assert result.read_bytes() == (data_dir / "tiny.fa").read_bytes()


def test_fetch_url_returns_the_processor_output(tmp_path: Path, data_dir: Path) -> None:
    dest = tmp_path / "sacCer3.fa.gz"
    dest.write_bytes((data_dir / "tiny.fa.gz").read_bytes())

    result = fetch_url(
        "https://example.org/sacCer3.fa.gz",
        tmp_path,
        known_hash=_sha256(dest),
        fname="sacCer3.fa.gz",
        processor=pooch.Decompress(method="gzip", name="sacCer3.fa"),
        progressbar=False,
    )

    assert result.resolve() == (tmp_path / "sacCer3.fa").resolve()
    assert result.read_text() == (data_dir / "tiny.fa").read_text()


def test_fetch_passes_arguments_to_the_fetch_step(fake_fetch: FakeFetch, tmp_path: Path) -> None:
    dl = Downloader(cache_dir=tmp_path)
    result = dl.fetch(
        "https://example.org/big.bed.gz",
        known_hash="md5:abc",
        fname="big.bed.gz",
    )

    assert result == tmp_path / "big.bed.gz"
    call = fake_fetch.last
    assert call.url == "https://example.org/big.bed.gz"
    assert call.known_hash == "md5:abc"
    assert call.fname == "big.bed.gz"
    assert call.dest_dir == tmp_path
    assert call.progressbar is True


def test_ucsc_fasta_url() -> None:
    assert (
        UCSCGenomeDownloader("hg38").fasta_url
        == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    )
    assert UCSCGenomeDownloader("mm39").fasta_url.endswith("mm39/bigZips/mm39.fa.gz")


def test_ucsc_assembly_url() -> None:
    assert (
        UCSCGenomeDownloader("hg38").assembly_url
        == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/"
    )


def test_validate_assembly_ok_hits_directory_url(head_recorder: _HeadRecorder) -> None:
    dl = UCSCGenomeDownloader("hg38")
    dl.validate_assembly()  # 200 by default → no raise
    assert head_recorder.calls[0]["url"] == dl.assembly_url


def test_validate_assembly_404_raises_value_error(head_recorder: _HeadRecorder) -> None:
    head_recorder.status_code = 404
    dl = UCSCGenomeDownloader("nope99")
    with pytest.raises(ValueError, match="Unknown UCSC assembly 'nope99'"):
        dl.validate_assembly()


def test_validate_assembly_other_status_raises_http_error(head_recorder: _HeadRecorder) -> None:
    head_recorder.status_code = 500
    with pytest.raises(requests.exceptions.HTTPError):
        UCSCGenomeDownloader("hg38").validate_assembly()


def test_fetch_fasta_validates_a_derived_url(
    fake_fetch: FakeFetch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    # "tiny" is in no table, so the URL is derived and the name is checked against UCSC.
    fake_fetch.serve("tiny.fa.gz")
    UCSCGenomeDownloader("tiny", cache_dir=tmp_path).fetch_fasta()
    assert len(head_recorder.calls) == 1


def test_fetch_fasta_aborts_before_download_on_bad_assembly(
    fake_fetch: FakeFetch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    head_recorder.status_code = 404

    with pytest.raises(ValueError, match="Unknown UCSC assembly"):
        UCSCGenomeDownloader("bad", cache_dir=tmp_path).fetch_fasta()

    assert fake_fetch.calls == []  # validation fails before anything is fetched


def test_fetch_fasta_decompresses_by_default(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)

    result = dl.fetch_fasta()

    assert fake_fetch.last.url == dl.fasta_url
    # the compressed download is unpacked to <assembly>.fa beside it
    assert result == tmp_path / "hg38.fa"
    assert result.read_text() == (data_dir / "tiny.fa").read_text()


def test_fetch_fasta_without_decompress_keeps_the_archive(
    fake_fetch: FakeFetch, tmp_path: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)

    result = dl.fetch_fasta(decompress=False)

    assert fake_fetch.last.processor is None
    assert result.name.endswith(".fa.gz")


def test_fetch_genome_runs_full_pipeline(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # fetch_genome chains fetch_fasta (network) and prepare_fasta (native tools);
    # the fetch step is replaced by the fixture, prepare_fasta is stubbed here, so the
    # test stays offline and binary-free and asserts the wiring.
    fake_fetch.serve("tiny.fa.gz")
    prepared: dict[str, object] = {}

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        prepared["fasta"] = fasta_path
        prepared["overwrite"] = overwrite
        fasta = Path(fasta_path)
        return GenomeFiles(
            fasta=fasta,
            fai=fasta.with_name(fasta.name + ".fai"),
            twobit=fasta.with_name("hg38.2bit"),
            chrom_sizes=fasta.with_name("hg38.chrom.sizes"),
        )

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    files = dl.fetch_genome()

    assert fake_fetch.last.url == dl.fasta_url
    # fetch_fasta's decompressed output is handed to prepare_fasta...
    assert prepared["fasta"] == tmp_path / "hg38.fa"
    assert prepared["overwrite"] is False  # default: caches reused
    # ...and every derived path is surfaced on the returned record.
    assert files.fasta == tmp_path / "hg38.fa"
    assert files.fai == tmp_path / "hg38.fa.fai"
    assert files.twobit == tmp_path / "hg38.2bit"
    assert files.chrom_sizes == tmp_path / "hg38.chrom.sizes"


def test_fetch_genome_forwards_overwrite(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    prepared: dict[str, object] = {}

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        prepared["overwrite"] = overwrite
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    dl.fetch_genome(overwrite=True)

    assert prepared["overwrite"] is True


def test_fetch_genome_forwards_known_hash_and_decompresses(
    fake_fetch: FakeFetch, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        assert overwrite is False  # default: no forced regeneration
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    dl.fetch_genome(known_hash="md5:abc")

    call = fake_fetch.last
    assert call.url == dl.fasta_url
    assert call.known_hash == "md5:abc"
    # the pipeline always decompresses, so a Decompress processor is selected.
    assert isinstance(call.processor, pooch.Decompress)


# --- a row's pinned source and checksum -------------------------------------


def test_a_row_that_pins_a_source_is_the_url_fetched(fake_fetch: FakeFetch, tmp_path: Path) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    dl.fetch_fasta()

    assert dl.fasta_url == _PINNED_URL
    assert fake_fetch.last.url == _PINNED_URL


def test_a_record_passed_in_beats_the_shipped_table(fake_fetch: FakeFetch, tmp_path: Path) -> None:
    # hg38 has a row of its own; an explicit record replaces it wholesale.
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    dl.fetch_fasta()

    assert fake_fetch.last.url == _PINNED_URL


def test_a_pinned_source_skips_the_ucsc_name_check(
    fake_fetch: FakeFetch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    # Validation is a property of the source, and a pinned URL *is* the source, so
    # there is nothing left to guess about the name (ADR-0003).
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    dl.fetch_fasta()

    assert head_recorder.calls == []


def test_an_assembly_with_no_row_still_uses_the_golden_path(
    fake_fetch: FakeFetch, tmp_path: Path, head_recorder: _HeadRecorder, no_native_prepare: None
) -> None:
    # The table is a cross-reference, not an allow-list: no row takes nothing away.
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    assert dl.metadata is None

    files = dl.fetch_genome()

    assert dl.fasta_url == "https://hgdownload.soe.ucsc.edu/goldenPath/tiny/bigZips/tiny.fa.gz"
    assert fake_fetch.last.url == dl.fasta_url
    assert head_recorder.calls[0]["url"] == dl.assembly_url  # still validated
    assert files.fasta == tmp_path / "tiny.fa"


def test_registering_accepts_a_fasta_matching_the_pinned_checksum(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader(
        "tiny",
        cache_dir=tmp_path,
        metadata=_row(source_url=_PINNED_URL, sha256=_TINY_FA_SHA256),
    )

    files = dl.fetch_genome()

    assert files.fasta == tmp_path / "tiny.fa"
    assert files.fasta.is_file()


def test_registering_rejects_a_fasta_that_is_not_the_pinned_one(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None
) -> None:
    wrong = "0" * 64
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader(
        "tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL, sha256=wrong)
    )

    with pytest.raises(ChecksumMismatchError) as excinfo:
        dl.fetch_genome()

    message = str(excinfo.value)
    assert wrong in message  # what the row expected...
    assert _TINY_FA_SHA256 in message  # ...and what actually arrived


def test_the_archives_own_digest_is_not_what_is_checked(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path, no_native_prepare: None
) -> None:
    # The whole point of hashing unpacked content: pinning the .fa.gz's digest — which
    # is what pooch's known_hash would check — fails against the FASTA inside it.
    archive_digest = sha256_file(data_dir / "tiny.fa.gz")
    assert archive_digest != _TINY_FA_SHA256
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader(
        "tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL, sha256=archive_digest)
    )

    with pytest.raises(ChecksumMismatchError):
        dl.fetch_genome()


def test_a_blank_checksum_registers_and_reports_the_computed_value(
    fake_fetch: FakeFetch, tmp_path: Path, no_native_prepare: None
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))

    files = dl.fetch_genome()

    assert files.fasta.is_file()  # nothing to compare against, so nothing to fail
    assert dl.verify_fasta(files.fasta) == _TINY_FA_SHA256


def test_verify_fasta_defaults_to_the_assemblys_own_fasta(
    fake_fetch: FakeFetch, tmp_path: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path, metadata=_row(source_url=_PINNED_URL))
    dl.fetch_fasta()

    assert dl.verify_fasta() == _TINY_FA_SHA256


def test_table_row_fills_in_the_computed_checksum(fake_fetch: FakeFetch, tmp_path: Path) -> None:
    fake_fetch.serve("tiny.fa.gz")

    row = assembly_table_row("hg38", cache_dir=tmp_path, progressbar=False)

    assert row["sha256"] == _TINY_FA_SHA256
    assert row["source_url"] == "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    assert row["ncbi_name"] == "GRCh38"  # the curated identifiers are carried through


def test_table_row_leaves_a_new_assemblys_identifiers_blank(
    fake_fetch: FakeFetch, tmp_path: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")

    row = assembly_table_row("newAsm", cache_dir=tmp_path, progressbar=False)

    assert row["assembly_name"] == "newAsm"
    assert row["species"] is None  # only a person can supply this
    assert str(row["source_url"]).endswith("/newAsm/bigZips/newAsm.fa.gz")
    assert row["sha256"] == _TINY_FA_SHA256


# --- seeding from a user-provided FASTA (path_or_url) -----------------------


@pytest.mark.parametrize(
    ("source", "expected"),
    [
        ("https://x/y.fa.gz", True),
        ("http://x/y.fa", True),
        ("ftp://x/y.fa", True),
        ("sftp://x/y.fa", True),
        ("ftps://x/y.fa", False),  # pooch ships no ftps downloader
        ("/data/ref.fa", False),
        ("ref.fa.gz", False),
        ("~/ref.fa", False),
        ("relative/path.fa", False),
    ],
)
def test_looks_like_url(source: str, expected: bool) -> None:
    assert download_mod._looks_like_url(source) is expected


def test_materialize_fasta_copies_local_plain(tmp_path: Path, data_dir: Path) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "cache")

    out = dl._materialize_fasta(data_dir / "tiny.fa")

    assert out == tmp_path / "cache" / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()


def test_materialize_fasta_decompresses_local_gz(tmp_path: Path, data_dir: Path) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "cache")

    out = dl._materialize_fasta(data_dir / "tiny.fa.gz")

    assert out == tmp_path / "cache" / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()
    # the compressed download is kept alongside the decompressed FASTA
    assert (tmp_path / "cache" / "tiny.fa.gz").is_file()


def test_materialize_fasta_missing_local_raises(tmp_path: Path) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    with pytest.raises(FileNotFoundError, match="local FASTA source not found"):
        dl._materialize_fasta(tmp_path / "nope.fa")


def test_materialize_fasta_unfetchable_scheme_raises(tmp_path: Path) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    with pytest.raises(ValueError, match="no downloader for the 'ftps' scheme"):
        dl._materialize_fasta("ftps://example.org/tiny.fa")


def test_materialize_fasta_reuses_existing_unless_overwrite(tmp_path: Path) -> None:
    src = tmp_path / "src.fa"
    src.write_text("ONE")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path / "cache")

    dl._materialize_fasta(src)
    src.write_text("TWO")

    # a fresh <assembly>.fa is reused without re-reading the source...
    assert dl._materialize_fasta(src).read_text() == "ONE"
    # ...unless overwrite forces a refresh.
    assert dl._materialize_fasta(src, overwrite=True).read_text() == "TWO"


def test_materialize_fasta_url_goes_through_the_fetch_step(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)

    out = dl._materialize_fasta("https://example.org/whatever.fa", progressbar=False)

    assert out == tmp_path / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()
    call = fake_fetch.last
    assert call.url == "https://example.org/whatever.fa"
    assert call.dest_dir == tmp_path
    assert call.fname == "tiny.fa"  # the download lands as <assembly>.fa
    assert call.progressbar is False


def test_materialize_fasta_url_gz_is_decompressed(
    fake_fetch: FakeFetch, tmp_path: Path, data_dir: Path
) -> None:
    fake_fetch.serve("tiny.fa.gz")
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)

    out = dl._materialize_fasta("https://example.org/whatever.fa.gz")

    assert fake_fetch.last.fname == "tiny.fa.gz"
    assert out == tmp_path / "tiny.fa"
    assert out.read_text() == (data_dir / "tiny.fa").read_text()


def test_materialize_fasta_url_overwrite_refetches(fake_fetch: FakeFetch, tmp_path: Path) -> None:
    dl = UCSCGenomeDownloader("tiny", cache_dir=tmp_path)
    url = "https://example.org/whatever.fa"

    dl._materialize_fasta(url)
    dl._materialize_fasta(url, overwrite=True)

    # overwrite discards the kept download, so the source is fetched again rather
    # than served from disk.
    assert len(fake_fetch.calls) == 2


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
        return tmp_path / "tiny.fa"

    def fake_prepare(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        seen["prepared"] = fasta_path
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

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
