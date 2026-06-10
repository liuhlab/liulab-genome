"""Tests for genome.io.download.

Network access is avoided by monkeypatching pooch.retrieve: we assert that our
wrappers wire the right arguments through, and that the UCSC subclass builds the
right URL and selects a decompression processor.
"""

from __future__ import annotations

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
    liulab_data_dir,
)
from genome.io.fasta import GenomeFiles


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


def test_fetch_passes_arguments_to_pooch(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    def fake_retrieve(**kwargs: object) -> str:
        captured.update(kwargs)
        return str(tmp_path / "downloaded.bin")

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)

    dl = Downloader(cache_dir=tmp_path)
    result = dl.fetch(
        "https://example.org/big.bed.gz",
        known_hash="md5:abc",
        fname="big.bed.gz",
    )

    assert result == tmp_path / "downloaded.bin"
    assert captured["url"] == "https://example.org/big.bed.gz"
    assert captured["known_hash"] == "md5:abc"
    assert captured["fname"] == "big.bed.gz"
    assert captured["path"] == tmp_path
    assert captured["progressbar"] is True


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


def test_fetch_fasta_validates_by_default(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    monkeypatch.setattr(pooch, "retrieve", lambda **_kwargs: str(tmp_path / "hg38.fa"))
    UCSCGenomeDownloader("hg38", cache_dir=tmp_path).fetch_fasta()
    assert len(head_recorder.calls) == 1


def test_fetch_fasta_validate_false_skips_check(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    monkeypatch.setattr(pooch, "retrieve", lambda **kwargs: str(tmp_path / "hg38.fa"))
    UCSCGenomeDownloader("hg38", cache_dir=tmp_path).fetch_fasta(validate=False)
    assert head_recorder.calls == []


def test_fetch_fasta_aborts_before_download_on_bad_assembly(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, head_recorder: _HeadRecorder
) -> None:
    head_recorder.status_code = 404

    def fail_retrieve(**kwargs: object) -> str:
        raise AssertionError("download must not start when validation fails")

    monkeypatch.setattr(pooch, "retrieve", fail_retrieve)
    with pytest.raises(ValueError, match="Unknown UCSC assembly"):
        UCSCGenomeDownloader("bad", cache_dir=tmp_path).fetch_fasta()


def test_fetch_fasta_decompresses_by_default(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: dict[str, object] = {}

    def fake_retrieve(**kwargs: object) -> str:
        captured.update(kwargs)
        return str(tmp_path / "hg38.fa")

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    result = dl.fetch_fasta()

    assert result == tmp_path / "hg38.fa"
    assert captured["url"] == dl.fasta_url
    processor = captured["processor"]
    assert isinstance(processor, pooch.Decompress)
    assert processor.name == "hg38.fa"


def test_fetch_fasta_without_decompress_has_no_processor(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: dict[str, object] = {}

    def fake_retrieve(**kwargs: object) -> str:
        captured.update(kwargs)
        return str(tmp_path / "hg38.fa.gz")

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    dl.fetch_fasta(decompress=False)

    assert captured["processor"] is None


def test_fetch_genome_runs_full_pipeline(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    # fetch_genome chains fetch_fasta (network) and prepare_fasta (native tools);
    # stub both so the test stays offline and binary-free, then assert the wiring.
    def fake_retrieve(**kwargs: object) -> str:
        assert kwargs["url"] == dl.fasta_url
        return str(tmp_path / "hg38.fa")

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

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)
    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    files = dl.fetch_genome()

    # fetch_fasta's decompressed output is handed to prepare_fasta...
    assert prepared["fasta"] == tmp_path / "hg38.fa"
    assert prepared["overwrite"] is False  # default: caches reused
    # ...and every derived path is surfaced on the returned record.
    assert files.fasta == tmp_path / "hg38.fa"
    assert files.fai == tmp_path / "hg38.fa.fai"
    assert files.twobit == tmp_path / "hg38.2bit"
    assert files.chrom_sizes == tmp_path / "hg38.chrom.sizes"


def test_fetch_genome_forwards_overwrite(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    prepared: dict[str, object] = {}

    def fake_retrieve(**kwargs: object) -> str:
        assert kwargs["url"] == dl.fasta_url
        return str(tmp_path / "hg38.fa")

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        prepared["overwrite"] = overwrite
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)
    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    dl.fetch_genome(overwrite=True)

    assert prepared["overwrite"] is True


def test_fetch_genome_forwards_known_hash_and_decompresses(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    captured: dict[str, object] = {}

    def fake_retrieve(**kwargs: object) -> str:
        captured.update(kwargs)
        return str(tmp_path / "hg38.fa")

    def fake_prepare_fasta(fasta_path: Path, *, overwrite: bool = False) -> GenomeFiles:
        assert overwrite is False  # default: no forced regeneration
        fasta = Path(fasta_path)
        return GenomeFiles(fasta=fasta, fai=fasta, twobit=fasta, chrom_sizes=fasta)

    monkeypatch.setattr(pooch, "retrieve", fake_retrieve)
    monkeypatch.setattr(download_mod, "prepare_fasta", fake_prepare_fasta)

    dl = UCSCGenomeDownloader("hg38", cache_dir=tmp_path)
    dl.fetch_genome(known_hash="md5:abc")

    assert captured["url"] == dl.fasta_url
    assert captured["known_hash"] == "md5:abc"
    # the pipeline always decompresses, so a Decompress processor is selected.
    assert isinstance(captured["processor"], pooch.Decompress)
