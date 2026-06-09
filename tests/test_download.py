"""Tests for genome.io.download.

Network access is avoided by monkeypatching pooch.retrieve: we assert that our
wrappers wire the right arguments through, and that the UCSC subclass builds the
right URL and selects a decompression processor.
"""

from __future__ import annotations

from pathlib import Path

import pooch
import pytest

from genome.io.download import Downloader, UCSCGenomeDownloader


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
