"""Tests for genome.io.registration — the seam, not the downloader.

What test_download already proves is that registering a *download* still works end to
end; nothing here repeats it. These are the assertions that only exist because the
registration steps are now reachable without one: an assembly finishes on disk with
nothing fetched and no metadata row in sight, and the downloader is still one of these
plus a fetch.
"""

from __future__ import annotations

from pathlib import Path

import pooch
import pytest

from genome.io.completion import UnfinishedRegistrationError, read_record
from genome.io.download import Downloader, UCSCGenomeDownloader
from genome.io.fasta import GenomeFiles
from genome.io.registration import (
    ANNOTATIONS_SUBDIR,
    INDEXES_SUBDIR,
    AssemblyDir,
    AssemblyRegistration,
    assembly_data_dir,
)


def _derive(fasta: Path) -> GenomeFiles:
    """Write what a real preparation derives from ``fasta`` and return the whole set."""
    files = GenomeFiles(
        fasta=fasta,
        fai=fasta.with_name(fasta.name + ".fai"),
        twobit=fasta.with_name(fasta.stem + ".2bit"),
        chrom_sizes=fasta.with_name(fasta.stem + ".chrom.sizes"),
    )
    for derived in (files.fai, files.twobit, files.chrom_sizes):
        derived.write_text(f"derived from {fasta.name}\n")
    return files


def test_a_registration_defaults_to_the_assemblys_own_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # Construction moved with the steps: where an assembly's files belong is answered
    # once, so a builder that never downloads cannot answer it differently.
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
    assert AssemblyRegistration("hg38").cache_dir == tmp_path / "genome" / "hg38"
    assert AssemblyRegistration("hg38", tmp_path / "elsewhere").cache_dir == tmp_path / "elsewhere"


def test_an_assembly_registers_with_nothing_fetched(tmp_path: Path) -> None:
    # The whole reason for the seam: stage in the working area, place the FASTA, write
    # the record — no URL, no pinned digest, no network. A source_url of None is a
    # legitimate finished registration, not a missing field.
    registration = AssemblyRegistration("tiny", tmp_path)
    registration._work_dir.mkdir(parents=True)
    (registration._work_dir / "built.fa").write_text(">chrI\nACGT\n")

    fasta = registration._place_fasta(registration._work_dir / "built.fa")
    registration._record_completion(_derive(fasta), source_url=None, sha256="abc")

    assert fasta == tmp_path / "tiny.fa"  # placed under the assembly's own name
    record = read_record(tmp_path)
    assert record is not None
    assert (record.kind, record.name, record.source_url, record.sha256) == (
        "genome",
        "tiny",
        None,
        "abc",
    )
    assert sorted(record.files) == ["tiny.2bit", "tiny.chrom.sizes", "tiny.fa", "tiny.fa.fai"]
    assert not registration._work_dir.exists()  # disposable once the record vouches


def test_the_record_is_what_answers_a_second_registration(tmp_path: Path) -> None:
    registration = AssemblyRegistration("tiny", tmp_path)
    (tmp_path / "tiny.fa").write_text(">chrI\nACGT\n")
    files = _derive(tmp_path / "tiny.fa")
    registration._record_completion(files, source_url=None, sha256=None)

    again = AssemblyRegistration("tiny", tmp_path)

    assert again._completed_genome(overwrite=False, repair=again._repair_command()) == files
    # ...unless the caller is repairing, which skips the question entirely.
    assert again._completed_genome(overwrite=True, repair=again._repair_command()) is None


def test_a_directory_with_no_record_raises_naming_the_repair(tmp_path: Path) -> None:
    (tmp_path / "tiny.fa").write_text(">chrI\nACGT\n")
    registration = AssemblyRegistration("tiny", tmp_path)

    with pytest.raises(UnfinishedRegistrationError, match="genome register tiny --force"):
        registration._completed_genome(overwrite=False, repair=registration._repair_command())


def test_another_contexts_subtree_does_not_make_an_assembly_look_broken(tmp_path: Path) -> None:
    # The layout constants travelled with the steps that read them: an annotation or an
    # index registered first is not an interrupted assembly.
    (tmp_path / ANNOTATIONS_SUBDIR).mkdir()
    (tmp_path / INDEXES_SUBDIR).mkdir()
    registration = AssemblyRegistration("tiny", tmp_path)

    assert registration._completed_genome(overwrite=False, repair="unused") is None


def test_the_downloader_is_a_registration_that_also_fetches(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("LIULAB_DATA", str(tmp_path))
    dl = UCSCGenomeDownloader("hg38")

    assert isinstance(dl, AssemblyRegistration)
    # And *only* a registration. It used to be a Downloader as well, which gave it a
    # second answer to "which directory?" that its constructor then had to overrule; the
    # assembly's own directory is the only one it ever wanted.
    assert not isinstance(dl, Downloader)
    assert dl.cache_dir == assembly_data_dir("hg38")
    assert dl.cache_dir != Path(pooch.os_cache("genome"))


def test_the_plain_downloader_still_caches_where_it_always_did(tmp_path: Path) -> None:
    # The sibling it stopped inheriting from is untouched: fetch_url bound to a cache
    # directory, defaulting to pooch's per-user one.
    assert Downloader().cache_dir == Path(pooch.os_cache("genome"))
    assert Downloader(cache_dir=tmp_path).cache_dir == tmp_path


def test_a_directory_answers_whether_it_is_finished_without_a_registration(
    tmp_path: Path,
) -> None:
    # What a verification asks: it holds a name and a directory and no registration, and
    # the question is the same one a registration asks itself.
    here = AssemblyDir.locate("tiny", tmp_path)
    assert here.completed_files(repair="unused") is None  # empty: a fresh registration

    (tmp_path / "tiny.fa").write_text(">chrI\nACGT\n")
    files = _derive(tmp_path / "tiny.fa")
    AssemblyRegistration("tiny", tmp_path)._record_completion(files, source_url=None, sha256=None)

    assert here.completed_files(repair="unused") == files
