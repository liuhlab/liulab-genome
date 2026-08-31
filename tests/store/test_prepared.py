"""Tests for genome.store.prepared — the pipeline every **Prepared set** is prepared by.

**The point of this file is the fourth source.** Three contexts each own one of these —
a motif set, an xref set, a homology set — and each used to carry the pipeline written out
again. So a fictitious fourth one is prepared here end to end, and everything it declares
is a URL, a checksum, a directory, a reader and an error class of its own. If preparing it
had needed a fetch, a working area, a digest or a marker written here, the abstraction
would be the wrong one and this file would be the proof.

Offline throughout: the ``fake_fetch`` fixture stands in for the package's one fetch step,
and what it serves is written into the test's own directory rather than committed — the
publisher is invented, so there are no real bytes to subsample.

The unit lane, unmarked: nothing here needs a binary.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

import pytest

from genome.homology.compara import homology_data_dir
from genome.store import fetch as fetch_mod
from genome.store import prepared as prepared_mod
from genome.store.completion import (
    RECORD_NAME,
    RegistrationMismatchError,
    UnfinishedRegistrationError,
    read_record,
)
from genome.store.prepared import (
    ARCHIVE,
    UNPACKED,
    Checksum,
    PreparedChecksumError,
    PreparedSetNotDownloadedError,
    PreparedSource,
    login_node_help,
    prepare,
    unpacked_digest,
    unpacked_lines,
    write_through,
)
from genome.tf.motif.jaspar import motif_data_dir
from genome.xref.xref import xref_data_dir

from .._sources import PACKAGE, imports_any, sources
from ..conftest import FakeFetch

# ---------------------------------------------------------------------------
# A fictitious fourth Prepared set: a URL, a checksum and a reader, and no more
# ---------------------------------------------------------------------------

#: What the invented publisher publishes: one row per sighting, several species in one
#: file, the shape every real source here has.
PUBLISHED = (
    "species\tstem\tsightings\n"
    "beast\tBEAST00000001\t3\n"
    "quokka\tQUOK00000001\t9\n"
    "beast\tBEAST00000002\t1\n"
)

#: The set's own columns, which are the publisher's minus the one it was sliced on.
SLICE_HEADER = "stem\tsightings\n"

#: Where the invented publisher publishes it, and what a caller runs to prepare it.
URL = "https://example.invalid/beasts/2026/sightings.tsv.gz"
PREPARE_COMMAND = "python -c \"from beasts import Sightings; Sightings('beast')\""


class SightingsNotDownloadedError(PreparedSetNotDownloadedError):
    """This fourth topic's own not-downloaded error, which is all it takes to have one."""


def read_sightings(lines: Iterator[str], staged: Path, *, origin: str) -> Mapping[str, Any]:
    """Slice the publisher's rows to one species — the whole of what a fourth source adds."""
    kept: list[str] = []
    header = next(lines, "")
    if header != "species\tstem\tsightings\n":
        raise ValueError(f"{origin} is not the sightings table")
    for line in lines:
        species, stem, sightings = line.rstrip("\n").split("\t")
        if species == "beast":
            kept.append(f"{stem}\t{sightings}\n")
    with gzip.open(staged, "wt", encoding="utf-8") as out:
        out.write(SLICE_HEADER)
        out.writelines(kept)
    return {"rows": len(kept)}


def source(directory: Path, *, checksum: Checksum | None) -> PreparedSource:
    """Declare the fourth set. Everything not here is the shared pipeline's."""
    return PreparedSource(
        url=URL,
        directory=directory,
        stored_name="beast.sightings.tsv.gz",
        kind="sightings",
        name="beasts/2026",
        prepare_command=PREPARE_COMMAND,
        description="the beasts 2026 sightings set",
        read=read_sightings,
        not_downloaded=SightingsNotDownloadedError,
        checksum=checksum,
        details={"publisher": "The Beast Survey", "release": "2026"},
    )


@pytest.fixture
def published(tmp_path: Path) -> Path:
    """The publisher's file, gzipped as it is served."""
    path = tmp_path / "published" / "sightings.tsv.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as out:
        out.write(PUBLISHED)
    return path


@pytest.fixture
def pinned() -> Checksum:
    """The pin the publisher would state: sha256 of its **unpacked** bytes (ADR-0006)."""
    return Checksum("sha256", hashlib.sha256(PUBLISHED.encode()).hexdigest())


@pytest.fixture
def serving(fake_fetch: FakeFetch, published: Path) -> FakeFetch:
    """A fetch step serving the invented publisher's file."""
    fake_fetch.serve(published)
    return fake_fetch


class TestAFourthSource:
    def test_a_url_a_checksum_and_a_reader_prepare_a_set_end_to_end(
        self, serving: FakeFetch, pinned: Checksum, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        result = prepare(source(directory, checksum=pinned), progressbar=False)

        assert serving.last.url == URL
        # The pin covers the unpacked bytes, so pooch is handed none: it would compare it
        # against the compressed file and reject every download.
        assert serving.last.known_hash is None
        assert result.path == directory / "beast.sightings.tsv.gz"
        with gzip.open(result.path, "rt", encoding="utf-8") as handle:
            assert handle.read() == (f"{SLICE_HEADER}BEAST00000001\t3\nBEAST00000002\t1\n")

        # The marker is the pipeline's, and it carries what the source declared and what
        # the reader measured — no code here wrote one.
        record = read_record(directory)
        assert record is not None
        assert record == result.record
        assert (record.kind, record.name) == ("sightings", "beasts/2026")
        assert record.source_url == URL
        assert record.details["publisher"] == "The Beast Survey"
        assert record.details["rows"] == 2
        assert record.sha256 == unpacked_digest(result.path)
        assert record.files == {result.path.name: result.path.stat().st_size}
        # Working state, never a claimed output: it goes once the marker is written.
        assert not (directory / ".work").exists()

    def test_a_second_prepare_reads_the_marker_and_fetches_nothing(
        self, serving: FakeFetch, pinned: Checksum, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        first = prepare(source(directory, checksum=pinned), progressbar=False)
        again = prepare(source(directory, checksum=pinned), progressbar=False)

        assert len(serving.calls) == 1
        assert again.path == first.path
        assert again.record == first.record

    def test_a_file_that_does_not_match_its_pin_is_refused_and_names_both_halves(
        self, serving: FakeFetch, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        wrong = Checksum("sha256", "0" * 64)

        with pytest.raises(PreparedChecksumError) as raised:
            prepare(source(directory, checksum=wrong), progressbar=False)

        message = str(raised.value)
        assert "unpacked bytes" in message
        assert f"rm -rf {directory}" in message
        assert PREPARE_COMMAND in message
        # Nothing was placed and nothing was recorded: the set is absent, not short.
        assert list(directory.glob("*.tsv.gz")) == []
        assert read_record(directory) is None

    def test_a_set_that_cannot_be_fetched_raises_its_own_error_and_names_the_login_node(
        self, monkeypatch: pytest.MonkeyPatch, pinned: Checksum, tmp_path: Path
    ) -> None:
        def no_internet(url: str, dest_dir: Path, **kwargs: Any) -> Path:
            raise ConnectionError("the compute node has no internet")

        monkeypatch.setattr(fetch_mod, "fetch_url", no_internet)
        directory = tmp_path / "beasts" / "2026"

        with pytest.raises(SightingsNotDownloadedError) as raised:
            prepare(source(directory, checksum=pinned), progressbar=False)

        message = str(raised.value)
        assert "the beasts 2026 sightings set" in message  # this set, not "a set"
        assert "login node" in message
        assert PREPARE_COMMAND in message

    def test_an_interrupted_run_reads_as_unfinished_and_the_repair_names_the_rebuild(
        self, serving: FakeFetch, pinned: Checksum, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        prepare(source(directory, checksum=pinned), progressbar=False)
        (directory / RECORD_NAME).unlink()

        with pytest.raises(UnfinishedRegistrationError) as raised:
            prepare(source(directory, checksum=pinned), progressbar=False)

        assert f"rm -rf {directory}" in str(raised.value)
        assert PREPARE_COMMAND in str(raised.value)

    def test_a_marker_that_disagrees_with_the_file_is_a_mismatch(
        self, serving: FakeFetch, pinned: Checksum, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        result = prepare(source(directory, checksum=pinned), progressbar=False)
        result.path.write_bytes(b"")

        with pytest.raises(RegistrationMismatchError, match=RECORD_NAME):
            prepare(source(directory, checksum=pinned), progressbar=False)


# ---------------------------------------------------------------------------
# The steps the fourth source did not have to write
# ---------------------------------------------------------------------------


class TestTheWorkingArea:
    def test_what_nothing_can_vouch_for_is_swept_before_the_fetch(
        self, serving: FakeFetch, pinned: Checksum, tmp_path: Path
    ) -> None:
        # pooch serves a file already sitting at the destination, so a leftover download
        # nothing pins would be adopted as a finished one. With no archive pin, it goes.
        directory = tmp_path / "beasts" / "2026"
        work = directory / ".work"
        work.mkdir(parents=True)
        stale = work / "sightings.tsv.gz"
        stale.write_bytes(b"not a gzip")

        prepare(source(directory, checksum=pinned), progressbar=False)
        assert not work.exists()

    def test_an_archive_pin_is_pooch_s_to_check_and_its_download_is_kept(
        self, serving: FakeFetch, published: Path, tmp_path: Path
    ) -> None:
        # The other half of the same rule: a pin over the bytes as served is handed to
        # pooch, which re-checks a leftover download and fetches again when it does not
        # match — so the working area survives an interrupted run and a repeat costs no
        # second download.
        directory = tmp_path / "beasts" / "2026"
        archive = Checksum("md5", hashlib.md5(published.read_bytes()).hexdigest(), covers=ARCHIVE)
        work = directory / ".work"
        work.mkdir(parents=True)
        (work / "kept-by-an-interrupted-run").write_text("still here")

        prepare(source(directory, checksum=archive), progressbar=False)

        assert serving.last.known_hash == f"md5:{archive.digest}"
        # Nothing was verified a second time here: what the pin covers decides where it is
        # checked, and the record still carries the stored slice's own digest.
        record = read_record(directory)
        assert record is not None
        assert record.sha256 != archive.digest


class TestDecompressWhileHashing:
    def test_one_implementation_reads_packed_and_plain_alike_and_digests_the_content(
        self, tmp_path: Path
    ) -> None:
        plain = tmp_path / "sightings.tsv"
        plain.write_text(PUBLISHED, encoding="utf-8")
        packed = tmp_path / "sightings.tsv.gz"
        with gzip.open(packed, "wt", encoding="utf-8") as out:
            out.write(PUBLISHED)

        assert list(unpacked_lines(packed)) == PUBLISHED.splitlines(keepends=True)
        assert list(unpacked_lines(plain)) == list(unpacked_lines(packed))

        content = hashlib.sha256(PUBLISHED.encode()).hexdigest()
        assert unpacked_digest(plain) == content
        assert unpacked_digest(packed) == content
        # The archive digests to something else entirely, which is the whole of ADR-0006:
        # recompressing at another level breaks a match that ought to hold.
        assert hashlib.sha256(packed.read_bytes()).hexdigest() != content

    def test_the_digest_is_fed_the_bytes_as_they_are_read(self, tmp_path: Path) -> None:
        packed = tmp_path / "sightings.tsv.gz"
        with gzip.open(packed, "wt", encoding="utf-8") as out:
            out.write(PUBLISHED)

        digest = hashlib.md5()
        read = "".join(unpacked_lines(packed, digest))
        assert read == PUBLISHED
        assert digest.hexdigest() == hashlib.md5(PUBLISHED.encode()).hexdigest()

    def test_write_through_stores_what_it_was_handed_byte_for_byte(self, tmp_path: Path) -> None:
        staged = tmp_path / "through" / "sightings.tsv"
        measured = write_through(iter(PUBLISHED.splitlines(keepends=True)), staged, origin="x")

        assert staged.read_text(encoding="utf-8") == PUBLISHED
        assert measured == {"lines": len(PUBLISHED.splitlines())}


class TestWhatASourceDeclares:
    def test_a_checksum_carries_its_algorithm_and_what_it_covers(self) -> None:
        assert Checksum.parse("md5:abc") == Checksum("md5", "abc", covers=UNPACKED)
        assert str(Checksum.parse("md5:abc", covers=ARCHIVE)) == "md5:abc"

        with pytest.raises(ValueError, match="algorithm"):
            Checksum.parse("abc")
        with pytest.raises(ValueError, match="unpacked or archive"):
            Checksum.parse("md5:abc", covers="the middle")

    def test_the_repair_is_both_halves_and_the_login_node_sentence_is_one_sentence(
        self, tmp_path: Path
    ) -> None:
        declared = source(tmp_path / "a b" / "2026", checksum=None)
        # Quoted, because a Data dir with a space in it is still a Data dir.
        assert declared.repair == f"rm -rf '{tmp_path / 'a b' / '2026'}' && {PREPARE_COMMAND}"
        assert login_node_help(PREPARE_COMMAND).endswith(
            f"`{PREPARE_COMMAND}`, after which "
            "it is read from the Data dir and "
            "shared by every project on the "
            "machine."
        )

    def test_a_source_that_pins_nothing_records_the_digest_of_what_it_stored(
        self, serving: FakeFetch, tmp_path: Path
    ) -> None:
        directory = tmp_path / "beasts" / "2026"
        result = prepare(source(directory, checksum=None), progressbar=False)

        assert serving.last.known_hash is None
        assert result.record.sha256 == unpacked_digest(result.path)


class TestEveryRootIsDeclaredWhereItIsRead:
    def test_the_three_prepared_set_roots_sit_beside_the_assembly_tree(
        self, liulab_data: Path
    ) -> None:
        assert motif_data_dir() == liulab_data / "motif"
        assert xref_data_dir() == liulab_data / "xref"
        assert homology_data_dir() == liulab_data / "homology"
        # Siblings of the assembly tree, never tenants of it.
        assert not any(
            root.is_relative_to(liulab_data / "genome")
            for root in (motif_data_dir(), xref_data_dir(), homology_data_dir())
        )

    def test_the_pipeline_declares_none_of_them(self) -> None:
        # Each root is imported above from the context that reads what lives under it —
        # `motif/` from the JASPAR loader, `xref/` from the xref set, `homology/` from the
        # Compara set — and the pipeline they share knows none of the three names. That is
        # what keeps a fourth prepared set from having to be declared here: it brings its
        # own root, in its own package, and this module gains nothing.
        assert not [name for name in vars(prepared_mod) if name.endswith("_data_dir")]


def test_no_module_in_the_store_imports_a_context() -> None:
    # What makes this a place every context can reach: nothing here reaches one back. The
    # store is the **Data dir** root, the fetch step, the **Completion marker**, the digest
    # and this pipeline — bytes, directories and hashes, and no idea what an assembly or an
    # annotation is. Were one of these five to import a context, the other contexts would
    # be importing that one through it, and the package that exists to belong to nobody
    # would belong to whichever context got there first. Read off the source at any depth,
    # since a deferred import is still an import.
    forbidden = ("genome.assembly", "genome.annotation")
    offenders = [path.name for path in sources(PACKAGE / "store") if imports_any(path, forbidden)]

    assert offenders == []


def test_the_shared_pipeline_reaches_the_fetch_step_through_its_module() -> None:
    # The offline seam the whole suite rests on: one `monkeypatch.setattr` on
    # `genome.store.fetch` must take every download in the package offline, which a module
    # that imported the name would defeat by holding a reference the rebinding never
    # reaches. Asserted on the source, because a module that got this wrong would still
    # pass every other test in this file — `fake_fetch` patches the same object either way
    # only while the import is spelled this way.
    source_text = Path(prepared_mod.__file__).read_text(encoding="utf-8")
    assert "from genome.store import fetch" in source_text
    assert "fetch.fetch_url(" in source_text
    assert "from genome.store.fetch import" not in source_text


def test_the_marker_is_written_last_after_the_file_it_claims(
    serving: FakeFetch, pinned: Checksum, tmp_path: Path
) -> None:
    # Not an ordering detail: the marker is the only thing that says a set is finished, so
    # it must not exist while the file it claims does not.
    directory = tmp_path / "beasts" / "2026"
    result = prepare(source(directory, checksum=pinned), progressbar=False)

    claimed = json.loads((directory / RECORD_NAME).read_text(encoding="utf-8"))["files"]
    assert claimed == {result.path.name: result.path.stat().st_size}
    assert (directory / RECORD_NAME).stat().st_mtime >= result.path.stat().st_mtime
