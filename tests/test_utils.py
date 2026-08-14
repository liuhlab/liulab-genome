"""Tests for genome.io.utils — checksums and decompression, and nothing that runs a tool.

None of these need the native binaries, because this module no longer reaches one.
Running an **External tool** — resolution, the two run flavours, the failure message and
the freshness rule — is tested once in test_external against ``genome.external``, and
the preparation steps that drive one are tested in test_fasta.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from genome.io import utils
from genome.io.utils import (
    ChecksumMismatchError,
    _gunzip,
    sha256_file,
)


def test_sha256_file_digests_a_real_fixture(data_dir: Path) -> None:
    # The committed sacCer3 subsample, hashed independently with `shasum -a 256`.
    assert (
        sha256_file(data_dir / "tiny.fa")
        == "9316629bab14f9298a043f8b92e1e04a573b12d6a367ccc07c8f8040e5a13981"
    )


def test_sha256_of_a_gzipped_file_differs_from_its_contents(data_dir: Path) -> None:
    # Why the metadata table records the unpacked digest: the archive's is a different
    # number entirely, and changes whenever the file is recompressed.
    assert sha256_file(data_dir / "tiny.fa.gz") != sha256_file(data_dir / "tiny.fa")


def test_sha256_file_spans_many_read_buffers(tmp_path: Path) -> None:
    # Several megabytes, so the digest is assembled from many chunks rather than one.
    payload = b"ACGTN" * 1_000_000
    path = tmp_path / "big.fa"
    path.write_bytes(payload)

    assert sha256_file(path) == hashlib.sha256(payload).hexdigest()


def test_sha256_file_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        sha256_file(tmp_path / "nope.fa")


@pytest.fixture(scope="session")
def digest_scratch(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A reusable file path for the property test (session-scoped: hypothesis-safe)."""
    return tmp_path_factory.mktemp("digest") / "payload.bin"


@given(payload=st.binary(max_size=1 << 16))
def test_chunking_is_invisible(digest_scratch: Path, payload: bytes) -> None:
    # Streaming the file must give the digest of its bytes, whatever they are.
    digest_scratch.write_bytes(payload)
    assert sha256_file(digest_scratch) == hashlib.sha256(payload).hexdigest()


def test_checksum_mismatch_names_the_file_and_both_values() -> None:
    err = ChecksumMismatchError(Path("/data/hg38.fa"), "1a2b3c", "9f8e7d")

    assert isinstance(err, ValueError)
    message = str(err)
    assert "/data/hg38.fa" in message
    assert "1a2b3c" in message  # expected
    assert "9f8e7d" in message  # actual
    assert (err.expected, err.actual) == ("1a2b3c", "9f8e7d")


def test_gunzip_round_trips(tmp_path: Path) -> None:
    payload = b"".join(b">seq%d\nACGTACGTNN\n" % i for i in range(1000))
    src = tmp_path / "data.gz"
    with gzip.open(src, "wb") as fh:
        fh.write(payload)

    dest = tmp_path / "data"
    result = _gunzip(src, dest)

    assert result == dest
    assert dest.read_bytes() == payload


def test_this_layer_carries_no_run_of_its_own() -> None:
    # The fold, asserted. `_run` and `_run_to` were two lines each over
    # `ExternalTool.run`/`run_to` — a name-addressed restatement that existed so a test
    # could patch a module global, and that kept a second copy of the freshness branch
    # in a module whose subject is checksums. `io/fasta` holds the tool now, so the copy
    # is gone; this keeps the layer from growing back.
    assert not hasattr(utils, "_run")
    assert not hasattr(utils, "_run_to")
