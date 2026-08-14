"""Tests for scripts/require_tools.sh — the guard that fronts the unit lane.

The guard is driven as the lane drives it: a real ``bash`` on the real script, with
``PATH`` holding stub binaries written into ``tmp_path``. Nothing here needs samtools or
either UCSC binary installed, which is the point — the cases worth testing are the ones a
developed machine cannot produce, a binary that is absent, one that resolves but will not
execute, and one that runs and refuses to say what it is.

``PATH`` is the stub directory and nothing else, plus a ``cat`` linked into it, because
the guard writes its message with a heredoc and a message that cannot print is not a
message.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from genome.external import REQUIRED_TOOLS

from .conftest import StubBinary

_GUARD = Path(__file__).resolve().parents[1] / "scripts" / "require_tools.sh"
_BASH = shutil.which("bash") or "/bin/bash"

# What the real binaries do, which is the behaviour the stubs stand in for: samtools
# answers `--version` with exit 0; faToTwoBit and twoBitInfo reject every flag and exit
# 255, so a bare run is the only thing that distinguishes them from an absent one.
_ANSWERS = {
    "samtools": "echo 'samtools 1.23.1'\nexit 0",
    "faToTwoBit": "echo 'faToTwoBit - Convert DNA from fasta to 2bit format' >&2\nexit 255",
    "twoBitInfo": "echo 'twoBitInfo - get information about sequences in a .2bit file' >&2\nexit 255",
}


@pytest.fixture
def stub_bin(tmp_path: Path) -> Path:
    """Return an otherwise empty directory that will be the guard's whole ``PATH``."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "cat").symlink_to(shutil.which("cat") or "/bin/cat")
    return bin_dir


@pytest.fixture
def answering_bin(stub_bin: Path, stub_binary: StubBinary) -> Path:
    """Return that directory with all three binaries installed and answering."""
    for name, body in _ANSWERS.items():
        stub_binary(stub_bin, name, body)
    return stub_bin


def _run_guard(bin_dir: Path) -> subprocess.CompletedProcess[str]:
    """Run the guard with ``bin_dir`` as its entire ``PATH``."""
    return subprocess.run(
        [_BASH, str(_GUARD)],
        env={"PATH": str(bin_dir)},
        capture_output=True,
        text=True,
        check=False,
    )


def test_the_guard_names_every_missing_binary_at_once(stub_bin: Path) -> None:
    result = _run_guard(stub_bin)

    assert result.returncode != 0
    for tool in REQUIRED_TOOLS:
        assert tool in result.stderr


def test_the_message_names_what_to_do_next(stub_bin: Path) -> None:
    result = _run_guard(stub_bin)

    assert "pixi install" in result.stderr


def test_the_guard_passes_when_every_binary_answers(answering_bin: Path) -> None:
    # The case that makes the probe what it is: two of the three exit 255 at everything
    # they are offered, so a guard testing for exit 0 would fail an installed machine.
    result = _run_guard(answering_bin)

    assert result.returncode == 0, result.stderr


def test_a_binary_that_resolves_but_cannot_execute_is_missing(
    answering_bin: Path, stub_binary: StubBinary
) -> None:
    # `which` would answer yes here. The guard asks the binary instead, and a file the
    # shell cannot execute answers 126 — which is what "answer, do not merely resolve"
    # means on a half-installed environment.
    stub_binary(answering_bin, "faToTwoBit", _ANSWERS["faToTwoBit"], executable=False)

    result = _run_guard(answering_bin)

    assert result.returncode != 0
    assert "faToTwoBit" in result.stderr
    assert "twoBitInfo" not in result.stderr


def test_a_samtools_that_will_not_report_its_version_is_missing(
    answering_bin: Path, stub_binary: StubBinary
) -> None:
    # A samtools whose shared libraries are gone runs and fails rather than resolving to
    # nothing, and it cannot index a FASTA either way.
    stub_binary(
        answering_bin, "samtools", "echo 'error while loading shared libraries' >&2\nexit 1"
    )

    result = _run_guard(answering_bin)

    assert result.returncode != 0
    assert "samtools" in result.stderr


def test_the_guard_probes_exactly_the_tools_the_package_requires() -> None:
    # The list is REQUIRED_TOOLS or it is a second copy of it, drifting quietly. Only the
    # code is read: the comment above it names STAR, chromap and bedtools to say why each
    # is deliberately not probed.
    text = _GUARD.read_text()
    body = text.split("set -uo pipefail", 1)[1].split('if [ -n "$missing" ]', 1)[0]
    code = "\n".join(line for line in body.splitlines() if not line.lstrip().startswith("#"))

    for tool in REQUIRED_TOOLS:
        assert tool in code
    for optional in ("STAR", "chromap", "bedtools"):
        assert optional not in code
