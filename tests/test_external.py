"""Tests for genome.external — the one module that shells out to a native binary.

Almost nothing here needs a tool installed: a stub on ``PATH`` is a real binary as far as
this module is concerned, so resolution, the two run flavours, the failure message and the
freshness cache are all driven against one — see the ``stub_binary`` fixture for what one
is and why it is built the way it is. The handful that do need the pixi env's own binaries
say so and skip cleanly outside it.

``RecordingTool`` is exercised as the adapter it is — the same base class, only the
execution replaced — because every aligner test in the suite is standing on it.
"""

from __future__ import annotations

import shutil
import sys
from collections.abc import Callable
from pathlib import Path

import pytest

import genome.external as external
from genome.external import (
    NO_VERSION_REPORTED,
    REQUIRED_TOOLS,
    ExternalTool,
    InstalledTool,
    RecordingTool,
    ToolNotFoundError,
    clear_version_cache,
    doctor,
    is_fresh,
)
from genome.io.fasta import PREPARATION_TOOLS

from .conftest import StubBinary

_BINARIES_PRESENT = all(shutil.which(t) is not None for t in REQUIRED_TOOLS)


@pytest.fixture
def on_path(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stub_binary: StubBinary
) -> Callable[[str, str], Path]:
    """Return a helper that puts a stub on ``PATH`` under a chosen tool name.

    ``PATH`` becomes exactly this directory, and the interpreter is pointed somewhere
    empty, so the tools the pixi env really has cannot answer for a test that never
    installed them — the resolver's second lookup would otherwise find them.
    """
    bin_dir = tmp_path / "bin"
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "nowhere" / "python"))

    def install(name: str, body: str) -> Path:
        return stub_binary(bin_dir, name, body)

    return install


# -- locating a binary ------------------------------------------------------


def test_a_tool_on_path_resolves_to_it_falls_back_to_the_interpreters_bin_and_is_cached(
    on_path: Callable[[str, str], Path],
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    stub_binary: StubBinary,
) -> None:
    written = on_path("faToTwoBit", "exit 0")
    tool = InstalledTool("faToTwoBit")

    assert tool.path == str(written)
    written.unlink()  # gone from disk; the answer was remembered, not looked up again

    assert tool.path == str(written)

    # Running the env's interpreter without PATH activated: the normal lookup misses,
    # but the tool sits beside sys.executable (the conda/pixi bin/), so resolution falls
    # back to that directory.
    bin_dir = tmp_path / "bin"
    fallback = stub_binary(bin_dir, "twoBitInfo", "exit 0")
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(bin_dir / "python"))

    assert InstalledTool("twoBitInfo").path == str(fallback)


def test_a_tool_that_is_nowhere_raises_naming_the_command_that_installs_it(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "bin" / "python"))

    with pytest.raises(ToolNotFoundError) as raised:
        _ = InstalledTool("faToTwoBit").path

    message = str(raised.value)
    # The bit that was wrong before this module knew: the conda package is not the
    # binary's own name, so `pixi add faToTwoBit` would have failed for whoever ran it.
    assert "pixi add ucsc-fatotwobit" in message
    assert "pixi shell" in message

    # The same instructions, asked for directly rather than raised: a tool nobody
    # catalogued still names a plausible command, and a catalogued one quotes its own
    # homepage too.
    assert "pixi add bowtie2" in RecordingTool("bowtie2").install_instructions()
    assert "https://github.com/alexdobin/STAR" in RecordingTool("STAR").install_instructions()


# -- asking a tool what it is -----------------------------------------------


def test_version_is_the_first_line_of_stdout_or_falls_back_to_stderr(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo 'samtools 1.21'; echo 'Using htslib 1.21'")
    assert InstalledTool("samtools").version == "samtools 1.21"

    # Different tools choose differently; chromap answers on stderr.
    on_path("chromap", "echo '0.3.2-r518' >&2")
    assert InstalledTool("chromap").version == "0.3.2-r518"


def test_a_tool_that_declines_to_identify_itself_reports_no_version_but_absence_raises(
    on_path: Callable[[str, str], Path], monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # What several UCSC binaries do: `faToTwoBit --version` exits non-zero saying the
    # option is not valid. That is a tool declining to identify itself, not a failure,
    # and the "" it reports is remembered exactly as a real answer would be.
    tally = tmp_path / "asked"
    on_path(
        "faToTwoBit", f"echo x >> {tally}; echo '--version is not a valid option' >&2; exit 255"
    )

    assert InstalledTool("faToTwoBit").version == ""
    assert InstalledTool("faToTwoBit").version == ""
    assert tally.read_text() == "x\n"

    # The two answers must stay apart: "" means it ran and would not say, so a tool that
    # is not there at all still raises rather than being folded into that same "".
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "nowhere" / "python"))
    with pytest.raises(ToolNotFoundError):
        _ = InstalledTool("twoBitInfo").version


def test_the_version_is_asked_for_once_per_binary(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    # Recording a build's provenance must not cost a subprocess per mention of it — not
    # for the same object asked twice, and not for a second object naming the same path,
    # which is the waste a build re-asking per step would otherwise pay.
    tally = tmp_path / "asked"
    on_path("samtools", f"echo x >> {tally}; echo 'samtools 1.21'")
    tool = InstalledTool("samtools")

    assert tool.version == tool.version == "samtools 1.21"
    assert InstalledTool("samtools").version == "samtools 1.21"
    assert tally.read_text() == "x\n"


def test_two_binaries_of_the_same_name_each_answer_for_themselves(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stub_binary: StubBinary
) -> None:
    # Why the answer is remembered against the path and not against the name: a name is
    # only as stable as PATH, and a caller that has just been pointed at another samtools
    # must get that one's answer rather than the one asked earlier in this process.
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "nowhere" / "python"))
    old, new = tmp_path / "old", tmp_path / "new"
    stub_binary(old, "samtools", "echo 'samtools 1.20'")
    stub_binary(new, "samtools", "echo 'samtools 1.21'")

    monkeypatch.setenv("PATH", str(old))
    assert InstalledTool("samtools").version == "samtools 1.20"

    monkeypatch.setenv("PATH", str(new))
    assert InstalledTool("samtools").version == "samtools 1.21"

    # Absence is never remembered, either. Only an answer a binary actually gave is
    # kept, so locating stays per object and a tool installed midway through a process
    # is found.
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))

    with pytest.raises(ToolNotFoundError):
        _ = InstalledTool("faToTwoBit").version

    stub_binary(empty, "faToTwoBit", "echo 'faToTwoBit v456'")

    assert InstalledTool("faToTwoBit").version == "faToTwoBit v456"


def test_clearing_the_cache_sends_the_next_ask_back_but_a_recording_tool_bypasses_it(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    # The documented way out, and what the suite's autouse fixture calls: without it a
    # test that stubs a tool could be answered from what another test learned.
    tally = tmp_path / "asked"
    on_path("samtools", f"echo x >> {tally}; echo 'samtools 1.21'")
    assert InstalledTool("samtools").version == "samtools 1.21"

    clear_version_cache()

    assert InstalledTool("samtools").version == "samtools 1.21"
    assert tally.read_text() == "x\nx\n"

    # The cache belongs to the installed adapter, where the subprocess is. Stand-ins of
    # one name share a path and report the versions they were each told to all the same —
    # which every aligner test in the suite is standing on.
    assert RecordingTool("STAR", version="2.7.11b").version == "2.7.11b"
    assert RecordingTool("STAR", version="2.7.10a").version == "2.7.10a"


# -- running a tool ---------------------------------------------------------


def test_run_returns_captured_stdout_with_arguments_and_working_directory_honoured(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    on_path("samtools", 'echo "ran: $@"')
    assert InstalledTool("samtools").run(["faidx", "hg38.fa"]) == "ran: faidx hg38.fa\n"

    on_path("STAR", "pwd")
    where = tmp_path / "index"
    where.mkdir()
    assert InstalledTool("STAR").run(["--runMode"], cwd=where).strip() == str(where)


def test_a_captured_failure_carries_stderr_and_a_missing_tool_raises_before_anything_else(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo boom >&2; exit 1")

    with pytest.raises(RuntimeError) as raised:
        InstalledTool("samtools").run(["faidx", "nope.fa"])

    message = str(raised.value)
    assert "samtools failed (exit 1)" in message
    assert "boom" in message
    assert "'nope.fa'" in message  # the args, so the failing call is identifiable

    with pytest.raises(ToolNotFoundError, match="pixi"):
        InstalledTool("definitely-not-a-real-tool-xyz").run([])


def test_run_with_no_capture_streams_instead_of_returning_and_still_raises_on_failure(
    on_path: Callable[[str, str], Path], capfd: pytest.CaptureFixture[str]
) -> None:
    on_path("STAR", "echo 'started ..... done'")
    assert InstalledTool("STAR").run(["--runMode"], capture=False) == ""
    assert "started" in capfd.readouterr().out

    # The long-build flavour: nothing is captured, so the tool's diagnostics reach the
    # console live and the message sends the reader there rather than repeating nothing.
    on_path("STAR", "echo 'EXITING because of FATAL ERROR' >&2; exit 1")
    with pytest.raises(RuntimeError, match="see its output above") as raised:
        InstalledTool("STAR").run(["--runMode", "genomeGenerate"], capture=False)

    assert "STAR failed (exit 1)" in str(raised.value)
    assert "FATAL ERROR" in capfd.readouterr().err


# -- the freshness rule, and the run built on it ----------------------------


def test_is_fresh_rules(tmp_path: Path, touch_newer_than: Callable[..., None]) -> None:
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    assert is_fresh(out, [src]) is False  # missing output
    out.write_text("")
    assert is_fresh(out, [src]) is False  # empty output
    out.write_text("y")
    touch_newer_than(out, src)
    assert is_fresh(out, [src]) is True  # non-empty and newer
    touch_newer_than(src, out)
    assert is_fresh(out, [src]) is False  # input regenerated -> stale

    # The caller validates the inputs it requires; a missing one never makes a built
    # output look stale, which would rebuild it forever.
    assert is_fresh(out, [tmp_path / "never-existed"]) is True


def test_run_to_runs_missing_skips_fresh_reruns_stale_and_overwrite_forces(
    tmp_path: Path, touch_newer_than: Callable[..., None], monkeypatch: pytest.MonkeyPatch
) -> None:
    tool = RecordingTool("samtools")
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    # missing output -> runs
    result = tool.run_to(["build", str(out)], output=out, inputs=[src])
    assert result == out
    assert [call.args for call in tool.calls] == [("build", str(out))]

    out.write_text("cached")
    touch_newer_than(out, src)
    # fresh -> skipped, the tool not run again
    assert tool.run_to(["build"], output=out, inputs=[src]) == out
    assert len(tool.calls) == 1

    touch_newer_than(src, out)  # input regenerated after the output
    # stale -> reruns
    assert tool.run_to(["build"], output=out, inputs=[src]) == out
    assert [call.args for call in tool.calls][-1] == ("build",)
    assert len(tool.calls) == 2

    touch_newer_than(out, src)
    # fresh again, but overwrite forces the run regardless
    tool.run_to(["build"], output=out, inputs=[src], overwrite=True)
    assert len(tool.calls) == 3

    # And freshness is decided before anything is located, so re-preparing an assembly
    # on a machine that has lost samtools still answers rather than raising.
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "nowhere" / "python"))
    real_src = tmp_path / "real-in"
    real_src.write_text("x")
    real_out = tmp_path / "real-out"
    real_out.write_text("cached")
    touch_newer_than(real_out, real_src)

    result = InstalledTool("samtools").run_to(["faidx"], output=real_out, inputs=[real_src])
    assert result == real_out


# -- the recording adapter --------------------------------------------------


def test_the_recording_tool_is_a_real_adapter_that_records_calls_back_and_fails_like_one(
    tmp_path: Path,
) -> None:
    tool = RecordingTool("STAR", version="2.7.11b")
    tool.run(["--runMode", "genomeGenerate"], cwd=tmp_path, capture=False)
    (call,) = tool.calls
    assert call.args == ("--runMode", "genomeGenerate")
    assert call.cwd == tmp_path
    assert call.capture is False
    assert tool.version == "2.7.11b"
    assert tool.path == "/fake/STAR"

    seen: list[Path | None] = []
    callback_tool = RecordingTool("STAR", on_run=lambda call: seen.append(call.cwd))
    callback_tool.run(["a"], cwd=tmp_path)
    callback_tool.run(["b"])
    assert seen == [tmp_path, None]

    # The point of it being a real adapter: the failure message is the base class's, so
    # a test driving the failure path is driving the code that ships.
    failing_tool = RecordingTool("chromap")
    failing_tool.exit_code = 1
    with pytest.raises(RuntimeError, match=r"chromap failed \(exit 1\)"):
        failing_tool.run(["--build-index"])

    assert issubclass(InstalledTool, ExternalTool)
    assert issubclass(RecordingTool, ExternalTool)


# -- doctor -----------------------------------------------------------------


def test_doctor_checks_exactly_the_tools_preparing_an_assembly_uses() -> None:
    # The defect this closes: `doctor` used to check bedtools, which the package never
    # shells out to, and to check neither of the two UCSC binaries that prepare_fasta
    # cannot run without. Pinned against fasta's own list so the two cannot drift.
    assert set(REQUIRED_TOOLS) == set(PREPARATION_TOOLS)
    assert "bedtools" not in REQUIRED_TOOLS


def test_doctor_raises_naming_a_missing_tool_and_reports_one_that_declines_to_identify(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo 'samtools 1.21'")

    with pytest.raises(ToolNotFoundError, match="faToTwoBit is not installed"):
        doctor()

    on_path("faToTwoBit", "echo 'nope' >&2; exit 255")
    on_path("twoBitInfo", "echo 'nope' >&2; exit 255")

    assert doctor() == {
        "samtools": "samtools 1.21",
        "faToTwoBit": NO_VERSION_REPORTED,
        "twoBitInfo": NO_VERSION_REPORTED,
    }


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="preparation tools not on PATH")
def test_doctor_and_version_answer_for_every_required_tool_against_the_real_binaries() -> None:
    # The claim `doctor` above cannot make: non-empty is not the same as *this* tool's,
    # and a resolver pointing at the wrong binary would satisfy non-empty happily.
    assert "samtools" in InstalledTool("samtools").version.lower()

    result = doctor()
    assert set(result) == set(REQUIRED_TOOLS)
    for name, reported in result.items():
        assert reported.strip() != "", f"{name} reported nothing at all"


def test_the_interpreter_is_a_tool_like_any_other() -> None:
    # A real binary, resolved and run and asked its version, with nothing stubbed. An
    # absolute path resolves as itself, which is what makes this need nothing on PATH.
    tool = InstalledTool(sys.executable)

    assert tool.version.startswith("Python")
    assert tool.run(["-c", "print('ok')"]) == "ok\n"
