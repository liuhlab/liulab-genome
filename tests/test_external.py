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


def test_a_tool_on_path_resolves_to_it(on_path: Callable[[str, str], Path]) -> None:
    written = on_path("faToTwoBit", "exit 0")

    assert InstalledTool("faToTwoBit").path == str(written)


def test_resolution_falls_back_to_the_interpreters_own_bin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stub_binary: StubBinary
) -> None:
    # Running the env's interpreter without PATH activated: the normal lookup misses,
    # but the tool sits beside sys.executable (the conda/pixi bin/), so resolution falls
    # back to that directory.
    bin_dir = tmp_path / "bin"
    written = stub_binary(bin_dir, "faToTwoBit", "exit 0")

    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(bin_dir / "python"))

    assert InstalledTool("faToTwoBit").path == str(written)


def test_a_tool_that_is_nowhere_raises_with_the_command_that_installs_it(
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


def test_a_tool_nobody_catalogued_still_names_a_plausible_command() -> None:
    assert "pixi add bowtie2" in RecordingTool("bowtie2").install_instructions()


def test_install_instructions_quote_the_tools_own_homepage() -> None:
    assert "https://github.com/alexdobin/STAR" in RecordingTool("STAR").install_instructions()


def test_the_path_is_located_once(on_path: Callable[[str, str], Path]) -> None:
    written = on_path("samtools", "exit 0")
    tool = InstalledTool("samtools")

    assert tool.path == str(written)
    written.unlink()  # gone from disk; the answer was remembered, not looked up again

    assert tool.path == str(written)


# -- asking a tool what it is -----------------------------------------------


def test_version_is_the_first_line_of_what_the_tool_says(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo 'samtools 1.21'; echo 'Using htslib 1.21'")

    assert InstalledTool("samtools").version == "samtools 1.21"


def test_version_falls_back_to_stderr(on_path: Callable[[str, str], Path]) -> None:
    # Different tools choose differently; chromap answers on stderr.
    on_path("chromap", "echo '0.3.2-r518' >&2")

    assert InstalledTool("chromap").version == "0.3.2-r518"


def test_a_tool_that_rejects_the_flag_reports_no_version_rather_than_raising(
    on_path: Callable[[str, str], Path],
) -> None:
    # What several UCSC binaries do: `faToTwoBit --version` exits non-zero saying the
    # option is not valid. That is a tool declining to identify itself, not a failure.
    on_path("faToTwoBit", "echo '--version is not a valid option' >&2; exit 255")

    assert InstalledTool("faToTwoBit").version == ""


def test_a_tool_that_is_absent_raises_rather_than_reporting_no_version(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # The two answers must stay apart: "" means it ran and would not say.
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "bin" / "python"))

    with pytest.raises(ToolNotFoundError):
        _ = InstalledTool("faToTwoBit").version


def test_the_version_is_asked_for_once(on_path: Callable[[str, str], Path], tmp_path: Path) -> None:
    # Recording a build's provenance must not cost a subprocess per mention of it.
    tally = tmp_path / "asked"
    on_path("samtools", f"echo x >> {tally}; echo 'samtools 1.21'")
    tool = InstalledTool("samtools")

    assert tool.version == tool.version == "samtools 1.21"
    assert tally.read_text() == "x\n"


def test_the_version_is_asked_for_once_per_binary_not_once_per_object(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    # The waste the per-object cache above never reached: a build constructs a fresh tool
    # for each step, so preparing one assembly asked samtools its version once per step
    # for an answer that cannot change under a running process.
    tally = tmp_path / "asked"
    on_path("samtools", f"echo x >> {tally}; echo 'samtools 1.21'")

    assert InstalledTool("samtools").version == "samtools 1.21"
    assert InstalledTool("samtools").version == "samtools 1.21"

    assert tally.read_text() == "x\n"


def test_a_tool_that_declines_to_identify_itself_is_remembered_as_declining(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    # "" is an answer rather than a missing one, and it costs the same subprocess to
    # learn. Remembering on truthiness instead of on presence would re-probe exactly the
    # two UCSC binaries every preparation runs, which is most of what there was to save.
    tally = tmp_path / "asked"
    on_path("faToTwoBit", f"echo x >> {tally}; echo 'nope' >&2; exit 255")

    assert InstalledTool("faToTwoBit").version == ""
    assert InstalledTool("faToTwoBit").version == ""

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


def test_a_tool_that_was_missing_is_found_once_it_is_installed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, stub_binary: StubBinary
) -> None:
    # Absence is never remembered. Only an answer a binary actually gave is kept, so
    # locating stays per object and a tool installed midway through a process is found.
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "nowhere" / "python"))

    with pytest.raises(ToolNotFoundError):
        _ = InstalledTool("faToTwoBit").version

    stub_binary(bin_dir, "faToTwoBit", "echo 'faToTwoBit v456'")

    assert InstalledTool("faToTwoBit").version == "faToTwoBit v456"


def test_clearing_the_cache_sends_the_next_ask_back_to_the_binary(
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


def test_a_recording_tool_answers_for_itself_rather_than_from_the_shared_cache() -> None:
    # The cache belongs to the installed adapter, where the subprocess is. Stand-ins of
    # one name share a path and report the versions they were each told to all the same —
    # which every aligner test in the suite is standing on.
    assert RecordingTool("STAR", version="2.7.11b").version == "2.7.11b"
    assert RecordingTool("STAR", version="2.7.10a").version == "2.7.10a"


# -- running a tool ---------------------------------------------------------


def test_run_returns_captured_stdout(on_path: Callable[[str, str], Path]) -> None:
    on_path("samtools", "echo hello")

    assert InstalledTool("samtools").run(["faidx"]) == "hello\n"


def test_run_passes_its_arguments_through(on_path: Callable[[str, str], Path]) -> None:
    on_path("samtools", 'echo "$@"')

    assert InstalledTool("samtools").run(["faidx", "hg38.fa"]) == "faidx hg38.fa\n"


def test_run_uses_the_working_directory_it_is_given(
    on_path: Callable[[str, str], Path], tmp_path: Path
) -> None:
    on_path("STAR", "pwd")
    where = tmp_path / "index"
    where.mkdir()

    assert InstalledTool("STAR").run(["--runMode"], cwd=where).strip() == str(where)


def test_a_captured_failure_carries_the_tools_own_stderr(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo boom >&2; exit 1")

    with pytest.raises(RuntimeError) as raised:
        InstalledTool("samtools").run(["faidx", "nope.fa"])

    message = str(raised.value)
    assert "samtools failed (exit 1)" in message
    assert "boom" in message
    assert "'nope.fa'" in message  # the args, so the failing call is identifiable


def test_an_inherited_failure_points_at_the_output_already_printed(
    on_path: Callable[[str, str], Path], capfd: pytest.CaptureFixture[str]
) -> None:
    # The long-build flavour: nothing is captured, so the tool's diagnostics reach the
    # console live and the message sends the reader there rather than repeating nothing.
    on_path("STAR", "echo 'EXITING because of FATAL ERROR' >&2; exit 1")

    with pytest.raises(RuntimeError, match="see its output above") as raised:
        InstalledTool("STAR").run(["--runMode", "genomeGenerate"], capture=False)

    assert "STAR failed (exit 1)" in str(raised.value)
    assert "FATAL ERROR" in capfd.readouterr().err


def test_run_with_no_capture_streams_rather_than_returning(
    on_path: Callable[[str, str], Path], capfd: pytest.CaptureFixture[str]
) -> None:
    on_path("STAR", "echo 'started ..... done'")

    assert InstalledTool("STAR").run(["--runMode"], capture=False) == ""
    assert "started" in capfd.readouterr().out


def test_running_a_tool_that_is_not_there_raises_before_anything_else() -> None:
    with pytest.raises(ToolNotFoundError, match="pixi"):
        InstalledTool("definitely-not-a-real-tool-xyz").run([])


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


def test_an_input_that_does_not_exist_is_ignored(tmp_path: Path) -> None:
    # The caller validates the inputs it requires; a missing one never makes a built
    # output look stale, which would rebuild it forever.
    out = tmp_path / "out"
    out.write_text("y")

    assert is_fresh(out, [tmp_path / "never-existed"]) is True


def test_run_to_runs_when_the_output_is_missing(tmp_path: Path) -> None:
    tool = RecordingTool("samtools")
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"

    result = tool.run_to(["build", str(out)], output=out, inputs=[src])

    assert result == out
    assert [call.args for call in tool.calls] == [("build", str(out))]


def test_run_to_skips_the_tool_when_the_output_is_fresh(
    tmp_path: Path, touch_newer_than: Callable[..., None]
) -> None:
    tool = RecordingTool("samtools")
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    assert tool.run_to(["build"], output=out, inputs=[src]) == out
    assert tool.calls == []


def test_run_to_reruns_when_an_input_is_newer(
    tmp_path: Path, touch_newer_than: Callable[..., None]
) -> None:
    # The stale case, driven through `run_to` rather than through `is_fresh` alone: an
    # output that exists is not thereby fresh, and a regenerated input must rebuild it.
    tool = RecordingTool("samtools")
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(src, out)  # input regenerated after the output

    assert tool.run_to(["build"], output=out, inputs=[src]) == out
    assert [call.args for call in tool.calls] == [("build",)]


def test_run_to_overwrite_forces_the_run(
    tmp_path: Path, touch_newer_than: Callable[..., None]
) -> None:
    tool = RecordingTool("samtools")
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    tool.run_to(["build"], output=out, inputs=[src], overwrite=True)

    assert len(tool.calls) == 1


def test_a_fresh_output_is_served_without_the_tool_being_installed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, touch_newer_than: Callable[..., None]
) -> None:
    # Freshness is decided before anything is located, so re-preparing an assembly on a
    # machine that has lost samtools still answers rather than raising.
    monkeypatch.setenv("PATH", "")
    monkeypatch.setattr(external.sys, "executable", str(tmp_path / "bin" / "python"))
    src = tmp_path / "in"
    src.write_text("x")
    out = tmp_path / "out"
    out.write_text("cached")
    touch_newer_than(out, src)

    assert InstalledTool("samtools").run_to(["faidx"], output=out, inputs=[src]) == out


# -- the recording adapter --------------------------------------------------


def test_the_recording_tool_records_everything_a_call_carries(tmp_path: Path) -> None:
    tool = RecordingTool("STAR", version="2.7.11b")

    tool.run(["--runMode", "genomeGenerate"], cwd=tmp_path, capture=False)

    (call,) = tool.calls
    assert call.args == ("--runMode", "genomeGenerate")
    assert call.cwd == tmp_path
    assert call.capture is False
    assert tool.version == "2.7.11b"
    assert tool.path == "/fake/STAR"


def test_the_recording_tool_fails_the_way_a_real_one_does() -> None:
    # The point of it being a real adapter: the failure message is the base class's, so a
    # test driving the failure path is driving the code that ships.
    tool = RecordingTool("chromap")
    tool.exit_code = 1

    with pytest.raises(RuntimeError, match=r"chromap failed \(exit 1\)"):
        tool.run(["--build-index"])


def test_the_recording_tool_calls_back_as_each_call_is_made(tmp_path: Path) -> None:
    seen: list[Path | None] = []
    tool = RecordingTool("STAR", on_run=lambda call: seen.append(call.cwd))

    tool.run(["a"], cwd=tmp_path)
    tool.run(["b"])

    assert seen == [tmp_path, None]


def test_both_adapters_are_the_same_interface() -> None:
    assert issubclass(InstalledTool, ExternalTool)
    assert issubclass(RecordingTool, ExternalTool)


# -- doctor -----------------------------------------------------------------


def test_doctor_checks_exactly_the_tools_preparing_an_assembly_uses() -> None:
    # The defect this closes: `doctor` used to check bedtools, which the package never
    # shells out to, and to check neither of the two UCSC binaries that prepare_fasta
    # cannot run without. Pinned against fasta's own list so the two cannot drift.
    assert set(REQUIRED_TOOLS) == set(PREPARATION_TOOLS)
    assert "bedtools" not in REQUIRED_TOOLS


def test_doctor_reports_a_tool_that_will_not_identify_itself_rather_than_hiding_it(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo 'samtools 1.21'")
    on_path("faToTwoBit", "echo 'nope' >&2; exit 255")
    on_path("twoBitInfo", "echo 'nope' >&2; exit 255")

    assert doctor() == {
        "samtools": "samtools 1.21",
        "faToTwoBit": NO_VERSION_REPORTED,
        "twoBitInfo": NO_VERSION_REPORTED,
    }


def test_doctor_raises_naming_the_tool_that_is_missing(
    on_path: Callable[[str, str], Path],
) -> None:
    on_path("samtools", "echo 'samtools 1.21'")

    with pytest.raises(ToolNotFoundError, match="faToTwoBit is not installed"):
        doctor()


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="preparation tools not on PATH")
def test_a_version_string_identifies_the_tool_it_came_from() -> None:
    # The claim `doctor` below cannot make: non-empty is not the same as *this* tool's,
    # and a resolver pointing at the wrong binary would satisfy non-empty happily.
    assert "samtools" in InstalledTool("samtools").version.lower()


@pytest.mark.skipif(not _BINARIES_PRESENT, reason="preparation tools not on PATH")
def test_doctor_answers_for_every_required_tool_against_the_real_binaries() -> None:
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
