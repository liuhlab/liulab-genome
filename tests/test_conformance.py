"""Negative tests for `scripts/conformance.py` — one tree per rule that VIOLATES it.

A rule that checks nothing passes the happy path, and a gate that fires on nothing looks exactly
like a gate that passes. That already happened once here: overlapping Vale globs left every ADR
checked by no rule at all, with zero alerts and green CI. So the claim these tests make is not
"the template conforms" — `pixi run check` already says that — but "each rule can still fail".

Every test builds a tmp_path tree that is correct in every respect but one, and asserts the run
marks THAT rule FAIL and names the offending file. The rule name alone would prove nothing: the
status table prints all fifteen names on every run, passing or failing, so the assertions read the
status word out of the table rather than grepping the output for a name.

The exit status is checked here too, all three of it: 0, 1 when a rule failed, and 2 when the
checker could not run at all. A checker that could not run has checked nothing, and it must not
report that as a rule this repo broke.

Conformance is invoked as a SUBPROCESS, so what is under test is the command an agent runs and
CI runs, not an internal function that could be refactored into something the command no longer
calls.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CONFORMANCE = REPO / "scripts" / "conformance.py"
INSTALLER = REPO / "skills" / "install.py"

# The one thing this file cannot assume in a derived repo: rule 7 DELEGATES to the installer, so
# a tree with no `skills/install.py` to copy in has nothing for that rule to run and it reports
# `--`. This repo ships no skills tree, so the two below are what keep the rest of the file
# running here and unchanged in a repo that does have one: with the installer present the
# verdict is `ok` and nothing skips, which is the template's own state.
SKILL_SYMLINKS_PASS = "ok" if INSTALLER.is_file() else "--"
needs_installer = pytest.mark.skipif(
    not INSTALLER.is_file(),
    reason="rule 7 delegates to skills/install.py, and this repo ships no skills tree",
)

# Spelled in two pieces for the same reason `conformance.py` spells it that way: `init-repo`
# substitutes the placeholder through the whole tree and the dogfood job then greps a rendered
# repo for it. A literal here would be rewritten, and would be found.
PLACEHOLDER = "new" + "pkg"

FRONT_MATTER = "---\nsearch:\n  exclude: true\n---\n\n"

PYPROJECT = """\
[project]
name = "example"

[tool.pixi.tasks]
build = "python -m build"
test = "pytest"

[tool.liulab.agent-docs]
"AGENTS.md" = "LengthDoc"
"CONTEXT.md" = "LengthDoc"
"docs/agents/" = "LengthDoc"
"skills/" = "LengthDoc"
"docs/adr/" = "LengthAdr"
"docs/research/" = false

[tool.liulab.waived]
"""

# Four sections, every Length rule named in every one, exactly as the real `.vale.ini` does and
# for the same reason: Vale accumulates per-rule settings across the sections that match a file.
VALE_INI = """\
StylesPath = styles
MinAlertLevel = error

[{README.md,CHANGELOG.md,docs/**/*.md}]
BasedOnStyles = Lab
Lab.LengthDoc = NO
Lab.LengthAdr = NO

[{AGENTS.md,CONTEXT.md,docs/agents/**/*.md,skills/**/SKILL.md}]
BasedOnStyles = Lab
Lab.LengthDoc = YES
Lab.LengthAdr = NO

[docs/adr/**/*.md]
BasedOnStyles = Lab
Lab.LengthDoc = NO
Lab.LengthAdr = YES

[docs/research/**/*.md]
BasedOnStyles = Lab
Lab.LengthDoc = NO
Lab.LengthAdr = NO
"""

# The nav is NESTED, which the template's is not. A flat walk would pass every one of these tests
# and let an agent-facing page into a derived repo's navbar.
MKDOCS = """\
site_name: example
nav:
  - Home: index.md
  - Reference:
      - API: api.md
"""

# The publishing workflow, spelled out here rather than imported from the checker: a test that
# borrowed the path would agree with rule 9 by construction and prove nothing about the
# convention. Rule 8 wants a CHANGELOG.md wherever this file is, so `publishes` writes both.
RELEASE_YML = ".github/workflows/release.yml"

_STATUS_RE = re.compile(
    r"^\s+(?:rule|warn) \d+\s+(?P<name>\S+)\s+(?P<status>ok|FAIL|--|warn|waived)(\s|$)"
)


def write(root: Path, rel: str, text: str) -> Path:
    """Put one file in the tree, making its parents."""
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def publishes(root: Path, on: str) -> None:
    """Give the tree a publishing workflow with one `on:` block, and the changelog rule 8 wants.

    The body is a real job with real steps, not a bare `on:` block: the accessor the rule reads
    parses the whole file, and a fixture that carried triggers alone would never exercise that.
    Its one command invokes `build`, which `PYPROJECT` declares, so rule 12 stays green here and
    rule 9 is the only thing these trees vary.
    """
    write(
        root,
        RELEASE_YML,
        f"name: release\n{on}"
        "jobs:\n"
        "  build:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v7\n"
        "      - run: pixi run build\n",
    )
    write(root, "CHANGELOG.md", "# Changelog\n\nNothing yet.\n")


def runs(root: Path, *commands: str, job: str = "check") -> None:
    """A `ci.yml` whose one job checks out an action and then runs each command given.

    Not `release.yml`: that path is rule 9's premise and rule 8's, and a rule 12 fixture must vary
    nothing but the commands. The leading `uses:` step is in every tree here on purpose — it is
    what proves an action step is passed over rather than merely absent.
    """
    steps = "".join(f"      - run: {command}\n" for command in commands)
    write(
        root,
        ".github/workflows/ci.yml",
        "name: ci\non:\n  pull_request:\njobs:\n"
        f"  {job}:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v7\n" + steps,
    )


def add_skill(root: Path, name: str) -> None:
    """A skill and both committed discovery symlinks, relative, as `install.py` writes them."""
    write(root, f"skills/{name}/SKILL.md", f"---\nname: {name}\n---\n\n# {name}\n")
    for discovery in (".claude/skills", ".agents/skills"):
        link = root / discovery / name
        link.parent.mkdir(parents=True, exist_ok=True)
        link.symlink_to(f"../../skills/{name}", target_is_directory=True)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A tree that passes every rule, for one test to break in exactly one way.

    It has no `src/`, no workflow at all and so no release workflow, so rule 8's two conditionals,
    rule 9 and rule 11 are all vacuous here; the tests that care add the premise, through
    `publishes` and `runs`. It has no `skills/init-repo/`, so rule 1 and warning 14 are both live
    — the state a derived repo is in, and the only state where they check anything.
    """
    root = tmp_path / "repo"
    write(root, "pyproject.toml", PYPROJECT)
    write(root, ".vale.ini", VALE_INI)
    write(root, "mkdocs.yml", MKDOCS)
    write(root, "AGENTS.md", "# example\n\nWhat this repo is and how to work in it.\n")
    write(root, "CONTEXT.md", "# Context\n\n## Glossary\n\n### Term\n\nWhat the word means here.\n")
    write(root, "docs/index.md", "# example\n\nThe front page.\n")
    write(root, "docs/api.md", "# API\n\nThe reference.\n")
    write(root, "docs/agents/writing.md", f"{FRONT_MATTER}# Writing rules\n\nBe concise.\n")
    write(root, "docs/adr/0001-a-decision.md", f"{FRONT_MATTER}# A decision\n\nWe chose this.\n")
    write(root, "docs/research/a-note.md", f"{FRONT_MATTER}# A note\n\nWhat was found.\n")
    # The real installer, not a stand-in: rule 7 delegates to this file, so a copy of it is what
    # makes the delegation the thing under test. Absent in a repo with no skills tree, where
    # rule 7 has nothing to delegate to and says so.
    (root / "skills").mkdir(parents=True, exist_ok=True)
    if INSTALLER.is_file():
        shutil.copy(INSTALLER, root / "skills" / "install.py")
    subprocess.run(["git", "init", "-q"], cwd=root, check=True)
    return root


def conformance(
    root: Path, *, stage: bool = True, ceiling: Path | None = None
) -> subprocess.CompletedProcess[str]:
    """Run the checker over one tree, the way a person and CI run it.

    `git add -A` first, because every rule reads `git ls-files`: until a fixture's edit is staged
    it does not exist as far as conformance is concerned. Nothing is committed — no rule reads
    history, and a commit would need an identity this environment may not have. `stage=False` is
    for the one tree that is not a repo at all, where there is nothing to stage into.

    `ceiling` is what makes such a tree genuinely repo-less. `git ls-files` walks UP from its
    working directory, so a scratch directory inherits whatever repository happens to be above it —
    and a test that relied on the temporary directory sitting outside one would pass or fail on
    where the machine puts its temporary files. `GIT_CEILING_DIRECTORIES` stops the walk there.
    """
    if stage:
        subprocess.run(["git", "add", "-A"], cwd=root, check=True)
    env = dict(os.environ)
    if ceiling is not None:
        env["GIT_CEILING_DIRECTORIES"] = str(ceiling)
    return subprocess.run(
        [sys.executable, str(CONFORMANCE), "--root", str(root)],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def statuses(proc: subprocess.CompletedProcess[str]) -> dict[str, str]:
    """The verdict the run gave each rule, read off its own status table.

    Asserting on the status word rather than on the presence of a rule NAME is the point: the
    table prints all fifteen names on every run, so `"repo-shape" in output` is true of a green run
    and would prove nothing at all.
    """
    found: dict[str, str] = {}
    for line in proc.stdout.splitlines():
        match = _STATUS_RE.match(line)
        if match:
            found[match["name"]] = match["status"]
    return found


def flat(text: str) -> str:
    """One long line, so an assertion is not hostage to where a note happened to wrap."""
    return " ".join(text.split())


def test_the_fixture_tree_passes_every_rule(repo: Path) -> None:
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    verdicts = statuses(proc)
    assert set(verdicts) == {
        "placeholder-rename",
        "agent-docs-unpublished",
        "vale-length-sections",
        "glossary-entry-length",
        "skill-file-location",
        "skill-name-prefix",
        "skill-symlinks",
        "repo-shape",
        "release-trigger",
        "single-toolchain",
        "python-version-agreement",
        "workflow-step-tasks",
        "nav-target-exists",
        "init-sentinel",
        "waivers",
    }
    assert verdicts["placeholder-rename"] == "ok"  # live, not gated off
    assert verdicts["agent-docs-unpublished"] == "ok"
    assert verdicts["vale-length-sections"] == "ok"


def test_a_tree_that_has_skills_passes_every_rule(repo: Path) -> None:
    # The worktree this was written in had no skills, so the rules that read `skills/` would have
    # been vacuous in every other test here. This is the tree they are written for.
    add_skill(repo, "template-dev")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    verdicts = statuses(proc)
    assert verdicts["skill-file-location"] == "ok"
    assert verdicts["skill-name-prefix"] == "ok"
    assert verdicts["skill-symlinks"] == SKILL_SYMLINKS_PASS
    assert "1 skill" in proc.stdout


# The exit status is a three-way answer and not a boolean, so all three are asserted here rather
# than left implied by the rule tests below. "A rule failed" means fix the repo; "could not run"
# means the checker never ran and the repo was never checked, and one status for both sends the
# reader of a CI log hunting a violation that does not exist. 2 is the same word
# `scripts/check.sh` uses for a usage error and for a capture directory it could not make.


def test_it_exits_0_when_it_ran_and_every_rule_passed(repo: Path) -> None:
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "All 13 rules pass." in proc.stdout


def test_it_exits_1_when_it_ran_and_a_rule_failed(repo: Path) -> None:
    write(repo, "README.md", f"# liulab-{PLACEHOLDER}\n\nA repo that was never renamed.\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert "1 of 13 rules failed" in proc.stderr


def test_it_exits_2_when_there_is_no_git_to_list_tracked_files_from(tmp_path: Path) -> None:
    # Every rule reads `git ls-files`, so without a repository NOTHING was checked — the state CI
    # would be in given a source archive instead of a checkout. The tree is otherwise complete, so
    # the missing `.git` is the only thing under test.
    plain = tmp_path / "plain"
    write(plain, "pyproject.toml", PYPROJECT)
    proc = conformance(plain, stage=False, ceiling=tmp_path)
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "cannot list tracked files" in proc.stderr
    assert "not a rule failure" in proc.stderr


def test_it_exits_2_when_there_is_no_pyproject_toml(repo: Path) -> None:
    # The declarations every rule reads about themselves live there — the agent-docs table, the
    # waivers, the task list. Without it the answer is unknown, not bad.
    (repo / "pyproject.toml").unlink()
    proc = conformance(repo)
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "no pyproject.toml" in proc.stderr


def test_it_exits_2_when_the_manifest_will_not_parse(repo: Path) -> None:
    # Same class as an absent manifest, and it must not be a traceback: an uncaught exception exits
    # 1, which would read as a failed rule with an unusually bad report.
    write(repo, "pyproject.toml", "[project\nname = 'example'\n")
    proc = conformance(repo)
    assert proc.returncode == 2, proc.stdout + proc.stderr
    assert "cannot read pyproject.toml" in proc.stderr
    assert "Traceback" not in proc.stderr


def test_rule_1_fires_on_a_tracked_placeholder(repo: Path) -> None:
    write(repo, "README.md", f"# liulab-{PLACEHOLDER}\n\nA repo that was never renamed.\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["placeholder-rename"] == "FAIL"
    assert f"README.md still names `{PLACEHOLDER}` in its contents" in proc.stderr


def test_rule_1_fires_on_a_placeholder_in_a_path_not_only_in_a_file(repo: Path) -> None:
    # A `py.typed` under the placeholder module holds no text to grep. The directory name is
    # the whole violation. (Spelled around, not out: rule 1 greps this file too.)
    write(repo, f"src/{PLACEHOLDER}/py.typed", "")
    write(repo, "tests/test_it.py", "def test_it() -> None:\n    assert True\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["placeholder-rename"] == "FAIL"
    assert f"src/{PLACEHOLDER}/py.typed still names `{PLACEHOLDER}` in its path" in proc.stderr


def test_rule_1_is_not_checked_while_the_init_skill_is_there(repo: Path) -> None:
    # The template ships the placeholder ON PURPOSE, and `skills/init-repo/` is what says the
    # rename has not happened yet. Same discriminator warning 14 uses.
    write(repo, "README.md", f"# liulab-{PLACEHOLDER}\n\nStill the template.\n")
    add_skill(repo, "init-repo")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["placeholder-rename"] == "--"
    assert "not checked: skills/init-repo/ is here" in flat(proc.stdout)


def test_rule_1_refuses_to_be_waived(repo: Path) -> None:
    write(repo, "README.md", f"# liulab-{PLACEHOLDER}\n\nA repo that was never renamed.\n")
    write(
        repo,
        "pyproject.toml",
        PYPROJECT + 'placeholder-rename = "we like the name"\n',
    )
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["placeholder-rename"] == "FAIL"
    assert "cannot be waived" in proc.stderr
    assert "REFUSED" in proc.stdout  # and the waiver table still prints it


def test_rule_2_fires_on_a_page_that_is_not_excluded_from_search(repo: Path) -> None:
    write(repo, "docs/agents/writing.md", "# Writing rules\n\nBe concise.\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["agent-docs-unpublished"] == "FAIL"
    assert "docs/agents/writing.md has no `search: exclude: true` front matter" in proc.stderr


def test_rule_2_fires_on_an_agent_page_in_the_nav(repo: Path) -> None:
    # Nav paths are relative to `docs_dir`, so `adr/0001-a-decision.md` here IS
    # `docs/adr/0001-a-decision.md` in `git ls-files`. Compare the two without that prefix and
    # the rule passes vacuously — which is the failure this whole file exists to catch.
    write(repo, "mkdocs.yml", MKDOCS + "      - ADR: adr/0001-a-decision.md\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["agent-docs-unpublished"] == "FAIL"
    assert (
        "docs/adr/0001-a-decision.md appears in mkdocs.yml nav as `adr/0001-a-decision.md`"
        in proc.stderr
    )


def test_rule_2_follows_a_docs_dir_the_site_moved(repo: Path) -> None:
    write(
        repo, "mkdocs.yml", "site_name: example\ndocs_dir: site-source\nnav:\n  - Home: index.md\n"
    )
    # The page the moved nav names, so rule 13 is green here and rule 2 is the only rule this
    # tree is about.
    write(repo, "site-source/index.md", "# example\n\nThe front page.\n")
    proc = conformance(repo)
    # Nothing under `site-source/` is declared agent-facing, so the rule has nothing to check —
    # and says so, rather than reporting a green it did not earn.
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["agent-docs-unpublished"] == "--"
    assert "not checked: no [tool.liulab.agent-docs] key is under site-source/" in flat(proc.stdout)


def test_rule_3_fires_when_a_section_turns_the_declared_cap_off(repo: Path) -> None:
    # The config this reproduces is not hypothetical: it is what left the ADRs checked by
    # nothing. Every rule is present, so asserting mere presence would call this green.
    write(repo, ".vale.ini", VALE_INI.replace("Lab.LengthDoc = YES", "Lab.LengthDoc = NO"))
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["vale-length-sections"] == "FAIL"
    assert "should set Lab.LengthDoc on for `AGENTS.md`, and says `NO`" in proc.stderr


def test_rule_3_fires_when_a_key_has_no_section_at_all(repo: Path) -> None:
    truncated = VALE_INI.split("[docs/adr/**/*.md]")[0]
    write(repo, ".vale.ini", truncated)
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["vale-length-sections"] == "FAIL"
    assert "`docs/adr/` is spelled in 0 .vale.ini section header(s), not 1" in proc.stderr


def test_rule_3_does_not_require_a_key_for_every_header_token(repo: Path) -> None:
    # `CHANGELOG.md` is named in section 1 and is not an agent-docs key. A rule written from the
    # headers inward rather than from the keys outward would fail this tree.
    write(repo, "CHANGELOG.md", "# Changelog\n\nNothing yet.\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["vale-length-sections"] == "ok"


def test_rule_4_fires_on_an_overlong_glossary_entry(repo: Path) -> None:
    body = " ".join(["word"] * 250)
    write(repo, "CONTEXT.md", f"# Context\n\n## Glossary\n\n### Long term\n\n{body}\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["glossary-entry-length"] == "FAIL"
    assert "CONTEXT.md glossary entry `Long term` is 252 words, over 200" in proc.stderr


def test_rule_4_caps_each_entry_and_not_the_file(repo: Path) -> None:
    # Four entries, 150 words each: 600 words of glossary, and no violation. The cap is per
    # entry so a large vocabulary is not punished for the size of its domain.
    body = " ".join(["word"] * 150)
    entries = "".join(f"### Term {n}\n\n{body}\n\n" for n in range(4))
    write(repo, "CONTEXT.md", f"# Context\n\n## Glossary\n\n{entries}")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["glossary-entry-length"] == "ok"


def test_rule_5_fires_on_a_skill_file_the_installer_cannot_see(repo: Path) -> None:
    write(repo, "docs/SKILL.md", "# A skill in the wrong place\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["skill-file-location"] == "FAIL"
    assert "docs/SKILL.md is a SKILL.md outside skills/<name>/" in proc.stderr


def test_rule_6_fires_on_a_skill_that_shadows_the_shared_plugin(repo: Path) -> None:
    # Correctly installed in both discovery paths, so rule 7 is green and the NAME is the only
    # thing wrong — which is the point: this one is invisible to every other check.
    add_skill(repo, "lab-hpc")
    proc = conformance(repo)
    assert proc.returncode == 1
    verdicts = statuses(proc)
    assert verdicts["skill-name-prefix"] == "FAIL"
    assert verdicts["skill-symlinks"] == SKILL_SYMLINKS_PASS
    assert "is a repo-local skill named `lab-hpc`" in proc.stderr


@needs_installer
def test_rule_7_delegates_to_the_installer(repo: Path) -> None:
    # A skill with no symlinks. The assertion is not just that rule 7 fails but that the
    # installer's OWN report comes through, fix line included — proof the rule delegates rather
    # than reimplementing three invariants that would then drift from the installer's.
    write(repo, "skills/demo/SKILL.md", "---\nname: demo\n---\n\n# demo\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["skill-symlinks"] == "FAIL"
    assert ".claude/skills/demo is missing" in proc.stderr
    assert "fix: python skills/install.py --target claude" in proc.stderr


def test_rule_8_fires_on_a_package_with_no_tests(repo: Path) -> None:
    write(repo, "src/example/__init__.py", "")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["repo-shape"] == "FAIL"
    assert "src/ is tracked and tests/ is not" in proc.stderr


def test_rule_8_fires_on_a_release_workflow_with_no_changelog(repo: Path) -> None:
    write(repo, ".github/workflows/release.yml", "on:\n  release:\n    types: [published]\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["repo-shape"] == "FAIL"
    assert ".github/workflows/release.yml is tracked and CHANGELOG.md is not" in proc.stderr


def test_rule_8_is_vacuous_and_not_failing_when_neither_premise_holds(repo: Path) -> None:
    # One lab repo has neither `src/` nor `tests/`. A repo with no package is not a repo that
    # failed to have tests, and the run has to say it checked nothing rather than report a green.
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["repo-shape"] == "--"
    assert "not checked: no tracked src/, so tests/ is not required" in flat(proc.stdout)
    assert "no .github/workflows/release.yml, so CHANGELOG.md is not required" in flat(proc.stdout)


def test_rule_8_is_satisfied_and_not_merely_skipped_when_both_premises_hold(repo: Path) -> None:
    write(repo, "src/example/__init__.py", "")
    write(repo, "tests/test_it.py", "def test_it() -> None:\n    assert True\n")
    write(repo, ".github/workflows/release.yml", "on:\n  release:\n    types: [published]\n")
    write(repo, "CHANGELOG.md", "# Changelog\n\nNothing yet.\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["repo-shape"] == "ok"


def test_rule_9_fires_on_a_publishing_workflow_triggered_by_a_tag_push(repo: Path) -> None:
    publishes(repo, "on:\n  push:\n    tags: ['v*']\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["release-trigger"] == "FAIL"
    assert f"{RELEASE_YML} can be triggered by a tag push: `on: push:` names `tags:`" in proc.stderr
    assert "cannot be undone" in flat(proc.stderr)  # the rule and the fix, not a rewritten assert
    assert "trigger it on `release:` with `types: [published]`" in flat(proc.stderr)


def test_rule_9_fires_on_a_push_with_no_branch_filter(repo: Path) -> None:
    # The hole a rule written only for `tags:` would leave wide open. An absent filter is not
    # "no refs" — GitHub reads it as every ref, so this fires on `v0.0.0` like any other.
    publishes(repo, "on:\n  push:\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["release-trigger"] == "FAIL"
    assert "`on: push:` names no branch filter, so it fires on every ref, tags too" in proc.stderr


def test_rule_9_fires_on_create_which_a_pushed_tag_also_fires(repo: Path) -> None:
    publishes(repo, "on:\n  create:\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["release-trigger"] == "FAIL"
    assert "`on: create:` fires when a tag comes into existence" in proc.stderr


def test_rule_9_passes_a_workflow_triggered_by_a_published_release(repo: Path) -> None:
    publishes(repo, "on:\n  release:\n    types: [published]\n  workflow_dispatch:\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["release-trigger"] == "ok"
    # The note is the anti-vacuity assertion, and it is the only one here that matters. YAML 1.1
    # resolves an unquoted `on` to the boolean true, so an accessor reading `document["on"]` finds
    # no triggers at all — and a rule handed no triggers passes every workflow ever written.
    assert f"{RELEASE_YML} is triggered by release, workflow_dispatch" in flat(proc.stdout)


def test_rule_9_allows_a_push_that_names_branches(repo: Path) -> None:
    # Naming any branch filter is what stops tags reaching a workflow, so this one is a branch
    # trigger and not a violation. A rule that failed every `push:` would fail `ci.yml` too.
    publishes(repo, "on:\n  push:\n    branches: [main]\n  release:\n    types: [published]\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["release-trigger"] == "ok"


def test_rule_9_is_vacuous_and_not_passing_when_the_repo_publishes_nothing(repo: Path) -> None:
    # `init-repo` deletes the workflow for a repo that does not publish. Nothing was checked, and
    # the run has to say so: a green here would read as "the trigger is safe", which is a claim
    # this run has no evidence for.
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["release-trigger"] == "--"
    assert f"not checked: no {RELEASE_YML}, so this repo publishes to no package index" in flat(
        proc.stdout
    )


def test_warning_14_warns_about_the_sentinel_and_does_not_fail(repo: Path) -> None:
    write(
        repo,
        "AGENTS.md",
        "# example\n\n> **Replace this file.** Run `/init`, then delete this line.\n",
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr  # a warning, never a failure
    assert statuses(proc)["init-sentinel"] == "warn"
    assert "AGENTS.md still carries the `/init` sentinel" in flat(proc.stdout)


def test_warning_14_is_silent_while_the_init_skill_is_the_nag(repo: Path) -> None:
    write(
        repo,
        "AGENTS.md",
        "# example\n\n> **Replace this file.** Run `/init`, then delete this line.\n",
    )
    add_skill(repo, "init-repo")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["init-sentinel"] == "--"


def test_a_waiver_suppresses_its_rule_and_is_still_printed(repo: Path) -> None:
    body = " ".join(["word"] * 250)
    write(repo, "CONTEXT.md", f"# Context\n\n## Glossary\n\n### Long term\n\n{body}\n")
    reason = "the domain needs one long entry, tracked in issue 41"
    write(repo, "pyproject.toml", PYPROJECT + f'glossary-entry-length = "{reason}"\n')
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["glossary-entry-length"] == "waived"
    assert f"1 problem(s) suppressed: {reason}" in flat(proc.stdout)
    # and warning 15 prints it a second time, in the waiver table
    assert f"glossary-entry-length: 1 problem(s) suppressed — {reason}" in flat(proc.stdout)


def test_every_waiver_is_printed_on_a_green_run_too(repo: Path) -> None:
    # An escape hatch that stops being visible becomes permanent, so a waiver suppressing nothing
    # is printed exactly as loudly as one suppressing a failure.
    reason = "this repo has no package"
    write(repo, "pyproject.toml", PYPROJECT + f'repo-shape = "{reason}"\n')
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert f"repo-shape: 0 problem(s) suppressed — {reason}" in flat(proc.stdout)


def test_a_waiver_naming_no_rule_says_so(repo: Path) -> None:
    write(repo, "pyproject.toml", PYPROJECT + 'no-such-rule = "a rule that was renamed"\n')
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "no-such-rule: matches no rule, so it waives nothing" in flat(proc.stdout)


def test_every_failure_is_reported_and_not_just_the_first(repo: Path) -> None:
    # `check.sh` collects every step and conformance collects every rule, for the same reason:
    # one run has to tell you everything to fix.
    write(repo, "docs/agents/writing.md", "# Writing rules\n\nBe concise.\n")
    write(repo, "docs/SKILL.md", "# A skill in the wrong place\n")
    write(repo, "src/example/__init__.py", "")
    proc = conformance(repo)
    assert proc.returncode == 1
    verdicts = statuses(proc)
    assert verdicts["agent-docs-unpublished"] == "FAIL"
    assert verdicts["skill-file-location"] == "FAIL"
    assert verdicts["repo-shape"] == "FAIL"
    assert "3 of 13 rules failed" in proc.stderr


def test_rule_10_fires_on_a_pre_commit_configuration(repo: Path) -> None:
    write(repo, ".pre-commit-config.yaml", "repos:\n  - repo: local\n    hooks: []\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["single-toolchain"] == "FAIL"
    assert ".pre-commit-config.yaml declares dependencies for pre-commit" in proc.stderr
    assert "pyproject.toml is the single source of truth" in flat(proc.stderr)


def test_rule_10_fires_on_a_requirements_file(repo: Path) -> None:
    # Not `requirements.txt`: the marker covers the family, and a rule written for the one
    # spelling one repo happened to use would pass this tree.
    write(repo, "requirements-dev.txt", "ruff==0.6.0\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["single-toolchain"] == "FAIL"
    assert "requirements-dev.txt declares dependencies for pip" in proc.stderr


def test_rule_10_fires_on_a_conda_environment_file(repo: Path) -> None:
    write(repo, "environment.yml", "name: example\ndependencies:\n  - python=3.13\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["single-toolchain"] == "FAIL"
    assert "environment.yml declares dependencies for conda" in proc.stderr


def test_rule_10_fires_on_a_second_build_backends_dependency_section(repo: Path) -> None:
    # No file of its own: the versions are a table in the same manifest, which is the shape a
    # rule that only looked at filenames would miss.
    write(repo, "pyproject.toml", PYPROJECT + '\n[tool.poetry.dependencies]\nrequests = "*"\n')
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["single-toolchain"] == "FAIL"
    assert "pyproject.toml has [tool.poetry], which configures poetry" in proc.stderr


def test_rule_10_fires_on_another_resolvers_lockfile(repo: Path) -> None:
    write(repo, "uv.lock", "version = 1\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["single-toolchain"] == "FAIL"
    assert "uv.lock declares dependencies for uv" in proc.stderr


def test_rule_10_leaves_the_pixi_lock_and_a_fixture_alone(repo: Path) -> None:
    # `pixi.lock` is the lock this rule protects, and a requirements file BELOW the root is some
    # test's input, not this repo's dependencies. Markers are named per tool and root-anchored so
    # that neither of these is a failure.
    write(repo, "pixi.lock", "version: 6\n")
    write(repo, "tests/fixtures/requirements.txt", "requests==2.0.0\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["single-toolchain"] == "ok"


def test_rule_10_says_which_second_toolchains_it_looked_for(repo: Path) -> None:
    # This rule can never be vacuous, so it never prints `--`. It still has to say what it knew
    # to look for, or its green means only that nothing on an unstated list was found.
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["single-toolchain"] == "ok"
    assert "8 second toolchain(s) looked for" in flat(proc.stdout)
    assert "pre-commit, pip, conda, pipenv, setuptools, poetry, uv, pdm" in flat(proc.stdout)


def pyproject_python(
    *,
    floor: str | None = None,
    pin: str | None = None,
    feature_pin: str | None = None,
    ruff: str | None = None,
    pyright: str | None = None,
) -> str:
    """The fixture's pyproject plus whichever Python declarations one test wants.

    A builder rather than five literals because rule 11 is about the COMBINATION: every test
    below differs only in which sites it fills in and what each of them says.
    """
    text = PYPROJECT
    if floor is not None:
        text = text.replace('name = "example"', f'name = "example"\nrequires-python = "{floor}"')
    if pin is not None:
        text += f'\n[tool.pixi.dependencies]\npython = "{pin}"\n'
    if feature_pin is not None:
        text += f'\n[tool.pixi.feature.floor.dependencies]\npython = "{feature_pin}"\n'
    if ruff is not None:
        text += f'\n[tool.ruff]\ntarget-version = "{ruff}"\n'
    if pyright is not None:
        text += f'\n[tool.pyright]\npythonVersion = "{pyright}"\n'
    return text


def test_rule_11_is_vacuous_when_the_repo_pins_no_interpreter(repo: Path) -> None:
    # The fixture declares no Python at all. With no pinned interpreter there is nothing for the
    # other declarations to agree WITH, so the run has to say it checked nothing.
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["python-version-agreement"] == "--"
    assert "not checked: no `python` in [tool.pixi.dependencies] or in any pixi feature" in flat(
        proc.stdout
    )


def test_rule_12_fires_on_a_step_that_runs_a_bare_command(repo: Path) -> None:
    runs(repo, "pytest -q")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    assert (
        ".github/workflows/ci.yml job `check` step 1 runs a command rather than `pixi run <task>`"
        in proc.stderr
    )
    assert "declare what it runs as a task in [tool.pixi.tasks]" in flat(proc.stderr)


def test_rule_12_passes_a_step_that_invokes_a_declared_task(repo: Path) -> None:
    # Two spellings, because the template writes both: a plain task, and one reached through an
    # environment flag. A rule that read `pixi run <word>` and stopped would fail the second.
    runs(repo, "pixi run build", "pixi run -e test test")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "ok"
    # The anti-vacuity assertion. A rule handed no steps passes every workflow ever written, so
    # the green above means something only once the run says how many steps it read.
    assert "2 step(s) that run a command, across 1 workflow(s), against 2 declared task(s)" in flat(
        proc.stdout
    )


def test_rule_11_passes_on_a_floor_and_a_pin_that_agree(repo: Path) -> None:
    # Two declarations, which is what the template ships: ruff takes its target from the floor
    # and pyright takes its version from the interpreter, so neither writes a level down.
    write(repo, "pyproject.toml", pyproject_python(floor=">=3.13", pin="3.13.*"))
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["python-version-agreement"] == "ok"
    assert "2 declaration(s)" in flat(proc.stdout)


def test_rule_11_passes_when_four_declarations_agree(repo: Path) -> None:
    # Writing the derived levels down is allowed; disagreeing with the floor is not. A rule that
    # counted declarations would fail this tree, and this tree is fine.
    write(
        repo,
        "pyproject.toml",
        pyproject_python(floor=">=3.13", pin="3.13.*", ruff="py313", pyright="3.13"),
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["python-version-agreement"] == "ok"
    assert "4 declaration(s)" in flat(proc.stdout)


def test_rule_11_fires_when_a_tool_holds_a_level_the_floor_does_not(repo: Path) -> None:
    write(
        repo,
        "pyproject.toml",
        pyproject_python(floor=">=3.13", pin="3.13.*", ruff="py312", pyright="3.13"),
    )
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["python-version-agreement"] == "FAIL"
    assert (
        "[tool.ruff] target-version `py312` says 3.12, and [project] requires-python `>=3.13` "
        "says 3.13" in proc.stderr
    )
    # Named the site, what it says, and what it should say instead.
    assert 'write `target-version = "py313"`' in flat(proc.stderr)


def test_rule_11_fires_on_a_pin_below_the_floor(repo: Path) -> None:
    # The other direction: the repo runs an interpreter its own metadata refuses to install on.
    write(repo, "pyproject.toml", pyproject_python(floor=">=3.13", pin="3.12.*"))
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["python-version-agreement"] == "FAIL"
    assert (
        "[tool.pixi.dependencies] python `3.12.*` pins 3.12, below the [project] requires-python "
        "floor `>=3.13`" in proc.stderr
    )


def test_rule_11_measures_a_tool_against_the_pin_when_there_is_no_floor(repo: Path) -> None:
    write(repo, "pyproject.toml", pyproject_python(pin="3.13.*", ruff="py312"))
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["python-version-agreement"] == "FAIL"
    assert "[tool.pixi.dependencies] python `3.13.*` says 3.13" in proc.stderr


def test_rule_11_passes_a_range_that_pixi_actually_resolves(repo: Path) -> None:
    # The case the rule exists to leave alone: a repo supporting 3.10 through 3.13 on purpose. The
    # floor is below the default pin, the tools are held to the floor, and a feature gives 3.10 an
    # environment — so the range is a claim pixi resolves rather than one nothing exercises.
    write(
        repo,
        "pyproject.toml",
        pyproject_python(
            floor=">=3.10", pin="3.13.*", feature_pin="3.10.*", ruff="py310", pyright="3.10"
        ),
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["python-version-agreement"] == "ok"
    assert "the floor 3.10 is itself pinned, so pixi resolves it" in flat(proc.stdout)


def test_rule_11_says_it_cannot_tell_a_range_from_a_floor_nobody_lowered(repo: Path) -> None:
    # Copied from a real Liu Lab repo, which carried exactly this: a >=3.12 floor, a 3.13 pin,
    # ruff and pyright at 3.12, and a lockfile holding only 3.13. The declarations AGREE — the
    # tools are at the floor, which is what a deliberate range looks like — so the rule does not
    # fail it. What it must not do is call that a clean pass in silence: nothing pins 3.12, so the
    # promise of 3.12 support is exercised by nothing, and the notes say so on every run.
    write(
        repo,
        "pyproject.toml",
        pyproject_python(floor=">=3.12", pin="3.13.*", ruff="py312", pyright="3.12"),
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["python-version-agreement"] == "ok"
    printed = flat(proc.stdout)
    # Every site and what each one says, so the split is readable at a glance.
    for said in (
        "[project] requires-python `>=3.12`",
        "[tool.pixi.dependencies] python `3.13.*`",
        "[tool.ruff] target-version `py312`",
        "[tool.pyright] pythonVersion `3.12`",
    ):
        assert said in printed
    assert "not checked: whether 3.12 is ever resolved or tested" in printed


def test_rule_12_reads_past_the_environment_flag_to_the_task(repo: Path) -> None:
    # `test` is BOTH an environment and a declared task here, which is the trap: a rule that
    # accepted any declared name among the tokens would call this green and let the inlined
    # `pytest` through — the exact drift the rule exists to catch, hidden behind a `-e`.
    runs(repo, "pixi run -e test pytest")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    assert "step 1 invokes `pytest`, which this repo declares no task by" in proc.stderr


def test_rule_12_fires_on_a_task_nobody_declared(repo: Path) -> None:
    # Not an inlined command: it follows the convention exactly and names a task that does not
    # exist. Nothing runs it until the workflow does, and this one is `ci.yml`; on `release.yml`
    # the first run is the release.
    runs(repo, "pixi run typecheck")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    assert "invokes `typecheck`, which this repo declares no task by" in proc.stderr
    assert "declare `typecheck` in pyproject.toml under [tool.pixi.tasks]" in flat(proc.stderr)


def test_rule_12_finds_a_task_declared_on_a_feature(repo: Path) -> None:
    # `docs-build` is on the `docs` FEATURE in the real manifest, not on the workspace. A rule
    # reading `[tool.pixi.tasks]` alone would fail the docs job of a workflow that is correct.
    write(
        repo, "pyproject.toml", PYPROJECT + '\n[tool.pixi.feature.docs.tasks]\ndocs-build = "z"\n'
    )
    runs(repo, "pixi run docs-build")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "ok"
    assert "against 3 declared task(s)" in flat(proc.stdout)


def test_rule_12_does_not_flag_a_step_that_uses_an_action(repo: Path) -> None:
    # Three `uses:` steps, one of them carrying a command line in its `with:` block. An action is
    # a dependency the workflow pulls in, not a step of this repo's build, so none of them is a
    # task that went missing — and the count in the note is what proves they were passed over
    # rather than merely absent.
    write(
        repo,
        ".github/workflows/ci.yml",
        "name: ci\non:\n  pull_request:\njobs:\n"
        "  check:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v7\n"
        "      - uses: prefix-dev/setup-pixi@v0.10.2\n"
        "        with:\n"
        "          environments: default\n"
        "      - run: pixi run build\n"
        "      - uses: actions/upload-artifact@v7\n"
        "        with:\n"
        "          args: rm -rf dist\n",
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "ok"
    assert "1 step(s) that run a command" in flat(proc.stdout)


def test_rule_12_fires_on_a_command_chained_onto_a_task(repo: Path) -> None:
    # It does invoke a declared task. What follows the task is the same drift one level down: a
    # laptop running `pixi run build` never deletes anything.
    runs(repo, "pixi run build && rm -rf dist")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    assert "step 1 runs a command rather than `pixi run <task>`" in proc.stderr


def test_rule_12_fires_on_a_multi_line_script(repo: Path) -> None:
    write(
        repo,
        ".github/workflows/ci.yml",
        "name: ci\non:\n  pull_request:\njobs:\n"
        "  check:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - name: Everything at once\n"
        "        run: |\n"
        "          pixi run build\n"
        "          pixi run -e test test\n",
    )
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    # Named, because this step has a `name:` and a reader should not have to count steps.
    assert "step 0 (Everything at once) runs a command" in proc.stderr


def test_rule_12_reports_a_task_named_by_an_expression_as_unresolved(repo: Path) -> None:
    # A matrix over tasks. GitHub substitutes the expression when the job runs, so the second step
    # names a task this cannot look up — reported, not failed, since a rule that cannot evaluate
    # its premise has checked nothing. The first step shows an expression somewhere OTHER than the
    # task name is resolved past and not counted.
    runs(repo, "pixi run -e ${{ matrix.environment }} test", "pixi run ${{ matrix.task }}")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "ok"
    assert "2 step(s) that run a command" in flat(proc.stdout)
    assert "1 step(s) name their task with a `${{}}` expression" in flat(proc.stdout)


def test_rule_12_is_vacuous_and_not_passing_when_no_workflow_runs_a_command(repo: Path) -> None:
    # A repo whose workflows are all actions, or that has none at all. Nothing was checked, and a
    # green here would read as "every step invokes a task", which this run has no evidence for.
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "--"
    assert "not checked: no tracked workflow has a step that runs a command" in flat(proc.stdout)


def test_the_workflow_reader_flattens_every_job_and_not_just_the_first(repo: Path) -> None:
    # The step the rule must find is the LAST step of the SECOND job. A reader that stopped at the
    # first job, or at a job's first step, would call this tree green.
    write(
        repo,
        ".github/workflows/ci.yml",
        "name: ci\non:\n  pull_request:\njobs:\n"
        "  check:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v7\n"
        "      - run: pixi run build\n"
        "  test:\n"
        "    runs-on: ubuntu-latest\n"
        "    steps:\n"
        "      - uses: actions/checkout@v7\n"
        "      - run: pixi run -e test test\n"
        "      - run: coverage report\n",
    )
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["workflow-step-tasks"] == "FAIL"
    assert ".github/workflows/ci.yml job `test` step 2 runs a command" in proc.stderr
    # One problem and not three: the other two steps of the two jobs were read and were fine.
    assert "workflow-step-tasks FAIL 1 problem(s)" in flat(proc.stdout)


def test_a_workflow_that_will_not_parse_is_read_as_empty_and_not_as_an_exception(
    repo: Path,
) -> None:
    # A derived repo's unusual workflow must not turn a conformance rule into a crash. The file
    # counts as a workflow and contributes no steps, so the rules that read workflows keep
    # reporting on the ones they could read.
    write(repo, ".github/workflows/broken.yml", "name: broken\non: [\njobs: }{\n")
    runs(repo, "pixi run build")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Traceback" not in proc.stderr
    assert statuses(proc)["workflow-step-tasks"] == "ok"
    assert "1 step(s) that run a command, across 2 workflow(s)" in flat(proc.stdout)


def test_rule_13_fires_on_a_nav_entry_naming_a_file_that_is_not_there(repo: Path) -> None:
    # The hole this rule exists for. Measured on the pinned site generator: it does not validate
    # the nav at all, so `--strict` reports "No issues found", exits 0, and publishes a menu item
    # whose href points at a page that was never rendered.
    write(repo, "mkdocs.yml", MKDOCS + "  - Ghost: ghost.md\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["nav-target-exists"] == "FAIL"
    assert (
        "mkdocs.yml nav entry `ghost.md` names docs/ghost.md, which is not tracked" in proc.stderr
    )
    assert "delete it from the `nav:` list in mkdocs.yml" in flat(proc.stderr)


def test_rule_13_passes_a_nav_whose_pages_are_all_there(repo: Path) -> None:
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["nav-target-exists"] == "ok"
    # The anti-vacuity assertion: a rule handed no entries passes every nav ever written, so the
    # green above means something only once the run says how many entries it resolved.
    assert "mkdocs.yml: 2 nav entry(s) against docs/" in flat(proc.stdout)


def test_rule_13_is_vacuous_and_not_passing_when_the_repo_publishes_no_site(repo: Path) -> None:
    # One lab repo ships no docs site at all. Nothing was checked, and a green here would read as
    # "every menu item resolves", which this run has no evidence for.
    (repo / "mkdocs.yml").unlink()
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["nav-target-exists"] == "--"
    assert "not checked: no tracked mkdocs.yml, so this repo publishes no site" in flat(proc.stdout)


def test_rule_13_is_vacuous_when_the_site_declares_no_nav(repo: Path) -> None:
    # No `nav:` is not an empty nav. The generator builds the navbar from the directory tree
    # instead, so there is no entry to resolve and nothing was checked.
    write(repo, "mkdocs.yml", "site_name: example\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["nav-target-exists"] == "--"
    assert "not checked: mkdocs.yml declares no `nav:`" in flat(proc.stdout)


def test_rule_13_does_not_flag_a_link_or_an_anchor(repo: Path) -> None:
    # A nav entry may point out of the site entirely. None of these three names a file, so none
    # of them is a file that can be missing — and the count proves the two real pages were still
    # the thing checked.
    write(
        repo,
        "mkdocs.yml",
        MKDOCS
        + "  - Issues: https://github.com/liuhlab/example/issues\n"
        + "  - Mail: mailto:lab@example.org\n"
        + "  - Top: '#top'\n",
    )
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["nav-target-exists"] == "ok"
    assert "mkdocs.yml: 2 nav entry(s) against docs/" in flat(proc.stdout)


def test_rule_13_fires_on_a_page_that_exists_but_sits_outside_the_site_source(repo: Path) -> None:
    # The second class, and the reason it is the same rule rather than a separate concern:
    # `README.md` really is there, so a check that only asked whether the file exists would call
    # this fine. The builder renders the site source and nothing else, so the menu item 404s
    # exactly as a missing page does — and the fix is a different one, so it says so.
    write(repo, "README.md", "# example\n\nThe repo.\n")
    write(repo, "mkdocs.yml", MKDOCS + "  - Readme: ../README.md\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["nav-target-exists"] == "FAIL"
    assert "nav entry `../README.md` resolves to README.md, outside docs/" in proc.stderr
    assert "move the page under docs/" in flat(proc.stderr)


def test_rule_13_fires_on_a_page_that_is_present_on_disk_but_untracked(repo: Path) -> None:
    # It builds on the laptop that wrote it and nowhere else, because CI builds from a checkout.
    # The problem says which of the two it is, since "it is right there" is what the person
    # reading the failure is about to say.
    write(repo, ".gitignore", "docs/draft.md\n")
    write(repo, "docs/draft.md", "# Draft\n\nNever committed.\n")
    write(repo, "mkdocs.yml", MKDOCS + "  - Draft: draft.md\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["nav-target-exists"] == "FAIL"
    assert "names docs/draft.md, which is not tracked" in proc.stderr
    assert "present on disk but untracked, and CI builds from a checkout" in flat(proc.stderr)


def test_rule_13_reads_a_second_configuration_that_inherits_the_first(repo: Path) -> None:
    # The shape this template itself has: a second config INHERITs `mkdocs.yml`, and its `nav:`
    # APPENDS to the inherited one rather than replacing it, so both files' entries are in the
    # navbar the build renders. A rule reading only `mkdocs.yml` would call this tree green.
    write(repo, "mkdocs.other.yml", "INHERIT: mkdocs.yml\nnav:\n  - Extra: extra.md\n")
    proc = conformance(repo)
    assert proc.returncode == 1
    assert statuses(proc)["nav-target-exists"] == "FAIL"
    assert (
        "mkdocs.other.yml nav entry `extra.md` names docs/extra.md, which is not tracked"
        in proc.stderr
    )


def test_rule_13_resolves_an_inherited_docs_dir_and_not_the_default(repo: Path) -> None:
    # The second configuration names no `docs_dir` of its own, so its nav paths resolve against
    # the one it inherits. A rule that assumed `docs/` here would report a page that is there —
    # a false failure on the exact two-file shape this template ships.
    write(
        repo, "mkdocs.yml", "site_name: example\ndocs_dir: site-source\nnav:\n  - Home: index.md\n"
    )
    write(repo, "site-source/index.md", "# example\n\nThe front page.\n")
    write(repo, "site-source/extra.md", "# Extra\n\nA second page.\n")
    write(repo, "mkdocs.other.yml", "INHERIT: mkdocs.yml\nnav:\n  - Extra: extra.md\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert statuses(proc)["nav-target-exists"] == "ok"
    assert "mkdocs.other.yml: 1 nav entry(s) against site-source/" in flat(proc.stdout)


def test_rule_13_reads_a_site_config_that_will_not_parse_as_empty(repo: Path) -> None:
    # A derived repo's unusual site config must not turn a conformance rule into a crash. It
    # contributes no nav, and the configuration that could be read is still reported on.
    write(repo, "mkdocs.broken.yml", "nav: [\n  - x: }{\n")
    proc = conformance(repo)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Traceback" not in proc.stderr
    assert statuses(proc)["nav-target-exists"] == "ok"
    assert "not checked: mkdocs.broken.yml declares no `nav:`" in flat(proc.stdout)
