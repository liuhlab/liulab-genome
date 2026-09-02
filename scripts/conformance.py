#!/usr/bin/env python3
"""The conformance check: the rules a Liu Lab repo keeps after it leaves the template.

    pixi run conformance
    python scripts/conformance.py [--root PATH]

It states the RULE, not the file contents — a pull model, so a repo that legitimately diverges
stays green while the shared conventions stay checked. Ten rules: eight fail, two only warn.

Three things it does on purpose:

- **It collects every failure.** One run tells you everything to fix, the same idiom
  `scripts/check.sh` uses. Nothing short-circuits.
- **It names the rule and a fix.** You should never have to read this file to act on its output,
  so nothing here prints a rewritten assert.
- **It says what it did NOT check.** A rule whose premise is absent is vacuous, not failing — a
  repo with no `src/` is not a repo that failed to have tests — but a vacuous rule looks exactly
  like a passing one, so each is reported and says why.

The declared agent-facing set is `[tool.liulab.agent-docs]` in `pyproject.toml`: keys are plain
path prefixes matched against `git ls-files`, values name the Vale Length rule that applies or
`false` where none does. Rules 2 and 3 both read it, so the two consumers cannot drift.

Waivers are `[tool.liulab.waived]`, keyed by RULE NAME with the reason as the value. Per-rule,
never per-file. No expiry, and every waiver is printed on every run — green ones included —
because an escape hatch that stops being visible becomes permanent. Rule 1 refuses to be waived.

This lives in `scripts/` and not `tests/` because one lab repo has neither `src/` nor `tests/`,
and pytest is in no default environment. Its negative tests are `tests/test_conformance.py`: one
tree per rule that violates it, because a rule that checks nothing passes the happy path and a
gate that fires on nothing looks identical to a gate that passes.
"""

from __future__ import annotations

import argparse
import posixpath
import re
import subprocess
import sys
import textwrap
import tomllib
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

#: The template's placeholder module name, SPELLED IN TWO PIECES on purpose. `init-repo`
#: substitutes that name throughout the tree and the dogfood job then greps a rendered repo for
#: it. A literal here would be rewritten — turning rule 1 into a check for the new module's own
#: name, which is everywhere — and would leave a copy for that grep to find. Two pieces, and
#: neither happens.
PLACEHOLDER = "new" + "pkg"

#: Rule 1 and warning 9 are both gated on this directory, and it is the only discriminator either
#: needs. While it is here, `init-repo` has not run: the placeholder is still the repo's own name,
#: and the auto-discovered skill is itself the nag. Once it is gone the repo claims to be its own,
#: and both rules start checking.
INIT_SKILL = "skills/init-repo"

#: Rule 1's name, needed by name because it is the one rule that refuses to be waived.
UNWAIVABLE = "placeholder-rename"

#: Where `mkdocs.yml` puts the site source when it does not say. Rule 2 is scoped to the pages the
#: site publishes, so it reads this rather than assuming.
DEFAULT_DOCS_DIR = "docs"

#: Vale's spelling of on and off. Both forms appear in the wild; the gate accepts either.
VALE_ON = {"YES", "TRUE", "ON"}
VALE_OFF = {"NO", "FALSE", "OFF"}

_SECTION_RE = re.compile(r"^\[(?P<header>.+)\]\s*$")
_SETTING_RE = re.compile(r"^(?P<key>[A-Za-z][^=]*?)\s*=\s*(?P<value>.*?)\s*$")
_SKILL_FILE_RE = re.compile(r"^skills/[^/]+/SKILL\.md$")
#: A skill directory, at the source of truth and at both committed discovery paths. `skills/`
#: needs the trailing slash so `skills/install.py` is not read as a skill called `install.py`;
#: the discovery paths do not, because `git ls-files` reports a symlink to a directory as a plain
#: path with nothing after it — and a symlink is exactly how a name gets shadowed.
_SKILL_DIR_RES = (
    re.compile(r"^skills/(?P<name>[^/]+)/"),
    re.compile(r"^\.(?:claude|agents)/skills/(?P<name>[^/]+)(?:/|$)"),
)


class _TolerantLoader(yaml.SafeLoader):
    """A SafeLoader that reads a config it does not fully understand.

    mkdocs configurations carry local tags — `!ENV`, `!!python/name:...` for an emoji index — and
    a plain `safe_load` raises on the first one. Rule 2 only wants `nav:` and `docs_dir:`, so an
    unknown tag becomes ``None`` rather than an exception: a derived repo that adds one must not
    turn a conformance rule into a crash.
    """


# `None` is PyYAML's catch-all prefix: it is consulted only for a tag no constructor claimed.
_TolerantLoader.add_multi_constructor(None, lambda _loader, _suffix, _node: None)


def _problem(what: str, why: str, fix: str) -> str:
    """Format one problem for whoever has to act on it.

    Three lines and always the same three — WHAT is wrong, WHY it matters, and how to repair it —
    which is the shape `skills/install.py --check` already prints, so the two commands an agent
    meets in this repo read the same way. The what line is left unwrapped so it stays one
    greppable line.
    """
    body = textwrap.fill(why, width=94, initial_indent="    ", subsequent_indent="    ")
    return f"{what}\n{body}\n    fix: {fix}"


@dataclass
class Result:
    """What one rule found when it ran.

    Attributes
    ----------
    problems
        One entry per violation, formatted by :func:`_problem`. Empty means the rule held.
    notes
        What the rule examined, or — when its premise was absent — why it examined nothing. A
        vacuous rule prints identically to a passing one unless it says so itself.
    """

    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


@dataclass
class Repo:
    """The tree under inspection and the declarations it makes about itself."""

    root: Path
    tracked: tuple[str, ...]
    agent_docs: dict[str, Any]
    waived: dict[str, str]

    def exists(self, rel: str) -> bool:
        """Whether a path is present on disk, tracked or not."""
        return (self.root / rel).exists()

    def read(self, rel: str) -> str | None:
        """Read a file, or return ``None`` when it is not there."""
        path = self.root / rel
        if not path.is_file():
            return None
        return path.read_text(encoding="utf-8", errors="replace")

    def under(self, prefix: str) -> list[str]:
        """Every tracked path under one `[tool.liulab.agent-docs]` key.

        Keys are plain path prefixes, deliberately not globs: reimplementing Vale's glob matcher
        would reproduce two known surprises. A key ending in `/` matches everything beneath that
        directory; any other key matches that one path exactly, so `AGENTS.md` cannot also claim
        `AGENTS.md.bak`.
        """
        if prefix.endswith("/"):
            return [p for p in self.tracked if p.startswith(prefix)]
        return [p for p in self.tracked if p == prefix]

    def tracked_under(self, prefix: str) -> bool:
        """Whether anything tracked lives under a directory prefix."""
        return any(p.startswith(prefix) for p in self.tracked)

    @cached_property
    def mkdocs(self) -> dict[str, Any]:
        """`mkdocs.yml`, or an empty mapping when the repo publishes no site."""
        text = self.read("mkdocs.yml")
        if text is None:
            return {}
        # Not `safe_load`: `_TolerantLoader` is SafeLoader plus a shrug at unknown tags.
        loaded = yaml.load(text, Loader=_TolerantLoader)
        return loaded if isinstance(loaded, dict) else {}

    @cached_property
    def docs_dir(self) -> str:
        """The site source directory, repo-relative and without a trailing slash."""
        declared = self.mkdocs.get("docs_dir", DEFAULT_DOCS_DIR)
        return posixpath.normpath(str(declared)).strip("/")

    @cached_property
    def nav(self) -> dict[str, str]:
        """Repo-relative path -> the `nav:` entry that names it.

        Nav paths are relative to `docs_dir`, so they are prefixed with it before they can be
        compared with anything `git ls-files` said. Getting that wrong is not a false negative
        that shows up later — it makes rule 2 pass vacuously, which is the failure this whole file
        exists to prevent. External links are skipped: a nav entry may be a URL.
        """
        entries: dict[str, str] = {}
        for raw in _walk_strings(self.mkdocs.get("nav")):
            if "://" in raw or raw.startswith("#"):
                continue
            entries[posixpath.normpath(f"{self.docs_dir}/{raw}")] = raw
        return entries


def _walk_strings(node: Any) -> Iterator[str]:
    """Every string value in a `nav:` tree, at any depth.

    The template's nav is a flat list of one-key mappings, but a derived repo will nest sections,
    and a nested entry that escaped this walk would be an agent-facing page silently allowed into
    the navbar.
    """
    if isinstance(node, str):
        yield node
    elif isinstance(node, list):
        for item in node:  # pyright: ignore[reportUnknownVariableType]
            yield from _walk_strings(item)
    elif isinstance(node, dict):
        for value in node.values():  # pyright: ignore[reportUnknownVariableType]
            yield from _walk_strings(value)


def _front_matter(text: str) -> dict[str, Any] | None:
    """Parse a markdown file's YAML front matter, or return ``None`` when it has none."""
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        return None
    for index, line in enumerate(lines[1:], start=1):
        if line.strip() in {"---", "..."}:
            loaded = yaml.safe_load("\n".join(lines[1:index]))
            return loaded if isinstance(loaded, dict) else {}
    return None


def _excluded_from_search(text: str) -> bool:
    """Whether a page carries `search: exclude: true`."""
    matter = _front_matter(text)
    if matter is None:
        return False
    search = matter.get("search")
    return isinstance(search, dict) and search.get("exclude") is True


def _vale_sections(text: str) -> list[tuple[str, dict[str, str]]]:
    """`.vale.ini` as (section header, settings) pairs, in file order.

    Hand-parsed rather than handed to `configparser`, which lowercases keys and rejects a
    duplicate header — neither of which is what this file means.
    """
    sections: list[tuple[str, dict[str, str]]] = []
    current: dict[str, str] | None = None
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", ";")):
            continue
        header = _SECTION_RE.match(stripped)
        if header:
            current = {}
            sections.append((header["header"], current))
            continue
        setting = _SETTING_RE.match(stripped)
        if setting and current is not None:
            current[setting["key"]] = setting["value"]
    return sections


def _glossary_entries(text: str) -> list[tuple[str, int]]:
    """Each `###` entry under `## Glossary` in `CONTEXT.md`, with its word count.

    The cap is per ENTRY and not per file, so a repo with a large vocabulary is not punished for
    the size of its domain. An entry runs from its own heading to the next heading at the same or
    a higher level.
    """
    entries: list[tuple[str, int]] = []
    in_glossary = False
    title: str | None = None
    words = 0
    for line in text.splitlines():
        if line.startswith("## "):
            if title is not None:
                entries.append((title, words))
                title = None
            in_glossary = line[3:].strip().lower() == "glossary"
            continue
        if not in_glossary:
            continue
        if line.startswith("### "):
            if title is not None:
                entries.append((title, words))
            title = line[4:].strip()
            words = len(title.split())
            continue
        if title is not None:
            words += len(line.split())
    if title is not None:
        entries.append((title, words))
    return entries


# --------------------------------------------------------------------------------------
# The rules. Eight fail, two warn. Each returns what it found; nothing here exits or prints.
# --------------------------------------------------------------------------------------


def rule_placeholder_rename(repo: Repo) -> Result:
    """1. No tracked file contains the placeholder name."""
    result = Result()
    # The refusal is reported even where the rule is vacuous: a waiver that silently does nothing
    # is the invisible escape hatch this design exists to prevent.
    if UNWAIVABLE in repo.waived:
        result.problems.append(
            _problem(
                f"[tool.liulab.waived] names `{UNWAIVABLE}`, which cannot be waived",
                "a half-renamed repo publishes package metadata and canonical URLs naming a "
                "package that does not exist, and one grep is what proves the rename complete. "
                "The waiver changes nothing; the rule ran anyway",
                f"delete the `{UNWAIVABLE}` line from [tool.liulab.waived]",
            )
        )
    if repo.exists(INIT_SKILL):
        result.notes.append(
            f"not checked: {INIT_SKILL}/ is here, so `init-repo` has not run and "
            f"`{PLACEHOLDER}` is still this repo's own name"
        )
        return result
    for rel in repo.tracked:
        path = repo.root / rel
        if PLACEHOLDER in rel:
            where = "its path"
        elif path.is_symlink():
            where = "its link target" if PLACEHOLDER in str(path.readlink()) else ""
        elif path.is_file() and PLACEHOLDER.encode() in path.read_bytes():
            where = "its contents"
        else:
            where = ""
        if where:
            result.problems.append(
                _problem(
                    f"{rel} still names `{PLACEHOLDER}` in {where}",
                    "one grep proves the rename complete. A leftover copy means an import, a "
                    "distribution name or a site_url still points at a package that does not "
                    "exist, and a redirect can hide that for months",
                    "run `init-repo`, or rename it by hand and re-run this check",
                )
            )
    result.notes.append(f"{len(repo.tracked)} tracked path(s) read for `{PLACEHOLDER}`")
    return result


def rule_agent_docs_unpublished(repo: Repo) -> Result:
    """2. Agent-facing pages the site would publish stay out of the navbar and out of search."""
    result = Result()
    prefix = f"{repo.docs_dir}/"
    # Scoped to the agent-docs keys UNDER THE SITE SOURCE. `AGENTS.md`, `CONTEXT.md` and
    # `skills/` are agent-facing too, but the site never publishes them: front matter there means
    # nothing, and on a SKILL.md it lands in the frontmatter a skill loader reads.
    keys = sorted(k for k in repo.agent_docs if k.startswith(prefix))
    if not keys:
        result.notes.append(f"not checked: no [tool.liulab.agent-docs] key is under {prefix}")
        return result
    pages = sorted({p for key in keys for p in repo.under(key) if p.endswith(".md")})
    for page in pages:
        text = repo.read(page)
        if text is not None and not _excluded_from_search(text):
            result.problems.append(
                _problem(
                    f"{page} has no `search: exclude: true` front matter",
                    "an unlisted page is still indexed, so it turns up in the site's search box "
                    "for a human who was looking for the documentation",
                    "add this to the top of the file:\n"
                    "               ---\n"
                    "               search:\n"
                    "                 exclude: true\n"
                    "               ---",
                )
            )
        if page in repo.nav:
            result.problems.append(
                _problem(
                    f"{page} appears in mkdocs.yml nav as `{repo.nav[page]}`",
                    "nav is the navbar, and this page is written for an agent. Front matter "
                    "does nothing about the navbar and nav does nothing about search, which is "
                    "why both are required",
                    f"delete `{repo.nav[page]}` from the nav: list in mkdocs.yml. The page stays "
                    "reachable by URL",
                )
            )
    result.notes.append(f"{len(pages)} markdown file(s) under {', '.join(keys)}")
    return result


def rule_vale_length_sections(repo: Repo) -> Result:
    """3. `.vale.ini` gives every agent-docs key a section that turns on its declared cap."""
    result = Result()
    text = repo.read(".vale.ini")
    if text is None:
        result.problems.append(
            _problem(
                ".vale.ini is missing",
                "the declared word caps are enforced by Vale and by nothing else, so without "
                "this file every agent-docs key is checked by nothing",
                "restore .vale.ini from the template",
            )
        )
        return result
    sections = _vale_sections(text)
    # The Length rules the repo actually declares, taken from the table rather than hardcoded, so
    # adding a third cap needs no edit here.
    caps = sorted({v for v in repo.agent_docs.values() if isinstance(v, str)})
    for key, declared in sorted(repo.agent_docs.items()):
        matched = [(header, settings) for header, settings in sections if key in header]
        if len(matched) != 1:
            result.problems.append(
                _problem(
                    f"`{key}` is spelled in {len(matched)} .vale.ini section header(s), not 1",
                    "conformance pairs a key with its section by looking for the key inside the "
                    "header, which is exact only while each key appears in exactly one. A key in "
                    "none is checked by nothing; a key in two cannot be paired",
                    f"give `{key}` one section header that spells it verbatim, as the comment at "
                    "the top of .vale.ini asks",
                )
            )
            continue
        header, settings = matched[0]
        for cap in caps:
            setting = f"Lab.{cap}"
            wanted = VALE_ON if cap == declared else VALE_OFF
            actual = settings.get(setting, "").upper()
            if actual in wanted:
                continue
            state = "on" if cap == declared else "off"
            saw = f"`{settings[setting]}`" if setting in settings else "nothing"
            result.problems.append(
                _problem(
                    f"[{header}] should set {setting} {state} for `{key}`, and says {saw}",
                    "Vale's `*` matches `/`, so more than one section matches the same file and "
                    "per-rule settings ACCUMULATE, later sections winning. An omitted rule is "
                    "not default-on: it keeps whatever an earlier section said, which can leave "
                    "these files checked by nothing, with zero alerts and green CI",
                    f"write `{setting} = {'YES' if state == 'on' else 'NO'}` in that section. "
                    "Every rule is named in every section for this reason",
                )
            )
    result.notes.append(
        f"{len(repo.agent_docs)} agent-docs key(s) against {len(sections)} .vale.ini section(s)"
    )
    return result


def rule_glossary_entry_length(repo: Repo) -> Result:
    """4. No glossary entry in `CONTEXT.md` runs past 200 words."""
    cap = 200
    result = Result()
    text = repo.read("CONTEXT.md")
    if text is None:
        result.notes.append("not checked: no CONTEXT.md")
        return result
    entries = _glossary_entries(text)
    if not entries:
        result.notes.append("not checked: CONTEXT.md has no entries under `## Glossary`")
        return result
    for title, words in entries:
        if words > cap:
            result.problems.append(
                _problem(
                    f"CONTEXT.md glossary entry `{title}` is {words} words, over {cap}",
                    "the glossary says what a word MEANS in this repo. An entry that long is an "
                    "explanation, and it belongs in a document the glossary can point at",
                    "cut it to a definition and move the rest to docs/, or to an ADR if it is a "
                    "decision",
                )
            )
    longest = max(words for _, words in entries)
    result.notes.append(f"{len(entries)} glossary entry(s), longest {longest} of {cap} words")
    return result


def rule_skill_file_location(repo: Repo) -> Result:
    """5. Every `SKILL.md` sits at `skills/<name>/SKILL.md`."""
    result = Result()
    found = [p for p in repo.tracked if PurePosixPath(p).name == "SKILL.md"]
    for rel in found:
        if not _SKILL_FILE_RE.match(rel):
            result.problems.append(
                _problem(
                    f"{rel} is a SKILL.md outside skills/<name>/",
                    "the installer discovers skills/*/SKILL.md and nothing else, so a skill "
                    "written anywhere else sits there unnoticed by every agent product",
                    "move it to skills/<name>/SKILL.md, then run "
                    "`python skills/install.py --target all`",
                )
            )
    result.notes.append(f"{len(found)} tracked SKILL.md")
    return result


def rule_skill_name_prefix(repo: Repo) -> Result:
    """6. No repo-local skill directory begins with `lab-`."""
    result = Result()
    names: dict[str, str] = {}
    for rel in repo.tracked:
        for pattern in _SKILL_DIR_RES:
            match = pattern.match(rel)
            if match:
                names.setdefault(match["name"], rel)
                break
    for name, rel in sorted(names.items()):
        if name.startswith("lab-"):
            result.problems.append(
                _problem(
                    f"{rel} is a repo-local skill named `{name}`",
                    "the lab's shared plugin owns the `lab-` prefix, and a repo-local skill of "
                    "the same name shadows it for everyone working in this repo — silently, "
                    "since both are discovered the same way",
                    "rename it to something without the `lab-` prefix, then run "
                    "`python skills/install.py --target all`",
                )
            )
    result.notes.append(f"{len(names)} skill director(ies) named across the discovery paths")
    return result


def rule_skill_symlinks(repo: Repo) -> Result:
    """7. `python skills/install.py --check` passes."""
    result = Result()
    installer = "skills/install.py"
    if not repo.exists(installer):
        result.notes.append(f"not checked: no {installer}")
        return result
    # DELEGATED, not reimplemented: the command that reports a broken link is the command that
    # repairs it, and its own three invariants are not restated here to drift from.
    proc = subprocess.run(
        [sys.executable, installer, "--check"],
        cwd=repo.root,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        detail = (proc.stdout + proc.stderr).strip() or f"exit {proc.returncode}, no output"
        result.problems.append(
            _problem(
                f"`python {installer} --check` failed:",
                "it holds the two committed discovery paths to three invariants — a link in "
                "both, nothing dangling, every link relative. Its own report follows",
                "run the command it names below",
            )
            + "\n"
            + textwrap.indent(detail, "      ")
        )
    else:
        result.notes.append((proc.stdout.strip() or "ok").splitlines()[-1])
    return result


def rule_repo_shape(repo: Repo) -> Result:
    """8. If there is a package there are tests; if there is a release workflow there is a log."""
    result = Result()
    release = ".github/workflows/release.yml"
    # Both halves are CONDITIONAL, and an absent premise is vacuous rather than failing: a repo
    # with no package is not a repo that failed to have tests. That is the whole shape of the
    # rule, and it is why the notes below say which premise was absent.
    if not repo.tracked_under("src/"):
        result.notes.append("not checked: no tracked src/, so tests/ is not required")
    elif not repo.tracked_under("tests/"):
        result.problems.append(
            _problem(
                "src/ is tracked and tests/ is not",
                "a package with no test directory has nothing pytest can find, so the `test` "
                "step of `pixi run check` passes having run nothing",
                "add tests/, or delete src/ if this repo has no package",
            )
        )
    else:
        result.notes.append("src/ is tracked and so is tests/")
    if release not in repo.tracked:
        result.notes.append(f"not checked: no {release}, so CHANGELOG.md is not required")
    elif "CHANGELOG.md" not in repo.tracked:
        result.problems.append(
            _problem(
                f"{release} is tracked and CHANGELOG.md is not",
                "this repo publishes releases, and the people who install it from git have "
                "nowhere to read what changed",
                "add CHANGELOG.md. `init-repo` never subtracts it, so it was deleted by hand",
            )
        )
    else:
        result.notes.append(f"{release} is tracked and so is CHANGELOG.md")
    return result


def warning_init_sentinel(repo: Repo) -> Result:
    """9. `AGENTS.md` still carries the line telling you to run `/init`."""
    result = Result()
    if repo.exists(INIT_SKILL):
        result.notes.append(f"not checked: {INIT_SKILL}/ is here, and that skill is the nag")
        return result
    text = repo.read("AGENTS.md")
    if text is None:
        result.notes.append("not checked: no AGENTS.md")
        return result
    # The sentinel is a blockquote line naming `/init`. It is phrased as an instruction and it
    # says to delete itself, so replacing the file is what makes this warning stop — "done" needs
    # no separate step. This WARNS and never fails: a new repo must not start red for a task its
    # owner is allowed to defer.
    if any(line.lstrip().startswith(">") and "/init" in line for line in text.splitlines()):
        result.problems.append(
            _problem(
                "AGENTS.md still carries the `/init` sentinel",
                "the file is the template's generic copy: it says how ANY Liu Lab repo works, "
                "not what this one does, and an agent that reads it will act on a description "
                "that was never true here",
                "run `/init`, then delete the sentinel line. This is a warning, not a failure",
            )
        )
    else:
        result.notes.append("AGENTS.md no longer carries it")
    return result


@dataclass(frozen=True)
class Rule:
    """One rule, its number, and the one-line statement printed when it fails."""

    number: int
    name: str
    statement: str
    check: Callable[[Repo], Result]
    warns_only: bool = False


RULES: tuple[Rule, ...] = (
    Rule(
        1,
        UNWAIVABLE,
        f"no tracked file contains `{PLACEHOLDER}` — the one grep that proves the rename "
        f"complete. Checked only once {INIT_SKILL}/ is gone, and never waivable",
        rule_placeholder_rename,
    ),
    Rule(
        2,
        "agent-docs-unpublished",
        "every markdown file under an agent-docs key below the site source carries "
        "`search: exclude: true` and appears in no nav: entry",
        rule_agent_docs_unpublished,
    ),
    Rule(
        3,
        "vale-length-sections",
        ".vale.ini has one section per agent-docs key, turning its declared Length rule on and "
        "every other Length rule off",
        rule_vale_length_sections,
    ),
    Rule(
        4,
        "glossary-entry-length",
        "no glossary entry in CONTEXT.md exceeds 200 words",
        rule_glossary_entry_length,
    ),
    Rule(
        5,
        "skill-file-location",
        "every SKILL.md lives at skills/<name>/SKILL.md, where the installer looks",
        rule_skill_file_location,
    ),
    Rule(
        6,
        "skill-name-prefix",
        "no repo-local skill directory begins with `lab-`, which the shared plugin owns",
        rule_skill_name_prefix,
    ),
    Rule(
        7,
        "skill-symlinks",
        "`python skills/install.py --check` passes: a link in both discovery paths, nothing "
        "dangling, every link relative",
        rule_skill_symlinks,
    ),
    Rule(
        8,
        "repo-shape",
        "if src/ is tracked then tests/ is too; if release.yml is tracked then CHANGELOG.md is",
        rule_repo_shape,
    ),
    Rule(
        9,
        "init-sentinel",
        "AGENTS.md is still the template's generic copy",
        warning_init_sentinel,
        warns_only=True,
    ),
)

#: Warning 10 is not a rule with a tree to inspect: it is the waiver table itself, printed on
#: every run. It lives in `report` because only the reporter knows what each waiver suppressed.
WAIVER_RULE = (10, "waivers")


def load(root: Path) -> Repo:
    """Read the tree and the declarations, or explain why that was not possible."""
    proc = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise SystemExit(f"conformance: cannot list tracked files in {root}\n{proc.stderr.strip()}")
    tracked = tuple(sorted(p for p in proc.stdout.split("\0") if p))
    pyproject = root / "pyproject.toml"
    if not pyproject.is_file():
        raise SystemExit(f"conformance: no pyproject.toml in {root}")
    tool = tomllib.loads(pyproject.read_text(encoding="utf-8")).get("tool", {})
    liulab = tool.get("liulab", {})
    return Repo(
        root=root,
        tracked=tracked,
        agent_docs=liulab.get("agent-docs", {}),
        waived=liulab.get("waived", {}),
    )


def _verdict(repo: Repo, rule: Rule, result: Result) -> tuple[str, list[str]]:
    """One rule's status word and the lines printed beside it.

    `ok` and `--` are deliberately different words. A rule whose premise was absent checked
    nothing, and a gate that fires on nothing looks exactly like one that passes unless the run
    says which happened.
    """
    waiver = repo.waived.get(rule.name) if rule.name != UNWAIVABLE else None
    if result.problems and waiver is not None:
        return "waived", [f"{len(result.problems)} problem(s) suppressed: {waiver}"]
    if result.problems and rule.warns_only:
        return "warn", [problem.splitlines()[0] for problem in result.problems]
    if result.problems:
        return "FAIL", [f"{len(result.problems)} problem(s), listed below"]
    if not result.notes:
        return "ok", ["nothing to check"]
    # `--` only when EVERY note says the premise was absent. Rule 8 has two independent halves,
    # and a run that really checked one of them has not checked nothing.
    if all(note.startswith("not checked") for note in result.notes):
        return "--", result.notes
    return "ok", result.notes


def _waiver_lines(repo: Repo, results: list[tuple[Rule, Result]]) -> list[str]:
    """Warning 10: every waiver, with what it actually suppressed.

    A waiver naming no rule is called out rather than ignored, and so is one naming rule 1, which
    refuses to be waived. A waiver that quietly does nothing is the same invisible escape hatch
    this table exists to keep visible.
    """
    known = {rule.name: result for rule, result in results}
    lines: list[str] = []
    for name, reason in sorted(repo.waived.items()):
        if name == UNWAIVABLE:
            lines.append(f"{name}: REFUSED, this rule cannot be waived — {reason}")
        elif name not in known:
            lines.append(f"{name}: matches no rule, so it waives nothing — {reason}")
        else:
            lines.append(f"{name}: {len(known[name].problems)} problem(s) suppressed — {reason}")
    return lines or ["none"]


def report(repo: Repo, results: list[tuple[Rule, Result]]) -> int:
    """Print every rule's verdict and every waiver, and return the exit status.

    The status table is printed whatever happened, failures included: what a run did NOT check is
    as much a part of the answer as what it did. Failures go to stderr, in the same
    what / why / fix shape `python skills/install.py --check` prints.
    """
    width = max(len(rule.name) for rule, _ in results)
    rows: list[tuple[str, str, str, list[str]]] = []
    failed: list[tuple[Rule, Result]] = []
    warned: list[tuple[Rule, Result]] = []
    for rule, result in results:
        status, notes = _verdict(repo, rule, result)
        rows.append(
            (f"{'warn' if rule.warns_only else 'rule'} {rule.number}", rule.name, status, notes)
        )
        if status == "FAIL":
            failed.append((rule, result))
        elif status == "warn":
            warned.append((rule, result))
    rows.append((f"warn {WAIVER_RULE[0]}", WAIVER_RULE[1], "--", _waiver_lines(repo, results)))

    print("conformance — the rules a Liu Lab repo keeps")
    for label, name, status, notes in rows:
        head = f"  {label:<8}{name:<{width}}  {status:<7}"
        for index, note in enumerate(notes):
            print(
                textwrap.fill(
                    note,
                    width=110,
                    initial_indent=head if index == 0 else " " * len(head),
                    subsequent_indent=" " * len(head),
                )
            )

    for rule, result in warned:
        print(f"\nwarning {rule.number}  {rule.name} — {rule.statement}")
        for problem in result.problems:
            print(textwrap.indent(problem, "  "))

    total = sum(1 for rule, _ in results if not rule.warns_only)
    if not failed:
        print(f"\nAll {total} rules pass.")
        return 0
    # Flushed before the first byte reaches stderr. `check.sh` captures both streams into one
    # file, where stdout is block-buffered and stderr is not, so without this the failures print
    # above the table that says which rules they came from.
    sys.stdout.flush()
    print(f"\nconformance: {len(failed)} of {total} rules failed.\n", file=sys.stderr)
    for rule, result in failed:
        statement = textwrap.fill(
            rule.statement, width=94, initial_indent="  ", subsequent_indent="  "
        )
        print(f"rule {rule.number}  {rule.name}\n{statement}", file=sys.stderr)
        for problem in result.problems:
            print(textwrap.indent(problem, "  "), file=sys.stderr)
        print("", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    """Run every rule against one tree and report."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="Repository to check. Defaults to the repo this script lives in.",
    )
    args = parser.parse_args(argv)
    repo = load(args.root.resolve())
    results = [(rule, rule.check(repo)) for rule in RULES]
    return report(repo, results)


if __name__ == "__main__":
    raise SystemExit(main())
