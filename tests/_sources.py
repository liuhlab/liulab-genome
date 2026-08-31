"""Reading the package's own source, for the guards that assert structure rather than behaviour.

Four test modules parse `src/genome` with `ast` and ask a structural question of the text:
which modules import a name they are banned from (three of them), and which strings a user
can be shown (one). They agreed on the technique and disagreed on everything around it —
five spellings of the package root, two of them repo-relative paths that answer for the
checkout rather than for the installed package; two readers, one of which took the
locale's encoding for the package's own sources; and three near-copies of *pull the module
names out of an import node*, each subtly its own.

The technique is what belongs here, so that a fifth guard inherits it instead of deriving
it again. **What each guard bans stays with the guard** — the ban is the test's subject and
reads where it is asserted, not in a list over here.

**Read, never imported.** A structural ban asked of `sys.modules` would answer for whatever
the suite happened to import first, and importing a module to inspect it drags its whole
package in behind it. So every question here is answered off the source text. That is also
what makes an import deferred into a function body, or hidden under ``TYPE_CHECKING``,
count the same as one at the top: a ban that only the laziest evasion defeats is not a
structural guarantee.

Not a conftest addition. `conftest.py` is shared by every module in the suite, and this is
one technique four of them happen to want; not a fixture either, since nothing here has
setup or teardown.
"""

from __future__ import annotations

import ast
from collections.abc import Iterable
from pathlib import Path

import genome

#: The package as installed, which is what a guard should be asking about — a path built
#: out of this file's own location answers for the checkout instead, and would keep
#: passing against a tree that no longer ships.
PACKAGE = Path(genome.__file__).parent


def sources(root: Path | None = None) -> list[Path]:
    """Return every ``.py`` file under ``root`` (the package by default), in path order.

    Empty is refused rather than returned. A guard that walked nothing would pass, and
    would go on passing — the one failure mode a structural test cannot afford.
    """
    found = sorted((root or PACKAGE).rglob("*.py"))
    assert found, f"no sources found under {root or PACKAGE}, so nothing would be checked"
    return found


def parse(path: Path) -> ast.Module:
    """Return the parsed source of ``path``, read as UTF-8.

    The encoding is named because these are *this package's own* sources: what they are
    written in is a fact about the repository, not about the machine reading them.
    """
    return ast.parse(path.read_text(encoding="utf-8"))


def imported_modules(path: Path) -> set[str]:
    """Return every module ``path`` names in an import, at any depth in the file.

    ``import a.b`` contributes ``a.b``; ``from a.b import c`` contributes ``a.b``. A
    relative import contributes ``""``, which no ban matches and none should — resolving
    it would mean knowing the package the file sits in, and every ban here is written
    against absolute names.
    """
    found: set[str] = set()
    for node in ast.walk(parse(path)):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            found.add(node.module or "")
    return found


def imports_any(path: Path, prefixes: Iterable[str]) -> bool:
    """Whether ``path`` imports any of ``prefixes``, or a module underneath one.

    The dot is load-bearing: ``genome.tf`` is matched by ``genome.tf`` and by
    ``genome.tf.gene.census``, and not by a module that merely starts with those letters.
    A bare `startswith` would ban a `genome.tfidf` nobody meant to ban, which is the kind
    of guard that gets deleted rather than fixed the first time it is wrong.
    """
    wanted = tuple(prefixes)
    return any(
        name == prefix or name.startswith(f"{prefix}.")
        for name in imported_modules(path)
        for prefix in wanted
    )
