"""The suite's two autouse guards reach the docstring examples, which sit outside ``tests/``.

A doctest item collected from ``src/`` is a test, and it gets the same two promises every
other test in the lane gets — it reaches no network, and it sees a **Data dir** of its own
— or it is the hole in them. Neither promise is free here: a conftest's fixtures reach the
directory it sits in and nothing above, so for as long as both guards lived under
``tests/`` they had no reach over ``src/`` at all. ``conftest.py`` at the repository root
is what gives them that reach, and this is what proves it.

**A structural assertion would not.** That the root conftest exists, or that it names
:mod:`tests._guards`, says nothing about whether a guard is *in force* where it matters.
So this runs a throwaway module in a pytest of its own — outside ``tests/``, collected only
because ``--doctest-modules`` is on, given the real root conftest to load — and asks the
run what happened. The run is a subprocess on purpose: in-process it would inherit the
guards this very test is running behind, and would pass just as happily with the wiring
removed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genome.store.data_dir import LIULAB_DATA_ENV

#: The wiring under test, read rather than restated. The throwaway run is handed the real
#: file, so a ``pytest_plugins`` line that stops naming the guards fails this test.
ROOT_CONFTEST = Path(__file__).resolve().parents[1] / "conftest.py"

#: The repository root, which is what has to be importable for the root conftest to find
#: :mod:`tests._guards` from inside a subprocess that starts somewhere else entirely.
REPO_ROOT = ROOT_CONFTEST.parent

#: Reserved by RFC 2606 and resolvable by nothing, so a guard that failed to fire costs a
#: failed DNS lookup rather than a download. What proves the guard fired is the message
#: below, which only the guard can produce — not the mere fact that the example failed.
UNREACHABLE = "https://example.invalid/hg38.fa.gz"

#: An example that would leave the machine. Written as if it succeeds, because what it is
#: here to show is the refusal.
REACHES_THE_NETWORK = f'''
"""A module standing in for one under ``src/``: outside ``tests/``, carrying an example.

Examples
--------
>>> import requests
>>> requests.get({UNREACHABLE!r})
<Response [200]>
"""
'''

#: An example that asks where the **Data dir** points. It passes when the guard re-pointed
#: it at a directory of this item's own, and raises ``KeyError`` when nothing did.
ASKS_FOR_THE_DATA_ROOT = f'''
"""A module standing in for one under ``src/``, whose example reads the data root.

Examples
--------
>>> import os
>>> from pathlib import Path
>>> Path(os.environ[{LIULAB_DATA_ENV!r}]).is_dir()
True
"""
'''


def test_a_docstring_example_outside_tests_runs_behind_both_autouse_guards(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch
) -> None:
    pytester.makeconftest(ROOT_CONFTEST.read_text(encoding="utf-8"))
    pytester.makepyfile(
        reaches_the_network=REACHES_THE_NETWORK,
        asks_for_the_data_root=ASKS_FOR_THE_DATA_ROOT,
    )
    # The subprocess starts with neither the repository on its path nor this machine's
    # data root in its environment: the first is what lets the root conftest resolve the
    # guards, and unsetting the second is what makes the data-root example prove the
    # fixture rather than inherit an answer from whoever ran the suite.
    monkeypatch.setenv("PYTHONPATH", str(REPO_ROOT))
    monkeypatch.delenv(LIULAB_DATA_ENV, raising=False)

    result = pytester.runpytest_subprocess("--doctest-modules")

    # One of each, and which is which is not in doubt: only the network example can fail
    # while the guards hold, and only it fails with a message the guard writes.
    result.assert_outcomes(passed=1, failed=1)
    result.stdout.fnmatch_lines([f"*blocked network call: requests GET {UNREACHABLE}*"])
