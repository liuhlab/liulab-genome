"""The suite's guards, time limits and slow lane, declared where they reach ``src/``.

A conftest's fixtures reach the directory it sits in and everything under it, and nothing
above — so ``tests/conftest.py`` reaches tests and only tests. The package's docstring
examples are collected from ``src/`` (``--doctest-modules``, see ``pyproject.toml``), which
is outside that tree, and a doctest item is a test: it runs behind both autouse guards or
it is the hole in them. This file is the one conftest above both trees, so declaring the
guards here is what gives them the reach the promise already claimed. The guards themselves
stay in :mod:`tests._guards`, which is also where the rest of the suite reaches for
``install_network_guard``.

``pytester`` is enabled alongside them because a plugin list is declared in the root
conftest or nowhere: ``tests/test_doctest_guards.py`` runs a throwaway module in a pytest
of its own to prove the reach rather than assert the shape of this file.

`suite_budget` in `[tool.pytest.ini_options]` is the most seconds a run may take, and
`timeout`, read by pytest-timeout, the most one test may take. A test marked `slow` is held
to neither, and the default run deselects it.
"""

from __future__ import annotations

import time

import pytest

pytest_plugins = ["pytester", "tests._guards"]

SLOW = "slow"


class SuiteBudget:
    """Times one run, less its `slow` tests, and fails it when that is over budget."""

    def __init__(self, budget: float) -> None:
        self.budget = budget
        self.started = time.perf_counter()
        self.slow_ids: set[str] = set()
        self.slow_seconds = 0.0
        self.over = False
        self.verdict = ""

    def pytest_sessionstart(self) -> None:
        """Start the clock where pytest starts its own."""
        self.started = time.perf_counter()

    def pytest_collection_modifyitems(self, items: list[pytest.Item]) -> None:
        """Note which tests carry the mark. Keywords would also match a directory named `slow`."""
        self.slow_ids = {item.nodeid for item in items if item.get_closest_marker(SLOW)}

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        """Count a `slow` test's setup, call and teardown, so the budget leaves them out."""
        if report.nodeid in self.slow_ids:
            self.slow_seconds += report.duration

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session: pytest.Session) -> None:
        """Fail a run that took longer than its budget."""
        took = time.perf_counter() - self.started - self.slow_seconds
        excluded = (
            f" ({self.slow_seconds:.1f} s of slow tests not counted)" if self.slow_seconds else ""
        )
        self.over = took > self.budget
        if not self.over:
            self.verdict = f"suite time: {took:.1f} s of its {self.budget:g} s budget{excluded}"
            return
        self.verdict = (
            f"suite time: {took:.1f} s, over its {self.budget:g} s budget{excluded}. "
            "Find the slowest tests with --durations=10, then make them faster or mark them slow."
        )
        if session.exitstatus == pytest.ExitCode.OK:
            session.exitstatus = pytest.ExitCode.TESTS_FAILED

    def pytest_terminal_summary(self, terminalreporter: pytest.TerminalReporter) -> None:
        """Print the time against the budget."""
        terminalreporter.write_line(self.verdict, red=self.over, bold=self.over)


def pytest_addoption(parser: pytest.Parser) -> None:
    """Register `suite_budget`."""
    parser.addini(
        "suite_budget",
        "Most seconds a run may take, slow tests excluded. 0 or unset: no budget.",
        type="float",
        default=0.0,
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register the `slow` marker, and the budget when one is set."""
    config.addinivalue_line(
        "markers",
        f"{SLOW}: exempt from timeout and suite_budget, and deselected unless -m {SLOW}",
    )
    budget = config.getini("suite_budget")
    if budget:
        config.pluginmanager.register(SuiteBudget(budget), "suite-budget")


def pytest_timeout_set_timer(item: pytest.Item) -> bool | None:
    """Set no timer for a `slow` test; pytest-timeout sets one for every other."""
    return True if item.get_closest_marker(SLOW) else None
