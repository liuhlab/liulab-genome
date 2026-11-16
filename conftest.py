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
to neither: the default run deselects it, and a run in which one ran is not timed.
"""

from __future__ import annotations

import time
from collections.abc import Generator

import pytest

pytest_plugins = ["pytester", "tests._guards"]

SLOW = "slow"

#: Rides on a `slow` test's reports, which pytest-xdist carries from a worker to the controller.
SLOW_TAG = (SLOW, True)


class SuiteBudget:
    """Times one run's tests, and fails the run when that is over budget, unless a `slow` test ran."""

    def __init__(self, budget: float) -> None:
        self.budget = budget
        self.started = time.perf_counter()
        self.ended: float | None = None
        self.slow_ran = False
        self.over = False
        self.verdict = ""

    def pytest_sessionstart(self) -> None:
        """Start the clock where pytest starts its own."""
        self.started = time.perf_counter()

    @pytest.hookimpl(wrapper=True)
    def pytest_runtestloop(self) -> Generator[None, object, object]:
        """Stop the clock where the tests end, so later session-finish work is not counted."""
        try:
            return (yield)
        finally:
            self.ended = time.perf_counter()

    def pytest_runtest_logreport(self, report: pytest.TestReport) -> None:
        """Note a `slow` test."""
        if SLOW_TAG in report.user_properties:
            self.slow_ran = True

    @pytest.hookimpl(trylast=True)
    def pytest_sessionfinish(self, session: pytest.Session) -> None:
        """Fail a run that took longer than its budget."""
        if self.slow_ran:
            self.verdict = "suite time: not checked, because slow tests ran"
            return
        took = (self.ended or time.perf_counter()) - self.started
        self.over = took > self.budget
        if not self.over:
            self.verdict = f"suite time: {took:.1f} s of its {self.budget:g} s budget"
            return
        self.verdict = (
            f"suite time: {took:.1f} s, over its {self.budget:g} s budget. "
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
        "Most seconds a run may take; a run in which a slow test ran is not timed. "
        "0 or unset: no budget.",
        type="float",
        default=0.0,
    )


def pytest_configure(config: pytest.Config) -> None:
    """Register the `slow` marker, and the budget when one is set, outside xdist's workers."""
    config.addinivalue_line(
        "markers",
        f"{SLOW}: exempt from timeout and suite_budget, and deselected unless -m {SLOW}",
    )
    budget = config.getini("suite_budget")
    if budget and not hasattr(config, "workerinput"):
        config.pluginmanager.register(SuiteBudget(budget), "suite-budget")


@pytest.hookimpl(wrapper=True)
def pytest_runtest_makereport(
    item: pytest.Item,
) -> Generator[None, pytest.TestReport, pytest.TestReport]:
    """Tag a `slow` test's report."""
    report = yield
    if item.get_closest_marker(SLOW):
        report.user_properties.append(SLOW_TAG)
    return report


def pytest_timeout_set_timer(item: pytest.Item) -> bool | None:
    """Set no timer for a `slow` test; pytest-timeout sets one for every other."""
    return True if item.get_closest_marker(SLOW) else None
