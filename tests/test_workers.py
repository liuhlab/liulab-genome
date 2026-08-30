"""Tests for genome.tf.motif.workers — how many processes a scan may use.

A pure function over the environment, so it is tested as one: **no process is started
anywhere in this file**, and one test proves that by making the ways of starting one raise.
That is the point of the split — the count is decided by arithmetic on environment
variables, and the arithmetic is what can be wrong.

The two things that can be wrong: defaulting to the machine's cores inside a two-CPU Slurm
allocation, which would start fourteen workers too many on a sixteen-core node; and
defaulting to anything but 1 in the library, which would have an imported package start a
process pool that re-imports its caller's script.

Slurm variables are deleted rather than assumed absent, so this passes on a cluster.

The unit lane, unmarked: nothing here needs a binary — or a subprocess.
"""

from __future__ import annotations

import concurrent.futures
import multiprocessing
import multiprocessing.process
import os

import pytest

from genome.tf.motif.workers import DEFAULT_WORKERS, SLURM_CPU_VARS, resolve_workers


@pytest.fixture(autouse=True)
def outside_an_allocation(monkeypatch: pytest.MonkeyPatch) -> None:
    """Clear every Slurm variable, so a test that wants one sets it itself."""
    for name in SLURM_CPU_VARS:
        monkeypatch.delenv(name, raising=False)


class TestTheLibraryDefault:
    def test_it_is_one_and_ignores_the_environment(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # An imported library that started a pool would make an unguarded caller script
        # re-execute itself under the spawn start method.
        assert DEFAULT_WORKERS == 1
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        assert resolve_workers(DEFAULT_WORKERS) == 1


class TestAnExplicitCountIsHonoured:
    def test_an_explicit_count_is_taken_as_given_and_beats_the_allocation(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        for workers in (1, 2, 4, 64):
            assert resolve_workers(workers) == workers
        # A caller who says 2 gets 2, on a laptop or in a sixteen-CPU allocation.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        assert resolve_workers(2) == 2

    def test_fewer_than_one_says_what_one_and_none_mean(self) -> None:
        for workers in (0, -1):
            with pytest.raises(ValueError, match="at least 1"):
                resolve_workers(workers)


class TestResolvingFromTheEnvironment:
    def test_slurm_variables_are_read_in_priority_order(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        assert SLURM_CPU_VARS == ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "3")
        assert resolve_workers(None) == 3
        # The task's own share is what this process may use; the node's is the job's.
        monkeypatch.setenv("SLURM_CPUS_ON_NODE", "16")
        assert resolve_workers(None) == 3
        monkeypatch.delenv("SLURM_CPUS_PER_TASK")
        assert resolve_workers(None) == 16  # falls back to the node's count alone

    def test_a_value_that_is_not_a_count_falls_through_rather_than_crashing(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        # Slurm writes "2(x3)" on a heterogeneous job. Falling through beats crashing, and
        # beats reading a 2 out of it that was never meant.
        for value in ["", "   ", "2(x3)", "many", "0", "-4"]:
            monkeypatch.setenv("SLURM_CPUS_PER_TASK", value)
            assert resolve_workers(None) == resolve_workers(None)
            assert resolve_workers(None) >= 1
        # And falls through to the next variable, not straight to 1.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2(x3)")
        monkeypatch.setenv("SLURM_CPUS_ON_NODE", "7")
        assert resolve_workers(None) == 7

    def test_outside_an_allocation_it_is_the_machine_and_at_least_one(self) -> None:
        found = resolve_workers(None)
        assert found >= 1
        assert found <= (os.cpu_count() or 1)


class TestNothingIsStarted:
    def test_resolving_starts_no_process(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # The criterion this module exists to be testable against: worker-count resolution
        # is arithmetic, and arithmetic does not fork.
        def refuse(*args: object, **kwargs: object) -> object:
            raise AssertionError("a process was started")

        monkeypatch.setattr(multiprocessing, "Process", refuse)
        monkeypatch.setattr(multiprocessing.process.BaseProcess, "start", refuse)
        monkeypatch.setattr(concurrent.futures, "ProcessPoolExecutor", refuse)
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "3")
        assert resolve_workers(None) == 3
        assert resolve_workers(2) == 2
