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
    def test_it_is_one(self) -> None:
        # An imported library that started a pool would make an unguarded caller script
        # re-execute itself under the spawn start method.
        assert DEFAULT_WORKERS == 1

    def test_one_resolves_to_one_without_reading_anything(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
        assert resolve_workers(DEFAULT_WORKERS) == 1


class TestAnExplicitCountIsHonoured:
    @pytest.mark.parametrize("workers", [1, 2, 4, 64])
    def test_a_number_is_taken_as_given(self, workers: int) -> None:
        assert resolve_workers(workers) == workers

    def test_it_beats_the_allocation(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # A caller who says 2 gets 2, on a laptop or in a sixteen-CPU allocation.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
        assert resolve_workers(2) == 2

    @pytest.mark.parametrize("workers", [0, -1])
    def test_fewer_than_one_says_what_one_and_none_mean(self, workers: int) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            resolve_workers(workers)


class TestResolvingFromTheEnvironment:
    def test_the_slurm_allocation_comes_first(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "3")
        assert resolve_workers(None) == 3

    def test_cpus_per_task_beats_cpus_on_node(self, monkeypatch: pytest.MonkeyPatch) -> None:
        # The task's own share is what this process may use; the node's is the job's.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
        monkeypatch.setenv("SLURM_CPUS_ON_NODE", "16")
        assert resolve_workers(None) == 2

    def test_cpus_on_node_is_the_fallback(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("SLURM_CPUS_ON_NODE", "5")
        assert resolve_workers(None) == 5

    @pytest.mark.parametrize("value", ["", "   ", "2(x3)", "many", "0", "-4"])
    def test_a_value_that_is_not_a_count_is_ignored(
        self, monkeypatch: pytest.MonkeyPatch, value: str
    ) -> None:
        # Slurm writes "2(x3)" on a heterogeneous job. Falling through beats crashing, and
        # beats reading a 2 out of it that was never meant.
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", value)
        assert resolve_workers(None) == resolve_workers(None)
        assert resolve_workers(None) >= 1

    def test_an_unreadable_value_falls_through_to_the_next_variable(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2(x3)")
        monkeypatch.setenv("SLURM_CPUS_ON_NODE", "7")
        assert resolve_workers(None) == 7

    def test_outside_an_allocation_it_is_the_machine_and_at_least_one(self) -> None:
        found = resolve_workers(None)
        assert found >= 1
        assert found <= (os.cpu_count() or 1)

    def test_the_variables_it_reads_are_the_two_named(self) -> None:
        assert SLURM_CPU_VARS == ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")


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
