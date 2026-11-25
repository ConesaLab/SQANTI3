import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.parallel import chunk_cpu_plan


class TestChunkCpuPlan:
    """Test suite for chunk_cpu_plan."""

    @pytest.mark.parametrize("cpus, n_chunks, available, expected", [
        (8, 8, 8, (8, 1)),    # one CPU per chunk
        (8, 3, 8, (3, 2)),    # fewer chunks than CPUs: the CPUs are shared
        (8, 2, 8, (2, 4)),
        (4, 8, 8, (4, 1)),    # more chunks than CPUs: the rest wait
        (16, 4, 8, (4, 2)),   # -t above the CPUs the job may use
        (1, 4, 8, (1, 1)),
        (8, 1, 8, (1, 8)),
    ])
    def test_plan(self, cpus, n_chunks, available, expected):
        assert chunk_cpu_plan(cpus, n_chunks, available) == expected

    @pytest.mark.parametrize("cpus", [1, 2, 4, 7, 8, 16, 64])
    @pytest.mark.parametrize("n_chunks", [1, 2, 3, 8, 20])
    @pytest.mark.parametrize("available", [1, 4, 8, 48])
    def test_total_stays_within_the_budget(self, cpus, n_chunks, available):
        """Chunks at once x CPUs of each never pass -t nor the CPUs of the job."""
        workers, chunk_cpus = chunk_cpu_plan(cpus, n_chunks, available)
        assert workers >= 1 and chunk_cpus >= 1
        assert workers <= n_chunks
        assert workers * chunk_cpus <= min(cpus, available)
