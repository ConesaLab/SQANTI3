import os
import sys

import pytest

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src import parallel
from src.parallel import memory_use_percent, read_cgroup_memory

GB = 1024 ** 3


def fake_cgroup(tmp_path, monkeypatch, version, limit, usage, inactive_file, active_file):
    """A job cgroup with a memory limit, as SLURM leaves it with --mem."""
    job = tmp_path / "slurm" / "job_1"
    job.mkdir(parents=True)
    if version == 1:
        controller, limit_name, usage_name, prefix = "memory", "memory.limit_in_bytes", "memory.usage_in_bytes", "total_"
        cgroup_line = "5:memory:/slurm/job_1\n"
    else:
        controller, limit_name, usage_name, prefix = "", "memory.max", "memory.current", ""
        cgroup_line = "0::/slurm/job_1\n"
    (job / limit_name).write_text(f"{limit}\n")
    (job / usage_name).write_text(f"{usage}\n")
    (job / "memory.stat").write_text(f"{prefix}inactive_file {inactive_file}\n{prefix}active_file {active_file}\n")
    proc_cgroup = tmp_path / "proc_cgroup"
    proc_cgroup.write_text(cgroup_line)
    layout = (str(tmp_path),) + parallel.CGROUP_MEMORY_LAYOUTS[controller][1:]
    monkeypatch.setattr(parallel, "CGROUP_MEMORY_LAYOUTS", {controller: layout})
    real_open = open
    monkeypatch.setattr("builtins.open", lambda path, *a, **k: real_open(
        str(proc_cgroup) if path == "/proc/self/cgroup" else path, *a, **k))
    monkeypatch.setattr(parallel.psutil, "virtual_memory",
                        lambda: type("vm", (), {"total": 256 * GB, "percent": 5.0})())


class TestCgroupMemory:
    """Test suite for the memory check of the chunk scheduler."""

    @pytest.mark.parametrize("version", [1, 2])
    def test_file_cache_is_not_counted_as_used(self, tmp_path, monkeypatch, version):
        """Job with 32 GB: 27 GB charged, 21 GB of them are files that were read (6 inactive, 15 active)."""
        fake_cgroup(tmp_path, monkeypatch, version, limit=32 * GB, usage=27 * GB,
                    inactive_file=6 * GB, active_file=15 * GB)
        assert read_cgroup_memory() == (6 * GB, 32 * GB)
        assert memory_use_percent() == pytest.approx(18.8)

    def test_memory_of_the_processes_is_counted(self, tmp_path, monkeypatch):
        fake_cgroup(tmp_path, monkeypatch, 2, limit=32 * GB, usage=30 * GB,
                    inactive_file=1 * GB, active_file=1 * GB)
        assert memory_use_percent() == pytest.approx(87.5)

    def test_cgroup_without_limit(self, tmp_path, monkeypatch):
        """Without a limit below the memory of the machine, only the machine is checked."""
        fake_cgroup(tmp_path, monkeypatch, 2, limit=1024 * GB, usage=27 * GB,
                    inactive_file=0, active_file=0)
        assert read_cgroup_memory() is None
        assert memory_use_percent() == 5.0
