"""Tests for the output-stall watchdog in _run_command_in_process_group."""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import pytest

from metainformant.rna.retrieval.ena_downloader import _run_command_in_process_group

STALLED_WRITER = (
    "import sys, time, pathlib\n"
    "out = pathlib.Path(sys.argv[1])\n"
    "out.mkdir(parents=True, exist_ok=True)\n"
    "(out / 'seed.txt').write_text('x')\n"
    "time.sleep(600)\n"
)
ACTIVE_WRITER = (
    "import time, pathlib, sys\n"
    "out = pathlib.Path(sys.argv[1])\n"
    "out.mkdir(parents=True, exist_ok=True)\n"
    "for i in range(6):\n"
    "    (out / f'chunk{i}.txt').write_text('x' * 1024)\n"
    "    time.sleep(1)\n"
)


def test_stalled_writer_terminates_early(tmp_path: Path) -> None:
    """A writer with frozen output is killed at the stall deadline, not the timeout."""
    out_dir = tmp_path / "quant" / "SRR_X"
    started = time.monotonic()
    with pytest.raises(subprocess.TimeoutExpired) as excinfo:
        _run_command_in_process_group(
            [sys.executable, "-c", STALLED_WRITER, str(out_dir)],
            timeout=600,
            output_dirs=[str(out_dir)],
            stall_timeout=5,
        )
    elapsed = time.monotonic() - started
    assert elapsed < 120, f"stall watchdog did not fire early: {elapsed:.1f}s"
    assert excinfo.value.output == b"stalled"


def test_active_writer_runs_to_completion(tmp_path: Path) -> None:
    """A growing output directory is never treated as stalled."""
    out_dir = tmp_path / "quant" / "SRR_OK"
    result = _run_command_in_process_group(
        [sys.executable, "-c", ACTIVE_WRITER, str(out_dir)],
        timeout=120,
        output_dirs=[str(out_dir)],
        stall_timeout=5,
    )
    assert result.returncode == 0


def test_watchdog_disabled_by_default(tmp_path: Path) -> None:
    """Without output_dirs/stall_timeout the plain path is unchanged."""
    out_dir = tmp_path / "quant" / "SRR_OFF"
    result = _run_command_in_process_group(
        [sys.executable, "-c", "print('ok')"],
        timeout=60,
    )
    assert result.returncode == 0
    assert "ok" in (result.stdout or "")
    assert out_dir.exists() is False
