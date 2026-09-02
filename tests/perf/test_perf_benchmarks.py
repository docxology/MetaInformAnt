"""Fast pytest wrappers over the perf benchmarks.

These are NOT timing assertions (machine-dependent flakiness); they are
smoke tests proving each benchmark script runs end to end on a small
frame, produces its RESULT-LINE, and reports CONSISTENCY: OK. The
full-scale runs and recorded numbers live in tests/perf/BASELINES.md and
the scripts themselves.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PERF_DIR = Path(__file__).resolve().parent


def _run_bench(script: str, *flags: str) -> str:
    result = subprocess.run(
        [sys.executable, str(PERF_DIR / script), *flags],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert result.returncode == 0, f"{script} failed: {result.stderr[-2000:]}"
    return result.stdout


def test_bench_tau_small_frame_consistent() -> None:
    out = _run_bench("bench_tau_expression.py", "--genes", "600", "--samples", "20")
    assert "RESULT-LINE:" in out
    assert "CONSISTENCY: OK" in out


def test_bench_ortholog_small_frame_consistent() -> None:
    out = _run_bench("bench_ortholog_classification.py", "--orthogroups", "2000", "--species", "8")
    assert "RESULT-LINE:" in out
    assert "CONSISTENCY: OK" in out


def test_bench_batches_small_frame_consistent() -> None:
    out = _run_bench("bench_batches_memory.py", "--rows", "20000")
    assert "RESULT-LINE:" in out
    assert "CONSISTENCY: OK" in out
