"""Benchmark: memory-hygiene helpers over a campaign-scale frame.

Times chunked vs whole-frame aggregation of a synthetic 500k-row
per-sample expression-like table, and reports deep memory. Real execution,
no mocks; asserts the chunked aggregation is numerically consistent.

    uv run python tests/perf/bench_batches_memory.py [--rows N]

Recorded baselines live in tests/perf/BASELINES.md.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from metainformant.core.utils.batches import frame_memory_mb, frame_memory_report, iter_frame_chunks

DEFAULT_ROWS = 500_000
SEED = 20260901


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, default=DEFAULT_ROWS)
    args = parser.parse_args()

    rng = np.random.default_rng(SEED)
    frame = pd.DataFrame(
        {
            "transcript": [f"XM_{i:07d}" for i in range(args.rows)],
            "sample": rng.integers(0, 120, size=args.rows),
            "count": rng.exponential(10.0, size=args.rows),
        }
    )
    report = frame_memory_report(frame)

    t0 = time.perf_counter()
    whole_total = float(frame["count"].sum())
    t_whole = time.perf_counter() - t0

    t0 = time.perf_counter()
    chunk_total = 0.0
    n_chunks = 0
    for chunk in iter_frame_chunks(frame, rows_per_chunk=50_000):
        chunk_total += float(chunk["count"].sum())
        n_chunks += 1
    t_chunked = time.perf_counter() - t0

    peak_ok = abs(whole_total - chunk_total) <= 1e-9 * max(1.0, abs(whole_total))
    print(
        f"frame: {report['n_rows']} rows x {report['n_cols']} cols, "
        f"{report['memory_mb']:.1f} MiB deep, {report['memory_per_row_bytes']:.0f} B/row"
    )
    print(f"whole-frame sum:   {t_whole:.4f}s")
    print(
        f"chunked sum:       {t_chunked:.4f}s  chunks={n_chunks} "
        f"(whole/chunked ratio={t_whole / max(t_chunked, 1e-9):.2f}x)"
    )
    print(f"memory_mb (deep):  {frame_memory_mb(frame):.1f}")
    print("CONSISTENCY: OK" if peak_ok else "CONSISTENCY: FAIL")
    print(
        "RESULT-LINE: "
        f"rows={args.rows} whole_s={t_whole:.4f} chunked_s={t_chunked:.4f} "
        f"chunks={n_chunks} mem_mib={report['memory_mb']:.1f}"
    )
    return 0 if peak_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
