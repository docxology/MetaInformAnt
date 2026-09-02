"""Campaign-scale benchmark: tau tissue-specificity over a large expression frame.

Synthetic data only (8,000 genes x 120 tissue/sample columns, the campaign's
upper per-species envelope). Times real pandas/numpy execution end to end.
Run directly:

    uv run python tests/perf/bench_tau_expression.py [--genes N] [--samples N]

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

from metainformant.core.utils.batches import frame_memory_report, iter_frame_chunks
from metainformant.rna.analysis.tissue_specificity import compute_tau, tau_summary

DEFAULT_GENES = 8_000
DEFAULT_SAMPLES = 120
SEED = 20260901


def synthetic_expression_frame(n_genes: int, n_samples: int) -> pd.DataFrame:
    """Build a deterministic TPM-like expression frame (real data in memory)."""
    rng = np.random.default_rng(SEED)
    log_means = rng.normal(3.0, 2.0, size=n_genes)
    tissue_bias = rng.normal(0.0, 1.0, size=(n_genes, n_samples))
    values = np.exp(log_means[:, None] + tissue_bias)
    values[::50] = 0.0  # include all-zero rows (tau-undefined path)
    index = pd.Index([f"GENE{i:05d}" for i in range(n_genes)], name="gene")
    columns = [f"tissue_{j:03d}" for j in range(n_samples)]
    return pd.DataFrame(values, index=index, columns=columns)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, default=DEFAULT_GENES)
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    args = parser.parse_args()

    t0 = time.perf_counter()
    frame = synthetic_expression_frame(args.genes, args.samples)
    t_gen = time.perf_counter() - t0

    t0 = time.perf_counter()
    tau = compute_tau(frame)
    t_tau = time.perf_counter() - t0

    t0 = time.perf_counter()
    summary = tau_summary(tau)
    t_summary = time.perf_counter() - t0

    # Chunked pass: identical tau sums must come out of row-block iteration.
    # PITFALL (measured): compute_tau applies its lowest-10% mean-expression
    # filter over the rows it is given, so naively chunking raw rows selects
    # different genes per chunk and shifts results (verified ~0.14% drift on
    # the default frame). The correct chunked pattern filters GLOBALLY first,
    # then chunks the filtered frame — pure tau arithmetic is row-local.
    from metainformant.rna.analysis.tissue_specificity import filter_low_expression

    t0 = time.perf_counter()
    filtered = filter_low_expression(frame, lowest_fraction=0.10)
    chunk_valid_total = 0.0
    chunk_rows = 0
    for chunk in iter_frame_chunks(filtered, rows_per_chunk=2_000):
        chunk_tau = compute_tau(chunk, lowest_fraction=0.0)
        chunk_valid_total += float(chunk_tau.dropna().sum())
        chunk_rows += len(chunk)
    t_chunked = time.perf_counter() - t0
    full_valid_total = float(tau.dropna().sum())

    report = frame_memory_report(frame)
    print(
        f"frame: {report['n_rows']} genes x {report['n_cols']} samples, "
        f"{report['memory_mb']:.1f} MiB deep, {report['memory_per_row_bytes']:.0f} B/row"
    )
    print(f"generate:       {t_gen:.3f}s")
    print(f"compute_tau:    {t_tau:.3f}s  ({args.genes * args.samples / max(t_tau, 1e-9) / 1e6:.1f} M cell/s)")
    print(f"tau_summary:    {t_summary:.4f}s  valid={summary['valid_count']} nan={summary['nan_count']}")
    print(f"chunked pass:   {t_chunked:.3f}s  rows={chunk_rows}")
    print(f"tau consistency: full={full_valid_total:.6f} chunked={chunk_valid_total:.6f}")
    if abs(full_valid_total - chunk_valid_total) > 1e-6 * max(1.0, abs(full_valid_total)):
        print("CONSISTENCY: FAIL")
        return 1
    print("CONSISTENCY: OK")
    print(
        "RESULT-LINE: "
        f"genes={args.genes} samples={args.samples} gen_s={t_gen:.3f} "
        f"tau_s={t_tau:.3f} summary_s={t_summary:.4f} chunked_s={t_chunked:.3f} "
        f"mem_mib={report['memory_mb']:.1f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
