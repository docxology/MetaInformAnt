"""Campaign-scale benchmark: ortholog bridge-table parsing and classification.

Builds a synthetic transcript_orthogroups.tsv-shaped bridge table (150,000
orthogroups x 27 species — the campaign roster scale) and times the
orthology-class stratification path, including a row-chunked variant using
``metainformant.core.utils.batches.chunked``. Real pandas execution, no test doubles.

    uv run python tests/perf/bench_ortholog_classification.py [--orthogroups N] [--species N]

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

from metainformant.core.utils.batches import chunked, frame_memory_mb
from metainformant.rna.analysis.tissue_specificity import classify_orthogroups, join_expression_with_orthology

DEFAULT_ORTHOGROUPS = 150_000
DEFAULT_SPECIES = 27
SEED = 20260901


def synthetic_bridge_table(n_orthogroups: int, n_species: int) -> pd.DataFrame:
    """Build a deterministic orthogroup x species table of comma-separated IDs."""
    rng = np.random.default_rng(SEED)
    species_cols = [f"species_{j:02d}" for j in range(n_species)]
    data: dict[str, list[str]] = {c: [] for c in species_cols}
    for _ in range(n_orthogroups):
        for col in species_cols:
            r = rng.random()
            if r < 0.25:
                data[col].append("")  # unmapped
            elif r < 0.90:
                data[col].append(f"XM_{rng.integers(1, 10**7):07d}")
            else:
                # multicopy cell: 2-3 comma-separated IDs
                k = int(rng.integers(2, 4))
                data[col].append(",".join(f"XM_{rng.integers(1, 10**7):07d}" for _ in range(k)))
    index = pd.Index([f"OG{n:06d}" for n in range(n_orthogroups)], name="orthogroup")
    return pd.DataFrame(data, index=index)


def synthetic_expression_frame(n_transcripts: int) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 1)
    index = pd.Index([f"XM_{i:07d}" for i in range(1, n_transcripts)], name="transcript")
    return pd.DataFrame(
        rng.exponential(5.0, size=(len(index), 4)),
        index=index,
        columns=["brain", "ovary", "midgut", "malpighian"],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", type=int, default=DEFAULT_ORTHOGROUPS)
    parser.add_argument("--species", type=int, default=DEFAULT_SPECIES)
    args = parser.parse_args()

    t0 = time.perf_counter()
    bridge = synthetic_bridge_table(args.orthogroups, args.species)
    t_gen = time.perf_counter() - t0

    t0 = time.perf_counter()
    classified = classify_orthogroups(bridge)
    t_classify = time.perf_counter() - t0

    value_counts = classified["orthology_class"].value_counts().to_dict()

    expr = synthetic_expression_frame(200_000)
    species_col = bridge.columns[0]

    t0 = time.perf_counter()
    joined = join_expression_with_orthology(expr, bridge, species=species_col)
    t_join = time.perf_counter() - t0

    # Chunked classification: row-chunk the bridge, classify per chunk,
    # then combine. Class counts must match the whole-table pass.
    t0 = time.perf_counter()
    counts_chunked: dict[str, int] = {}
    for chunk_rows in chunked(range(len(bridge)), 25_000):
        part = classify_orthogroups(bridge.iloc[chunk_rows])
        for cls, cnt in part["orthology_class"].value_counts().items():
            counts_chunked[cls] = counts_chunked.get(cls, 0) + int(cnt)
    t_chunked = time.perf_counter() - t0

    print(f"bridge: {args.orthogroups} OG x {args.species} species, {frame_memory_mb(bridge):.1f} MiB deep")
    print(f"generate:    {t_gen:.3f}s")
    print(f"classify:    {t_classify:.3f}s  ({args.orthogroups / max(t_classify, 1e-9):.0f} OG/s)")
    print(f"             classes={ {k: int(v) for k, v in value_counts.items()} }")
    print(f"join expr:   {t_join:.3f}s  joined_rows={len(joined)}")
    print(f"chunked:     {t_chunked:.3f}s  counts={counts_chunked}")
    ok = all(counts_chunked.get(k, 0) == int(v) for k, v in value_counts.items())
    print("CONSISTENCY: OK" if ok else "CONSISTENCY: FAIL")
    print(
        "RESULT-LINE: "
        f"ogs={args.orthogroups} species={args.species} gen_s={t_gen:.3f} "
        f"classify_s={t_classify:.3f} join_s={t_join:.3f} chunked_s={t_chunked:.3f} "
        f"mem_mib={frame_memory_mb(bridge):.1f}"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
