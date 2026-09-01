"""Descriptive multi-omics summaries.

Provides per-omics descriptive statistics and cross-omics pairwise block
correlations over shared samples. Purely descriptive: no inferential testing,
no p-value claims — the cross-omics correlations are point estimates over
shared samples only, suitable for campaign dashboards and atlas-style
summarization regardless of evidence-manifest freeze state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


@dataclass
class OmicsBlockSummary:
    """Descriptive summary of one omics data block.

    Attributes:
        name: Block name (e.g. "rna", "protein").
        n_samples: Number of samples (rows).
        n_features: Number of features (columns).
        mean: Per-feature mean (NaN for all-NaN columns).
        std: Per-feature sample standard deviation (ddof=1).
        missing_fraction: Fraction of missing values across the block.
        quantiles: Requested quantiles of all non-NaN values, keyed by q.
    """

    name: str
    n_samples: int
    n_features: int
    mean: Any  # np.ndarray (n_features,)
    std: Any  # np.ndarray (n_features,)
    missing_fraction: float
    quantiles: Dict[str, float]


@dataclass
class BlockCorrelationSummary:
    """Cross-omics pairwise correlation summary for one block pair.

    Attributes:
        block_a / block_b: Block names.
        n_shared_samples: Samples present in both blocks.
        corr_mean / corr_median: Mean/median of per-feature absolute Spearman
            correlations between matched features (union by name).
        n_feature_pairs: Number of matched feature pairs used.
        top_pairs: Highest-|rho| feature pairs as (feature, rho) tuples,
            descending, capped at ``top_n``.
    """

    block_a: str
    block_b: str
    n_shared_samples: int
    corr_mean: float
    corr_median: float
    n_feature_pairs: int
    top_pairs: List[Tuple[str, float]]


@dataclass
class MultiOmicsDescriptiveSummary:
    """Aggregated descriptive summary across omics blocks."""

    blocks: List[OmicsBlockSummary]
    block_correlations: List[BlockCorrelationSummary]

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dict representation for dashboards/manifests."""
        return {
            "blocks": [
                {
                    "name": b.name,
                    "n_samples": b.n_samples,
                    "n_features": b.n_features,
                    "missing_fraction": b.missing_fraction,
                    "quantiles": b.quantiles,
                    "mean_head": {
                        f: (None if not np.isfinite(m) else round(float(m), 6))
                        for f, m in list(zip(_feature_names_safe(b), b.mean))[:10]
                    },
                }
                for b in self.blocks
            ],
            "block_correlations": [
                {
                    "pair": [c.block_a, c.block_b],
                    "n_shared_samples": c.n_shared_samples,
                    "n_feature_pairs": c.n_feature_pairs,
                    "abs_spearman_mean": round(c.corr_mean, 6),
                    "abs_spearman_median": round(c.corr_median, 6),
                    "top_pairs": [[f, round(r, 6)] for f, r in c.top_pairs],
                }
                for c in self.block_correlations
            ],
        }


def _feature_names_safe(block: OmicsBlockSummary) -> List[str]:
    return [f"feature_{i}" for i in range(len(block.mean))]


def summarize_omics_block(
    data: pd.DataFrame,
    name: str,
    *,
    quantiles: Tuple[float, ...] = (0.25, 0.5, 0.75),
) -> OmicsBlockSummary:
    """Compute descriptive statistics for one omics block.

    Args:
        data: Samples x features DataFrame.
        name: Block label.
        quantiles: Quantiles of the pooled non-NaN values to report.

    Returns:
        OmicsBlockSummary with per-feature mean/std and pooled quantiles.
    """
    if not isinstance(data, pd.DataFrame):
        raise TypeError("summarize_omics_block expects a pandas DataFrame")
    arr = data.to_numpy(dtype=np.float64, na_value=np.nan)
    mean = np.nanmean(arr, axis=0) if arr.size else np.array([])
    std = np.nanstd(arr, axis=0, ddof=1) if arr.size else np.array([])
    flat = arr.ravel()
    flat = flat[~np.isnan(flat)]
    pooled = {str(q): float(np.quantile(flat, q)) for q in quantiles if flat.size}
    n_missing = int(np.isnan(arr).sum())
    total = int(arr.size) if arr.size else 1
    logger.info("Block summary %s: %d x %d", name, data.shape[0], data.shape[1])
    return OmicsBlockSummary(
        name=name,
        n_samples=int(data.shape[0]),
        n_features=int(data.shape[1]),
        mean=mean,
        std=std,
        missing_fraction=n_missing / total,
        quantiles=pooled,
    )


def cross_omics_block_correlation(
    data_a: pd.DataFrame,
    data_b: pd.DataFrame,
    name_a: str,
    name_b: str,
    *,
    top_n: int = 10,
) -> BlockCorrelationSummary:
    """Compute descriptive cross-omics feature-pair correlations (Spearman).

    Features are matched by name (inner join on columns); correlations are
    computed per feature pair across shared samples. Descriptive point
    estimates only — no hypothesis tests are performed.

    Args:
        data_a / data_b: Samples x features DataFrames with named features.
        name_a / name_b: Block labels for the summary.
        top_n: Number of top |rho| pairs to retain.

    Returns:
        BlockCorrelationSummary; n_feature_pairs=0 when no features overlap.
    """
    shared = data_a.index.intersection(data_b.index)
    if len(shared) < 3:
        return BlockCorrelationSummary(name_a, name_b, int(len(shared)), float("nan"), float("nan"), 0, [])
    a = data_a.loc[shared]
    b = data_b.loc[shared]
    common_features = a.columns.intersection(b.columns)
    if len(common_features) == 0:
        logger.info("No shared features between %s and %s", name_a, name_b)
        return BlockCorrelationSummary(name_a, name_b, int(len(shared)), float("nan"), float("nan"), 0, [])

    rhos: List[Tuple[str, float]] = []
    for feature in common_features:
        x = a[feature].to_numpy(dtype=np.float64)
        y = b[feature].to_numpy(dtype=np.float64)
        mask = ~(np.isnan(x) | np.isnan(y))
        if mask.sum() < 3:
            continue
        xs, ys = x[mask], y[mask]
        if np.std(xs) == 0 or np.std(ys) == 0:
            continue
        rx = _rankdata(xs)
        ry = _rankdata(ys)
        rho = float(np.corrcoef(rx, ry)[0, 1])
        if np.isfinite(rho):
            rhos.append((str(feature), rho))

    if not rhos:
        return BlockCorrelationSummary(name_a, name_b, int(len(shared)), float("nan"), float("nan"), 0, [])

    abs_vals = np.array([abs(r) for _, r in rhos])
    top = sorted(rhos, key=lambda t: abs(t[1]), reverse=True)[: max(1, top_n)]
    logger.info("Cross-omics %s~%s: %d feature pairs", name_a, name_b, len(rhos))
    return BlockCorrelationSummary(
        block_a=name_a,
        block_b=name_b,
        n_shared_samples=int(len(shared)),
        corr_mean=float(abs_vals.mean()),
        corr_median=float(np.median(abs_vals)),
        n_feature_pairs=len(rhos),
        top_pairs=top,
    )


def _rankdata(values: Any) -> Any:
    """Average-rank data (ascending, ties averaged), matching scipy.stats.rankdata."""
    arr = np.asarray(values, dtype=np.float64)
    order = np.argsort(arr, kind="mergesort")
    ranks = np.empty(arr.size, dtype=np.float64)
    sorted_vals = arr[order]
    i = 0
    while i < arr.size:
        j = i
        while j + 1 < arr.size and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        ranks[order[i : j + 1]] = avg_rank
        i = j + 1
    return ranks


def multi_omics_descriptive_summary(
    omics_data: Dict[str, pd.DataFrame],
    *,
    quantiles: Tuple[float, ...] = (0.25, 0.5, 0.75),
    top_n: int = 10,
) -> MultiOmicsDescriptiveSummary:
    """Build a full descriptive summary across omics blocks.

    Args:
        omics_data: Mapping of block name to samples x features DataFrame.
        quantiles: Quantiles for per-block pooled summaries.
        top_n: Top |rho| feature pairs retained per block pair.

    Returns:
        MultiOmicsDescriptiveSummary over blocks (alphabetical order) and all
        block pairs.
    """
    if not omics_data:
        raise ValueError("omics_data must contain at least one block")

    blocks = [summarize_omics_block(data, name, quantiles=quantiles) for name, data in sorted(omics_data.items())]

    names = sorted(omics_data.keys())
    correlations: List[BlockCorrelationSummary] = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            correlations.append(
                cross_omics_block_correlation(
                    omics_data[names[i]],
                    omics_data[names[j]],
                    names[i],
                    names[j],
                    top_n=top_n,
                )
            )

    return MultiOmicsDescriptiveSummary(blocks=blocks, block_correlations=correlations)
