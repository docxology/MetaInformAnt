"""Tests for descriptive multi-omics summaries.

Real synthetic omics DataFrames with shared samples and known planted
correlations; descriptive-only assertions (no inferential claims).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.multiomics.analysis.descriptive_summary import (
    MultiOmicsDescriptiveSummary,
    cross_omics_block_correlation,
    multi_omics_descriptive_summary,
    summarize_omics_block,
)


def _rna_block(rng: np.random.Generator, samples: list[str]) -> pd.DataFrame:
    latent = rng.normal(size=len(samples))
    return pd.DataFrame(
        {
            "GENE1": latent + rng.normal(scale=0.1, size=len(samples)),
            "GENE2": -latent + rng.normal(scale=0.1, size=len(samples)),
            "GENE3": rng.normal(size=len(samples)),
        },
        index=samples,
    )


def _protein_block(rng: np.random.Generator, samples: list[str]) -> pd.DataFrame:
    latent = None  # reconstructed below from GENE1-like latent via fresh draw
    del latent
    shared = rng.normal(size=len(samples))
    return pd.DataFrame(
        {
            "PROT1": shared + rng.normal(scale=0.1, size=len(samples)),
            "PROT2": rng.normal(size=len(samples)),
        },
        index=samples,
    )


def test_summarize_omics_block_statistics():
    rng = np.random.default_rng(5)
    samples = [f"S{i}" for i in range(20)]
    data = _rna_block(rng, samples)
    summary = summarize_omics_block(data, "rna")

    assert summary.name == "rna"
    assert summary.n_samples == 20
    assert summary.n_features == 3
    assert summary.missing_fraction == 0.0
    np.testing.assert_allclose(summary.mean, data.mean().to_numpy())
    np.testing.assert_allclose(summary.std, data.std(ddof=1).to_numpy())
    assert set(summary.quantiles) == {"0.25", "0.5", "0.75"}


def test_summarize_block_missing_values_and_type_check():
    data = pd.DataFrame({"A": [1.0, np.nan, 3.0], "B": [4.0, 5.0, np.nan]})
    summary = summarize_omics_block(data, "meth")
    assert summary.missing_fraction == pytest.approx(2 / 6)
    assert summary.n_samples == 3
    with pytest.raises(TypeError):
        summarize_omics_block([[1, 2], [3, 4]], "bad")  # type: ignore[arg-type]


def test_cross_omics_correlation_planted_signal():
    rng = np.random.default_rng(9)
    samples = [f"S{i}" for i in range(30)]
    latent = rng.normal(size=len(samples))
    # Same gene symbol in both blocks: mRNA level vs protein level of GENE1.
    rna = pd.DataFrame({"GENE1": latent}, index=samples)
    prot = pd.DataFrame(
        {"GENE1": latent + rng.normal(scale=0.05, size=len(samples))},
        index=samples,
    )
    corr = cross_omics_block_correlation(rna, prot, "rna", "protein", top_n=5)

    assert corr.n_shared_samples == 30
    assert corr.n_feature_pairs == 1
    assert corr.corr_mean > 0.99
    assert corr.top_pairs[0][0] == "GENE1"
    assert corr.top_pairs[0][1] > 0.99


def test_cross_omics_no_overlap_returns_empty_summary():
    rng = np.random.default_rng(1)
    samples = [f"S{i}" for i in range(10)]
    a = pd.DataFrame({"X": rng.normal(size=10)}, index=samples)
    b = pd.DataFrame({"Y": rng.normal(size=10)}, index=samples)
    corr = cross_omics_block_correlation(a, b, "a", "b")
    assert corr.n_feature_pairs == 0
    assert np.isnan(corr.corr_mean)

    # fewer than 3 shared samples -> empty
    b2 = b.copy()
    b2.index = [f"T{i}" for i in range(10)]
    corr2 = cross_omics_block_correlation(a, b2, "a", "b")
    assert corr2.n_shared_samples == 0
    assert corr2.n_feature_pairs == 0


def test_full_summary_structure_and_dict_roundtrip():
    rng = np.random.default_rng(13)
    samples = [f"S{i}" for i in range(25)]
    rna = _rna_block(rng, samples)
    prot = _protein_block(rng, samples)

    summary = multi_omics_descriptive_summary({"rna": rna, "protein": prot})
    assert isinstance(summary, MultiOmicsDescriptiveSummary)
    assert [b.name for b in summary.blocks] == ["protein", "rna"]
    assert len(summary.block_correlations) == 1
    assert summary.block_correlations[0].block_a == "protein"
    assert summary.block_correlations[0].block_b == "rna"

    payload = summary.to_dict()
    assert set(payload) == {"blocks", "block_correlations"}
    assert payload["blocks"][0]["name"] == "protein"
    pair = payload["block_correlations"][0]
    assert pair["pair"] == ["protein", "rna"]
    assert pair["n_shared_samples"] == 25


def test_empty_dict_raises():
    with pytest.raises(ValueError):
        multi_omics_descriptive_summary({})
