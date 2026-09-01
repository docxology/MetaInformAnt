"""Tests for tissue_specificity: tau, orthology classes, coupling summaries.

Zero-mocks: all tests use small synthetic expression frames and real
computation. Covers edge cases required by the round-3 lane brief:
all-zero gene, one-tissue, tie values, n=2 tissues.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.tissue_specificity import (
    classify_orthogroups,
    compute_tau,
    duplication_specificity_summary,
    filter_low_expression,
    join_expression_with_orthology,
    tau_summary,
    wilcoxon_duplication_specificity,
)


def _frame(data: dict, tissues=("brain", "ovary", "gut")) -> pd.DataFrame:
    return pd.DataFrame.from_dict(data, orient="index", columns=list(tissues))


class TestFilterLowExpression:
    def test_drops_all_zero_gene(self) -> None:
        df = _frame(
            {
                "g1": [10.0, 10.0, 10.0],
                "gzero": [0.0, 0.0, 0.0],
            }
        )
        out = filter_low_expression(df)
        assert "gzero" not in out.index
        assert "g1" in out.index

    def test_drops_lowest_ten_percent_by_mean(self) -> None:
        # 10 nonzero genes: means 1..10; lowest 10% of those (gene with
        # mean 1) must be dropped.
        data = {f"g{i}": [float(i)] * 3 for i in range(1, 11)}
        out = filter_low_expression(_frame(data))
        assert "g1" not in out.index
        assert len(out) == 9

    def test_boundary_ties_kept(self) -> None:
        # 10 genes where the boundary mean (cutoff) is tied across genes:
        # means 2,2,3,3,...  quantile(0.1) of means = 2 -> both mean-2 kept.
        data = {f"g{i}": [float(i)] * 3 for i in range(1, 11)}
        data["gA"] = [2.0] * 3
        data["gB"] = [2.0] * 3
        out = filter_low_expression(_frame(data))
        assert "gA" in out.index and "gB" in out.index

    def test_invalid_fraction_raises(self) -> None:
        df = _frame({"g1": [1.0, 1.0, 1.0]})
        with pytest.raises(ValueError):
            filter_low_expression(df, lowest_fraction=1.0)

    def test_empty_frame(self) -> None:
        out = filter_low_expression(pd.DataFrame())
        assert out.empty


class TestComputeTau:
    def test_perfectly_specific_gene(self) -> None:
        # Expressed only in one of 3 tissues -> tau = 1. Bypass the 10%
        # filter: with a 2-gene frame it would legitimately drop a gene.
        df = _frame({"g_specific": [0.0, 100.0, 0.0], "g_house": [50.0, 50.0, 50.0]})
        tau = compute_tau(df, lowest_fraction=0.0)
        assert tau["g_specific"] == pytest.approx(1.0)
        assert tau["g_house"] == pytest.approx(0.0)

    def test_manual_reference_two_tissues(self) -> None:
        # n=2: tau = (1 - x1/max) + (1 - x2/max)) / 1. With values 4 and 8,
        # log2(5)/log2(9) and log2(9)/log2(9).
        df = _frame({"g": [4.0, 8.0]}, tissues=("t1", "t2"))
        tau = compute_tau(df)
        expected = (1 - np.log2(5) / np.log2(9)) + (1 - 1.0)
        assert tau["g"] == pytest.approx(expected)

    def test_single_tissue_is_nan(self) -> None:
        df = pd.DataFrame({"g1": [10.0]}, index=["t1"])
        tau = compute_tau(df.T)  # genes as rows, one tissue column
        assert tau.shape[0] == 1
        assert np.isnan(tau.iloc[0])

    def test_all_zero_genes_excluded(self) -> None:
        df = _frame({"gzero": [0.0, 0.0, 0.0], "g1": [1.0, 0.0, 0.0]})
        tau = compute_tau(df)
        assert "gzero" not in tau.index
        assert tau["g1"] == pytest.approx(1.0)

    def test_tie_values_give_zero_tau(self) -> None:
        df = _frame({"g": [7.0, 7.0, 7.0]})
        tau = compute_tau(df)
        assert tau["g"] == pytest.approx(0.0)

    def test_nan_propagates(self) -> None:
        df = _frame({"gnan": [np.nan, 5.0, 5.0]})
        tau = compute_tau(df, lowest_fraction=0.0)
        assert np.isnan(tau["gnan"])

    def test_no_log2_changes_value(self) -> None:
        df = _frame({"g": [4.0, 8.0]}, tissues=("t1", "t2"))
        tau_lin = compute_tau(df, log2=False)
        assert tau_lin["g"] == pytest.approx(0.5)

    def test_empty_input(self) -> None:
        tau = compute_tau(pd.DataFrame())
        assert tau.empty


class TestTauSummary:
    def test_summary_statistics(self) -> None:
        tau = pd.Series([0.1, 0.5, 0.9, np.nan])
        s = tau_summary(tau)
        assert s["count"] == 4
        assert s["nan_count"] == 1
        assert s["valid_count"] == 3
        assert s["median"] == pytest.approx(0.5)

    def test_all_nan(self) -> None:
        s = tau_summary(pd.Series([np.nan, np.nan]))
        assert s["median"] is None
        assert s["nan_count"] == 2


def _bridge() -> pd.DataFrame:
    """Orthogroup x species bridge table in build_ortholog_bridge.py schema."""
    return pd.DataFrame(
        {
            "Bombus_terrestris": ["BT1", "BT2,BT3", "BT4", "", "BT5"],
            "Apis_mellifera": ["AM1", "AM2", "AM3,AM4", "", ""],
            "Nasonia_vitripennis": ["NV1", "NV2", "NV3", "NV4", ""],
        },
        index=["OG_single", "OG_multi", "OG_multi2", "OG_one_species", "OG_two_species"],
    )


class TestClassifyOrthogroups:
    def test_single_copy_vs_multicopy(self) -> None:
        classes = classify_orthogroups(_bridge())
        assert classes.loc["OG_single", "orthology_class"] == "single_copy"
        assert classes.loc["OG_multi", "orthology_class"] == "multicopy"
        assert classes.loc["OG_multi2", "orthology_class"] == "multicopy"
        assert classes.loc["OG_one_species", "orthology_class"] == "insufficient"
        # OG_two_species maps in only one species (BT5), so min_species=2
        # default makes it insufficient despite its name.
        assert classes.loc["OG_two_species", "orthology_class"] == "insufficient"

    def test_counts(self) -> None:
        classes = classify_orthogroups(_bridge())
        assert classes.loc["OG_multi", "n_species_mapped"] == 3
        assert classes.loc["OG_multi", "max_copies"] == 2
        assert classes.loc["OG_one_species", "n_species_mapped"] == 1

    def test_min_species_threshold(self) -> None:
        classes = classify_orthogroups(_bridge(), min_species=3)
        assert classes.loc["OG_one_species", "orthology_class"] == "insufficient"
        assert classes.loc["OG_two_species", "orthology_class"] == "insufficient"

    def test_empty_table(self) -> None:
        classes = classify_orthogroups(pd.DataFrame())
        assert classes.empty

    def test_invalid_min_species(self) -> None:
        with pytest.raises(ValueError):
            classify_orthogroups(_bridge(), min_species=0)


class TestJoinExpressionWithOrthology:
    def test_join_labels_and_passthrough(self) -> None:
        expr = _frame(
            {
                "BT1": [10.0, 1.0, 1.0],
                "BT2": [1.0, 9.0, 1.0],
                "BT_unmapped": [5.0, 5.0, 5.0],
            }
        )
        joined = join_expression_with_orthology(expr, _bridge(), "Bombus_terrestris")
        assert joined.loc["BT1", "orthology_class"] == "single_copy"
        assert joined.loc["BT2", "orthology_class"] == "multicopy"
        assert joined.loc["BT_unmapped", "orthology_class"] == "unmapped"
        assert joined.loc["BT1", "orthogroup"] == "OG_single"
        # Expression values pass through untouched.
        assert joined.loc["BT1", "brain"] == pytest.approx(10.0)

    def test_unknown_species_raises(self) -> None:
        with pytest.raises(ValueError):
            join_expression_with_orthology(_frame({"BT1": [1.0, 1.0, 1.0]}), _bridge(), "Nope")


class TestDuplicationSpecificitySummary:
    def test_descriptive_only(self) -> None:
        tau = pd.Series([0.2, 0.8, 0.4, 0.9], index=["a", "b", "c", "d"])
        cls = pd.Series(
            ["single_copy", "multicopy", "single_copy", "multicopy"],
            index=["a", "b", "c", "d"],
        )
        summary = duplication_specificity_summary(tau, cls)
        assert summary.loc["single_copy", "median"] == pytest.approx(0.3)
        assert summary.loc["multicopy", "median"] == pytest.approx(0.85)
        # No test-statistic columns in descriptive output.
        assert not any("p_value" in str(col) for col in summary.columns)

    def test_misaligned_index_raises(self) -> None:
        tau = pd.Series([0.2, 0.8])
        cls = pd.Series(["single_copy", "multicopy"], index=["x", "y"])
        with pytest.raises(ValueError):
            duplication_specificity_summary(tau, cls)


class TestGatedWilcoxon:
    def test_refuses_when_unfrozen(self) -> None:
        tau = pd.Series([0.2, 0.8, 0.4, 0.9])
        cls = pd.Series(["single_copy", "multicopy", "single_copy", "multicopy"])
        with pytest.raises(RuntimeError, match="post-freeze"):
            wilcoxon_duplication_specificity(tau, cls)

    def test_runs_when_gated_open(self) -> None:
        tau = pd.Series([0.1, 0.2, 0.8, 0.9])
        cls = pd.Series(["single_copy", "single_copy", "multicopy", "multicopy"])
        result = wilcoxon_duplication_specificity(tau, cls, evidence_manifest_frozen=True)
        assert result["gate"] == "post-freeze"
        assert 0.0 <= result["p_value"] <= 1.0

    def test_refuses_when_too_few_observations(self) -> None:
        tau = pd.Series([0.2, 0.8])
        cls = pd.Series(["single_copy", "multicopy"])
        with pytest.raises(RuntimeError, match="Insufficient"):
            wilcoxon_duplication_specificity(tau, cls, evidence_manifest_frozen=True)
