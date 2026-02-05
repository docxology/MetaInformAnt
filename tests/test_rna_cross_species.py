"""Comprehensive tests for the cross-species gene expression comparison module.

Tests all 7 public functions in metainformant.rna.analysis.cross_species:
    - build_ortholog_map
    - map_expression_to_orthologs
    - compute_expression_conservation
    - identify_divergent_genes
    - compare_expression_across_species
    - compute_expression_divergence_matrix
    - phylogenetic_expression_profile
    - cross_species_pca

All tests use real numpy/pandas/scipy data. No mocking. No external fixtures.
"""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.cross_species import (
    build_ortholog_map,
    compare_expression_across_species,
    compute_expression_conservation,
    compute_expression_divergence_matrix,
    cross_species_pca,
    identify_divergent_genes,
    map_expression_to_orthologs,
    phylogenetic_expression_profile,
)


# =============================================================================
# Fixtures: small synthetic data builders (inline, no external files)
# =============================================================================


def _make_ortholog_df(
    pairs: List[Tuple[str, str]],
    source_col: str = "human",
    target_col: str = "mouse",
) -> pd.DataFrame:
    """Build a small ortholog DataFrame from a list of (source, target) pairs."""
    return pd.DataFrame({source_col: [p[0] for p in pairs], target_col: [p[1] for p in pairs]})


def _make_expression_df(
    data: Dict[str, List[float]],
    gene_ids: List[str],
) -> pd.DataFrame:
    """Build an expression DataFrame (genes x samples) from dict and gene list."""
    return pd.DataFrame(data, index=gene_ids)


def _make_species_tree_simple() -> Dict[str, Any]:
    """A simple 2-leaf phylogenetic tree."""
    return {
        "name": "root",
        "distance": 0.0,
        "children": [
            {"name": "human", "distance": 0.1},
            {"name": "mouse", "distance": 0.15},
        ],
    }


def _make_species_tree_three() -> Dict[str, Any]:
    """A 3-leaf tree with one internal node."""
    return {
        "name": "root",
        "distance": 0.0,
        "children": [
            {
                "name": "mammals",
                "distance": 0.3,
                "children": [
                    {"name": "human", "distance": 0.1},
                    {"name": "mouse", "distance": 0.15},
                ],
            },
            {"name": "zebrafish", "distance": 0.5},
        ],
    }


# =============================================================================
# Tests: build_ortholog_map
# =============================================================================


class TestBuildOrthologMap:
    """Tests for build_ortholog_map."""

    def test_basic_one_to_one(self) -> None:
        """One-to-one mappings produce single-element lists."""
        df = _make_ortholog_df([("TP53", "Trp53"), ("BRCA1", "Brca1"), ("EGFR", "Egfr")])
        result = build_ortholog_map(df, "human", "mouse")

        assert isinstance(result, dict)
        assert len(result) == 3
        assert result["TP53"] == ["Trp53"]
        assert result["BRCA1"] == ["Brca1"]
        assert result["EGFR"] == ["Egfr"]

    def test_one_to_many(self) -> None:
        """One source gene maps to multiple target genes."""
        df = _make_ortholog_df([("EGFR", "Egfr"), ("EGFR", "Erbb1"), ("TP53", "Trp53")])
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 2
        assert sorted(result["EGFR"]) == ["Egfr", "Erbb1"]
        assert result["TP53"] == ["Trp53"]

    def test_many_to_one(self) -> None:
        """Multiple source genes can map to the same target gene independently."""
        df = _make_ortholog_df([("GeneA", "Ortho1"), ("GeneB", "Ortho1")])
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 2
        assert result["GeneA"] == ["Ortho1"]
        assert result["GeneB"] == ["Ortho1"]

    def test_deduplication(self) -> None:
        """Duplicate source-target pairs are deduplicated."""
        df = _make_ortholog_df([("TP53", "Trp53"), ("TP53", "Trp53"), ("TP53", "Trp53")])
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 1
        assert result["TP53"] == ["Trp53"]

    def test_nan_rows_dropped(self) -> None:
        """Rows with NaN values in key columns are dropped."""
        df = pd.DataFrame({"human": ["TP53", None, "EGFR"], "mouse": ["Trp53", "Brca1", None]})
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 1
        assert "TP53" in result

    def test_empty_df_raises(self) -> None:
        """Empty DataFrame raises ValueError."""
        df = pd.DataFrame({"human": [], "mouse": []})
        with pytest.raises(ValueError, match="cannot be empty"):
            build_ortholog_map(df, "human", "mouse")

    def test_missing_source_col_raises(self) -> None:
        """Missing source column raises ValueError."""
        df = _make_ortholog_df([("TP53", "Trp53")])
        with pytest.raises(ValueError, match="Source column.*not found"):
            build_ortholog_map(df, "nonexistent", "mouse")

    def test_missing_target_col_raises(self) -> None:
        """Missing target column raises ValueError."""
        df = _make_ortholog_df([("TP53", "Trp53")])
        with pytest.raises(ValueError, match="Target column.*not found"):
            build_ortholog_map(df, "human", "nonexistent")

    def test_all_nan_returns_empty(self) -> None:
        """DataFrame with all NaN values in key columns returns empty dict."""
        df = pd.DataFrame({"human": [None, None], "mouse": [None, None]})
        # The df is not empty (has rows) but after dropna it is
        # Still raises because the original df is not empty but valid becomes empty
        result = build_ortholog_map(df, "human", "mouse")
        assert result == {}

    def test_single_row(self) -> None:
        """Single-row ortholog table works correctly."""
        df = _make_ortholog_df([("GeneA", "OrthoA")])
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 1
        assert result["GeneA"] == ["OrthoA"]

    def test_extra_columns_ignored(self) -> None:
        """Additional columns in the DataFrame are ignored."""
        df = pd.DataFrame({
            "human": ["TP53", "BRCA1"],
            "mouse": ["Trp53", "Brca1"],
            "score": [0.95, 0.88],
            "type": ["1:1", "1:1"],
        })
        result = build_ortholog_map(df, "human", "mouse")

        assert len(result) == 2
        assert result["TP53"] == ["Trp53"]

    def test_integer_gene_ids_converted_to_strings(self) -> None:
        """Numeric gene IDs are converted to strings."""
        df = pd.DataFrame({"human": [100, 200], "mouse": [300, 400]})
        result = build_ortholog_map(df, "human", "mouse")

        assert "100" in result
        assert result["100"] == ["300"]


# =============================================================================
# Tests: map_expression_to_orthologs
# =============================================================================


class TestMapExpressionToOrthologs:
    """Tests for map_expression_to_orthologs."""

    def test_basic_one_to_one_mapping(self) -> None:
        """One-to-one ortholog mapping preserves values."""
        expr = _make_expression_df(
            {"s1": [10.0, 20.0, 30.0], "s2": [15.0, 25.0, 35.0]},
            gene_ids=["GeneA", "GeneB", "GeneC"],
        )
        orth_map = {"GeneA": ["Ortho1"], "GeneB": ["Ortho2"], "GeneC": ["Ortho3"]}

        result = map_expression_to_orthologs(expr, orth_map, aggregation="mean")

        assert isinstance(result, pd.DataFrame)
        assert set(result.index) == {"Ortho1", "Ortho2", "Ortho3"}
        assert list(result.columns) == ["s1", "s2"]
        assert result.loc["Ortho1", "s1"] == 10.0
        assert result.loc["Ortho2", "s2"] == 25.0

    def test_many_to_one_mean_aggregation(self) -> None:
        """Multiple source genes mapping to same target are averaged with mean."""
        expr = _make_expression_df(
            {"s1": [10.0, 20.0, 30.0], "s2": [15.0, 25.0, 35.0]},
            gene_ids=["GeneA", "GeneB", "GeneC"],
        )
        # GeneA and GeneB both map to Ortho1
        orth_map = {"GeneA": ["Ortho1"], "GeneB": ["Ortho1"], "GeneC": ["Ortho2"]}

        result = map_expression_to_orthologs(expr, orth_map, aggregation="mean")

        assert len(result) == 2
        # Ortho1 = mean(GeneA, GeneB)
        assert result.loc["Ortho1", "s1"] == pytest.approx(15.0)  # (10+20)/2
        assert result.loc["Ortho1", "s2"] == pytest.approx(20.0)  # (15+25)/2
        assert result.loc["Ortho2", "s1"] == pytest.approx(30.0)

    def test_many_to_one_max_aggregation(self) -> None:
        """Max aggregation takes the maximum across source genes."""
        expr = _make_expression_df(
            {"s1": [10.0, 20.0], "s2": [15.0, 25.0]},
            gene_ids=["GeneA", "GeneB"],
        )
        orth_map = {"GeneA": ["Ortho1"], "GeneB": ["Ortho1"]}

        result = map_expression_to_orthologs(expr, orth_map, aggregation="max")

        assert result.loc["Ortho1", "s1"] == pytest.approx(20.0)
        assert result.loc["Ortho1", "s2"] == pytest.approx(25.0)

    def test_many_to_one_sum_aggregation(self) -> None:
        """Sum aggregation adds expression across source genes."""
        expr = _make_expression_df(
            {"s1": [10.0, 20.0], "s2": [15.0, 25.0]},
            gene_ids=["GeneA", "GeneB"],
        )
        orth_map = {"GeneA": ["Ortho1"], "GeneB": ["Ortho1"]}

        result = map_expression_to_orthologs(expr, orth_map, aggregation="sum")

        assert result.loc["Ortho1", "s1"] == pytest.approx(30.0)
        assert result.loc["Ortho1", "s2"] == pytest.approx(40.0)

    def test_genes_not_in_map_excluded(self) -> None:
        """Source genes not present in ortholog_map are excluded."""
        expr = _make_expression_df(
            {"s1": [10.0, 20.0, 30.0]},
            gene_ids=["GeneA", "GeneB", "GeneC"],
        )
        orth_map = {"GeneA": ["Ortho1"]}  # Only GeneA has a mapping

        result = map_expression_to_orthologs(expr, orth_map)

        assert len(result) == 1
        assert "Ortho1" in result.index

    def test_map_genes_not_in_expression_excluded(self) -> None:
        """Ortholog map entries for genes not in expression_df are skipped."""
        expr = _make_expression_df({"s1": [10.0]}, gene_ids=["GeneA"])
        orth_map = {"GeneA": ["Ortho1"], "GeneZ": ["Ortho2"]}

        result = map_expression_to_orthologs(expr, orth_map)

        assert len(result) == 1
        assert "Ortho1" in result.index
        assert "Ortho2" not in result.index

    def test_no_matching_genes_returns_empty(self) -> None:
        """No matching genes between expression and map returns empty DataFrame."""
        expr = _make_expression_df({"s1": [10.0]}, gene_ids=["GeneA"])
        orth_map = {"GeneZ": ["Ortho1"]}

        result = map_expression_to_orthologs(expr, orth_map)

        assert result.empty
        assert list(result.columns) == ["s1"]

    def test_empty_expression_raises(self) -> None:
        """Empty expression DataFrame raises ValueError."""
        expr = pd.DataFrame()
        with pytest.raises(ValueError, match="cannot be empty"):
            map_expression_to_orthologs(expr, {"A": ["B"]})

    def test_unknown_aggregation_raises(self) -> None:
        """Unknown aggregation method raises ValueError."""
        expr = _make_expression_df({"s1": [10.0]}, gene_ids=["GeneA"])
        with pytest.raises(ValueError, match="Unknown aggregation"):
            map_expression_to_orthologs(expr, {"GeneA": ["Ortho1"]}, aggregation="median")  # type: ignore[arg-type]

    def test_one_to_many_ortholog_produces_multiple_targets(self) -> None:
        """A single source gene mapping to multiple targets produces multiple rows."""
        expr = _make_expression_df({"s1": [42.0]}, gene_ids=["GeneA"])
        orth_map = {"GeneA": ["OrthoX", "OrthoY"]}

        result = map_expression_to_orthologs(expr, orth_map)

        assert len(result) == 2
        assert set(result.index) == {"OrthoX", "OrthoY"}
        # Each target gets the same value (single source, no aggregation needed)
        assert result.loc["OrthoX", "s1"] == pytest.approx(42.0)
        assert result.loc["OrthoY", "s1"] == pytest.approx(42.0)

    def test_integer_expression_values(self) -> None:
        """Integer expression values are handled correctly."""
        expr = pd.DataFrame({"s1": [10, 20], "s2": [30, 40]}, index=["GeneA", "GeneB"])
        orth_map = {"GeneA": ["Ortho1"], "GeneB": ["Ortho2"]}

        result = map_expression_to_orthologs(expr, orth_map)

        assert result.loc["Ortho1", "s1"] == 10
        assert result.loc["Ortho2", "s2"] == 40

    def test_result_index_name(self) -> None:
        """Result DataFrame index is named 'gene'."""
        expr = _make_expression_df({"s1": [10.0]}, gene_ids=["GeneA"])
        orth_map = {"GeneA": ["Ortho1"]}

        result = map_expression_to_orthologs(expr, orth_map)

        assert result.index.name == "gene"


# =============================================================================
# Tests: compute_expression_conservation
# =============================================================================


class TestComputeExpressionConservation:
    """Tests for compute_expression_conservation."""

    def test_perfectly_correlated_spearman(self) -> None:
        """Perfectly correlated expression yields correlation ~1.0."""
        expr_a = pd.DataFrame(
            {"c1": [1.0, 5.0], "c2": [2.0, 6.0], "c3": [3.0, 7.0], "c4": [4.0, 8.0]},
            index=["gene1", "gene2"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10.0, 50.0], "c2": [20.0, 60.0], "c3": [30.0, 70.0], "c4": [40.0, 80.0]},
            index=["gene1", "gene2"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="spearman")

        assert isinstance(result, pd.DataFrame)
        assert set(result.columns) == {"gene_id", "correlation", "p_value", "conserved"}
        assert len(result) == 2

        for _, row in result.iterrows():
            assert row["correlation"] == pytest.approx(1.0, abs=1e-6)
            assert bool(row["conserved"]) is True

    def test_perfectly_correlated_pearson(self) -> None:
        """Perfectly linearly correlated data yields pearson ~1.0."""
        expr_a = pd.DataFrame(
            {"c1": [1.0], "c2": [2.0], "c3": [3.0], "c4": [4.0]},
            index=["gene1"],
        )
        expr_b = pd.DataFrame(
            {"c1": [2.0], "c2": [4.0], "c3": [6.0], "c4": [8.0]},
            index=["gene1"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="pearson")

        assert result.iloc[0]["correlation"] == pytest.approx(1.0, abs=1e-6)

    def test_anticorrelated_spearman(self) -> None:
        """Anti-correlated expression yields negative correlation."""
        expr_a = pd.DataFrame(
            {"c1": [1.0], "c2": [2.0], "c3": [3.0], "c4": [4.0]},
            index=["gene1"],
        )
        expr_b = pd.DataFrame(
            {"c1": [4.0], "c2": [3.0], "c3": [2.0], "c4": [1.0]},
            index=["gene1"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="spearman")

        assert result.iloc[0]["correlation"] == pytest.approx(-1.0, abs=1e-6)
        assert bool(result.iloc[0]["conserved"]) is False

    def test_euclidean_identical_normalized_values(self) -> None:
        """Identical normalized values produce euclidean similarity near 1.0."""
        expr_a = pd.DataFrame(
            {"c1": [0.0, 10.0], "c2": [5.0, 20.0], "c3": [10.0, 30.0]},
            index=["gene1", "gene2"],
        )
        # Same normalized pattern
        expr_b = pd.DataFrame(
            {"c1": [0.0, 100.0], "c2": [5.0, 200.0], "c3": [10.0, 300.0]},
            index=["gene1", "gene2"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="euclidean")

        assert len(result) == 2
        for _, row in result.iterrows():
            # Both genes have identical normalized patterns
            assert row["correlation"] == pytest.approx(1.0, abs=0.05)

    def test_no_shared_genes_raises(self) -> None:
        """No shared genes between the two DataFrames raises ValueError."""
        expr_a = pd.DataFrame({"c1": [1.0]}, index=["geneA"])
        expr_b = pd.DataFrame({"c1": [2.0]}, index=["geneB"])

        with pytest.raises(ValueError, match="No shared genes"):
            compute_expression_conservation(expr_a, expr_b)

    def test_unknown_method_raises(self) -> None:
        """Unknown correlation method raises ValueError."""
        expr_a = pd.DataFrame({"c1": [1.0]}, index=["gene1"])
        expr_b = pd.DataFrame({"c1": [2.0]}, index=["gene1"])

        with pytest.raises(ValueError, match="Unknown method"):
            compute_expression_conservation(expr_a, expr_b, method="kendall")  # type: ignore[arg-type]

    def test_single_sample_per_gene(self) -> None:
        """Single sample (< 2 data points) produces NaN correlation."""
        expr_a = pd.DataFrame({"c1": [1.0, 2.0]}, index=["gene1", "gene2"])
        expr_b = pd.DataFrame({"c1": [3.0, 4.0]}, index=["gene1", "gene2"])

        result = compute_expression_conservation(expr_a, expr_b, method="spearman")

        assert len(result) == 2
        for _, row in result.iterrows():
            assert np.isnan(row["correlation"])
            assert bool(row["conserved"]) is False

    def test_mismatched_sample_counts_truncated(self) -> None:
        """Different sample counts are handled by truncating to minimum."""
        expr_a = pd.DataFrame(
            {"c1": [1.0], "c2": [2.0], "c3": [3.0], "c4": [4.0]},
            index=["gene1"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10.0], "c2": [20.0]},
            index=["gene1"],
        )

        # Should not raise; uses min_len = 2
        result = compute_expression_conservation(expr_a, expr_b, method="pearson")

        assert len(result) == 1
        # With only 2 points [1,2] vs [10,20] Pearson = 1.0
        assert result.iloc[0]["correlation"] == pytest.approx(1.0, abs=1e-6)

    def test_partial_overlap_genes(self) -> None:
        """Only shared genes are compared, non-shared genes are excluded."""
        expr_a = pd.DataFrame(
            {"c1": [1.0, 2.0, 3.0], "c2": [4.0, 5.0, 6.0]},
            index=["gene1", "gene2", "gene3"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10.0, 20.0], "c2": [40.0, 50.0]},
            index=["gene2", "gene3"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="spearman")

        assert len(result) == 2
        gene_ids = set(result["gene_id"])
        assert gene_ids == {"gene2", "gene3"}

    def test_constant_expression_handling(self) -> None:
        """Constant expression across samples produces NaN correlation handled gracefully."""
        import warnings

        expr_a = pd.DataFrame(
            {"c1": [5.0], "c2": [5.0], "c3": [5.0]},
            index=["gene1"],
        )
        expr_b = pd.DataFrame(
            {"c1": [5.0], "c2": [5.0], "c3": [5.0]},
            index=["gene1"],
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = compute_expression_conservation(expr_a, expr_b, method="spearman")

        # Spearman of constant vectors is undefined; code handles NaN -> 0.0
        assert len(result) == 1
        # The module replaces NaN correlation with 0.0
        assert result.iloc[0]["correlation"] == pytest.approx(0.0, abs=1e-6)

    def test_result_types(self) -> None:
        """Result DataFrame has correct column dtypes."""
        expr_a = pd.DataFrame(
            {"c1": [1.0, 5.0], "c2": [2.0, 6.0], "c3": [3.0, 7.0]},
            index=["gene1", "gene2"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10.0, 50.0], "c2": [20.0, 60.0], "c3": [30.0, 70.0]},
            index=["gene1", "gene2"],
        )

        result = compute_expression_conservation(expr_a, expr_b)

        assert pd.api.types.is_string_dtype(result["gene_id"])
        assert pd.api.types.is_float_dtype(result["correlation"])
        assert pd.api.types.is_float_dtype(result["p_value"])
        assert pd.api.types.is_bool_dtype(result["conserved"])


# =============================================================================
# Tests: identify_divergent_genes
# =============================================================================


class TestIdentifyDivergentGenes:
    """Tests for identify_divergent_genes."""

    def test_basic_filtering(self) -> None:
        """Genes below threshold are identified as divergent."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3", "gene4"],
            "correlation": [0.9, 0.5, 0.2, -0.3],
            "p_value": [0.001, 0.01, 0.1, 0.5],
            "conserved": [True, True, False, False],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert len(result) == 2
        assert set(result["gene_id"]) == {"gene3", "gene4"}
        # Sorted ascending: most divergent first
        assert result.iloc[0]["gene_id"] == "gene4"
        assert result.iloc[1]["gene_id"] == "gene3"

    def test_all_conserved(self) -> None:
        """No divergent genes when all are above threshold."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2"],
            "correlation": [0.8, 0.9],
            "p_value": [0.001, 0.001],
            "conserved": [True, True],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert len(result) == 0

    def test_all_divergent(self) -> None:
        """All genes divergent when all are below threshold."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2"],
            "correlation": [0.1, -0.2],
            "p_value": [0.5, 0.8],
            "conserved": [False, False],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert len(result) == 2

    def test_nan_correlation_treated_as_divergent(self) -> None:
        """NaN correlations are treated as divergent (filled with -1.0)."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2"],
            "correlation": [np.nan, 0.9],
            "p_value": [np.nan, 0.001],
            "conserved": [False, True],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert len(result) == 1
        assert result.iloc[0]["gene_id"] == "gene1"

    def test_threshold_at_boundary(self) -> None:
        """Genes exactly at threshold are included (<=)."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2"],
            "correlation": [0.3, 0.31],
            "p_value": [0.05, 0.05],
            "conserved": [False, False],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert len(result) == 1
        assert result.iloc[0]["gene_id"] == "gene1"

    def test_empty_conservation_df(self) -> None:
        """Empty input DataFrame returns empty DataFrame."""
        conservation_df = pd.DataFrame(columns=["gene_id", "correlation", "p_value", "conserved"])

        result = identify_divergent_genes(conservation_df)

        assert result.empty

    def test_missing_required_columns_raises(self) -> None:
        """Missing required columns raises ValueError."""
        df = pd.DataFrame({"gene_id": ["gene1"], "p_value": [0.05]})
        with pytest.raises(ValueError, match="missing required columns"):
            identify_divergent_genes(df)

    def test_custom_threshold(self) -> None:
        """Custom threshold values work correctly."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3"],
            "correlation": [0.1, 0.5, 0.8],
        })

        # High threshold captures more genes
        result_high = identify_divergent_genes(conservation_df, threshold=0.6)
        assert len(result_high) == 2

        # Low threshold captures fewer genes
        result_low = identify_divergent_genes(conservation_df, threshold=0.05)
        assert len(result_low) == 0

    def test_result_sorted_ascending(self) -> None:
        """Result is sorted by correlation ascending."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3"],
            "correlation": [0.2, -0.5, 0.0],
            "p_value": [0.1, 0.5, 0.3],
            "conserved": [False, False, False],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert list(result["gene_id"]) == ["gene2", "gene3", "gene1"]

    def test_result_index_reset(self) -> None:
        """Result index is reset to 0-based sequential."""
        conservation_df = pd.DataFrame({
            "gene_id": ["gene1", "gene2", "gene3"],
            "correlation": [0.1, 0.2, 0.0],
        })

        result = identify_divergent_genes(conservation_df, threshold=0.3)

        assert list(result.index) == [0, 1, 2]


# =============================================================================
# Tests: compute_expression_divergence_matrix
# =============================================================================


class TestComputeExpressionDivergenceMatrix:
    """Tests for compute_expression_divergence_matrix."""

    def test_basic_two_species(self) -> None:
        """Two species produce a 2x2 symmetric matrix with zero diagonal."""
        # Perfectly correlated: divergence should be low (~0)
        species_data = {
            "human": pd.DataFrame(
                {"s1": [1.0, 2.0, 3.0], "s2": [4.0, 5.0, 6.0]},
                index=["gene1", "gene2", "gene3"],
            ),
            "mouse": pd.DataFrame(
                {"s1": [10.0, 20.0, 30.0], "s2": [40.0, 50.0, 60.0]},
                index=["gene1", "gene2", "gene3"],
            ),
        }

        result = compute_expression_divergence_matrix(species_data)

        assert isinstance(result, pd.DataFrame)
        assert result.shape == (2, 2)
        assert sorted(result.index) == ["human", "mouse"]
        assert sorted(result.columns) == ["human", "mouse"]
        # Diagonal should be 0
        assert result.loc["human", "human"] == pytest.approx(0.0)
        assert result.loc["mouse", "mouse"] == pytest.approx(0.0)
        # Symmetric
        assert result.loc["human", "mouse"] == pytest.approx(result.loc["mouse", "human"])

    def test_three_species_symmetry(self) -> None:
        """Three species produce a 3x3 symmetric matrix."""
        rng = np.random.default_rng(123)
        genes = [f"gene{i}" for i in range(10)]

        species_data = {
            "sp_a": pd.DataFrame(rng.random((10, 3)), index=genes, columns=["s1", "s2", "s3"]),
            "sp_b": pd.DataFrame(rng.random((10, 3)), index=genes, columns=["s1", "s2", "s3"]),
            "sp_c": pd.DataFrame(rng.random((10, 3)), index=genes, columns=["s1", "s2", "s3"]),
        }

        result = compute_expression_divergence_matrix(species_data)

        assert result.shape == (3, 3)
        # Check symmetry
        for i in result.index:
            for j in result.columns:
                assert result.loc[i, j] == pytest.approx(result.loc[j, i], abs=1e-10)
        # Diagonal is 0
        for sp in result.index:
            assert result.loc[sp, sp] == pytest.approx(0.0)

    def test_identical_species_zero_divergence(self) -> None:
        """Identical expression data yields near-zero divergence."""
        expr = pd.DataFrame(
            {"s1": [1.0, 2.0, 3.0], "s2": [4.0, 5.0, 6.0], "s3": [7.0, 8.0, 9.0]},
            index=["gene1", "gene2", "gene3"],
        )
        species_data = {"sp_a": expr.copy(), "sp_b": expr.copy()}

        result = compute_expression_divergence_matrix(species_data)

        # 1 - median(1.0) = 0.0
        assert result.loc["sp_a", "sp_b"] == pytest.approx(0.0, abs=0.05)

    def test_single_sample_uses_gene_level_corr(self) -> None:
        """Single sample per species falls back to gene-level correlation."""
        species_data = {
            "sp_a": pd.DataFrame({"s1": [1.0, 2.0, 3.0, 4.0]}, index=["g1", "g2", "g3", "g4"]),
            "sp_b": pd.DataFrame({"s1": [10.0, 20.0, 30.0, 40.0]}, index=["g1", "g2", "g3", "g4"]),
        }

        result = compute_expression_divergence_matrix(species_data)

        # Gene-level spearman of [1,2,3,4] vs [10,20,30,40] = 1.0
        # divergence = 1 - 1.0 = 0.0
        assert result.loc["sp_a", "sp_b"] == pytest.approx(0.0, abs=0.05)

    def test_no_shared_genes_max_divergence(self) -> None:
        """No shared genes yields maximum divergence of 2.0."""
        species_data = {
            "sp_a": pd.DataFrame({"s1": [1.0]}, index=["geneA"]),
            "sp_b": pd.DataFrame({"s1": [2.0]}, index=["geneB"]),
        }

        result = compute_expression_divergence_matrix(species_data)

        assert result.loc["sp_a", "sp_b"] == pytest.approx(2.0)

    def test_fewer_than_two_species_raises(self) -> None:
        """Fewer than 2 species raises ValueError."""
        with pytest.raises(ValueError, match="at least 2 species"):
            compute_expression_divergence_matrix({"only_one": pd.DataFrame({"s1": [1.0]}, index=["g1"])})

    def test_divergence_range(self) -> None:
        """Divergence values are in expected range [0, 2]."""
        rng = np.random.default_rng(42)
        genes = [f"gene{i}" for i in range(20)]

        species_data = {
            "sp_a": pd.DataFrame(rng.random((20, 4)), index=genes, columns=["s1", "s2", "s3", "s4"]),
            "sp_b": pd.DataFrame(rng.random((20, 4)), index=genes, columns=["s1", "s2", "s3", "s4"]),
        }

        result = compute_expression_divergence_matrix(species_data)

        off_diag = result.loc["sp_a", "sp_b"]
        assert 0.0 <= off_diag <= 2.0


# =============================================================================
# Tests: phylogenetic_expression_profile
# =============================================================================


class TestPhylogeneticExpressionProfile:
    """Tests for phylogenetic_expression_profile."""

    def test_basic_two_leaf_tree(self) -> None:
        """Simple 2-leaf tree produces correct node statistics."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame(
            {"human": [10.0, 20.0], "mouse": [12.0, 18.0]},
            index=["gene1", "gene2"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        assert isinstance(result, pd.DataFrame)
        assert set(result.columns) == {
            "gene_id", "node", "mean_expression", "variance",
            "n_leaves", "branch_distance", "expression_change",
        }

        # Should have entries for: human, mouse, root (for each gene)
        nodes = set(result["node"])
        assert "human" in nodes
        assert "mouse" in nodes
        assert "root" in nodes

        # Leaf nodes have n_leaves=1
        human_rows = result[result["node"] == "human"]
        assert all(human_rows["n_leaves"] == 1)

        # Root node has n_leaves=2
        root_rows = result[result["node"] == "root"]
        assert all(root_rows["n_leaves"] == 2)

    def test_root_mean_is_average_of_leaves(self) -> None:
        """Root mean expression is the average of leaf values."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame(
            {"human": [10.0], "mouse": [20.0]},
            index=["gene1"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        root_row = result[(result["node"] == "root") & (result["gene_id"] == "gene1")]
        assert len(root_row) == 1
        assert root_row.iloc[0]["mean_expression"] == pytest.approx(15.0)  # (10+20)/2

    def test_root_variance(self) -> None:
        """Root variance is population variance across leaf values."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame(
            {"human": [10.0], "mouse": [20.0]},
            index=["gene1"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        root_row = result[(result["node"] == "root") & (result["gene_id"] == "gene1")]
        # Population variance of [10, 20]: mean=15, var = ((10-15)^2 + (20-15)^2)/2 = 25
        assert root_row.iloc[0]["variance"] == pytest.approx(25.0)

    def test_leaf_variance_is_zero(self) -> None:
        """Leaf nodes have zero variance (single value)."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame(
            {"human": [10.0], "mouse": [20.0]},
            index=["gene1"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        leaf_rows = result[result["node"].isin(["human", "mouse"])]
        assert all(leaf_rows["variance"] == 0.0)

    def test_branch_distances(self) -> None:
        """Branch distances from tree are correctly recorded."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame({"human": [10.0], "mouse": [20.0]}, index=["gene1"])

        result = phylogenetic_expression_profile(expr, tree)

        human_row = result[(result["node"] == "human") & (result["gene_id"] == "gene1")]
        mouse_row = result[(result["node"] == "mouse") & (result["gene_id"] == "gene1")]
        root_row = result[(result["node"] == "root") & (result["gene_id"] == "gene1")]

        assert human_row.iloc[0]["branch_distance"] == pytest.approx(0.1)
        assert mouse_row.iloc[0]["branch_distance"] == pytest.approx(0.15)
        assert root_row.iloc[0]["branch_distance"] == pytest.approx(0.0)

    def test_three_leaf_tree(self) -> None:
        """Three-leaf tree with internal node produces correct structure."""
        tree = _make_species_tree_three()
        expr = pd.DataFrame(
            {"human": [10.0], "mouse": [20.0], "zebrafish": [30.0]},
            index=["gene1"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        nodes = set(result["node"])
        assert "human" in nodes
        assert "mouse" in nodes
        assert "zebrafish" in nodes
        assert "mammals" in nodes
        assert "root" in nodes

        # Mammals internal node should aggregate human + mouse
        mammal_row = result[(result["node"] == "mammals") & (result["gene_id"] == "gene1")]
        assert mammal_row.iloc[0]["mean_expression"] == pytest.approx(15.0)
        assert mammal_row.iloc[0]["n_leaves"] == 2

        # Root should aggregate all three
        root_row = result[(result["node"] == "root") & (result["gene_id"] == "gene1")]
        assert root_row.iloc[0]["mean_expression"] == pytest.approx(20.0)  # (10+20+30)/3
        assert root_row.iloc[0]["n_leaves"] == 3

    def test_empty_expression_raises(self) -> None:
        """Empty expression DataFrame raises ValueError."""
        tree = _make_species_tree_simple()
        with pytest.raises(ValueError, match="cannot be empty"):
            phylogenetic_expression_profile(pd.DataFrame(), tree)

    def test_malformed_tree_raises(self) -> None:
        """Tree missing 'name' key raises ValueError."""
        expr = pd.DataFrame({"human": [10.0]}, index=["gene1"])
        with pytest.raises(ValueError, match="must contain a 'name' key"):
            phylogenetic_expression_profile(expr, {"children": []})

    def test_missing_species_in_expression(self) -> None:
        """Leaf nodes not found in expression data are skipped gracefully."""
        tree = _make_species_tree_simple()
        # Only provide data for 'human', not 'mouse'
        expr = pd.DataFrame({"human": [10.0]}, index=["gene1"])

        result = phylogenetic_expression_profile(expr, tree)

        # Should still have results for human and root
        nodes = set(result["node"])
        assert "human" in nodes
        assert "root" in nodes

    def test_multiple_genes(self) -> None:
        """Multiple genes produce distinct per-gene statistics."""
        tree = _make_species_tree_simple()
        expr = pd.DataFrame(
            {"human": [10.0, 100.0], "mouse": [20.0, 200.0]},
            index=["gene1", "gene2"],
        )

        result = phylogenetic_expression_profile(expr, tree)

        gene1_root = result[(result["node"] == "root") & (result["gene_id"] == "gene1")]
        gene2_root = result[(result["node"] == "root") & (result["gene_id"] == "gene2")]

        assert gene1_root.iloc[0]["mean_expression"] == pytest.approx(15.0)
        assert gene2_root.iloc[0]["mean_expression"] == pytest.approx(150.0)


# =============================================================================
# Tests: cross_species_pca
# =============================================================================


class TestCrossSpeciesPCA:
    """Tests for cross_species_pca."""

    def test_basic_two_species(self) -> None:
        """Basic PCA with two species produces expected output structure."""
        rng = np.random.default_rng(42)
        genes = [f"gene{i}" for i in range(10)]

        species_data = {
            "human": pd.DataFrame(
                rng.random((10, 3)), index=genes, columns=["h1", "h2", "h3"]
            ),
            "mouse": pd.DataFrame(
                rng.random((10, 3)), index=genes, columns=["m1", "m2", "m3"]
            ),
        }

        result = cross_species_pca(species_data, n_components=2)

        assert isinstance(result, dict)
        assert set(result.keys()) == {
            "coordinates", "explained_variance", "loadings",
            "species_labels", "shared_genes",
        }

        coords = result["coordinates"]
        assert isinstance(coords, pd.DataFrame)
        assert "PC1" in coords.columns
        assert "PC2" in coords.columns
        assert "species" in coords.columns
        # 3 human + 3 mouse = 6 samples
        assert len(coords) == 6

        loadings = result["loadings"]
        assert isinstance(loadings, pd.DataFrame)
        assert loadings.shape == (10, 2)  # 10 genes, 2 components

        ev = result["explained_variance"]
        assert len(ev) == 2
        assert all(v >= 0 for v in ev)

        assert result["shared_genes"] == genes
        assert len(result["species_labels"]) == 6

    def test_species_labels_correct(self) -> None:
        """Species labels correctly identify which samples belong to which species."""
        rng = np.random.default_rng(101)
        genes = ["gene1", "gene2", "gene3"]
        species_data = {
            "human": pd.DataFrame(
                rng.random((3, 2)), index=genes, columns=["h1", "h2"]
            ),
            "mouse": pd.DataFrame(
                rng.random((3, 2)), index=genes, columns=["m1", "m2"]
            ),
        }

        result = cross_species_pca(species_data, n_components=2)

        species_col = result["coordinates"]["species"].tolist()
        assert species_col.count("human") == 2
        assert species_col.count("mouse") == 2

    def test_sample_names_prefixed(self) -> None:
        """Sample names are prefixed with species name to avoid collisions."""
        rng = np.random.default_rng(202)
        genes = ["gene1", "gene2"]
        species_data = {
            "human": pd.DataFrame(rng.random((2, 2)), index=genes, columns=["s1", "s2"]),
            "mouse": pd.DataFrame(rng.random((2, 2)), index=genes, columns=["s1", "s2"]),
        }

        result = cross_species_pca(species_data, n_components=1)

        sample_names = result["coordinates"].index.tolist()
        assert "human:s1" in sample_names
        assert "mouse:s1" in sample_names

    def test_n_components_reduced_if_too_large(self) -> None:
        """n_components is clamped to min(n_samples, n_features)."""
        genes = ["gene1", "gene2"]
        species_data = {
            "sp_a": pd.DataFrame(np.array([[1.0], [2.0]]), index=genes, columns=["s1"]),
            "sp_b": pd.DataFrame(np.array([[3.0], [4.0]]), index=genes, columns=["s1"]),
        }

        # Request 10 components but only 2 samples and 2 features
        result = cross_species_pca(species_data, n_components=10)

        coords = result["coordinates"]
        # Should be clamped: min(2 samples, 2 genes) = 2 components
        pc_cols = [c for c in coords.columns if c.startswith("PC")]
        assert len(pc_cols) == 2

    def test_no_shared_genes_raises(self) -> None:
        """No shared genes raises ValueError."""
        species_data = {
            "sp_a": pd.DataFrame({"s1": [1.0]}, index=["geneA"]),
            "sp_b": pd.DataFrame({"s1": [1.0]}, index=["geneB"]),
        }

        with pytest.raises(ValueError, match="No genes shared"):
            cross_species_pca(species_data)

    def test_fewer_than_two_species_raises(self) -> None:
        """Fewer than 2 species raises ValueError."""
        with pytest.raises(ValueError, match="at least 2 species"):
            cross_species_pca({"only_one": pd.DataFrame({"s1": [1.0]}, index=["g1"])})

    def test_explained_variance_sums_to_at_most_one(self) -> None:
        """Sum of explained variance ratios should not exceed 1.0."""
        rng = np.random.default_rng(99)
        genes = [f"gene{i}" for i in range(20)]

        species_data = {
            "sp_a": pd.DataFrame(rng.random((20, 5)), index=genes, columns=[f"a{i}" for i in range(5)]),
            "sp_b": pd.DataFrame(rng.random((20, 5)), index=genes, columns=[f"b{i}" for i in range(5)]),
        }

        result = cross_species_pca(species_data, n_components=5)

        ev = result["explained_variance"]
        assert np.sum(ev) <= 1.0 + 1e-6

    def test_three_species(self) -> None:
        """Three species PCA works correctly."""
        rng = np.random.default_rng(77)
        genes = [f"gene{i}" for i in range(8)]

        species_data = {
            "sp_a": pd.DataFrame(rng.random((8, 3)), index=genes, columns=["a1", "a2", "a3"]),
            "sp_b": pd.DataFrame(rng.random((8, 3)), index=genes, columns=["b1", "b2", "b3"]),
            "sp_c": pd.DataFrame(rng.random((8, 3)), index=genes, columns=["c1", "c2", "c3"]),
        }

        result = cross_species_pca(species_data, n_components=3)

        coords = result["coordinates"]
        assert len(coords) == 9  # 3 + 3 + 3 samples
        assert set(result["species_labels"]) == {"sp_a", "sp_b", "sp_c"}

    def test_nan_handling(self) -> None:
        """NaN values in expression are imputed with column means."""
        genes = ["gene1", "gene2", "gene3"]
        species_data = {
            "sp_a": pd.DataFrame(
                {"s1": [1.0, 2.0, np.nan], "s2": [4.0, np.nan, 6.0]},
                index=genes,
            ),
            "sp_b": pd.DataFrame(
                {"s1": [7.0, 8.0, 9.0], "s2": [10.0, 11.0, 12.0]},
                index=genes,
            ),
        }

        # Should not raise
        result = cross_species_pca(species_data, n_components=2)

        assert not result["coordinates"].empty
        assert len(result["coordinates"]) == 4


# =============================================================================
# Tests: compare_expression_across_species
# =============================================================================


class TestCompareExpressionAcrossSpecies:
    """Tests for compare_expression_across_species."""

    def test_basic_two_species(self) -> None:
        """Basic comparison with two species and a simple ortholog map."""
        # Species A expression (in species A gene space)
        expr_a = pd.DataFrame(
            {
                "s1": [1.0, 2.0, 3.0],
                "s2": [4.0, 5.0, 6.0],
                "s3": [7.0, 8.0, 9.0],
                "s4": [10.0, 11.0, 12.0],
            },
            index=["gA1", "gA2", "gA3"],
        )
        # Species B expression (in ortholog gene space)
        expr_b = pd.DataFrame(
            {
                "s1": [10.0, 20.0, 30.0],
                "s2": [40.0, 50.0, 60.0],
                "s3": [70.0, 80.0, 90.0],
                "s4": [100.0, 110.0, 120.0],
            },
            index=["oG1", "oG2", "oG3"],
        )

        orth_maps: Dict[Tuple[str, str], Dict[str, List[str]]] = {
            ("sp_a", "sp_b"): {"gA1": ["oG1"], "gA2": ["oG2"], "gA3": ["oG3"]},
        }

        species_expressions = {"sp_a": expr_a, "sp_b": expr_b}

        result = compare_expression_across_species(species_expressions, orth_maps)

        assert isinstance(result, pd.DataFrame)
        assert "gene_id" in result.columns
        assert "mean_conservation" in result.columns
        assert "min_conservation" in result.columns
        assert "max_conservation" in result.columns
        assert "n_species_compared" in result.columns
        assert "conserved_in_all" in result.columns

    def test_fewer_than_two_species_raises(self) -> None:
        """Fewer than 2 species raises ValueError."""
        with pytest.raises(ValueError, match="at least 2 species"):
            compare_expression_across_species(
                {"only_one": pd.DataFrame({"s1": [1.0]}, index=["g1"])},
                {},
            )

    def test_missing_species_pair_skipped(self) -> None:
        """Ortholog maps referencing missing species are skipped gracefully."""
        expr_a = pd.DataFrame(
            {"s1": [1.0], "s2": [2.0], "s3": [3.0]},
            index=["gA1"],
        )
        expr_b = pd.DataFrame(
            {"s1": [10.0], "s2": [20.0], "s3": [30.0]},
            index=["oG1"],
        )

        # Reference a species "sp_c" not in species_expressions
        orth_maps: Dict[Tuple[str, str], Dict[str, List[str]]] = {
            ("sp_a", "sp_c"): {"gA1": ["oG1"]},
            ("sp_a", "sp_b"): {"gA1": ["oG1"]},
        }

        species_expressions = {"sp_a": expr_a, "sp_b": expr_b}

        # Should not raise; sp_c pair is skipped
        result = compare_expression_across_species(species_expressions, orth_maps)

        assert isinstance(result, pd.DataFrame)

    def test_no_matching_genes_returns_empty(self) -> None:
        """No matching ortholog genes returns empty DataFrame."""
        expr_a = pd.DataFrame(
            {"s1": [1.0], "s2": [2.0], "s3": [3.0]},
            index=["gA1"],
        )
        expr_b = pd.DataFrame(
            {"s1": [10.0], "s2": [20.0], "s3": [30.0]},
            index=["oG1"],
        )

        # Map to a gene not in expr_b
        orth_maps: Dict[Tuple[str, str], Dict[str, List[str]]] = {
            ("sp_a", "sp_b"): {"gA1": ["not_in_b"]},
        }

        result = compare_expression_across_species(
            {"sp_a": expr_a, "sp_b": expr_b}, orth_maps
        )

        assert result.empty

    def test_result_sorted_by_conservation(self) -> None:
        """Results are sorted by mean_conservation descending."""
        expr_a = pd.DataFrame(
            {
                "s1": [1.0, 100.0],
                "s2": [2.0, 50.0],
                "s3": [3.0, 25.0],
                "s4": [4.0, 12.5],
            },
            index=["gA1", "gA2"],
        )
        # Perfect correlation for gA1->oG1, weaker for gA2->oG2
        expr_b = pd.DataFrame(
            {
                "s1": [10.0, 1.0],
                "s2": [20.0, 2.0],
                "s3": [30.0, 3.0],
                "s4": [40.0, 4.0],
            },
            index=["oG1", "oG2"],
        )

        orth_maps: Dict[Tuple[str, str], Dict[str, List[str]]] = {
            ("sp_a", "sp_b"): {"gA1": ["oG1"], "gA2": ["oG2"]},
        }

        result = compare_expression_across_species({"sp_a": expr_a, "sp_b": expr_b}, orth_maps)

        if len(result) >= 2:
            assert result.iloc[0]["mean_conservation"] >= result.iloc[1]["mean_conservation"]

    def test_n_species_compared_count(self) -> None:
        """n_species_compared reflects the number of pairwise comparisons per gene."""
        genes = ["gA1", "gA2"]
        expr_a = pd.DataFrame(
            {"s1": [1.0, 2.0], "s2": [3.0, 4.0], "s3": [5.0, 6.0]},
            index=genes,
        )
        expr_b = pd.DataFrame(
            {"s1": [10.0, 20.0], "s2": [30.0, 40.0], "s3": [50.0, 60.0]},
            index=["oG1", "oG2"],
        )

        orth_maps: Dict[Tuple[str, str], Dict[str, List[str]]] = {
            ("sp_a", "sp_b"): {"gA1": ["oG1"], "gA2": ["oG2"]},
        }

        result = compare_expression_across_species({"sp_a": expr_a, "sp_b": expr_b}, orth_maps)

        for _, row in result.iterrows():
            assert row["n_species_compared"] >= 1


# =============================================================================
# Integration tests combining multiple functions
# =============================================================================


class TestIntegrationWorkflows:
    """Integration tests combining multiple cross-species functions."""

    def test_full_pipeline_build_map_then_compare(self) -> None:
        """End-to-end: build ortholog map, map expression, compute conservation."""
        # Step 1: Build ortholog map
        orthologs_df = pd.DataFrame({
            "human_gene": ["TP53", "BRCA1", "EGFR", "MYC"],
            "mouse_gene": ["Trp53", "Brca1", "Egfr", "Myc"],
        })
        orth_map = build_ortholog_map(orthologs_df, "human_gene", "mouse_gene")

        assert len(orth_map) == 4

        # Step 2: Create expression data
        human_expr = pd.DataFrame(
            {
                "cond1": [10.0, 20.0, 30.0, 40.0],
                "cond2": [15.0, 25.0, 35.0, 45.0],
                "cond3": [20.0, 30.0, 40.0, 50.0],
            },
            index=["TP53", "BRCA1", "EGFR", "MYC"],
        )
        mouse_expr = pd.DataFrame(
            {
                "cond1": [100.0, 200.0, 300.0, 400.0],
                "cond2": [150.0, 250.0, 350.0, 450.0],
                "cond3": [200.0, 300.0, 400.0, 500.0],
            },
            index=["Trp53", "Brca1", "Egfr", "Myc"],
        )

        # Step 3: Map human expression to ortholog space
        mapped_human = map_expression_to_orthologs(human_expr, orth_map)

        assert len(mapped_human) == 4
        assert set(mapped_human.index) == {"Trp53", "Brca1", "Egfr", "Myc"}

        # Step 4: Compute conservation
        conservation = compute_expression_conservation(mapped_human, mouse_expr, method="spearman")

        assert len(conservation) == 4
        assert "gene_id" in conservation.columns
        assert "correlation" in conservation.columns

        # Step 5: Identify divergent genes
        divergent = identify_divergent_genes(conservation, threshold=0.3)

        assert isinstance(divergent, pd.DataFrame)

    def test_divergence_matrix_then_pca(self) -> None:
        """Build divergence matrix and PCA from same species data."""
        rng = np.random.default_rng(42)
        genes = [f"gene{i}" for i in range(15)]

        species_data = {
            "sp_a": pd.DataFrame(
                rng.random((15, 4)), index=genes, columns=["a1", "a2", "a3", "a4"]
            ),
            "sp_b": pd.DataFrame(
                rng.random((15, 4)), index=genes, columns=["b1", "b2", "b3", "b4"]
            ),
            "sp_c": pd.DataFrame(
                rng.random((15, 4)), index=genes, columns=["c1", "c2", "c3", "c4"]
            ),
        }

        # Divergence matrix
        div_matrix = compute_expression_divergence_matrix(species_data)
        assert div_matrix.shape == (3, 3)

        # PCA
        pca_result = cross_species_pca(species_data, n_components=3)
        assert len(pca_result["coordinates"]) == 12  # 4+4+4

    def test_conservation_to_divergent_pipeline(self) -> None:
        """compute_expression_conservation output feeds directly into identify_divergent_genes."""
        # Create data where some genes are conserved and some divergent
        rng = np.random.default_rng(2024)

        # Gene1: perfectly correlated (conserved) - both increasing
        # Gene2: anti-correlated (divergent) - A increasing, B decreasing
        # Gene3: random/uncorrelated
        expr_a = pd.DataFrame(
            {
                "c1": [1.0, 1.0, rng.random()],
                "c2": [2.0, 2.0, rng.random()],
                "c3": [3.0, 3.0, rng.random()],
                "c4": [4.0, 4.0, rng.random()],
            },
            index=["gene1", "gene2", "gene3"],
        )
        # gene1: B also increasing (correlated), gene2: B decreasing (anti-correlated)
        expr_b = pd.DataFrame(
            {
                "c1": [10.0, 40.0, rng.random()],
                "c2": [20.0, 30.0, rng.random()],
                "c3": [30.0, 20.0, rng.random()],
                "c4": [40.0, 10.0, rng.random()],
            },
            index=["gene1", "gene2", "gene3"],
        )

        conservation = compute_expression_conservation(expr_a, expr_b, method="spearman")
        divergent = identify_divergent_genes(conservation, threshold=0.3)

        # gene1 should be conserved (corr ~ 1.0): A=[1,2,3,4] vs B=[10,20,30,40]
        gene1_corr = conservation[conservation["gene_id"] == "gene1"]["correlation"].values[0]
        assert gene1_corr > 0.5

        # gene2 should be anti-correlated: A=[1,2,3,4] vs B=[40,30,20,10]
        gene2_corr = conservation[conservation["gene_id"] == "gene2"]["correlation"].values[0]
        assert gene2_corr < 0.0

        # gene2 should appear in divergent genes
        assert "gene2" in divergent["gene_id"].values


# =============================================================================
# Edge case and data type tests
# =============================================================================


class TestEdgeCases:
    """Edge case tests for boundary conditions."""

    def test_single_gene_conservation(self) -> None:
        """Single gene in both expression matrices works correctly."""
        expr_a = pd.DataFrame(
            {"c1": [1.0], "c2": [2.0], "c3": [3.0]},
            index=["gene1"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10.0], "c2": [20.0], "c3": [30.0]},
            index=["gene1"],
        )

        result = compute_expression_conservation(expr_a, expr_b)

        assert len(result) == 1
        assert result.iloc[0]["gene_id"] == "gene1"

    def test_single_gene_divergence_matrix(self) -> None:
        """Divergence matrix works with only 1 shared gene."""
        species_data = {
            "sp_a": pd.DataFrame(
                {"s1": [1.0], "s2": [2.0], "s3": [3.0]},
                index=["gene1"],
            ),
            "sp_b": pd.DataFrame(
                {"s1": [10.0], "s2": [20.0], "s3": [30.0]},
                index=["gene1"],
            ),
        }

        result = compute_expression_divergence_matrix(species_data)

        assert result.shape == (2, 2)

    def test_large_gene_set_pca(self) -> None:
        """PCA with many genes is numerically stable."""
        rng = np.random.default_rng(555)
        genes = [f"gene{i}" for i in range(100)]

        species_data = {
            "sp_a": pd.DataFrame(
                rng.random((100, 5)), index=genes, columns=[f"a{i}" for i in range(5)]
            ),
            "sp_b": pd.DataFrame(
                rng.random((100, 5)), index=genes, columns=[f"b{i}" for i in range(5)]
            ),
        }

        result = cross_species_pca(species_data, n_components=5)

        assert len(result["coordinates"]) == 10
        assert not np.any(np.isnan(result["explained_variance"]))

    def test_float_and_int_mixed_expression(self) -> None:
        """Mixed int/float expression values are handled correctly."""
        expr_a = pd.DataFrame(
            {"c1": [1, 2, 3], "c2": [4.0, 5.0, 6.0]},
            index=["gene1", "gene2", "gene3"],
        )
        expr_b = pd.DataFrame(
            {"c1": [10, 20, 30], "c2": [40.0, 50.0, 60.0]},
            index=["gene1", "gene2", "gene3"],
        )

        # Should not raise
        result = compute_expression_conservation(expr_a, expr_b)

        assert len(result) == 3

    def test_ortholog_map_with_empty_map_returns_empty(self) -> None:
        """Empty ortholog map produces empty result."""
        expr = _make_expression_df({"s1": [10.0]}, gene_ids=["GeneA"])

        result = map_expression_to_orthologs(expr, {})

        assert result.empty

    def test_pca_with_partial_gene_overlap(self) -> None:
        """PCA uses only the intersection of genes across species."""
        species_data = {
            "sp_a": pd.DataFrame(
                {"s1": [1.0, 2.0, 3.0], "s2": [4.0, 5.0, 6.0]},
                index=["gene1", "gene2", "gene3"],
            ),
            "sp_b": pd.DataFrame(
                {"s1": [10.0, 20.0], "s2": [40.0, 50.0]},
                index=["gene1", "gene2"],
            ),
        }

        result = cross_species_pca(species_data, n_components=2)

        # Only gene1 and gene2 are shared
        assert result["shared_genes"] == ["gene1", "gene2"]
        assert result["loadings"].shape[0] == 2

    def test_euclidean_conservation_range(self) -> None:
        """Euclidean conservation scores are in expected range [0, 1]."""
        rng = np.random.default_rng(42)
        expr_a = pd.DataFrame(
            rng.random((5, 4)),
            index=[f"gene{i}" for i in range(5)],
            columns=["c1", "c2", "c3", "c4"],
        )
        expr_b = pd.DataFrame(
            rng.random((5, 4)),
            index=[f"gene{i}" for i in range(5)],
            columns=["c1", "c2", "c3", "c4"],
        )

        result = compute_expression_conservation(expr_a, expr_b, method="euclidean")

        for _, row in result.iterrows():
            assert 0.0 <= row["correlation"] <= 1.0 + 1e-6

    def test_phylogenetic_profile_single_gene_single_leaf(self) -> None:
        """Phylogenetic profile with a single-leaf tree (degenerate case)."""
        tree = {"name": "root", "distance": 0.0, "children": [{"name": "sp_a", "distance": 0.1}]}
        expr = pd.DataFrame({"sp_a": [42.0]}, index=["gene1"])

        result = phylogenetic_expression_profile(expr, tree)

        assert not result.empty
        leaf_row = result[result["node"] == "sp_a"]
        assert leaf_row.iloc[0]["mean_expression"] == pytest.approx(42.0)
