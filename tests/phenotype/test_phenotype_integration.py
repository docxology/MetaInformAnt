"""Tests for phenotype cross-omic integration: association, correlation, GxE.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import pytest

from metainformant.phenotype.integration.cross_omic import (
    multi_phenotype_integration,
    phenotype_environment_interaction,
    phenotype_genotype_association,
    trait_expression_correlation,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_phenotypes(n=20):
    return {f"sample_{i}": float(i * 2 + 5) for i in range(n)}


def _make_genotypes(n=20, n_variants=3):
    import random

    random.seed(42)
    genos = {}
    for v in range(n_variants):
        genos[f"rs{v + 1000}"] = [random.choice([0, 1, 2]) for _ in range(n)]
    return genos


def _make_expression(n=20, n_genes=5):
    import random

    random.seed(42)
    expr = {}
    sample_ids = [f"sample_{i}" for i in range(n)]
    for g in range(n_genes):
        gene_id = f"gene_{g}"
        expr[gene_id] = {s: random.gauss(10, 3) for s in sample_ids}
    return expr


# ---------------------------------------------------------------------------
# phenotype_genotype_association
# ---------------------------------------------------------------------------


class TestPhenotypeGenotypeAssociation:
    def test_basic_association(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        result = phenotype_genotype_association(phenotypes, genotypes)
        assert result["n_samples"] == 20
        assert result["n_variants"] == 3
        assert len(result["associations"]) > 0

    def test_linear_regression_method(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        result = phenotype_genotype_association(phenotypes, genotypes, method="linear_regression")
        for variant_id, stats in result["associations"].items():
            assert "beta" in stats
            assert "r_squared" in stats
            assert "p_value" in stats

    def test_correlation_method(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        result = phenotype_genotype_association(phenotypes, genotypes, method="correlation")
        for variant_id, stats in result["associations"].items():
            assert "correlation" in stats
            assert "p_value" in stats

    def test_empty_input(self):
        result = phenotype_genotype_association({}, {})
        assert "error" in result

    def test_p_values_in_range(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        result = phenotype_genotype_association(phenotypes, genotypes)
        for variant_id, stats in result["associations"].items():
            assert 0.0 <= stats["p_value"] <= 1.0


# ---------------------------------------------------------------------------
# trait_expression_correlation
# ---------------------------------------------------------------------------


class TestTraitExpressionCorrelation:
    def test_basic_correlation(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression)
        assert result["n_samples"] == 20
        assert result["n_genes_tested"] == 5
        assert len(result["correlations"]) > 0

    def test_pearson_method(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression, method="pearson")
        assert result["method"] == "pearson"

    def test_spearman_method(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression, method="spearman")
        assert result["method"] == "spearman"

    def test_gene_subset(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression, gene_list=["gene_0", "gene_1"])
        assert result["n_genes_tested"] == 2

    def test_top_genes_returned(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression)
        assert "top_genes" in result

    def test_empty_input(self):
        result = trait_expression_correlation({}, {})
        assert "error" in result

    def test_correlation_in_range(self):
        traits = _make_phenotypes()
        expression = _make_expression()
        result = trait_expression_correlation(traits, expression)
        for gene_id, stats in result["correlations"].items():
            assert -1.0 <= stats["correlation"] <= 1.0


# ---------------------------------------------------------------------------
# multi_phenotype_integration
# ---------------------------------------------------------------------------


class TestMultiPhenotypeIntegration:
    def test_basic_integration(self):
        phenotypes = {
            "height": {f"s{i}": float(i * 2) for i in range(10)},
            "weight": {f"s{i}": float(i * 3 + 5) for i in range(10)},
            "bmi": {f"s{i}": float(i * 0.5 + 20) for i in range(10)},
        }
        result = multi_phenotype_integration(phenotypes)
        assert result["n_phenotypes"] == 3
        assert len(result["correlation_matrix"]) == 3

    def test_strong_correlation_detected(self):
        phenotypes = {
            "x": {f"s{i}": float(i) for i in range(20)},
            "y": {f"s{i}": float(i * 2 + 1) for i in range(20)},
        }
        result = multi_phenotype_integration(phenotypes)
        assert len(result["strong_pairs"]) >= 1

    def test_insufficient_phenotypes(self):
        phenotypes = {"only_one": {f"s{i}": float(i) for i in range(10)}}
        result = multi_phenotype_integration(phenotypes)
        assert "error" in result

    def test_diagonal_is_one(self):
        phenotypes = {
            "a": {f"s{i}": float(i) for i in range(10)},
            "b": {f"s{i}": float(i * 2) for i in range(10)},
        }
        result = multi_phenotype_integration(phenotypes)
        for i in range(result["n_phenotypes"]):
            assert result["correlation_matrix"][i][i] == 1.0


# ---------------------------------------------------------------------------
# phenotype_environment_interaction
# ---------------------------------------------------------------------------


class TestPhenotypeEnvironmentInteraction:
    def test_basic_gxe(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        environment = {f"sample_{i}": float(i * 0.5) for i in range(20)}
        result = phenotype_environment_interaction(phenotypes, genotypes, environment)
        assert result["n_samples"] == 20
        assert "interactions" in result

    def test_multiplicative_model(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        environment = {f"sample_{i}": float(i * 0.5) for i in range(20)}
        result = phenotype_environment_interaction(
            phenotypes, genotypes, environment, interaction_model="multiplicative"
        )
        assert result["model"] == "multiplicative"

    def test_additive_model(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        environment = {f"sample_{i}": float(i * 0.5) for i in range(20)}
        result = phenotype_environment_interaction(
            phenotypes, genotypes, environment, interaction_model="additive"
        )
        assert result["model"] == "additive"

    def test_empty_input(self):
        result = phenotype_environment_interaction({}, {}, {})
        assert "error" in result

    def test_significant_interactions_sorted(self):
        phenotypes = _make_phenotypes()
        genotypes = _make_genotypes()
        environment = {f"sample_{i}": float(i * 0.5) for i in range(20)}
        result = phenotype_environment_interaction(phenotypes, genotypes, environment)
        sig = result.get("significant_interactions", [])
        for i in range(len(sig) - 1):
            assert sig[i]["p_value"] <= sig[i + 1]["p_value"]
