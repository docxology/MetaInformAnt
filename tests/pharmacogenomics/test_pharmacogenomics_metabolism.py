"""Tests for pharmacogenomics metabolism: metabolizer classification and dose adjustment.

Tests metabolizer classification from activity scores, activity score computation
from diplotype strings, the built-in allele function table, dose adjustment
recommendations, and the predict_metabolizer alias.

NO MOCKING -- all tests use real implementations.
"""

from __future__ import annotations

import pytest

# Also test the alias from the package __init__
from metainformant.pharmacogenomics import predict_metabolizer
from metainformant.pharmacogenomics.metabolism.metabolizer_status import (
    classify_metabolizer,
    compute_activity_score,
    default_allele_function_table,
    dose_adjustment,
    predict_metabolizer_status,
)


class TestDefaultAlleleFunctionTable:
    """Tests for the built-in allele function table."""

    def test_returns_dict(self) -> None:
        table = default_allele_function_table()
        assert isinstance(table, dict)

    def test_contains_key_genes(self) -> None:
        table = default_allele_function_table()
        assert "CYP2D6" in table
        assert "CYP2C19" in table
        assert "CYP2C9" in table
        assert "DPYD" in table
        assert "TPMT" in table

    def test_cyp2d6_allele_values(self) -> None:
        table = default_allele_function_table()
        cyp2d6 = table["CYP2D6"]
        assert cyp2d6["*1"] == 1.0
        assert cyp2d6["*4"] == 0.0
        assert cyp2d6["*9"] == 0.5

    def test_cyp2c19_star17_increased(self) -> None:
        table = default_allele_function_table()
        assert table["CYP2C19"]["*17"] == 1.5


class TestComputeActivityScore:
    """Tests for activity score computation from diplotype strings."""

    def test_normal_homozygous(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        score = compute_activity_score("*1/*1", table)
        assert score == 2.0

    def test_poor_metabolizer(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        score = compute_activity_score("*4/*4", table)
        assert score == 0.0

    def test_heterozygous(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        score = compute_activity_score("*1/*4", table)
        assert score == 1.0

    def test_decreased_function(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        score = compute_activity_score("*1/*9", table)
        assert score == 1.5

    def test_cyp2c19_rapid(self) -> None:
        table = default_allele_function_table()["CYP2C19"]
        score = compute_activity_score("*1/*17", table)
        assert score == 2.5

    def test_invalid_format_raises(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        with pytest.raises(ValueError, match="Invalid diplotype"):
            compute_activity_score("*1", table)

    def test_unknown_allele_raises(self) -> None:
        table = default_allele_function_table()["CYP2D6"]
        with pytest.raises(ValueError, match="not found"):
            compute_activity_score("*1/*999", table)


class TestClassifyMetabolizer:
    """Tests for metabolizer classification from activity score."""

    def test_cyp2d6_poor(self) -> None:
        result = classify_metabolizer(0.0, "CYP2D6")
        assert result == "poor"

    def test_cyp2d6_intermediate(self) -> None:
        result = classify_metabolizer(0.5, "CYP2D6")
        assert result == "intermediate"

    def test_cyp2d6_normal(self) -> None:
        result = classify_metabolizer(2.0, "CYP2D6")
        assert result == "normal"

    def test_cyp2d6_ultrarapid(self) -> None:
        result = classify_metabolizer(3.0, "CYP2D6")
        assert result == "ultrarapid"

    def test_cyp2c19_poor(self) -> None:
        result = classify_metabolizer(0.0, "CYP2C19")
        assert result == "poor"

    def test_cyp2c19_normal(self) -> None:
        result = classify_metabolizer(2.0, "CYP2C19")
        assert result == "normal"

    def test_case_insensitive_gene(self) -> None:
        result = classify_metabolizer(0.0, "cyp2d6")
        assert result == "poor"


class TestDoseAdjustment:
    """Tests for dose adjustment recommendations."""

    def test_codeine_poor_avoid(self) -> None:
        result = dose_adjustment("poor", "codeine")
        assert result["adjusted_dose_fraction"] == 0.0
        assert result["evidence_level"] == "A"
        assert "avoid" in result["recommendation"].lower()

    def test_codeine_normal_standard(self) -> None:
        result = dose_adjustment("normal", "codeine")
        assert result["adjusted_dose_fraction"] == 1.0

    def test_warfarin_poor_reduced(self) -> None:
        result = dose_adjustment("poor", "warfarin")
        assert result["adjusted_dose_fraction"] < 1.0

    def test_unknown_drug(self) -> None:
        result = dose_adjustment("normal", "unknown_drug")
        assert result["adjusted_dose_fraction"] == 1.0
        assert result["evidence_level"] == "D"

    def test_unknown_metabolizer_status(self) -> None:
        result = dose_adjustment("superrapid", "codeine")
        assert result["evidence_level"] == "D"

    def test_fluorouracil_poor(self) -> None:
        result = dose_adjustment("poor", "fluorouracil")
        assert result["adjusted_dose_fraction"] == 0.0
        assert "avoid" in result["recommendation"].lower()

    def test_azathioprine_intermediate(self) -> None:
        result = dose_adjustment("intermediate", "azathioprine")
        assert result["adjusted_dose_fraction"] == 0.5

    def test_omeprazole_ultrarapid(self) -> None:
        result = dose_adjustment("ultrarapid", "omeprazole")
        assert result["adjusted_dose_fraction"] > 1.0


class TestPredictMetabolizerStatus:
    """Tests for predict_metabolizer_status with full genotype dict."""

    def test_predict_from_diplotype(self) -> None:
        genotype = {"diplotype": "*1/*4"}
        result = predict_metabolizer_status(genotype, "CYP2D6")
        assert "phenotype" in result
        assert "activity_score" in result
        assert result["activity_score"] == 1.0
        assert result["gene"] == "CYP2D6"

    def test_predict_from_alleles(self) -> None:
        genotype = {"alleles": ["*4", "*4"]}
        result = predict_metabolizer_status(genotype, "CYP2D6")
        assert result["phenotype"] == "poor"
        assert result["activity_score"] == 0.0

    def test_predict_includes_thresholds(self) -> None:
        genotype = {"diplotype": "*1/*1"}
        result = predict_metabolizer_status(genotype, "CYP2D6")
        assert "thresholds_used" in result
        assert "poor" in result["thresholds_used"]

    def test_predict_missing_fields_raises(self) -> None:
        with pytest.raises(ValueError, match="diplotype.*alleles"):
            predict_metabolizer_status({}, "CYP2D6")

    def test_predict_unknown_gene_raises(self) -> None:
        with pytest.raises(ValueError, match="No allele function table"):
            predict_metabolizer_status({"diplotype": "*1/*1"}, "INVALID_GENE")

    def test_predict_cyp2c19_rapid(self) -> None:
        genotype = {"diplotype": "*1/*17"}
        result = predict_metabolizer_status(genotype, "CYP2C19")
        assert result["activity_score"] == 2.5
        assert result["phenotype"] == "rapid"


class TestPredictMetabolizerAlias:
    """Tests that predict_metabolizer alias works the same as predict_metabolizer_status."""

    def test_alias_returns_same_result(self) -> None:
        genotype = {"diplotype": "*1/*4"}
        result_alias = predict_metabolizer(genotype, "CYP2D6")
        result_orig = predict_metabolizer_status(genotype, "CYP2D6")
        assert result_alias == result_orig

    def test_alias_is_callable(self) -> None:
        assert callable(predict_metabolizer)
        genotype = {"diplotype": "*4/*4"}
        result = predict_metabolizer(genotype, "CYP2D6")
        assert result["phenotype"] == "poor"
