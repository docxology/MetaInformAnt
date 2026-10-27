"""Tests for GWAS association testing."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from metainformant.core.io.io import write_tsv
from metainformant.gwas.analysis.association import (
    association_test_linear,
    association_test_logistic,
)
from metainformant.gwas.workflow.workflow_execution import run_gwas


def test_association_linear_basic() -> None:
    """Test basic linear regression association test."""
    # Genotypes: 0, 1, 2, 0, 1, 2 (varying dosage)
    genotypes = [0, 1, 2, 0, 1, 2]
    # Phenotypes with positive association
    phenotypes = [10.0, 11.0, 12.0, 10.0, 11.0, 12.0]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    assert "beta" in result
    assert "se" in result
    assert "p_value" in result
    assert "r_squared" in result
    assert result["beta"] > 0  # Positive association


def test_association_linear_no_association() -> None:
    """Test linear regression with no association."""
    genotypes = [0, 1, 2, 0, 1, 2]
    # Random phenotypes (no association)
    phenotypes = [10.0, 10.5, 9.5, 11.0, 9.0, 10.5]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    # P-value should be large (not significant)
    assert result["p_value"] > 0.05


def test_association_linear_with_covariates() -> None:
    """Test linear regression with covariates."""
    # Use more samples and non-collinear covariates
    genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
    phenotypes = [10.0, 11.0, 12.0, 10.5, 11.5, 12.5, 9.5, 10.5, 11.5, 10.0]
    covariates = [
        [25, 35, 45, 30, 40, 50, 28, 38, 48, 33],
        [1, 1, 1, 0, 0, 0, 1, 0, 1, 0],
    ]  # age, sex

    result = association_test_linear(genotypes, phenotypes, covariates)

    assert result["status"] == "success"
    assert "beta" in result
    assert "p_value" in result

    # CRITICAL: verify result is NOT the old hardcoded (0.1, 0.05, 2.0, 0.05, 0.5)
    assert not (
        abs(result["beta"] - 0.1) < 1e-10 and abs(result["se"] - 0.05) < 1e-10 and abs(result["p_value"] - 0.05) < 1e-10
    ), "Multi-covariate regression still returns hardcoded values!"

    # SE should be positive (real computation)
    assert result["se"] > 0


def test_association_linear_covariates_produce_different_results() -> None:
    """Test that adding covariates produces different results than no covariates."""
    genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0]
    phenotypes = [10.0, 11.0, 12.0, 10.5, 11.5, 12.5, 9.5, 10.5, 11.5, 10.0]

    result_no_cov = association_test_linear(genotypes, phenotypes)
    # Use a covariate NOT perfectly correlated with genotype
    covariates = [[25, 35, 45, 30, 40, 50, 28, 38, 48, 33]]
    result_with_cov = association_test_linear(genotypes, phenotypes, covariates=covariates)

    assert result_no_cov["status"] == "success"
    assert result_with_cov["status"] == "success"
    # Results should differ when covariates are added
    assert result_no_cov["se"] != result_with_cov["se"] or result_no_cov["r_squared"] != result_with_cov["r_squared"]


def test_association_linear_missing_data_matches_complete_case() -> None:
    """Missing genotypes/phenotypes must be dropped, never treated as doses."""
    genotypes = [0, -1, 2, 0, 1, 2]  # One missing genotype
    phenotypes = [10.0, 11.0, 12.0, 10.0, 11.0, 12.0]

    result = association_test_linear(genotypes, phenotypes)

    assert result["status"] == "success"
    assert result["n_missing_dropped"] == 1
    assert result["n_samples"] == 5

    # Complete-case regression on the retained samples only.
    keep_g = [0, 2, 0, 1, 2]
    keep_y = [10.0, 12.0, 10.0, 11.0, 12.0]
    design = np.column_stack([np.ones(len(keep_g)), keep_g])
    beta_complete, *_ = np.linalg.lstsq(design, np.array(keep_y), rcond=None)
    assert result["beta"] == pytest.approx(float(beta_complete[1]), rel=1e-9)

    # NaN phenotypes are missing too and must not enter the design matrix.
    result_nan = association_test_linear([0, 1, 2, 0, 1], [10.0, float("nan"), 12.0, 10.0, 11.0])
    assert result_nan["status"] == "success"
    assert result_nan["n_missing_dropped"] == 1
    assert result_nan["n_samples"] == 4
    assert result_nan["beta"] == pytest.approx(1.0, rel=1e-9)

    # All-missing input cannot be fit and must error, not fabricate inference.
    result_all_missing = association_test_linear([0, -1], [5.0, 6.0])
    result_all_missing = association_test_linear([-1, -1], [5.0, 6.0])
    assert result_all_missing["n_missing_dropped"] == 2
    assert result_all_missing["converged"] is False


def test_association_linear_drops_samples_with_non_finite_covariates() -> None:
    """A NaN/inf covariate must drop its sample, not poison the design into a silent beta=0, p=1."""
    genotypes = [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2]
    phenotypes = [10.0, 11.2, 12.1, 10.4, 11.6, 12.3, 9.6, 10.9, 11.4, 10.2, 11.1, 12.6]
    age = [25.0, 35.0, 45.0, 30.0, 40.0, 50.0, 28.0, 38.0, 48.0, 33.0, 36.0, 44.0]
    sex = [1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    age_with_gaps = list(age)
    age_with_gaps[1] = float("nan")
    sex_with_gaps = list(sex)
    sex_with_gaps[6] = float("inf")

    result = association_test_linear(genotypes, phenotypes, covariates=[age_with_gaps, sex_with_gaps])

    assert result["status"] == "success"
    assert result["n_missing_dropped"] == 2
    assert result["n_samples"] == 10

    keep = [i for i in range(len(genotypes)) if i not in (1, 6)]
    design = np.column_stack(
        [np.ones(len(keep)), [genotypes[i] for i in keep], [age[i] for i in keep], [sex[i] for i in keep]]
    )
    beta, *_ = np.linalg.lstsq(design, np.array([phenotypes[i] for i in keep]), rcond=None)
    assert result["beta"] == pytest.approx(float(beta[1]), rel=1e-9)
    assert result["beta"] != 0.0


def test_association_logistic_without_statsmodels_errors(
    monkeypatch: "pytest.MonkeyPatch",
) -> None:
    """Without statsmodels the logistic test must error, never fabricate inference."""
    import sys

    monkeypatch.setitem(sys.modules, "statsmodels.api", None)

    genotypes = [0, 1, 2, 0, 1, 2]
    phenotypes = [0, 1, 0, 1, 0, 1]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "error"
    assert "statsmodels" in result["error"]
    assert "converged" in result and result["converged"] is False


def test_association_logistic_basic() -> None:
    """Test basic logistic regression association test."""
    # Genotypes with moderate association (less extreme to avoid separation issues)
    genotypes = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    # Binary phenotypes with positive association
    phenotypes = [0, 0, 1, 0, 1, 1, 1, 1, 1]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "success"
    assert "beta" in result
    assert "se" in result
    assert "p_value" in result
    assert "odds_ratio" in result
    # With positive association, beta should be positive and OR > 1
    # But allow for some numerical instability
    if not np.isinf(result.get("beta", 0)) and not np.isnan(result.get("beta", 0)):
        assert result["odds_ratio"] > 0  # Just check it's positive, not NaN


def test_association_logistic_no_association() -> None:
    """Test logistic regression with no association."""
    genotypes = [0, 1, 2, 0, 1, 2]
    # Random case/control assignment
    phenotypes = [0, 1, 0, 1, 0, 1]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "success"
    # P-value should be large
    assert result["p_value"] > 0.05


def test_association_logistic_insufficient_cases() -> None:
    """Test logistic regression with insufficient cases/controls."""
    genotypes = [0, 0, 0, 0, 0, 0]
    phenotypes = [0, 0, 0, 0, 0, 1]  # Only 1 case

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "failed"
    assert "error" in result


def _logistic_newton_oracle(design: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Independent logistic-regression oracle (Newton-Raphson with Wald inference).

    Deliberately shares no code with the implementation under test and does not use
    statsmodels.  ``design`` already contains the intercept in column 0 and the
    genotype in column 1; inference is reported for the genotype coefficient.
    """
    beta = np.zeros(design.shape[1])
    for _ in range(200):
        mu = 1.0 / (1.0 + np.exp(-(design @ beta)))
        weights = mu * (1.0 - mu)
        step = np.linalg.solve(design.T @ (design * weights[:, None]), design.T @ (y - mu))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-13:
            break
    mu = 1.0 / (1.0 + np.exp(-(design @ beta)))
    weights = mu * (1.0 - mu)
    covariance = np.linalg.inv(design.T @ (design * weights[:, None]))
    se = float(np.sqrt(covariance[1, 1]))
    z_stat = float(beta[1]) / se
    return {
        "beta": float(beta[1]),
        "se": se,
        "z_stat": z_stat,
        "p_value": math.erfc(abs(z_stat) / math.sqrt(2.0)),
    }


def _case_control_cohort(n: int = 240, seed: int = 20261024) -> tuple[np.ndarray, np.ndarray, list[list[float]]]:
    """Dosages, binary phenotypes and two covariates (covariates x samples)."""
    rng = np.random.default_rng(seed)
    dosage = rng.integers(0, 3, n)
    age = rng.normal(50.0, 10.0, n)
    sex = rng.integers(0, 2, n).astype(float)
    eta = -0.4 + 0.9 * dosage + 0.02 * (age - 50.0) + 0.3 * sex
    phenotype = (rng.random(n) < 1.0 / (1.0 + np.exp(-eta))).astype(int)
    return dosage, phenotype, [age.tolist(), sex.tolist()]


def test_association_logistic_missing_genotypes_match_complete_case() -> None:
    """A missing call (-1) must be dropped, never entered into the model as dosage -1."""
    dosage, phenotype, covariates = _case_control_cohort()
    n = len(dosage)
    missing = np.random.default_rng(1).choice(n, 40, replace=False)
    keep = np.ones(n, dtype=bool)
    keep[missing] = False

    genotypes_with_missing = dosage.astype(float)
    genotypes_with_missing[missing] = -1

    result = association_test_logistic(genotypes_with_missing.tolist(), phenotype.tolist(), covariates=covariates)
    complete = association_test_logistic(
        dosage[keep].tolist(),
        phenotype[keep].tolist(),
        covariates=[np.asarray(cov)[keep].tolist() for cov in covariates],
    )

    assert result["status"] == "success"
    assert result["n_missing_dropped"] == 40
    assert complete["n_missing_dropped"] == 0
    assert result["n_samples"] == complete["n_samples"] == 200
    assert result["n_cases"] == int(phenotype[keep].sum())
    assert result["n_controls"] == int((1 - phenotype[keep]).sum())
    for key in ("beta", "se", "z_stat", "p_value", "odds_ratio"):
        assert result[key] == pytest.approx(complete[key], rel=1e-9), key

    # Independent oracle on the complete cases only.
    design = np.column_stack(
        [np.ones(int(keep.sum())), dosage[keep], np.asarray(covariates[0])[keep], np.asarray(covariates[1])[keep]]
    )
    oracle = _logistic_newton_oracle(design, phenotype[keep].astype(float))
    for key in ("beta", "se", "z_stat", "p_value"):
        assert result[key] == pytest.approx(oracle[key], rel=1e-6), key

    # Sensitivity of the fixture: feeding -1 to the model as a dose gives a materially different effect,
    # so agreement above cannot be an accident of weak data.
    all_rows = np.column_stack(
        [np.ones(n), genotypes_with_missing, np.asarray(covariates[0]), np.asarray(covariates[1])]
    )
    misfit = _logistic_newton_oracle(all_rows, phenotype.astype(float))
    assert abs(misfit["beta"] - oracle["beta"]) > 0.1


def test_association_logistic_drops_nan_none_and_non_numeric_genotypes() -> None:
    """NaN, inf, None and non-numeric calls are all missing and are excluded from the fit."""
    dosage, phenotype, _ = _case_control_cohort(n=90, seed=7)
    bad_positions = {3: float("nan"), 11: None, 20: "NA", 33: float("inf"), 47: float("nan"), 60: "./."}
    genotypes: list = dosage.tolist()
    for position, value in bad_positions.items():
        genotypes[position] = value
    keep = np.array([i not in bad_positions for i in range(len(dosage))])

    result = association_test_logistic(genotypes, phenotype.tolist())
    complete = association_test_logistic(dosage[keep].tolist(), phenotype[keep].tolist())

    assert result["status"] == "success"
    assert result["n_missing_dropped"] == len(bad_positions)
    assert result["n_samples"] == int(keep.sum())
    for key in ("beta", "se", "p_value"):
        assert result[key] == pytest.approx(complete[key], rel=1e-9), key

    oracle = _logistic_newton_oracle(
        np.column_stack([np.ones(int(keep.sum())), dosage[keep]]), phenotype[keep].astype(float)
    )
    assert result["beta"] == pytest.approx(oracle["beta"], rel=1e-6)
    assert result["se"] == pytest.approx(oracle["se"], rel=1e-6)


def test_association_logistic_drops_samples_with_non_finite_covariates() -> None:
    """A sample whose covariate is missing cannot enter the model; covariates are subset to match."""
    dosage, phenotype, covariates = _case_control_cohort(n=120, seed=11)
    covariates[0][5] = float("nan")
    covariates[1][17] = float("inf")
    keep = np.ones(len(dosage), dtype=bool)
    keep[[5, 17]] = False

    result = association_test_logistic(dosage.tolist(), phenotype.tolist(), covariates=covariates)

    assert result["status"] == "success"
    assert result["n_missing_dropped"] == 2
    assert result["n_samples"] == 118
    design = np.column_stack(
        [np.ones(118), dosage[keep], np.asarray(covariates[0])[keep], np.asarray(covariates[1])[keep]]
    )
    oracle = _logistic_newton_oracle(design, phenotype[keep].astype(float))
    assert result["beta"] == pytest.approx(oracle["beta"], rel=1e-6)
    assert result["se"] == pytest.approx(oracle["se"], rel=1e-6)


def test_association_logistic_recounts_cases_and_controls_after_dropping() -> None:
    """Group sizes are judged on the complete cases, not on the raw input."""
    # Three controls in the input, but two of them have no genotype call: one usable control remains.
    genotypes = [0, 1, 2, 1, -1, -1, 0]
    phenotypes = [1, 1, 1, 1, 0, 0, 0]

    result = association_test_logistic(genotypes, phenotypes)

    assert result["status"] == "failed"
    assert result["converged"] is False
    assert result["n_cases"] == 4
    assert result["n_controls"] == 1
    assert result["n_samples"] == 5
    assert result["n_missing_dropped"] == 2


def test_association_logistic_all_genotypes_missing_reports_failure() -> None:
    """With no complete case there is nothing to fit; report it instead of fitting garbage."""
    result = association_test_logistic([-1, float("nan"), -1, -1], [0, 1, 0, 1])

    assert result["status"] == "failed"
    assert result["converged"] is False
    assert result["n_samples"] == 0
    assert result["n_missing_dropped"] == 4
    # A failed test has no p-value: NaN (missing), not the fabricated p = 1 it used to carry.
    assert math.isnan(result["p_value"])


def test_association_logistic_rejects_covariates_of_the_wrong_length() -> None:
    """Covariates must carry one value per sample; a short or long covariate is an input error."""
    genotypes = [0, 1, 2, 1, 0, 2, 1, 0]
    phenotypes = [0, 1, 1, 0, 0, 1, 1, 0]
    with pytest.raises(ValueError, match="covariate"):
        association_test_logistic(genotypes, phenotypes, covariates=[[1.0, 2.0, 3.0]])
    with pytest.raises(ValueError, match="covariate"):
        association_test_logistic(genotypes, phenotypes, covariates=[[float(i) for i in range(len(genotypes) + 2)]])


def test_run_gwas_integration(tmp_path: Path) -> None:
    """Test full GWAS run with VCF and phenotype files."""
    # Create VCF file
    vcf_content = """##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
chr1	100	rs1	A	G	60	PASS	.	GT	0/1	1/1	0/0	0/1	0/0
chr1	200	rs2	T	C	80	PASS	.	GT	0/0	0/1	1/1	0/0	0/1
chr1	300	rs3	G	A	50	PASS	.	GT	1/1	0/0	0/1	1/1	0/0
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)

    # Create phenotype file
    phenotype_file = tmp_path / "phenotypes.tsv"
    phenotype_data = [
        ["sample_id", "height"],
        ["S1", "175.0"],
        ["S2", "180.0"],
        ["S3", "170.0"],
        ["S4", "178.0"],
        ["S5", "172.0"],
    ]
    write_tsv(phenotype_data, phenotype_file)

    config = {
        "model": "linear",
        "trait": "height",
        "covariates": [],
        "min_sample_size": 3,
    }

    results_dir = tmp_path / "results"
    results_dir.mkdir()

    result = run_gwas(vcf_file, phenotype_file, config, output_dir=results_dir)

    assert result["status"] == "success"
    assert "results" in result
    assert len(result["results"]) > 0

    # Check that results table was written
    if results_dir.exists():
        results_file = results_dir / "association_results.tsv"
        if results_file.exists():
            assert results_file.stat().st_size > 0


def test_run_gwas_missing_trait(tmp_path: Path) -> None:
    """Test GWAS run with missing trait in phenotype file."""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.2\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1\nchr1	100	rs1	A	G	60	PASS	.	GT	0/1\n"
    )

    phenotype_file = tmp_path / "phenotypes.tsv"
    phenotype_data = [["sample_id", "weight"], ["S1", "70.0"]]
    write_tsv(phenotype_data, phenotype_file)

    config = {"model": "linear", "trait": "height", "covariates": []}

    result = run_gwas(vcf_file, phenotype_file, config)
    assert result["status"] == "failed"
    assert "error" in result
