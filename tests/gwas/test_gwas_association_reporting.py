"""Association tests report degenerate fits, convergence and p-values honestly.

Zero-mocks: scipy's ``linregress`` and statsmodels' OLS are the independent
oracles for the p-values; the bcftools test runs the real tool when installed.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import statsmodels.api as sm
from scipy import stats

from metainformant.gwas.analysis.association import association_test_linear, association_test_logistic
from metainformant.gwas.analysis.calling import _bcftools_filter_command, filter_variants_bcftools

GENOTYPES = [0, 1, 2, 0, 1, 2, 1, 0, 2, 1]
PHENOTYPES = [1.1, 2.3, 2.9, 0.7, 2.0, 3.4, 1.6, 1.2, 2.8, 2.4]


def test_single_predictor_p_value_uses_students_t() -> None:
    result = association_test_linear(GENOTYPES, PHENOTYPES)
    oracle = stats.linregress(GENOTYPES, PHENOTYPES)

    assert result["status"] == "success"
    assert result["beta"] == pytest.approx(oracle.slope, rel=1e-10)
    assert result["p_value"] == pytest.approx(oracle.pvalue, rel=1e-8)


def test_covariate_adjusted_p_value_matches_ols() -> None:
    covariate = [0.3, -1.2, 0.8, 1.5, -0.4, 0.1, 0.9, -0.7, 0.2, -1.1]
    result = association_test_linear(GENOTYPES, PHENOTYPES, covariates=[covariate])
    design = sm.add_constant(np.column_stack([GENOTYPES, covariate]))
    oracle = sm.OLS(PHENOTYPES, design).fit()

    assert result["beta"] == pytest.approx(oracle.params[1], rel=1e-8)
    assert result["p_value"] == pytest.approx(oracle.pvalues[1], rel=1e-6)


@pytest.mark.parametrize(
    ("genotypes", "covariates", "reason"),
    [
        (GENOTYPES, [[2 * g for g in GENOTYPES]], "singular"),
        ([1] * 10, None, "constant"),
        ([0, 1], None, "at least 3"),
    ],
)
def test_degenerate_linear_fits_are_errors_not_null_results(genotypes, covariates, reason) -> None:
    result = association_test_linear(genotypes, PHENOTYPES[: len(genotypes)], covariates=covariates)

    assert result["status"] == "error"
    assert reason in result["error"]


def test_logistic_reports_whether_the_fit_converged() -> None:
    genotypes = [0] * 10 + [1] * 10 + [2] * 10
    separated = [0] * 15 + [1] * 15  # genotype 0 all controls, genotype 2 all cases
    mixed = [0, 1] * 15

    assert association_test_logistic(genotypes, separated)["converged"] is False
    assert association_test_logistic(genotypes, mixed)["converged"] is True


def test_bcftools_filter_includes_sites_that_pass_every_threshold() -> None:
    cmd = _bcftools_filter_command("in.vcf", "out.vcf", {"min_qual": 20, "max_missing": 0.1, "min_depth": 5})

    assert cmd[cmd.index("-i") + 1] == "(QUAL>=20) && (F_MISSING<=0.1) && (DP>=5)"
    assert "-e" not in cmd
    assert cmd[-1] == "in.vcf"
    assert "-i" not in _bcftools_filter_command("in.vcf", "out.vcf", {})


@pytest.mark.external_tool
@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools is not installed")
def test_bcftools_filter_keeps_high_quality_sites(tmp_path: Path) -> None:
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tlow\tA\tG\t5\tPASS\t.\n"
        "chr1\t200\thigh\tC\tT\t60\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    filter_variants_bcftools(vcf, out, {"min_qual": 20})

    kept = [line.split("\t")[2] for line in out.read_text().splitlines() if not line.startswith("#")]
    assert kept == ["high"]
    assert subprocess.run(["bcftools", "--version"], capture_output=True).returncode == 0
