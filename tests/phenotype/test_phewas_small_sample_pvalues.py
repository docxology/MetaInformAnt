"""PheWAS p-values come from Student's t at every sample size.

Zero-mocks: ``scipy.stats.linregress`` on the same data is the oracle.
"""

from __future__ import annotations

import pytest
from scipy import stats

from metainformant.phenotype.gwas_integration.phewas import run_phewas


@pytest.mark.parametrize("n_copies", [1, 2, 5])
def test_p_values_match_linear_regression(n_copies: int) -> None:
    genotype = [0, 1, 2, 1, 0, 2, 1, 0, 2, 1] * n_copies
    trait = [0, 1, 0, 0, 1, 0, 1, 0, 0, 1] * n_copies
    result = run_phewas(genotype, {"t": trait})[0]

    assert result["p_value"] == pytest.approx(stats.linregress(genotype, trait).pvalue, rel=1e-6)
