"""Conditional eQTL analysis accepts integer genotype dosages.

Zero-mocks: the real conditional analysis runs on a small hand-built dataset.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from metainformant.gwas.finemapping.eqtl import conditional_eqtl


def _data() -> tuple[np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(3)
    lead = rng.integers(0, 3, size=60)
    second = rng.integers(0, 3, size=60)
    null = rng.integers(0, 3, size=60)
    expression = 2.0 * lead + 1.0 * second + rng.normal(0.0, 0.5, size=60)
    genotypes = pd.DataFrame([lead, second, null], index=["lead", "second", "null"])  # int64 dosages
    return expression, genotypes


def test_integer_dosages_are_conditioned_like_floats() -> None:
    expression, genotypes = _data()
    assert genotypes.dtypes.unique().tolist() == [np.dtype("int64")]

    from_ints = conditional_eqtl(expression, genotypes, "lead")
    from_floats = conditional_eqtl(expression, genotypes.astype(float), "lead")

    pd.testing.assert_frame_equal(from_ints, from_floats)
    second = from_ints.set_index("variant_id").loc["second"]
    assert second["pvalue"] < 1e-6
