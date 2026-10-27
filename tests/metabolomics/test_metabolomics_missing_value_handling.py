"""Missing values (NaN) stay honest through the analysis steps that follow ``read_csv``.

``read_csv`` now returns NaN for a missing cell (it used to return 0.0). Two steps of
``metainformant.metabolomics.analysis.identification`` would have turned that NaN into silently wrong numbers:

* ``normalize_intensities`` with ``total_ion_count`` or ``median``: the NaN total of the sample that contained it
  was replaced by 1.0 (``np.where(total > 0, total, 1.0)``), so that sample was multiplied by the median total
  instead of being scaled, and, because the replaced 1.0 entered the median, every other sample got a wrong factor;
* ``fold_change``: a NaN group mean was floored at 1e-10 like a genuine zero, a fabricated log2 fold change of about
  -40 for what is really unknown.

Now the scale-based normalisations refuse NaN (impute first) and the fold change of such a metabolite is NaN.

Real computation; no test doubles.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from metainformant.metabolomics.analysis.identification import (
    fold_change,
    missing_value_imputation,
    normalize_intensities,
)
from metainformant.metabolomics.io.formats import read_csv

# metabolite x sample; metabolite 1 has a missing value in sample 1 (group A below)
WITH_GAP = np.array(
    [
        [100.0, 200.0, 300.0, 400.0],
        [50.0, np.nan, 150.0, 200.0],
        [10.0, 20.0, 30.0, 40.0],
    ]
)


class TestNormalizationWithMissingValues:
    @pytest.mark.parametrize("method", ["total_ion_count", "median"])
    def test_scale_based_methods_refuse_data_with_nan(self, method: str) -> None:
        with pytest.raises(ValueError, match=r"needs complete data.*missing_value_imputation"):
            normalize_intensities(WITH_GAP, method=method)

    def test_log2_keeps_a_nan_where_it_is(self) -> None:
        result = normalize_intensities(WITH_GAP, method="log2")

        assert np.isnan(result[1, 1])
        assert np.isnan(result).sum() == 1
        observed = ~np.isnan(WITH_GAP)
        np.testing.assert_allclose(result[observed], np.log2(WITH_GAP[observed] + 1))

    def test_complete_data_is_scaled_exactly_as_before(self) -> None:
        data = np.array([[100.0, 200.0, 300.0], [50.0, 100.0, 150.0]])

        # column totals 150, 300, 450 have median 300, so every column ends up with total 300
        np.testing.assert_allclose(
            normalize_intensities(data, "total_ion_count"), [[200.0, 200.0, 200.0], [100.0, 100.0, 100.0]]
        )

    def test_the_documented_order_read_impute_normalise_gives_correct_totals(self, tmp_path: Path) -> None:
        path = tmp_path / "pipeline.csv"
        path.write_text("met,S1,S2,S3,S4\nA,100,200,300,400\nB,50,NA,150,200\nC,10,20,30,40\n")

        imputed = missing_value_imputation(read_csv(path).intensities, method="min_half")
        normalised = normalize_intensities(imputed, method="total_ion_count")

        # after imputation B/S2 = min(50, 150, 200) / 2 = 25, so the sample totals are 160, 245, 480, 640 and
        # their median is (245 + 480) / 2 = 362.5, which is the total every normalised sample must reach
        np.testing.assert_allclose(imputed[1], [50.0, 25.0, 150.0, 200.0])
        np.testing.assert_allclose(normalised.sum(axis=0), [362.5] * 4)


class TestFoldChangeWithMissingValues:
    def test_a_metabolite_with_a_missing_value_has_a_nan_fold_change_not_a_floored_one(self) -> None:
        result = fold_change(WITH_GAP, group_a=[0, 1], group_b=[2, 3])

        # group means: metabolite 0: 150 vs 350; metabolite 2: 15 vs 35; metabolite 1 has no defined group-A mean
        assert result[0] == pytest.approx(math.log2(150 / 350))
        assert math.isnan(result[1])
        assert result[2] == pytest.approx(math.log2(15 / 35))

    def test_a_missing_value_in_the_second_group_also_gives_nan(self) -> None:
        result = fold_change(WITH_GAP, group_a=[2, 3], group_b=[0, 1])

        assert result[0] == pytest.approx(math.log2(350 / 150))
        assert math.isnan(result[1])

    def test_a_group_without_the_missing_value_is_unaffected(self) -> None:
        # the NaN is in sample 1; comparing samples 2,3 with sample 0 alone never touches it
        result = fold_change(WITH_GAP, group_a=[2, 3], group_b=[0])

        assert result[1] == pytest.approx(math.log2(175 / 50))

    def test_a_genuine_zero_mean_is_still_floored_to_keep_the_ratio_finite(self) -> None:
        data = np.array([[0.0, 0.0, 100.0, 100.0]])

        assert fold_change(data, [0, 1], [2, 3])[0] == pytest.approx(math.log2(1e-10 / 100.0))

    def test_complete_data_gives_the_documented_log2_ratio(self) -> None:
        data = np.array([[200.0, 200.0, 100.0, 100.0]])

        assert fold_change(data, [0, 1], [2, 3])[0] == pytest.approx(1.0)
