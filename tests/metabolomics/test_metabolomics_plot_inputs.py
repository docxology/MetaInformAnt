"""``prepare_spectrum_plot``, ``volcano_plot_data`` and ``pca_metabolomics`` neither reinterpret nor clip bad arguments.

Regression tests for ``metainformant.metabolomics.visualization.plots``:

* ``prepare_spectrum_plot(top_n=0)`` sliced ``argsort(...)[-0:]``, which is the whole array, so it returned *every*
  peak; a negative ``top_n`` silently dropped the weakest peaks instead of selecting any;
* an empty spectrum reported a base peak at m/z ``0.0`` with intensity ``0.0``, a peak that does not exist;
* ``volcano_plot_data`` clipped every p-value into ``[1e-300, 1]``, so a negative p-value (a bug upstream) was
  plotted as the most significant metabolite and called significant, and a p-value above 1 was plotted as 1.
* ``pca_metabolomics(n_components=k)`` returned fewer than ``k`` components, without a word, whenever ``k`` exceeded
  the number of samples or metabolites (and, for data without variance, an ``explained_variance_ratio`` sized by the
  request instead of by the scores), and a NaN from ``read_csv`` failed with an obscure "SVD did not converge".

Real computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

from metainformant.metabolomics.visualization.plots import pca_metabolomics, prepare_spectrum_plot, volcano_plot_data


class TestTopN:
    MZ = np.array([100.0, 110.0, 120.0, 130.0, 140.0, 150.0])
    INTENSITY = np.array([10.0, 60.0, 20.0, 100.0, 40.0, 80.0])

    @pytest.mark.parametrize("bad", [0, -1, -20, 2.5, "3", None, float("nan")])
    def test_top_n_must_be_a_positive_integer(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"top_n must be a positive integer"):
            prepare_spectrum_plot(self.MZ, self.INTENSITY, top_n=bad)

    def test_the_n_most_intense_peaks_come_back_in_mz_order(self) -> None:
        result = prepare_spectrum_plot(self.MZ, self.INTENSITY, top_n=3)

        assert result.mz_values.tolist() == [110.0, 130.0, 150.0]  # 60, 100 and 80 are the three largest
        assert result.intensities.tolist() == [60.0, 100.0, 80.0]  # relative to the base peak, 100

    def test_top_n_one_keeps_only_the_base_peak(self) -> None:
        result = prepare_spectrum_plot(self.MZ, self.INTENSITY, top_n=1)

        assert result.mz_values.tolist() == [130.0] and result.base_peak_mz == 130.0

    def test_a_top_n_larger_than_the_spectrum_keeps_every_peak(self) -> None:
        result = prepare_spectrum_plot(self.MZ, self.INTENSITY, top_n=100)

        assert result.mz_values.tolist() == self.MZ.tolist()

    def test_equally_intense_peaks_are_chosen_by_a_stable_sort(self) -> None:
        mz = np.arange(100.0, 108.0)
        intensity = np.array([5.0, 9.0, 5.0, 9.0, 5.0, 9.0, 5.0, 9.0])

        picked = prepare_spectrum_plot(mz, intensity, top_n=6).mz_values.tolist()

        # the four 9s are kept; of the four tied 5s the two with the higher m/z are (never "whichever the sort
        # algorithm of this CPU happens to leave last")
        assert picked == [101.0, 103.0, 104.0, 105.0, 106.0, 107.0]

    def test_ties_at_the_cut_keep_the_higher_mz_peaks(self) -> None:
        mz = np.arange(100.0, 106.0)
        intensity = np.array([9.0, 9.0, 9.0, 9.0, 9.0, 9.0])  # all tie

        assert prepare_spectrum_plot(mz, intensity, top_n=2).mz_values.tolist() == [104.0, 105.0]

    def test_numpy_integer_top_n_is_accepted(self) -> None:
        result = prepare_spectrum_plot(self.MZ, self.INTENSITY, top_n=np.int64(2))

        assert result.mz_values.tolist() == [130.0, 150.0]


class TestShapesAndEmptySpectra:
    def test_arrays_of_different_lengths_raise(self) -> None:
        with pytest.raises(ValueError, match=r"same length"):
            prepare_spectrum_plot(np.array([100.0, 200.0]), np.array([1.0]))

    def test_a_two_dimensional_array_raises(self) -> None:
        with pytest.raises(ValueError, match=r"1D arrays"):
            prepare_spectrum_plot(np.ones((2, 2)), np.ones((2, 2)))

    def test_lists_are_accepted_like_arrays(self) -> None:
        result = prepare_spectrum_plot([100.0, 200.0], [1, 4])

        assert result.base_peak_mz == 200.0 and result.intensities.tolist() == [25.0, 100.0]

    def test_an_empty_spectrum_has_no_base_peak(self) -> None:
        result = prepare_spectrum_plot(np.array([]), np.array([]))

        assert result.mz_values.size == 0 and result.intensities.size == 0 and result.annotations == {}
        assert math.isnan(result.base_peak_mz) and math.isnan(result.base_peak_intensity)

    def test_a_spectrum_of_zero_intensities_keeps_its_zero_base_peak(self) -> None:
        # a real base peak that happens to be silent is not "no base peak"
        result = prepare_spectrum_plot(np.array([100.0, 200.0]), np.array([0.0, 0.0]))

        assert result.base_peak_intensity == 0.0 and result.base_peak_mz == 100.0
        assert result.intensities.tolist() == [0.0, 0.0]


class TestVolcanoPValues:
    NAMES = ["A", "B", "C"]
    FC = np.array([2.0, 2.0, 2.0])

    @pytest.mark.parametrize("bad", [-0.01, -1.0, 1.0000001, 7.0, float("-inf"), float("inf")])
    def test_a_p_value_outside_zero_one_raises_and_names_the_metabolite(self, bad: float) -> None:
        with pytest.raises(ValueError, match=r"p_values must lie in \[0, 1\]; 'B' has"):
            volcano_plot_data(self.NAMES, self.FC, np.array([0.01, bad, 0.5]))

    def test_a_negative_p_value_is_not_the_most_significant_point(self) -> None:
        with pytest.raises(ValueError):
            volcano_plot_data(["A"], np.array([3.0]), np.array([-0.5]))

    def test_the_boundary_values_are_valid(self) -> None:
        points = volcano_plot_data(self.NAMES, self.FC, np.array([0.0, 1.0, 0.05]))

        assert [p.neg_log10_p for p in points[:2]] == [300.0, 0.0]
        assert points[2].neg_log10_p == pytest.approx(-math.log10(0.05))
        assert [p.significant for p in points] == [True, False, True]  # 0.05 <= p_threshold

    def test_a_nan_p_value_is_unknown_and_never_significant(self) -> None:
        points = volcano_plot_data(self.NAMES, self.FC, np.array([0.01, np.nan, 0.5]))

        assert math.isnan(points[1].neg_log10_p)
        assert [p.significant for p in points] == [True, False, False]

    def test_significant_is_a_plain_bool(self) -> None:
        (point,) = volcano_plot_data(["A"], np.array([2.0]), np.array([0.01]))

        assert type(point.significant) is bool

    def test_names_and_arrays_must_have_the_same_length(self) -> None:
        with pytest.raises(ValueError, match=r"same length"):
            volcano_plot_data(["A", "B"], np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.2, 0.3]))
        with pytest.raises(ValueError, match=r"same length"):
            volcano_plot_data(["A", "B"], np.array([1.0, 2.0]), np.array([0.1]))

    def test_lists_are_accepted_like_arrays(self) -> None:
        points = volcano_plot_data(["A"], [1.5], [0.01])

        assert points[0].neg_log10_p == pytest.approx(2.0)


class TestPcaArguments:
    DATA = np.random.default_rng(3).random((6, 4))  # 6 metabolites x 4 samples

    @pytest.mark.parametrize("n_components", [0, -1, 5, 100])
    def test_n_components_outside_the_possible_range_raises(self, n_components: int) -> None:
        with pytest.raises(ValueError, match=r"n_components must be between 1 and min\(n_samples, n_metabolites\) = 4"):
            pca_metabolomics(self.DATA, n_components=n_components)

    def test_the_largest_possible_number_of_components_is_allowed_and_explains_all_variance(self) -> None:
        scores, explained = pca_metabolomics(self.DATA, n_components=4)

        assert scores.shape == (4, 4) and explained.shape == (4,)
        assert explained.sum() == pytest.approx(1.0)

    def test_constant_data_keeps_matching_shapes(self) -> None:
        scores, explained = pca_metabolomics(np.ones((5, 3)), n_components=3)

        assert scores.shape == (3, 3) and explained.shape == (3,)  # explained used to be sized by the request
        assert explained.tolist() == [0.0, 0.0, 0.0]

    @pytest.mark.parametrize("bad", [np.nan, np.inf])
    def test_non_finite_values_say_how_to_fix_it(self, bad: float) -> None:
        data = self.DATA.copy()
        data[2, 1] = bad

        with pytest.raises(ValueError, match=r"must be finite.*missing_value_imputation"):
            pca_metabolomics(data, n_components=2)

    def test_scores_still_reproduce_the_centred_data(self) -> None:
        scores, explained = pca_metabolomics(self.DATA, n_components=4)
        centred = self.DATA.T - self.DATA.T.mean(axis=0)

        # full-rank PCA: the squared score norms are the singular values squared, i.e. the total variation
        assert np.sum(scores**2) == pytest.approx(np.sum(centred**2))
