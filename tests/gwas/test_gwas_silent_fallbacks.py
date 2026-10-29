"""Regression tests: GWAS code that used to substitute a value, or silently pick another input, now says so.

Each class covers one defect found while hunting for fabricated values and silent fallbacks:

* ``load_traits`` and ``load_traits_with_ids`` ignored an explicitly requested column that is not in the file: the
  first returned *another* column (the last non-ID one), the second skipped every row and returned two empty lists;
* ``create_synthetic_data`` recorded an unknown ``calibration_mode`` as given while generating the LOW dataset,
  reported one causal variant in the pure-null mode (which has none), and reseeded the module-level ``random`` state;
* ``visualization.genomic.ld.compute_ld_decay`` reseeded the module-level ``random`` state to 42 whenever it
  subsampled pairs, with no way to choose the seed;
* ``transition_transversion_plot`` reported Ts/Tv = 0.0 for a set without transversions (the ratio is undefined),
  counted an ambiguity code or a lower-case base as a transversion, and averaged 0 for windows without transversions;
* ``regional_plot`` drew a variant that had no LD value as ``r2 = 0`` and placed a result without a position at 0.

Oracles are hand-counted values and the points read straight from the figures' axes.
"""

from __future__ import annotations

import json
import logging
import math
import random
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")  # optional dependency: skip, do not fail, without it
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from metainformant.gwas.data import traits  # noqa: E402
from metainformant.gwas.simulation import synthetic  # noqa: E402
from metainformant.gwas.visualization import general  # noqa: E402
from metainformant.gwas.visualization.genomic import ld, variants  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# traits: an explicit column that does not exist
# ---------------------------------------------------------------------------


class TestTraitColumns:
    PHENOTYPES = "sample_id,body_weight,aggression\nS1,10.5,1.5\nS2,11.0,2.5\nS3,9.25,bad\n"

    def test_load_traits_rejects_a_missing_trait_column_instead_of_loading_another(self, tmp_path: Path) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text(self.PHENOTYPES)

        with pytest.raises(ValueError, match=r"trait_column 'height' is not a column.*\['sample_id', 'body_weight'"):
            traits.load_traits(path, trait_column="height")

    def test_load_traits_still_picks_the_last_non_id_column_when_none_is_requested(self, tmp_path: Path) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text(self.PHENOTYPES)

        assert traits.load_traits(path) == pytest.approx([1.5, 2.5])

    def test_load_traits_with_ids_rejects_a_missing_trait_column(self, tmp_path: Path) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text(self.PHENOTYPES)

        with pytest.raises(ValueError, match="trait_column 'height' is not a column"):
            traits.load_traits_with_ids(path, trait_column="height")

    def test_load_traits_with_ids_rejects_a_missing_id_column_instead_of_returning_nothing(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text(self.PHENOTYPES)

        with pytest.raises(ValueError, match="id_column 'nonexistent' is not a column"):
            traits.load_traits_with_ids(path, id_column="nonexistent")

    def test_requested_columns_that_exist_are_used(self, tmp_path: Path) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text(self.PHENOTYPES)

        ids, values = traits.load_traits_with_ids(path, id_column="sample_id", trait_column="body_weight")

        assert ids == ["S1", "S2", "S3"] and values == pytest.approx([10.5, 11.0, 9.25])

    def test_a_file_with_no_id_column_still_gives_empty_lists_when_none_was_requested(self, tmp_path: Path) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text("a,b\n1,2\n")

        assert traits.load_traits_with_ids(path) == ([], [])

    def test_skipped_rows_are_counted_in_a_warning(self, tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
        path = tmp_path / "pheno.csv"
        path.write_text("sample_id,aggression\nS1,1.5\nS2,bad\nS3,\nS4,4.5\n")

        with caplog.at_level(logging.WARNING):
            values = traits.load_traits(path)
            ids, with_ids = traits.load_traits_with_ids(path)

        assert values == pytest.approx([1.5, 4.5]) and ids == ["S1", "S4"] and with_ids == pytest.approx([1.5, 4.5])
        assert sum("Skipped 2 row(s)" in record.getMessage() for record in caplog.records) == 2


# ---------------------------------------------------------------------------
# the synthetic dataset generator
# ---------------------------------------------------------------------------


class TestSyntheticGenerator:
    def test_an_unknown_mode_is_an_error_not_a_low_dataset_with_the_wrong_label(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match=r"Unknown calibration_mode 'HIGHEST'.*\['HIGH', 'LOW', 'NULL'\]"):
            synthetic.create_synthetic_data(tmp_path / "out", n_samples=6, n_variants=30, calibration_mode="HIGHEST")

        assert not (tmp_path / "out" / "synthetic_params.json").exists()

    def test_a_lower_case_mode_name_is_unknown_too(self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="Unknown calibration_mode"):
            synthetic.create_synthetic_data(tmp_path, n_samples=6, n_variants=30, calibration_mode="null")

    def test_the_pure_null_mode_has_no_causal_variant(self, tmp_path: Path) -> None:
        params = synthetic.create_synthetic_data(tmp_path, n_samples=10, n_variants=60, calibration_mode="NULL", seed=3)

        assert params["n_causal"] == 0
        assert params["mean_causal_power"] == 0.0
        assert json.loads((tmp_path / "synthetic_params.json").read_text())["n_causal"] == 0

    @pytest.mark.parametrize("mode, minimum", [("LOW", 1), ("HIGH", 1)])
    def test_the_other_modes_still_have_causal_variants(self, tmp_path: Path, mode: str, minimum: int) -> None:
        params = synthetic.create_synthetic_data(tmp_path, n_samples=10, n_variants=200, calibration_mode=mode, seed=3)

        assert params["n_causal"] >= minimum

    def test_the_global_random_state_is_left_alone(self, tmp_path: Path) -> None:
        random.seed(123)
        expected = [random.random() for _ in range(3)]

        random.seed(123)
        synthetic.create_synthetic_data(tmp_path, n_samples=8, n_variants=40, seed=9)
        after = [random.random() for _ in range(3)]

        assert after == expected  # it used to restart the stream from its own seed

    def test_the_output_does_not_depend_on_the_global_random_state(self, tmp_path: Path) -> None:
        random.seed(1)
        synthetic.create_synthetic_data(tmp_path / "a", n_samples=8, n_variants=40, seed=9)
        random.seed(2)
        synthetic.create_synthetic_data(tmp_path / "b", n_samples=8, n_variants=40, seed=9)

        for name in ("synthetic.vcf", "phenotype.csv"):
            assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes()

    def test_a_different_seed_gives_a_different_dataset(self, tmp_path: Path) -> None:
        synthetic.create_synthetic_data(tmp_path / "a", n_samples=8, n_variants=40, seed=1)
        synthetic.create_synthetic_data(tmp_path / "b", n_samples=8, n_variants=40, seed=2)

        assert (tmp_path / "a" / "phenotype.csv").read_bytes() != (tmp_path / "b" / "phenotype.csv").read_bytes()


# ---------------------------------------------------------------------------
# visualization LD decay: the pair subsample has its own seeded generator
# ---------------------------------------------------------------------------


def _ld_inputs(n_variants: int = 448, n_samples: int = 8) -> tuple[List[List[int]], List[int]]:
    """448 variants give 100 128 pairs, just over the 100 000 that trigger the subsample."""
    rng = np.random.default_rng(4)
    genotypes = rng.integers(0, 3, size=(n_variants, n_samples)).tolist()
    positions = [1000 * (i + 1) for i in range(n_variants)]
    return genotypes, positions


class TestVisualizationLdDecay:
    def test_the_same_seed_gives_the_same_curve_whatever_the_global_state(self) -> None:
        genotypes, positions = _ld_inputs()

        random.seed(1)
        first = ld.compute_ld_decay(genotypes, positions, max_distance=500_000, n_bins=5, seed=7)
        random.seed(2)
        second = ld.compute_ld_decay(genotypes, positions, max_distance=500_000, n_bins=5, seed=7)

        assert first["total_pairs"] == 100_000
        assert first == second

    def test_a_different_seed_draws_other_pairs(self) -> None:
        genotypes, positions = _ld_inputs()

        one = ld.compute_ld_decay(genotypes, positions, max_distance=500_000, n_bins=5, seed=1)
        two = ld.compute_ld_decay(genotypes, positions, max_distance=500_000, n_bins=5, seed=2)

        assert one["n_pairs"] != two["n_pairs"] or one["mean_r2"] != two["mean_r2"]

    def test_the_global_random_state_is_left_alone(self) -> None:
        genotypes, positions = _ld_inputs()
        random.seed(5)
        expected = random.random()

        random.seed(5)
        ld.compute_ld_decay(genotypes, positions, max_distance=500_000, n_bins=5)
        after = random.random()

        assert after == expected  # it used to be reseeded to 42


# ---------------------------------------------------------------------------
# transition/transversion plot
# ---------------------------------------------------------------------------


def _records(*pairs: tuple[str, str]) -> List[Dict[str, str]]:
    return [{"REF": ref, "ALT": alt} for ref, alt in pairs]


class TestTransitionTransversion:
    def test_without_a_transversion_the_ratio_is_undefined_not_zero(self, tmp_path: Path) -> None:
        out = tmp_path / "tstv.png"

        result = variants.transition_transversion_plot(_records(("A", "G"), ("C", "T"), ("G", "A")), out)

        assert result["status"] == "success"
        assert result["ts_tv_ratio"] is None
        assert (result["transitions"], result["transversions"]) == (3, 0)
        assert out.stat().st_size > 0

    def test_the_ratio_is_transitions_over_transversions(self) -> None:
        result = variants.transition_transversion_plot(
            _records(("A", "G"), ("C", "T"), ("G", "A"), ("A", "C"), ("T", "G"))
        )

        assert result["ts_tv_ratio"] == pytest.approx(1.5)
        assert (result["transitions"], result["transversions"]) == (3, 2)

    def test_lower_case_bases_are_classified_like_upper_case(self) -> None:
        result = variants.transition_transversion_plot(_records(("a", "g"), ("c", "t"), ("a", "c")))

        assert (result["transitions"], result["transversions"]) == (2, 1)

    def test_what_is_not_a_substitution_between_bases_is_unclassified(self) -> None:
        result = variants.transition_transversion_plot(
            _records(
                ("A", "G"),
                ("A", "N"),
                ("A", "A"),
                ("A", "."),
                ("A", "AT"),
                ("AT", "A"),
                ("A", "G,T"),
                ("A", ""),
                ("", "G"),
            )
        )

        assert (result["transitions"], result["transversions"]) == (1, 0)  # N, '.', equal bases and indels are not Tv
        assert result["n_unclassified"] == 8
        assert result["n_variants"] == 9

    def test_windows_without_a_transversion_are_left_out_of_the_mean(self, tmp_path: Path) -> None:
        data = {
            "transitions": [6, 4, 3, 9],
            "transversions": [3, 0, 2, 3],
            "positions": [100, 200, 300, 400],
        }

        result = variants.transition_transversion_plot(data, tmp_path / "w.png")

        assert result["ts_tv_ratio"] == pytest.approx((2.0 + 1.5 + 3.0) / 3)  # not (2 + 0 + 1.5 + 3) / 4
        assert result["n_windows_without_transversions"] == 1

    def test_no_window_with_a_transversion_gives_no_mean(self) -> None:
        result = variants.transition_transversion_plot({"transitions": [3, 4], "transversions": [0, 0]})

        assert result["ts_tv_ratio"] is None
        assert result["n_windows_without_transversions"] == 2


# ---------------------------------------------------------------------------
# regional_plot: missing LD values and unplaced results
# ---------------------------------------------------------------------------

RESULTS = [
    {"chrom": "1", "pos": 100, "p_value": 1e-8, "variant_id": "rs1"},
    {"chrom": "1", "pos": 200, "p_value": 1e-6, "variant_id": "rs2"},
    {"chrom": "1", "pos": 300, "p_value": 1e-4, "variant_id": "rs3"},
]


def _scatters(fig: Any) -> List[PathCollection]:
    return [c for c in fig.axes[0].collections if isinstance(c, PathCollection)]


class TestRegionalPlotLd:
    def test_a_variant_without_an_ld_value_is_grey_not_r2_zero(self) -> None:
        fig = general.regional_plot(RESULTS, "1", 0, 1000, ld_values={"rs1": 1.0, "rs2": 0.6})

        grey, coloured = _scatters(fig)[0], _scatters(fig)[1]
        assert np.asarray(grey.get_offsets())[:, 0].tolist() == [300]  # rs3 has no LD value
        assert np.asarray(coloured.get_offsets())[:, 0].tolist() == [100, 200]
        assert np.asarray(coloured.get_array()).tolist() == pytest.approx([1.0, 0.6])  # nothing was set to 0
        assert coloured.get_clim() == (0.0, 1.0)

    def test_every_variant_has_an_ld_value(self) -> None:
        fig = general.regional_plot(RESULTS, "1", 0, 1000, ld_values={"rs1": 1.0, "rs2": 0.6, "rs3": 0.0})

        scatters = _scatters(fig)
        assert len(scatters) == 1
        assert np.asarray(scatters[0].get_array()).tolist() == pytest.approx([1.0, 0.6, 0.0])  # a real 0 stays 0

    def test_ld_values_can_be_keyed_by_chromosome_and_position(self) -> None:
        fig = general.regional_plot(RESULTS, "1", 0, 1000, ld_values={"1:100": 0.9})

        grey, coloured = _scatters(fig)[0], _scatters(fig)[1]
        assert np.asarray(grey.get_offsets())[:, 0].tolist() == [200, 300]
        assert np.asarray(coloured.get_array()).tolist() == pytest.approx([0.9])

    def test_non_numeric_ld_values_count_as_missing(self) -> None:
        fig = general.regional_plot(RESULTS, "1", 0, 1000, ld_values={"rs1": "NA", "rs2": None, "rs3": 0.2})

        grey, coloured = _scatters(fig)[0], _scatters(fig)[1]
        assert np.asarray(grey.get_offsets())[:, 0].tolist() == [100, 200]
        assert np.asarray(coloured.get_array()).tolist() == pytest.approx([0.2])

    def test_without_ld_values_the_plot_is_a_single_plain_scatter(self) -> None:
        fig = general.regional_plot(RESULTS, "1", 0, 1000)

        assert len(_scatters(fig)) == 1

    def test_a_result_without_a_position_is_an_error_not_a_variant_at_position_zero(self) -> None:
        results = RESULTS + [{"chrom": "1", "p_value": 1e-9}]

        with pytest.raises(
            ValueError, match=r"regional_plot: result 3 has no position \(accepted keys: 'position', 'pos'"
        ):
            general.regional_plot(results, "1", 0, 1000)

    def test_a_result_without_a_chromosome_is_an_error_not_a_silent_omission(self) -> None:
        results = RESULTS + [{"pos": 150, "p_value": 1e-9}]

        with pytest.raises(
            ValueError, match=r"regional_plot: result 3 has no chromosome \(accepted keys: 'chromosome', 'chrom'"
        ):
            general.regional_plot(results, "1", 0, 1000)

    def test_the_chromosome_may_be_given_as_a_number(self) -> None:
        fig = general.regional_plot([{"CHR": 1, "BP": 100, "P": 1e-5}], 1.0, 0, 1000)

        assert fig is not None and len(_scatters(fig)) == 1
        assert math.isclose(float(np.asarray(_scatters(fig)[0].get_offsets())[0, 1]), 5.0)
