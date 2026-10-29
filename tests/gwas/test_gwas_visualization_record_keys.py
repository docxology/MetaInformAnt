"""Regression tests: GWAS plots read every spelling of a record, and refuse (not invent) what a record lacks.

Defects, each reproduced from the code before it was fixed:

* the regional plots raised a bare ``KeyError: 'CHR'`` for results with no chromosome key (and ``'BP'`` / ``'P'``
  for results with no position or p-value); they now raise ``ValueError`` naming the accepted keys;
* ``start`` and ``end`` of the regional plots were ignored whenever ``chrom`` was not given, although the
  docstring says an omitted end is open and a given one bounds the region;
* ``recombination_rate`` of ``regional_plot`` was silently ignored unless it was a DataFrame, and raised
  ``KeyError`` when its ``BP`` / ``RATE`` columns were missing;
* ``genome.manhattan_plot`` accepted every chromosome spelling but read the position only from ``POS`` / ``pos``
  (``BP`` and ``position`` gave position 0) and the p-value only from ``p_value`` / ``P`` (``pvalue`` and ``pval``
  gave p = 1, drawn at -log10 p = 0); a record with no p-value was drawn as p = 1, one with no chromosome on
  chromosome 1, one with no position at position 0;
* ``effect_direction_plot`` had the same position and chromosome defaults and drew a missing beta as an effect of 0;
* ``circular_manhattan_plot`` and ``chromosome_ideogram`` defaulted the same way (and the circular plot rejected a
  DataFrame whose columns were not literally ``CHR`` / ``BP`` / ``P``);
* ``generate_all_plots`` built its records with ``iterrows()``, which upcasts an integer chromosome next to a float
  p-value to ``1.0``, so the regional plots of a table with integer chromosomes were named ``regional_plot_chr1.0_...``.

Oracles: the points of a figure are read straight from its axes, with hand-counted values.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")  # optional dependency: skip, do not fail, without it
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402

from metainformant.gwas.visualization import utils as viz_utils  # noqa: E402
from metainformant.gwas.visualization.genomic import genome, regional, variants  # noqa: E402
from metainformant.gwas.visualization.interactive import suite  # noqa: E402
from metainformant.gwas.visualization.statistical import effects  # noqa: E402

NAN = float("nan")
CHROMOSOME_KEYS = ("chromosome", "chrom", "chr", "CHROM", "CHR")
POSITION_KEYS = ("position", "pos", "POS", "BP")
P_VALUE_KEYS = ("p_value", "pvalue", "pval", "P")


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def _scatters(ax) -> List[PathCollection]:
    return [c for c in ax.collections if isinstance(c, PathCollection)]


def _points(ax) -> np.ndarray:
    """All scatter points of an axes as an (n, 2) array, in drawing order."""
    return np.vstack([np.asarray(c.get_offsets(), dtype=float) for c in _scatters(ax)])


def _keys_named(message: str, keys) -> bool:
    return all(f"'{key}'" in message for key in keys)


# ---------------------------------------------------------------------------
# 1. Regional plots: ValueError naming the accepted keys, not KeyError
# ---------------------------------------------------------------------------

REGIONAL_PLOTS = [
    pytest.param(regional.regional_plot, id="regional_plot"),
    pytest.param(regional.gene_annotation_plot, id="gene_annotation_plot"),
]


@pytest.mark.parametrize("plot", REGIONAL_PLOTS)
class TestRegionalPlotsNameTheAcceptedKeys:
    def test_records_without_a_chromosome_key(self, plot):
        records = [{"pos": 1000, "p_value": 1e-8}, {"pos": 1500, "p_value": 1e-6}]

        with pytest.raises(ValueError) as excinfo:
            plot(records, chrom="1", start=500, end=2500)

        message = str(excinfo.value)
        assert "chromosome" in message
        assert _keys_named(message, CHROMOSOME_KEYS)

    def test_dataframe_without_a_chromosome_column(self, plot):
        frame = pd.DataFrame({"BP": [1000, 1500], "P": [1e-8, 1e-6]})

        with pytest.raises(ValueError) as excinfo:
            plot(frame, chrom=1)

        assert _keys_named(str(excinfo.value), CHROMOSOME_KEYS)

    def test_the_message_lists_the_keys_that_are_present(self, plot):
        with pytest.raises(ValueError, match=r"keys present are \['p_value', 'pos'\]"):
            plot([{"pos": 1, "p_value": 0.1}], chrom="1")

    def test_it_is_not_a_key_error(self, plot):
        with pytest.raises(Exception) as excinfo:
            plot([{"pos": 1000, "p_value": 1e-8}], chrom="1")

        assert not isinstance(excinfo.value, KeyError)

    def test_records_without_a_position_key(self, plot):
        with pytest.raises(ValueError) as excinfo:
            plot([{"chrom": "1", "p_value": 1e-8}], chrom="1")

        message = str(excinfo.value)
        assert "position" in message and _keys_named(message, POSITION_KEYS)

    def test_records_without_a_p_value_key(self, plot):
        with pytest.raises(ValueError) as excinfo:
            plot([{"chrom": "1", "pos": 1000}], chrom="1")

        message = str(excinfo.value)
        assert "p-value" in message and _keys_named(message, P_VALUE_KEYS)

    def test_every_missing_column_is_named_at_once(self, plot):
        with pytest.raises(ValueError) as excinfo:
            plot([{"id": "rs1"}], chrom="1")

        message = str(excinfo.value)
        assert _keys_named(message, CHROMOSOME_KEYS + POSITION_KEYS + P_VALUE_KEYS)

    def test_no_chromosome_is_needed_when_none_is_asked_for(self, plot):
        result = plot([{"pos": 1000, "p_value": 1e-8}, {"pos": 1500, "p_value": 1e-6}])

        assert result["status"] == "success" and result["n_variants"] == 2

    def test_an_empty_input_is_an_empty_region_not_a_missing_key(self, plot):
        result = plot([], chrom="1", start=0, end=10)

        assert result["status"] == "success" and result["n_variants"] == 0

    def test_start_and_end_bound_the_region_without_a_chromosome(self, plot):
        records = [{"pos": position, "p_value": 1e-3} for position in (100, 1000, 1500, 9000)]

        assert plot(records, start=500, end=2500)["n_variants"] == 2
        assert plot(records, start=1200)["n_variants"] == 2
        assert plot(records, end=1200)["n_variants"] == 2

    def test_results_without_a_p_value_are_neither_drawn_nor_counted(self, plot):
        records = [
            {"chrom": "1", "pos": 100, "p_value": 1e-8},
            {"chrom": "1", "pos": 200, "p_value": NAN},
            {"chrom": "1", "pos": 300, "p_value": None},
            {"chrom": "1", "pos": 400, "p_value": "NA"},
            {"chrom": "1", "pos": 500, "p_value": 0.2},
        ]

        result = plot(records, chrom="1")

        assert result["n_variants"] == 2


class TestRecombinationRateIsValidated:
    FRAME = pd.DataFrame({"CHR": [1, 1], "BP": [1000, 1500], "P": [1e-8, 1e-6]})

    def test_a_table_without_rate_is_an_error_that_names_the_columns(self):
        with pytest.raises(ValueError, match=r"recombination_rate has no \['RATE'\]"):
            regional.regional_plot(self.FRAME, chrom=1, recombination_rate=pd.DataFrame({"BP": [1000]}))

    def test_a_table_without_position_is_an_error(self):
        with pytest.raises(ValueError, match=r"recombination_rate has no \['BP'\]"):
            regional.regional_plot(self.FRAME, chrom=1, recombination_rate=pd.DataFrame({"RATE": [1.0]}))

    @pytest.mark.parametrize("not_a_table", [[(1000, 2.0)], {"BP": [1000], "RATE": [2.0]}, "rates.tsv"])
    def test_a_value_that_is_not_a_table_is_not_silently_ignored(self, not_a_table):
        with pytest.raises(ValueError, match="recombination_rate must be a DataFrame"):
            regional.regional_plot(self.FRAME, chrom=1, recombination_rate=not_a_table)

    def test_a_valid_table_still_draws_the_track(self):
        rates = pd.DataFrame({"CHR": [1, 1], "BP": [1000, 1500], "RATE": [2.0, 3.0]})

        result = regional.regional_plot(self.FRAME, chrom=1, recombination_rate=rates)

        assert result["status"] == "success"


# ---------------------------------------------------------------------------
# 2. genome.manhattan_plot: every spelling, no invented values
# ---------------------------------------------------------------------------


def _manhattan(rows) -> Dict[str, Any]:
    result = genome.manhattan_plot(rows)
    return result


@pytest.mark.parametrize("key", POSITION_KEYS)
def test_manhattan_reads_every_position_spelling(key):
    rows = [
        {"CHROM": "1", key: 100, "p_value": 1e-3},
        {"CHROM": "1", key: 300, "p_value": 1e-4},
        {"CHROM": "2", key: 50, "p_value": 1e-5},
    ]

    result = _manhattan(rows)

    assert result["status"] == "success", result
    scatters = _scatters(plt.gcf().axes[0])
    assert np.asarray(scatters[0].get_offsets())[:, 0].tolist() == [100, 300]  # not both at position 0
    assert np.asarray(scatters[1].get_offsets())[:, 0].tolist() == [351]  # chromosome 2 starts after chromosome 1's end


@pytest.mark.parametrize("key", P_VALUE_KEYS)
def test_manhattan_reads_every_p_value_spelling(key):
    rows = [{"CHROM": "1", "POS": 100, key: 1e-7}, {"CHROM": "1", "POS": 200, key: 1e-3}]

    result = _manhattan(rows)

    assert result["status"] == "success", result
    assert _points(plt.gcf().axes[0])[:, 1] == pytest.approx([7.0, 3.0])  # not flat at -log10 p = 0


def test_manhattan_reads_mixed_spellings_in_one_list():
    rows = [
        {"chrom": "1", "pos": 100, "pvalue": 1e-2},
        {"CHR": "1", "BP": 200, "P": 1e-4},
        {"chromosome": "1", "position": 300, "pval": 1e-6},
    ]

    result = _manhattan(rows)

    assert result["n_variants"] == 3
    assert _points(plt.gcf().axes[0]).tolist() == [[100, 2.0], [200, 4.0], [300, 6.0]]


def test_manhattan_does_not_draw_a_result_without_a_p_value_as_p_equal_one():
    rows = [
        {"CHROM": "1", "POS": 100, "p_value": 1e-6},
        {"CHROM": "1", "POS": 200, "p_value": NAN},
        {"CHROM": "1", "POS": 300, "p_value": None},
        {"CHROM": "1", "POS": 400, "p_value": "NA"},
        {"CHROM": "1", "POS": 500},
        {"CHROM": "1", "POS": 600, "p_value": 1e-2},
    ]

    result = _manhattan(rows)

    assert result["status"] == "success"
    assert result["n_variants"] == 2 and result["n_without_p_value"] == 4
    assert _points(plt.gcf().axes[0]).tolist() == [[100, 6.0], [600, 2.0]]  # nothing at -log10 p = 0


def test_manhattan_with_no_p_value_at_all_fails_instead_of_drawing_an_empty_plot():
    result = _manhattan([{"CHROM": "1", "POS": 100}, {"CHROM": "1", "POS": 200, "p_value": NAN}])

    assert result["status"] == "failed"
    assert "p-value" in result["error"]


def test_manhattan_names_the_accepted_chromosome_keys_for_a_result_without_one():
    result = _manhattan([{"CHROM": "1", "POS": 100, "p_value": 1e-3}, {"POS": 200, "p_value": 1e-4}])

    assert result["status"] == "failed"
    assert "result 1 has no chromosome" in result["error"]
    assert _keys_named(result["error"], CHROMOSOME_KEYS)


def test_manhattan_names_the_accepted_position_keys_for_a_result_without_one():
    result = _manhattan([{"CHROM": "1", "p_value": 1e-3}])

    assert result["status"] == "failed"
    assert "result 0 has no position" in result["error"]
    assert _keys_named(result["error"], POSITION_KEYS)


def test_manhattan_does_not_move_an_unplaced_result_to_chromosome_1_position_0():
    """The old code returned success for this, with the variant on chromosome 1 at position 0."""
    result = _manhattan([{"p_value": 1e-8}])

    assert result["status"] == "failed"


def test_manhattan_chromosome_labels_are_still_normalised():
    rows = [{"CHROM": 1.0, "POS": 5, "p_value": 0.1}, {"CHROM": 2, "POS": 5, "p_value": 0.1}]

    assert _manhattan(rows)["chromosomes"] == ["1", "2"]


# ---------------------------------------------------------------------------
# 3. effect_direction_plot
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key", POSITION_KEYS)
def test_effect_direction_reads_every_position_spelling(key):
    rows = [
        {"CHROM": "1", key: 100, "BETA": 0.5},
        {"CHROM": "1", key: 300, "BETA": -0.25},
        {"CHROM": "2", key: 50, "BETA": 0.1},
    ]

    result = effects.effect_direction_plot(rows)

    assert result["status"] == "success", result
    assert _points(plt.gcf().axes[0]).tolist() == [[100, 0.5], [300, -0.25], [351, 0.1]]


@pytest.mark.parametrize("key", ["BETA", "beta"])
def test_effect_direction_reads_both_beta_spellings(key):
    result = effects.effect_direction_plot([{"CHROM": "1", "POS": 1, key: -0.5}, {"CHROM": "1", "POS": 2, key: 0.5}])

    assert (result["positive_effects"], result["negative_effects"]) == (1, 1)


def test_effect_direction_does_not_draw_a_missing_beta_as_an_effect_of_zero():
    rows = [
        {"CHROM": "1", "POS": 100, "BETA": 0.5},
        {"CHROM": "1", "POS": 200},
        {"CHROM": "1", "POS": 300, "beta": None},
        {"CHROM": "1", "POS": 400, "BETA": NAN},
        {"CHROM": "1", "POS": 500, "BETA": "NA"},
    ]

    result = effects.effect_direction_plot(rows)

    assert result["n_variants"] == 1 and result["n_without_effect"] == 4
    assert _points(plt.gcf().axes[0]).tolist() == [[100, 0.5]]


def test_effect_direction_names_the_accepted_keys_for_a_result_that_cannot_be_placed():
    result = effects.effect_direction_plot([{"BETA": 0.5}])

    assert result["status"] == "failed"
    assert _keys_named(result["error"], CHROMOSOME_KEYS) and _keys_named(result["error"], POSITION_KEYS)


def test_effect_direction_with_no_beta_at_all_fails():
    assert effects.effect_direction_plot([{"CHROM": "1", "POS": 1}])["status"] == "failed"


# ---------------------------------------------------------------------------
# 4. circular_manhattan_plot and chromosome_ideogram
# ---------------------------------------------------------------------------


def test_circular_manhattan_reads_every_spelling_and_skips_results_without_a_p_value():
    rows = [
        {"chrom": "1", "pos": 100, "pvalue": 1e-6},
        {"CHR": "1", "BP": 200, "P": 1e-3},
        {"CHROM": "2", "position": 50, "pval": 1e-4},
        {"CHROM": "2", "POS": 70, "p_value": NAN},
        {"CHROM": "2", "POS": 90},
    ]

    result = genome.circular_manhattan_plot(rows)

    assert result["status"] == "success"
    assert result["n_variants"] == 3  # the two without a p-value are not drawn (they were drawn as p = 1)
    assert result["n_chromosomes"] == 2


def test_circular_manhattan_accepts_a_dataframe_with_the_pipelines_columns():
    frame = pd.DataFrame({"chrom": ["1", "1", "2"], "pos": [100, 200, 50], "p_value": [1e-6, 1e-3, 1e-4]})

    result = genome.circular_manhattan_plot(frame)

    assert result["status"] == "success" and result["n_variants"] == 3


def test_circular_manhattan_dataframe_error_names_the_accepted_keys():
    result = genome.circular_manhattan_plot(pd.DataFrame({"pos": [1], "p_value": [0.1]}))

    assert result["status"] == "failed"
    assert _keys_named(result["error"], CHROMOSOME_KEYS)


def test_circular_manhattan_does_not_place_an_unplaced_result_on_chromosome_1():
    result = genome.circular_manhattan_plot([{"POS": 5, "p_value": 1e-3}])

    assert result["status"] == "failed"
    assert _keys_named(result["error"], CHROMOSOME_KEYS)


def test_circular_manhattan_with_no_p_value_at_all_fails():
    assert genome.circular_manhattan_plot([{"CHROM": "1", "POS": 5}])["status"] == "failed"


def test_ideogram_infers_lengths_from_every_position_spelling():
    rows = [
        {"CHROM": "1", "position": 3_000_000},
        {"chrom": "1", "BP": 1_000_000},
        {"CHR": "2", "pos": 2_000_000},
    ]

    result = genome.chromosome_ideogram(rows)

    assert result["status"] == "success"
    assert result["chromosomes"] == ["1", "2"]
    lengths = [text.get_text() for text in plt.gcf().axes[0].texts if text.get_text().endswith(" Mb")]
    assert lengths == ["3.0 Mb", "2.0 Mb"]  # the ``position`` spelling used to read as 0


def test_ideogram_does_not_invent_a_chromosome_or_a_length():
    for row in ({"POS": 5}, {"CHROM": "1"}):
        result = genome.chromosome_ideogram([row])
        assert result["status"] == "failed"
        assert "accepted keys" in result["error"]


# ---------------------------------------------------------------------------
# 5. The shared helpers
# ---------------------------------------------------------------------------


class TestSharedHelpers:
    def test_get_p_value_reads_every_spelling(self):
        for key in P_VALUE_KEYS:
            assert viz_utils.get_p_value({key: 0.25}) == 0.25

    def test_get_p_value_priority_follows_p_value_keys(self):
        assert viz_utils.get_p_value({"P": 0.5, "p_value": 0.1}) == 0.1

    def test_a_missing_or_unusable_p_value_is_nan_never_one(self):
        for record in ({}, {"p_value": None}, {"p_value": "NA"}, {"p_value": NAN}, {"p_value": math.inf}, {"P": True}):
            assert math.isnan(viz_utils.get_p_value(record)), record

    def test_the_first_spelling_that_holds_a_value_decides(self):
        """A failed test with ``p_value`` NaN stays missing; it does not fall through to another spelling."""
        assert math.isnan(viz_utils.get_p_value({"p_value": NAN, "P": 0.5}))

    def test_get_beta(self):
        assert viz_utils.get_beta({"beta": -1.5}) == -1.5
        assert viz_utils.get_beta({"BETA": "0.25"}) == 0.25
        assert math.isnan(viz_utils.get_beta({}))
        assert math.isnan(viz_utils.get_beta({"beta": "NA"}))

    def test_record_placement(self):
        assert viz_utils.record_placement({"CHR": 1.0, "BP": "250"}, 0, plot="p") == ("1", 250.0)

    @pytest.mark.parametrize(
        "record", [{"pos": 1}, {"chrom": None, "pos": 1}, {"chrom": NAN, "pos": 1}, {"chrom": " ", "pos": 1}]
    )
    def test_record_placement_without_a_chromosome(self, record):
        with pytest.raises(ValueError, match="plot_name: result 7 has no chromosome"):
            viz_utils.record_placement(record, 7, plot="plot_name")

    @pytest.mark.parametrize("record", [{"chrom": "1"}, {"chrom": "1", "pos": None}, {"chrom": "1", "pos": NAN}])
    def test_record_placement_without_a_position(self, record):
        with pytest.raises(ValueError, match="has no position"):
            viz_utils.record_placement(record, 0, plot="plot_name")

    def test_record_placement_with_a_position_that_is_not_a_number(self):
        with pytest.raises(ValueError, match="not a finite number: 'ten'"):
            viz_utils.record_placement({"chrom": "1", "pos": "ten"}, 0, plot="plot_name")

    def test_require_result_columns_passes_when_every_column_is_there(self):
        viz_utils.require_result_columns(
            pd.DataFrame({"CHR": [1], "BP": [1], "P": [0.1]}), ("CHR", "BP", "P"), plot="x"
        )

    def test_require_result_columns_names_only_what_is_missing(self):
        with pytest.raises(ValueError) as excinfo:
            viz_utils.require_result_columns(pd.DataFrame({"CHR": [1], "BP": [1]}), ("CHR", "BP", "P"), plot="x")

        message = str(excinfo.value)
        assert _keys_named(message, P_VALUE_KEYS)
        assert "'chromosome'" not in message and "'POS'" not in message


# ---------------------------------------------------------------------------
# 6. generate_all_plots: regional file names from chromosome labels, not upcast floats
# ---------------------------------------------------------------------------


def _association_table(path: Path, chromosomes: List[Any], beta: bool = False) -> Path:
    """A tab-separated results table: three chromosomes with five variants each, p-values from 1e-1 to 1e-15."""
    header = ["CHROM", "POS", "p_value"] + (["BETA"] if beta else [])
    lines = ["\t".join(header)]
    for c, chromosome in enumerate(chromosomes):
        for i in range(5):
            cells = [str(chromosome), str(1000 * (i + 1) + c), f"{10.0 ** -(1 + 3 * c + i):.3g}"]
            lines.append("\t".join(cells + (["0.5"] if beta else [])))
    path.write_text("\n".join(lines) + "\n")
    return path


def _regional_names(plots: Path) -> List[str]:
    return sorted(p.name for p in plots.glob("regional_plot_*.png"))


class TestRegionalFileNames:
    def test_integer_chromosomes_are_not_named_with_a_decimal_point(self, tmp_path):
        table = _association_table(tmp_path / "assoc.tsv", [1, 2, 3])

        result = suite.generate_all_plots(table, tmp_path / "plots")

        assert result["status"] == "completed"
        names = _regional_names(tmp_path / "plots")
        assert len(names) == 5
        assert not any(".0_" in name for name in names), names
        assert all(name.startswith(("regional_plot_chr3_", "regional_plot_chr2_")) for name in names), names
        assert "regional_plot_chr3_5002.png" in names  # the best hit overall: chromosome 3, position 5002, p = 1e-11

    def test_the_regional_plots_are_keyed_in_the_result_by_rank(self, tmp_path):
        table = _association_table(tmp_path / "assoc.tsv", [1, 2, 3])

        result = suite.generate_all_plots(table, tmp_path / "plots")

        assert [key for key in result["plots"] if key.startswith("regional_")] == [f"regional_{i}" for i in range(1, 6)]
        assert all(
            Path(path).stat().st_size > 0 for key, path in result["plots"].items() if key.startswith("regional_")
        )

    def test_a_chr_prefix_is_not_doubled(self, tmp_path):
        table = _association_table(tmp_path / "assoc.tsv", ["chr1", "chr2", "chr3"])

        suite.generate_all_plots(table, tmp_path / "plots")

        names = _regional_names(tmp_path / "plots")
        assert names and not any("chrchr" in name for name in names), names

    def test_a_contig_name_cannot_reach_outside_the_output_directory(self, tmp_path):
        table = _association_table(tmp_path / "assoc.tsv", ["../escape", "a/b", "c|d"])

        suite.generate_all_plots(table, tmp_path / "plots")

        assert not (tmp_path / "escape_1000.png").exists()
        written = [p for p in tmp_path.rglob("regional_plot_*.png")]
        assert written and all(p.parent == tmp_path / "plots" for p in written)

    def test_a_table_with_a_beta_column_is_named_the_same_way(self, tmp_path):
        table = _association_table(tmp_path / "assoc.tsv", [1, 2, 3], beta=True)

        result = suite.generate_all_plots(table, tmp_path / "plots")

        assert result["status"] == "completed"
        assert not any(".0_" in name for name in _regional_names(tmp_path / "plots"))


class TestChromosomeStem:
    @pytest.mark.parametrize(
        "label, stem",
        [("1", "1"), ("chr1", "1"), ("Chr12", "12"), ("X", "X"), ("chr", "chr"), ("NC_037638.1", "NC_037638.1")],
    )
    def test_stem(self, label, stem):
        assert suite._chromosome_stem(label) == stem

    @pytest.mark.parametrize("label", ["a/b", "..\\x", "lcl|NC_1", "with space", "tab\t"])
    def test_path_characters_are_replaced(self, label):
        stem = suite._chromosome_stem(label)

        assert all(ch.isalnum() or ch in "._-" for ch in stem)
        assert "/" not in stem and "\\" not in stem


# ---------------------------------------------------------------------------
# variant_density_plot
# ---------------------------------------------------------------------------


class TestVariantDensity:
    def test_every_spelling_is_read(self):
        rows = [
            {"CHROM": "1", "POS": 100},
            {"chrom": 1.0, "pos": 200},
            {"chromosome": "2", "position": 300},
            {"CHR": "2", "BP": 400},
        ]

        result = variants.variant_density_plot(rows)

        assert result == {"status": "success", "n_variants": 4, "n_chromosomes": 2}

    def test_a_variant_without_a_chromosome_is_not_moved_to_chromosome_1(self):
        result = variants.variant_density_plot([{"CHROM": "1", "POS": 100}, {"POS": 200}])

        assert result["status"] == "failed"
        assert "result 1 has no chromosome" in result["error"]
        assert _keys_named(result["error"], CHROMOSOME_KEYS)

    def test_a_variant_without_a_position_is_not_moved_to_position_0(self):
        result = variants.variant_density_plot([{"CHROM": "1"}])

        assert result["status"] == "failed"
        assert "result 0 has no position" in result["error"]
        assert _keys_named(result["error"], POSITION_KEYS)

    def test_a_plain_list_of_positions_is_unchanged(self):
        assert variants.variant_density_plot([10, 20, 30])["n_chromosomes"] == 1
