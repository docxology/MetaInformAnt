"""Regression tests: the interactive Manhattan plot must honour each variant's chromosome.

Defect: ``interactive_manhattan`` looked the chromosome up under ``chromosome``,
``chr`` and ``CHROM`` only and silently defaulted to ``"1"``. Results from
``run_*_gwas``, VCF variant records and ``generate_all_plots`` carry it as
``chrom``, so every variant was drawn on chromosome 1 and the cumulative
position offsets collapsed to zero.

The data preparation shared by the Plotly and fallback paths is tested directly,
so these tests hold whether or not Plotly is installed. The end-to-end checks read
the fallback HTML table and are skipped when Plotly (which renders a figure
instead) is available.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest

from metainformant.gwas.visualization.interactive import interactive as interactive_module

HAS_PLOTLY = interactive_module.HAS_PLOTLY


def _prepare(results: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str], int]:
    return interactive_module._prepare_manhattan_records(results)


def _fallback_rows(html_path: Path) -> list[tuple[str, ...]]:
    """(variant, chromosome, position) of every row of the fallback HTML table."""
    html = html_path.read_text(encoding="utf-8")
    return re.findall(r"<tr><td>(.*?)</td><td>(.*?)</td><td>(.*?)</td>", html)


class TestChromosomeAssignment:
    """The data-preparation step behind both plot paths."""

    def test_chrom_key_assigns_each_variant_to_its_chromosome(self) -> None:
        results = [
            {"chrom": "1", "pos": 100, "p_value": 1e-3, "variant_id": "a"},
            {"chrom": "2", "pos": 100, "p_value": 1e-4, "variant_id": "b"},
            {"chrom": "X", "pos": 100, "p_value": 1e-9, "variant_id": "c"},
        ]
        records, chromosomes, n_skipped = _prepare(results)

        assert {r["variant_id"]: r["chromosome"] for r in records} == {"a": "1", "b": "2", "c": "X"}
        assert chromosomes == ["1", "2", "X"]
        assert n_skipped == 0

    @pytest.mark.parametrize("key", ["chromosome", "chrom", "chr", "CHROM", "CHR"])
    @pytest.mark.parametrize("position_key", ["position", "pos", "POS"])
    def test_every_documented_key_spelling_is_honoured(self, key: str, position_key: str) -> None:
        results = [
            {key: "3", position_key: 10, "p_value": 0.1, "variant_id": "v3"},
            {key: "11", position_key: 10, "p_value": 0.2, "variant_id": "v11"},
        ]
        records, chromosomes, _ = _prepare(results)

        assert [r["chromosome"] for r in records] == ["3", "11"]
        assert chromosomes == ["3", "11"]  # numeric order, not alphabetical

    def test_cumulative_positions_lay_chromosomes_end_to_end(self) -> None:
        results = [
            {"chrom": "1", "pos": 100, "p_value": 0.1},
            {"chrom": "1", "pos": 900, "p_value": 0.1},
            {"chrom": "2", "pos": 50, "p_value": 0.1},
            {"chrom": "X", "pos": 10, "p_value": 0.1},
        ]
        records, chromosomes, _ = _prepare(results)

        # chr1 spans 0..900 (offset 0), chr2 starts at 901 and spans 51 more, chrX starts at 952.
        assert [r["cum_pos"] for r in records] == [100, 900, 951, 962]
        assert chromosomes == ["1", "2", "X"]

    def test_default_variant_id_names_the_assigned_chromosome(self) -> None:
        records, _, _ = _prepare([{"chrom": "7", "pos": 5, "p_value": 0.2}])
        assert records[0]["variant_id"] == "7:5"

    def test_integral_floats_group_with_their_string_form(self) -> None:
        """A numeric chromosome column with gaps arrives as floats (1.0)."""
        results = [
            {"chrom": 1, "pos": 1, "p_value": 0.1},
            {"chrom": 1.0, "pos": 2, "p_value": 0.1},
            {"chrom": "1", "pos": 3, "p_value": 0.1},
            {"chrom": 2.0, "pos": 4, "p_value": 0.1},
            {"chrom": " chr3 ", "pos": 5, "p_value": 0.1},
        ]
        records, chromosomes, _ = _prepare(results)

        assert [r["chromosome"] for r in records] == ["1", "1", "1", "2", "chr3"]
        assert chromosomes == ["1", "2", "chr3"]

    def test_unplaceable_results_are_skipped_not_put_on_chromosome_one(self) -> None:
        results = [
            {"pos": 1, "p_value": 0.1},  # no chromosome
            {"chrom": "1", "p_value": 0.1},  # no position
            {"chrom": "1", "pos": 1},  # no p-value
            {"chrom": "1", "pos": 2, "p_value": 0.0},  # p = 0 cannot be placed on a -log10 axis
            {"chrom": "1", "pos": 3, "p_value": float("nan")},
            {"chrom": "", "pos": 4, "p_value": 0.1},  # blank chromosome
            {"chrom": float("nan"), "pos": 5, "p_value": 0.1},
            {"chrom": "7", "pos": 6, "p_value": 0.2, "variant_id": "kept"},
        ]
        records, chromosomes, n_skipped = _prepare(results)

        assert [r["variant_id"] for r in records] == ["kept"]
        assert chromosomes == ["7"]
        assert n_skipped == 7

    def test_prepared_records_carry_neg_log_p_and_beta(self) -> None:
        records, _, _ = _prepare([{"chrom": "1", "pos": 1, "pvalue": 1e-8, "effect_size": 0.25}])

        assert records[0]["neg_log_p"] == pytest.approx(8.0)
        assert records[0]["beta"] == 0.25


class TestInteractiveManhattanEndToEnd:
    def test_result_reports_skipped_results(self, tmp_path: Path) -> None:
        results = [
            {"chrom": "1", "pos": 10, "p_value": 1e-4, "variant_id": "a"},
            {"pos": 20, "p_value": 1e-5, "variant_id": "no_chromosome"},
        ]
        result = interactive_module.interactive_manhattan(results, output_file=tmp_path / "m.html")

        assert result["status"] == "success"
        assert result["n_variants"] == 1
        assert result["n_skipped"] == 1

    def test_nothing_plottable_is_reported_as_skipped(self, tmp_path: Path) -> None:
        out = tmp_path / "m.html"
        result = interactive_module.interactive_manhattan([{"pos": 1, "p_value": 0.1}], output_file=out)

        assert result["status"] == "skipped"
        assert result["n_variants"] == 0
        assert result["n_skipped"] == 1
        assert not out.exists()

    @pytest.mark.skipif(HAS_PLOTLY, reason="asserts on the fallback table; the Plotly path uses the same records")
    def test_fallback_table_shows_the_assigned_chromosomes(self, tmp_path: Path) -> None:
        """Regression: every row used to read chromosome 1."""
        results = [
            {"chrom": "1", "pos": 100, "p_value": 1e-3, "variant_id": "a"},
            {"chrom": "2", "pos": 100, "p_value": 1e-4, "variant_id": "b"},
            {"chrom": "X", "pos": 100, "p_value": 1e-9, "variant_id": "c"},
        ]
        out = tmp_path / "m.html"
        interactive_module.interactive_manhattan(results, output_file=out)

        assert {variant: chrom for variant, chrom, _ in _fallback_rows(out)} == {"a": "1", "b": "2", "c": "X"}

    @pytest.mark.skipif(HAS_PLOTLY, reason="asserts on the fallback table; the Plotly path uses the same records")
    def test_records_shaped_like_generate_all_plots_keep_their_chromosomes(self, tmp_path: Path) -> None:
        """``suite.generate_all_plots`` hands over ``chrom`` / ``pos`` / ``p_value`` records."""
        records = [
            {"chrom": str(1 + i % 3), "pos": 1000 * (i + 1), "p_value": 10 ** -(1 + i % 5), "beta": 0.1, "se": 0.05}
            for i in range(30)
        ]
        out = tmp_path / "m.html"
        result = interactive_module.interactive_manhattan(records, output_file=out)

        assert result["n_variants"] == 30
        rows = _fallback_rows(out)  # sorted by p-value; the table shows at most 50 rows
        assert sorted({chrom for _, chrom, _ in rows}) == ["1", "2", "3"]
        assert {chrom for _, chrom, _ in rows} != {"1"}


class TestInteractiveVolcanoChromosome:
    """The volcano plot had the same gap: its chromosome column ignored ``chrom``."""

    @pytest.mark.skipif(HAS_PLOTLY, reason="the Plotly hover text does not show the chromosome")
    def test_fallback_table_shows_chromosome_from_chrom_key(self, tmp_path: Path) -> None:
        results = [
            {"chrom": "4", "p_value": 1e-9, "beta": 0.9, "variant_id": "a"},
            {"chrom": "9", "p_value": 1e-3, "beta": -0.2, "variant_id": "b"},
            {"p_value": 0.5, "beta": 0.1, "variant_id": "c"},
        ]
        out = tmp_path / "v.html"
        interactive_module.interactive_volcano(results, output_file=out)

        rows = {variant: chrom for variant, chrom, _ in _fallback_rows(out)}
        assert rows == {"a": "4", "b": "9", "c": "N/A"}  # a missing chromosome stays visibly missing
