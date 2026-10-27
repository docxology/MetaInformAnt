"""Regression tests for the GWAS plotting defects fixed in round 3.

Plots that return only a status dictionary leave their figure open, so the tests read the artists back
through ``plt.gcf()``; plots that close their figure (``genomic.regional``) are checked through the
metadata they return.  Every assertion is on data that reached the figure or the result, never on
"something was returned".
"""

from __future__ import annotations

import gzip
import json
import math
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402

from metainformant.gwas.visualization import utils as viz_utils  # noqa: E402
from metainformant.gwas.visualization.genomic import genome, regional  # noqa: E402
from metainformant.gwas.visualization.interactive import suite  # noqa: E402
from metainformant.gwas.visualization.statistical import comparison, effects, statistical  # noqa: E402


def _scatters(ax) -> list[PathCollection]:
    return [collection for collection in ax.collections if isinstance(collection, PathCollection)]


def _offsets(collection) -> np.ndarray:
    return np.asarray(collection.get_offsets(), dtype=float)


# ---------------------------------------------------------------------------
# Shared chromosome helpers
# ---------------------------------------------------------------------------


class TestChromosomeHelpers:
    def test_aliases_match_the_interactive_manhattan(self):
        from metainformant.gwas.visualization.interactive import interactive

        assert viz_utils.CHROMOSOME_KEYS == interactive._CHROMOSOME_KEYS

    def test_get_chromosome_finds_every_alias(self):
        for key in ("chromosome", "chrom", "chr", "CHROM", "CHR"):
            assert viz_utils.get_chromosome({key: "7", "POS": 1}) == "7"

    def test_get_chromosome_prefers_the_earlier_alias_and_skips_none(self):
        assert viz_utils.get_chromosome({"chrom": "3", "CHROM": "9"}) == "3"
        assert viz_utils.get_chromosome({"chromosome": None, "CHR": "5"}) == "5"

    def test_get_chromosome_default(self):
        assert viz_utils.get_chromosome({"POS": 1}) is None
        assert viz_utils.get_chromosome({"POS": 1}, "?") == "?"

    @pytest.mark.parametrize(
        ("value", "label"),
        [(1, "1"), (1.0, "1"), (np.int64(4), "4"), (np.float64(2.0), "2"), ("chr3 ", "chr3"), (2.5, "2.5"), ("X", "X")],
    )
    def test_chromosome_label_makes_int_float_and_text_agree(self, value, label):
        assert viz_utils.chromosome_label(value) == label

    def test_chromosome_label_of_nan_never_matches_a_real_chromosome(self):
        assert viz_utils.chromosome_label(float("nan")) not in {"1", "2", "X"}

    def test_chromosome_sort_key_is_natural(self):
        labels = ["10", "2", "X", "1", "chr3", "MT", "Y"]
        assert sorted(labels, key=viz_utils.chromosome_sort_key) == ["1", "2", "chr3", "10", "MT", "X", "Y"]

    def test_cumulative_offsets_lay_chromosomes_end_to_end(self):
        offsets = viz_utils.cumulative_chromosome_offsets({"1": 300, "2": 200, "3": 50}, ["1", "2", "3"])
        assert offsets == {"1": 0, "2": 301, "3": 502}

    def test_cumulative_offsets_follow_the_requested_order(self):
        offsets = viz_utils.cumulative_chromosome_offsets({"1": 300, "2": 200}, ["2", "1"])
        assert offsets == {"2": 0, "1": 201}


# ---------------------------------------------------------------------------
# genomic.genome.manhattan_plot: cumulative offsets and one scatter per chromosome
# ---------------------------------------------------------------------------

MANHATTAN_ROWS = [
    {"CHROM": "1", "POS": 100, "p_value": 1e-3},
    {"CHROM": "1", "POS": 200, "p_value": 1e-4},
    {"CHROM": "1", "POS": 300, "p_value": 1e-5},
    {"CHROM": "2", "POS": 50, "p_value": 1e-2},
    {"CHROM": "2", "POS": 150, "p_value": 1e-6},
]


def _draw_genome_manhattan(rows, tmp_path: Path):
    result = genome.manhattan_plot(rows, tmp_path / "manhattan.png")
    assert result["status"] == "success", result
    return result, plt.gcf().axes[0]


class TestGenomeManhattan:
    def test_positions_are_cumulative_chromosome_offsets(self, tmp_path):
        _, ax = _draw_genome_manhattan(MANHATTAN_ROWS, tmp_path)

        scatters = _scatters(ax)
        assert len(scatters) == 2
        assert _offsets(scatters[0])[:, 0].tolist() == [100, 200, 300]
        # chromosome 2 starts one base after the largest position on chromosome 1
        assert _offsets(scatters[1])[:, 0].tolist() == [351, 451]

    def test_heights_are_minus_log10_p(self, tmp_path):
        _, ax = _draw_genome_manhattan(MANHATTAN_ROWS, tmp_path)

        scatters = _scatters(ax)
        assert _offsets(scatters[0])[:, 1] == pytest.approx([3.0, 4.0, 5.0])
        assert _offsets(scatters[1])[:, 1] == pytest.approx([2.0, 6.0])

    def test_chromosomes_are_ordered_naturally(self, tmp_path):
        rows = (
            [{"CHROM": "10", "POS": 100 * i, "p_value": 0.1} for i in range(1, 4)]
            + [{"CHROM": "2", "POS": 100 * i, "p_value": 0.1} for i in range(1, 3)]
            + [{"CHROM": "1", "POS": 100, "p_value": 0.1}]
        )

        result, ax = _draw_genome_manhattan(rows, tmp_path)

        assert result["chromosomes"] == ["1", "2", "10"]
        scatters = _scatters(ax)
        assert [len(_offsets(s)) for s in scatters] == [1, 2, 3]
        assert _offsets(scatters[0])[:, 0].tolist() == [100]
        assert _offsets(scatters[1])[:, 0].tolist() == [201, 301]
        assert _offsets(scatters[2])[:, 0].tolist() == [402, 502, 602]

    def test_many_variants_use_one_scatter_per_chromosome(self, tmp_path):
        rng = np.random.default_rng(5)
        n = 3000
        rows = [
            {"CHROM": str(c), "POS": int(pos), "p_value": float(p)}
            for c, pos, p in zip(rng.integers(1, 4, n), rng.integers(1, 10**6, n), rng.uniform(1e-8, 1, n))
        ]

        result, ax = _draw_genome_manhattan(rows, tmp_path)

        scatters = _scatters(ax)
        assert len(scatters) == 3
        assert sum(len(_offsets(s)) for s in scatters) == n
        assert result["n_variants"] == n

    def test_x_positions_agree_with_the_interactive_manhattan(self, tmp_path):
        from metainformant.gwas.visualization.interactive.interactive import _prepare_manhattan_records

        rng = np.random.default_rng(11)
        rows = [
            {"CHROM": str(c), "POS": int(pos), "p_value": float(p)}
            for c, pos, p in zip(
                rng.choice(["1", "2", "10", "X"], 200), rng.integers(1, 5_000_000, 200), rng.uniform(1e-9, 1, 200)
            )
        ]

        result, ax = _draw_genome_manhattan(rows, tmp_path)

        records, chrom_order, _ = _prepare_manhattan_records(
            [{"chrom": r["CHROM"], "pos": r["POS"], "p_value": r["p_value"]} for r in rows]
        )
        assert result["chromosomes"] == chrom_order == ["1", "2", "10", "X"]
        for chrom, scatter in zip(chrom_order, _scatters(ax)):
            expected = sorted(r["cum_pos"] for r in records if r["chromosome"] == chrom)
            assert sorted(_offsets(scatter)[:, 0].tolist()) == expected

    def test_chromosomes_alternate_between_two_colours(self, tmp_path):
        rows = [{"CHROM": str(c), "POS": 10 * i, "p_value": 0.1} for c in (1, 2, 3) for i in (1, 2)]

        _, ax = _draw_genome_manhattan(rows, tmp_path)

        colours = [tuple(np.round(s.get_facecolor()[0][:3], 3)) for s in _scatters(ax)]
        blue, orange = tuple(np.round(to_rgb("#1f77b4"), 3)), tuple(np.round(to_rgb("#ff7f0e"), 3))
        assert colours == [blue, orange, blue]

    def test_chromosome_labels_sit_at_the_mean_position_of_each_chromosome(self, tmp_path):
        _, ax = _draw_genome_manhattan(MANHATTAN_ROWS, tmp_path)

        labels = {text.get_text(): text.get_position()[0] for text in ax.texts}
        assert labels == {"1": pytest.approx(200.0), "2": pytest.approx(401.0)}

    def test_metadata_is_unchanged(self, tmp_path):
        result, _ = _draw_genome_manhattan(MANHATTAN_ROWS, tmp_path)

        assert result["n_variants"] == 5
        assert result["n_chromosomes"] == 2
        assert result["chromosomes"] == ["1", "2"]


# ---------------------------------------------------------------------------
# statistical.effects.effect_direction_plot: cumulative offsets (and the labels that never appeared)
# ---------------------------------------------------------------------------

EFFECT_ROWS = [
    {"CHROM": "1", "POS": 100, "BETA": 0.1},
    {"CHROM": "1", "POS": 200, "BETA": -0.2},
    {"CHROM": "1", "POS": 300, "BETA": 0.3},
    {"CHROM": "2", "POS": 50, "BETA": -0.1},
    {"CHROM": "2", "POS": 150, "BETA": 0.2},
]


class TestEffectDirectionPlot:
    def test_positions_are_cumulative_chromosome_offsets_not_a_running_sum(self):
        result = effects.effect_direction_plot(EFFECT_ROWS)

        assert result["status"] == "success", result
        ax = plt.gcf().axes[0]
        assert _offsets(ax.collections[0])[:, 0].tolist() == [100, 200, 300, 351, 451]
        assert _offsets(ax.collections[0])[:, 1] == pytest.approx([0.1, -0.2, 0.3, -0.1, 0.2])

    def test_chromosome_labels_are_drawn_at_the_chromosome_centres(self):
        effects.effect_direction_plot(EFFECT_ROWS)

        ax = plt.gcf().axes[0]
        labels = {text.get_text(): text.get_position()[0] for text in ax.texts}
        assert labels == {"1": pytest.approx(200.0), "2": pytest.approx(401.0)}

    def test_chromosomes_are_ordered_naturally(self):
        rows = [
            {"CHROM": "10", "POS": 100, "BETA": 0.1},
            {"CHROM": "2", "POS": 100, "BETA": 0.2},
            {"CHROM": "1", "POS": 100, "BETA": 0.3},
        ]

        effects.effect_direction_plot(rows)

        ax = plt.gcf().axes[0]
        # points keep the input order; offsets follow the chromosome order 1 (0), 2 (101), 10 (202)
        assert _offsets(ax.collections[0])[:, 0].tolist() == [302, 201, 100]

    def test_metadata_counts_are_unchanged(self):
        result = effects.effect_direction_plot(EFFECT_ROWS)

        assert (result["n_variants"], result["n_chromosomes"]) == (5, 2)
        assert (result["positive_effects"], result["negative_effects"]) == (3, 2)


# ---------------------------------------------------------------------------
# Chromosome key aliases: chromosome / chrom / chr / CHROM / CHR
# ---------------------------------------------------------------------------

ALIASES = ["chromosome", "chrom", "chr", "CHROM", "CHR"]


class TestChromosomeKeyAliases:
    @pytest.mark.parametrize("key", ALIASES)
    def test_genome_manhattan(self, key):
        result = genome.manhattan_plot(
            [{key: "7", "POS": 100, "p_value": 1e-3}, {key: "8", "POS": 100, "p_value": 1e-4}]
        )

        assert result["chromosomes"] == ["7", "8"]

    @pytest.mark.parametrize("key", ALIASES)
    def test_genome_circular_manhattan(self, key):
        result = genome.circular_manhattan_plot(
            [{key: "7", "POS": 100, "p_value": 1e-3}, {key: "8", "POS": 100, "p_value": 1e-4}]
        )

        assert result["n_chromosomes"] == 2

    @pytest.mark.parametrize("key", ALIASES)
    def test_genome_ideogram(self, key):
        result = genome.chromosome_ideogram([{key: "7", "POS": 100}, {key: "8", "POS": 5000}])

        assert result["chromosomes"] == ["7", "8"]

    @pytest.mark.parametrize("key", ALIASES)
    def test_effect_direction_plot(self, key):
        result = effects.effect_direction_plot(
            [{key: "7", "POS": 100, "BETA": 0.1}, {key: "8", "POS": 100, "BETA": -0.1}]
        )

        assert result["n_chromosomes"] == 2

    @pytest.mark.parametrize("key", ALIASES)
    def test_forest_plot_labels_use_the_real_chromosome(self, key):
        rows = [
            {key: "7", "POS": 100, "p_value": 1e-3, "BETA": 0.1},
            {key: "7", "POS": 200, "p_value": 1e-4, "BETA": -0.2},
        ]

        result = effects.effect_size_forest_plot(rows)

        assert result["status"] == "success", result
        labels = [text.get_text() for text in plt.gcf().axes[0].texts]
        assert "chr7:100" in labels and "chr7:200" in labels
        assert not any("?" in label for label in labels)

    @pytest.mark.parametrize("key", ALIASES)
    def test_allelic_series_gene_region_filter(self, key):
        rows = [
            {key: "7", "POS": 100, "BETA": 0.1, "p_value": 1e-3},
            {key: "7", "POS": 200, "BETA": -0.1, "p_value": 1e-2},
            {key: "8", "POS": 150, "BETA": 0.3, "p_value": 1e-4},
        ]

        result = effects.allelic_series_plot(rows, gene_region=("7", 0, 1000))

        assert result["status"] == "success", result
        assert result["n_variants"] == 2


# ---------------------------------------------------------------------------
# interactive.suite.generate_all_plots
# ---------------------------------------------------------------------------


def _write_association_tsv(path: Path, n: int = 30, seed: int = 1) -> Path:
    rng = np.random.default_rng(seed)
    lines = ["CHROM\tPOS\tp_value\tBETA\tSE"]
    for i in range(n):
        lines.append(f"{1 + i % 2}\t{1000 * (i + 1)}\t{rng.uniform(1e-6, 1):.6g}\t{rng.normal():.4f}\t0.1")
    path.write_text("\n".join(lines) + "\n")
    return path


def _population_arrays(n_samples: int = 12, seed: int = 4) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    pcs = rng.normal(size=(n_samples, 3))
    kinship = np.corrcoef(rng.normal(size=(n_samples, 40)))
    return pcs, kinship


class TestPopulationInputFiles:
    """PCA / kinship inputs in the formats the package writes (and TSV) used to be skipped with a warning."""

    def test_pca_tsv_without_header(self, tmp_path):
        pcs, _ = _population_arrays()
        path = tmp_path / "pca.tsv"
        np.savetxt(path, pcs, delimiter="\t")

        loaded = suite._load_pca_file(path)

        np.testing.assert_allclose(loaded["pcs"], pcs)

    def test_pca_tsv_with_header_and_sample_labels(self, tmp_path):
        pcs, _ = _population_arrays(5)
        path = tmp_path / "pca_labelled.tsv"
        rows = ["sample\tPC1\tPC2\tPC3"] + [
            f"s{i}\t" + "\t".join(repr(float(v)) for v in row) for i, row in enumerate(pcs)
        ]
        path.write_text("\n".join(rows) + "\n")

        loaded = suite._load_pca_file(path)

        np.testing.assert_allclose(loaded["pcs"], pcs)

    def test_pca_csv(self, tmp_path):
        pcs, _ = _population_arrays(6)
        path = tmp_path / "pca.csv"
        np.savetxt(path, pcs, delimiter=",")

        np.testing.assert_allclose(suite._load_pca_file(path)["pcs"], pcs)

    def test_pca_npy_is_samples_by_components(self, tmp_path):
        pcs, _ = _population_arrays()
        path = tmp_path / "pca.npy"
        np.save(path, pcs)

        np.testing.assert_allclose(suite._load_pca_file(path)["pcs"], pcs)

    def test_pca_npz_prefers_the_pcs_array(self, tmp_path):
        pcs, _ = _population_arrays()
        path = tmp_path / "pca.npz"
        np.savez(path, other=np.zeros(3), pcs=pcs, explained_variance_ratio=np.array([0.5, 0.3, 0.1]))

        loaded = suite._load_pca_file(path)

        np.testing.assert_allclose(loaded["pcs"], pcs)
        assert list(loaded["explained_variance_ratio"]) == pytest.approx([0.5, 0.3, 0.1])

    def test_pca_json_result_dict(self, tmp_path):
        pcs, _ = _population_arrays()
        path = tmp_path / "pca_results.json"
        path.write_text(json.dumps({"pcs": pcs.tolist(), "explained_variance_ratio": [0.4, 0.2, 0.1]}))

        loaded = suite._load_pca_file(path)

        np.testing.assert_allclose(loaded["pcs"], pcs)
        assert list(loaded["explained_variance_ratio"]) == pytest.approx([0.4, 0.2, 0.1])

    def test_pca_json_workflow_format_is_transposed_to_samples_by_components(self, tmp_path):
        pcs, _ = _population_arrays()
        path = tmp_path / "pca_workflow.json"
        payload = {"components": pcs.T.tolist(), "variance": [], "explained_variance": [0.4, 0.2, 0.1], "loadings": []}
        path.write_text(json.dumps(payload))

        loaded = suite._load_pca_file(path)

        np.testing.assert_allclose(loaded["pcs"], pcs)
        assert list(loaded["explained_variance_ratio"]) == pytest.approx([0.4, 0.2, 0.1])

    @pytest.mark.parametrize("suffix", [".tsv", ".csv", ".npy", ".json"])
    def test_kinship_formats(self, tmp_path, suffix):
        _, kinship = _population_arrays()
        path = tmp_path / f"kinship{suffix}"
        if suffix == ".npy":
            np.save(path, kinship)
        elif suffix == ".json":
            path.write_text(json.dumps(kinship.tolist()))
        else:
            np.savetxt(path, kinship, delimiter="\t" if suffix == ".tsv" else ",")

        np.testing.assert_allclose(suite._load_kinship_file(path), kinship)

    def test_kinship_json_result_dict(self, tmp_path):
        _, kinship = _population_arrays()
        for key in ("kinship_matrix", "matrix"):
            path = tmp_path / f"kinship_{key}.json"
            path.write_text(json.dumps({key: kinship.tolist(), "status": "success"}))

            np.testing.assert_allclose(suite._load_kinship_file(path), kinship)

    def test_non_square_kinship_is_rejected(self, tmp_path):
        path = tmp_path / "kinship.tsv"
        np.savetxt(path, np.zeros((3, 4)), delimiter="\t")

        with pytest.raises(ValueError, match="square"):
            suite._load_kinship_file(path)

    def test_unsupported_format_is_rejected(self, tmp_path):
        path = tmp_path / "pca.xlsx"
        path.write_bytes(b"not a table")

        with pytest.raises(ValueError, match="Unsupported"):
            suite._load_pca_file(path)

    def test_windows_line_endings_and_a_blank_header_corner(self, tmp_path):
        crlf = tmp_path / "crlf.tsv"
        crlf.write_bytes(b"1.0\t2.0\r\n3.0\t4.0\r\n")
        pandas_style = tmp_path / "pandas_style.tsv"
        pandas_style.write_bytes(b"\tPC1\tPC2\ns1\t0.5\t0.25\ns2\t-1\t1e-3\n")

        np.testing.assert_allclose(suite._load_pca_file(crlf)["pcs"], [[1.0, 2.0], [3.0, 4.0]])
        np.testing.assert_allclose(suite._load_pca_file(pandas_style)["pcs"], [[0.5, 0.25], [-1.0, 0.001]])

    @pytest.mark.parametrize(
        ("content", "reason"),
        [
            ("PC1\tPC2\n", "no data rows"),
            ("", "no data rows"),
            ("1\t2\n3\n", "rectangular"),
            ("1\t2\n3\tx\n", "non-numeric"),
        ],
    )
    def test_malformed_tables_are_rejected_with_the_reason(self, tmp_path, content, reason):
        path = tmp_path / "bad.tsv"
        path.write_text(content)

        with pytest.raises(ValueError, match=reason):
            suite._load_pca_file(path)

    def test_unreadable_pca_file_does_not_stop_the_other_plots(self, tmp_path):
        bad = tmp_path / "pca.tsv"
        bad.write_text("1\t2\n3\n")  # not rectangular
        association = _write_association_tsv(tmp_path / "assoc.tsv")

        result = suite.generate_all_plots(association, tmp_path / "plots", pca_file=bad)

        assert result["status"] == "completed"
        assert "pca" not in result["plots"]
        assert "manhattan" in result["plots"] and "summary_panel" in result["plots"]

    @pytest.mark.parametrize("fmt", ["npy", "tsv", "json"])
    def test_generate_all_plots_plots_pca_and_kinship_files(self, tmp_path, fmt):
        pcs, kinship = _population_arrays()
        association = _write_association_tsv(tmp_path / "assoc.tsv")
        if fmt == "npy":
            np.save(tmp_path / "pca.npy", pcs)
            np.save(tmp_path / "kinship.npy", kinship)
        elif fmt == "tsv":
            np.savetxt(tmp_path / "pca.tsv", pcs, delimiter="\t")
            np.savetxt(tmp_path / "kinship.tsv", kinship, delimiter="\t")
        else:
            (tmp_path / "pca.json").write_text(json.dumps({"pcs": pcs.tolist()}))
            (tmp_path / "kinship.json").write_text(json.dumps(kinship.tolist()))

        result = suite.generate_all_plots(
            association,
            tmp_path / "plots",
            pca_file=tmp_path / f"pca.{fmt}",
            kinship_file=tmp_path / f"kinship.{fmt}",
        )

        assert result["status"] == "completed"
        for name in ("pca", "kinship"):
            assert name in result["plots"], result["plots"]
            assert Path(result["plots"][name]).stat().st_size > 0
        # the files feed every plot that needs PCA / kinship data, not only the two standalone ones
        # (the interactive PCA rejects a NumPy array, so the loaded components must stay nested lists)
        for name in ("pca_interactive", "population_structure_panel", "kinship_dendrogram", "kinship_clustermap"):
            assert name in result["plots"], result["plots"]


class TestSummaryPanelReceivesPopulationData:
    @staticmethod
    def _summary_png(directory: Path, **kwargs) -> bytes:
        result = suite.generate_all_plots(
            _write_association_tsv(directory / "assoc.tsv"), directory / "plots", **kwargs
        )
        assert result["status"] == "completed"
        return Path(result["plots"]["summary_panel"]).read_bytes()

    def test_in_memory_pca_and_kinship_change_the_summary_panel(self, tmp_path):
        pcs, kinship = _population_arrays()
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()
        (tmp_path / "c").mkdir()

        without_1 = self._summary_png(tmp_path / "a")
        without_2 = self._summary_png(tmp_path / "b")
        with_data = self._summary_png(tmp_path / "c", pca_data={"pcs": pcs.tolist()}, kinship_matrix=kinship)

        assert without_1 == without_2  # control: rendering is deterministic
        assert with_data != without_1  # the PCA and kinship panels are no longer "No data available"

    def test_files_change_the_summary_panel_too(self, tmp_path):
        pcs, kinship = _population_arrays()
        np.savetxt(tmp_path / "pca.tsv", pcs, delimiter="\t")
        np.savetxt(tmp_path / "kinship.tsv", kinship, delimiter="\t")
        (tmp_path / "a").mkdir()
        (tmp_path / "b").mkdir()

        without = self._summary_png(tmp_path / "a")
        with_files = self._summary_png(
            tmp_path / "b", pca_file=tmp_path / "pca.tsv", kinship_file=tmp_path / "kinship.tsv"
        )

        assert with_files != without

    def test_in_memory_data_wins_over_files(self, tmp_path):
        pcs, kinship = _population_arrays()
        other_pcs, other_kinship = _population_arrays(seed=99)
        np.savetxt(tmp_path / "pca.tsv", other_pcs, delimiter="\t")
        np.savetxt(tmp_path / "kinship.tsv", other_kinship, delimiter="\t")
        for name in ("a", "b"):
            (tmp_path / name).mkdir()

        from_memory = self._summary_png(tmp_path / "a", pca_data={"pcs": pcs.tolist()}, kinship_matrix=kinship)
        both = self._summary_png(
            tmp_path / "b",
            pca_data={"pcs": pcs.tolist()},
            kinship_matrix=kinship,
            pca_file=tmp_path / "pca.tsv",
            kinship_file=tmp_path / "kinship.tsv",
        )

        assert both == from_memory


VCF_TEXT = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
    "1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n"
    "1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1\n"
    "1\t300\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/0\n"
    "2\t100\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/1\n"
    "2\t200\t.\tT\tG\t.\t.\t.\tGT\t0/1\t0/1\n"
    "2\t300\t.\tAT\tA\t.\t.\t.\tGT\t0/1\t0/0\n"
    "2\t400\t.\tA\tG,T\t.\t.\t.\tGT\t0/1\t1/2\n"
)


class TestVcfFile:
    """``vcf_file`` was documented ("Optional VCF file for variant plots") but never read."""

    def test_vcf_produces_density_and_transition_transversion_plots(self, tmp_path):
        vcf = tmp_path / "calls.vcf"
        vcf.write_text(VCF_TEXT)

        result = suite.generate_all_plots(
            _write_association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots", vcf_file=vcf
        )

        assert result["status"] == "completed"
        for name in ("variant_density", "transition_transversion"):
            assert Path(result["plots"][name]).stat().st_size > 0
        # 7 sites; 3 transitions (A>G, C>T, G>A) and 2 transversions (A>C, T>G); the indel and the
        # multi-allelic site are not counted as either
        assert result["statistics"]["vcf_variants"] == 7
        assert result["statistics"]["ts_tv_ratio"] == pytest.approx(1.5)

    def test_gzipped_vcf(self, tmp_path):
        vcf = tmp_path / "calls.vcf.gz"
        with gzip.open(vcf, "wt") as handle:
            handle.write(VCF_TEXT)

        result = suite.generate_all_plots(
            _write_association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots", vcf_file=vcf
        )

        assert result["statistics"]["vcf_variants"] == 7
        assert "variant_density" in result["plots"]

    def test_missing_vcf_is_reported_not_fatal(self, tmp_path):
        result = suite.generate_all_plots(
            _write_association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots", vcf_file=tmp_path / "absent.vcf"
        )

        assert result["status"] == "completed"
        assert "variant_density" not in result["plots"]
        assert "vcf_variants" not in result["statistics"]

    def test_no_vcf_plots_without_a_vcf(self, tmp_path):
        result = suite.generate_all_plots(_write_association_tsv(tmp_path / "assoc.tsv"), tmp_path / "plots")

        assert "variant_density" not in result["plots"]
        assert "transition_transversion" not in result["plots"]


# ---------------------------------------------------------------------------
# genomic.regional
# ---------------------------------------------------------------------------


def _regional_frame(chr_values) -> pd.DataFrame:
    return pd.DataFrame({"CHR": chr_values, "BP": [1000, 1500, 2000, 1000], "P": [1e-8, 1e-6, 0.01, 0.5]})


CHR_COMBINATIONS = [
    pytest.param([1, 1, 1, 2], "1", id="int-column/str-query"),
    pytest.param(["1", "1", "1", "2"], 1, id="str-column/int-query"),
    pytest.param([1.0, 1.0, 1.0, 2.0], "1", id="float-column/str-query"),
    pytest.param(["1", "1", "1", "2"], "1", id="str-column/str-query"),
]


class TestRegionalChromosomeMatching:
    @pytest.mark.parametrize(("chr_values", "chrom"), CHR_COMBINATIONS)
    def test_regional_plot(self, tmp_path, chr_values, chrom):
        out = tmp_path / "regional.png"

        result = regional.regional_plot(_regional_frame(chr_values), out, chrom=chrom, start=500, end=2500)

        assert result["status"] == "success"
        assert result["n_variants"] == 3
        assert out.stat().st_size > 0

    @pytest.mark.parametrize(("chr_values", "chrom"), CHR_COMBINATIONS)
    def test_gene_annotation_plot(self, tmp_path, chr_values, chrom):
        out = tmp_path / "genes.png"

        result = regional.gene_annotation_plot(_regional_frame(chr_values), out, chrom=chrom, start=500, end=2500)

        assert result["status"] == "success"
        assert result["n_variants"] == 3


class TestRegionalOpenEndedRegions:
    """A missing ``start``/``end`` raised ``TypeError`` while formatting the title and the recombination slice."""

    FRAME = _regional_frame([1, 1, 1, 2])

    def test_no_start_and_no_end(self, tmp_path):
        result = regional.regional_plot(self.FRAME, tmp_path / "r.png", chrom=1)

        assert result["status"] == "success"
        assert result["n_variants"] == 3

    def test_only_start(self):
        result = regional.regional_plot(self.FRAME, chrom=1, start=1200)

        assert result["n_variants"] == 2

    def test_only_end(self):
        result = regional.regional_plot(self.FRAME, chrom=1, end=1500)

        assert result["n_variants"] == 2

    def test_no_chromosome_and_no_region_plots_everything(self):
        result = regional.regional_plot(self.FRAME)

        assert result["n_variants"] == 4

    def test_recombination_data_without_bounds(self):
        recombination = pd.DataFrame({"CHR": [1, 1], "BP": [1200, 5000], "RATE": [2.0, 3.0]})

        result = regional.regional_plot(self.FRAME, chrom=1, recombination_rate=recombination)

        assert result["status"] == "success"
        assert result["n_variants"] == 3


class TestGeneTrackHonesty:
    @staticmethod
    def _track(gene_data):
        _, ax = plt.subplots()
        regional.plot_gene_track(ax, "1", 0, 100, gene_data)
        return ax

    def test_missing_gene_annotations_are_stated(self):
        ax = self._track({})

        assert len(ax.patches) == 0
        assert [text.get_text() for text in ax.texts] == ["Gene annotations not provided"]

    def test_empty_gene_list_is_stated(self):
        ax = self._track({"genes": []})

        assert len(ax.patches) == 0
        assert [text.get_text() for text in ax.texts] == ["Gene annotations not provided"]

    def test_genes_outside_the_region_are_stated(self):
        ax = self._track({"genes": [{"name": "far", "start": 500, "end": 900}]})

        assert len(ax.patches) == 0
        assert [text.get_text() for text in ax.texts] == ["No annotated genes overlap this region"]

    def test_gene_without_coordinates_is_not_drawn_as_a_full_width_box(self):
        ax = self._track({"genes": [{"name": "no_coordinates"}]})

        assert len(ax.patches) == 0
        (message,) = [text.get_text() for text in ax.texts]
        assert message == "No annotated genes overlap this region (1 gene without coordinates not drawn)"

    def test_gene_with_coordinates_is_drawn_where_it_belongs(self):
        ax = self._track({"genes": [{"name": "G1", "start": 20, "end": 60}, {"name": "no_end", "start": 5}]})

        (box,) = ax.patches
        assert (box.get_x(), box.get_width()) == (20, 40)
        assert [text.get_text() for text in ax.texts] == ["G1"]  # no "not provided" note beside a drawn gene

    def test_genes_are_clipped_to_the_region(self):
        ax = self._track({"genes": [{"name": "edge", "start": -50, "end": 30}]})

        (box,) = ax.patches
        assert (box.get_x(), box.get_width()) == (0, 30)


class TestRecombinationHonesty:
    RESULTS = [
        {"CHROM": "chr1", "POS": 1000, "p_value": 1e-8},
        {"CHROM": "chr1", "POS": 2000, "p_value": 1e-6},
    ]

    def test_without_a_recombination_map_nothing_is_plotted_and_the_result_says_so(self, tmp_path):
        out = tmp_path / "recombination.png"

        result = regional.recombination_rate_plot(self.RESULTS, out, chrom="chr1", start=500, end=2500)

        assert result["status"] == "skipped"
        assert "recombination" in result["reason"].lower()
        assert result["output_path"] is None
        assert not out.exists()

    def test_one_of_rates_or_positions_is_not_enough(self, tmp_path):
        result = regional.recombination_rate_plot(
            self.RESULTS, tmp_path / "r.png", recombination_rates={"chr1": [1.0, 2.0]}
        )

        assert result["status"] == "skipped"
        assert "recombination_rates" in result["reason"] and "positions" in result["reason"]
        assert result["output_path"] is None

    def test_explicit_rates_and_positions_are_plotted(self, tmp_path):
        out = tmp_path / "recombination.png"

        result = regional.recombination_rate_plot(
            self.RESULTS,
            out,
            recombination_rates={"chr1": [1.0, 2.0, 1.5]},
            positions={"chr1": [1000, 2000, 3000]},
        )

        assert result["status"] == "success"
        assert result["n_chromosomes"] == 1
        assert out.stat().st_size > 0


# ---------------------------------------------------------------------------
# statistical: Q-Q plots
# ---------------------------------------------------------------------------


class TestQQPlotExpectedValues:
    P_VALUES = [1e-8, 0.01, 0.2, 0.5]

    def test_expected_axis_uses_uniform_quantiles(self, tmp_path):
        result = statistical.qq_plot(self.P_VALUES, tmp_path / "qq.png")

        assert result["status"] == "success", result
        points = _offsets(_scatters(plt.gcf().axes[0])[0])
        n = len(self.P_VALUES)
        assert points[:, 0] == pytest.approx([-math.log10((i - 0.5) / n) for i in range(1, n + 1)])
        assert points[:, 1] == pytest.approx([-math.log10(p) for p in sorted(self.P_VALUES)])

    def test_median_expected_in_the_result_follows_the_same_quantiles(self, tmp_path):
        result = statistical.qq_plot(self.P_VALUES, tmp_path / "qq.png")

        n = len(self.P_VALUES)
        expected = [-math.log10((i - 0.5) / n) for i in range(1, n + 1)]
        assert result["median_expected"] == pytest.approx(float(np.median(expected)))

    def test_expected_values_do_not_depend_on_input_order(self, tmp_path):
        statistical.qq_plot(list(reversed(self.P_VALUES)), tmp_path / "qq.png")

        points = _offsets(_scatters(plt.gcf().axes[0])[0])
        assert points[:, 0] == pytest.approx([-math.log10((i - 0.5) / 4) for i in range(1, 5)])

    def test_p_value_of_exactly_zero_is_clamped_not_fatal(self, tmp_path):
        out = tmp_path / "qq_zero.png"

        result = statistical.qq_plot([0.0, 1e-8, 0.01, 0.5], out)

        assert result["status"] == "success", result
        assert out.stat().st_size > 0
        points = _offsets(_scatters(plt.gcf().axes[0])[0])
        assert points[0, 1] == pytest.approx(-math.log10(sys.float_info.min))
        assert points[1:, 1] == pytest.approx([-math.log10(p) for p in (1e-8, 0.01, 0.5)])
        assert math.isfinite(result["lambda_gc"])

    def test_negative_p_values_are_rejected_rather_than_clamped(self, tmp_path):
        result = statistical.qq_plot([-0.1, 0.5, 0.2], tmp_path / "qq_negative.png")

        assert result["status"] == "failed"
        assert "non-negative" in result["error"]


class TestStratifiedQQ:
    def test_expected_axis_uses_uniform_quantiles_and_zero_is_clamped(self, tmp_path):
        p_values = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9]
        strata = ["a", "a", "a", "b", "b", "b"]

        result = statistical.qq_plot_stratified(p_values, tmp_path / "strat.png", strata=strata)

        assert result["status"] == "success", result
        fig = plt.gcf()
        points_a = _offsets(fig.axes[0].collections[0])
        points_b = _offsets(fig.axes[1].collections[0])
        expected_x = [-math.log10((i - 0.5) / 3) for i in (1, 2, 3)]
        assert points_a[:, 0] == pytest.approx(expected_x)
        assert points_b[:, 0] == pytest.approx(expected_x)
        assert points_a[:, 1] == pytest.approx([-math.log10(sys.float_info.min), -math.log10(0.2), -math.log10(0.4)])
        assert points_b[:, 1] == pytest.approx([-math.log10(0.6), -math.log10(0.8), -math.log10(0.9)])

    def test_stratum_summaries_are_reported(self, tmp_path):
        result = statistical.qq_plot_stratified(
            [0.1, 0.2, 0.3, 0.4], tmp_path / "strat.png", strata=["a", "a", "b", "b"]
        )

        assert sorted(result["strata_results"]) == ["a", "b"]
        assert result["strata_results"]["a"]["n_p_values"] == 2


class TestLambdaGCQQ:
    @staticmethod
    def _rows():
        rows = []
        for chrom, upper in (("1", 0.5), ("2", 0.9), ("3", 0.99)):
            rows += [{"p_value": float(p), "CHROM": chrom} for p in np.linspace(0.01, upper, 20)]
        return rows

    def test_expected_lambdas_are_normal_quantiles_and_reproducible(self):
        first = statistical.lambda_gc_plot(self._rows())
        points_1 = _offsets(plt.gcf().axes[1].collections[0]).copy()
        plt.close("all")
        second = statistical.lambda_gc_plot(self._rows())
        points_2 = _offsets(plt.gcf().axes[1].collections[0]).copy()

        assert first["status"] == second["status"] == "success"
        np.testing.assert_array_equal(points_1, points_2)
        reference = NormalDist(1.0, 0.1)
        assert points_1[:, 0] == pytest.approx([reference.inv_cdf((i - 0.5) / 3) for i in (1, 2, 3)])
        assert points_1[:, 1] == pytest.approx(sorted(first["lambda_gc_by_chrom"].values()))


# ---------------------------------------------------------------------------
# statistical.comparison: method power
# ---------------------------------------------------------------------------


def _p_frame(n_significant: int, n_total: int, index=None) -> pd.DataFrame:
    p_values = [1e-9] * n_significant + list(np.linspace(0.05, 0.95, n_total - n_significant))
    return pd.DataFrame({"P": p_values}, index=index)


class TestMethodPower:
    def test_every_method_with_a_p_column_is_reported(self):
        out = comparison.calculate_method_power(
            {
                "Linear": {"results": _p_frame(2, 5)},
                "NoP": {"results": pd.DataFrame({"BETA": [1.0]})},
                "NotAFrame": {"results": {"P": [0.1]}},
            }
        )

        assert list(out) == ["Linear"]
        assert out["Linear"]["significant_snps"] == 2
        assert out["Linear"]["total_snps"] == 5
        assert out["Linear"]["power"] == pytest.approx(0.4)

    def test_without_ground_truth_the_value_is_labelled_as_the_significant_fraction(self):
        out = comparison.calculate_method_power({"Linear": {"results": _p_frame(2, 5)}})

        assert out["Linear"]["basis"] == "significant_fraction"
        assert "causal_snps" not in out["Linear"]

    def test_threshold_is_configurable(self):
        frame = pd.DataFrame({"P": [1e-9, 1e-5, 1e-3, 0.5]})

        out = comparison.calculate_method_power({"M": {"results": frame}}, significance_threshold=1e-4)

        assert out["M"]["significant_snps"] == 2
        assert out["M"]["power"] == pytest.approx(0.5)

    def test_missing_p_values_are_not_counted_as_tests(self):
        frame = pd.DataFrame({"P": [1e-9, np.nan, 0.5, np.nan]})

        out = comparison.calculate_method_power({"M": {"results": frame}})

        assert out["M"]["total_snps"] == 2
        assert out["M"]["power"] == pytest.approx(0.5)

    def test_ground_truth_gives_real_power(self):
        frame = pd.DataFrame({"P": [1e-9, 1e-3, 0.5, 0.7, 1e-10]}, index=["s1", "s2", "s3", "s4", "s5"])

        out = comparison.calculate_method_power({"M": {"results": frame, "causal_snps": ["s1", "s2", "s5", "absent"]}})

        assert out["M"]["basis"] == "causal_snps"
        assert out["M"]["causal_snps"] == 3  # "absent" was never tested
        assert out["M"]["causal_detected"] == 2
        assert out["M"]["power"] == pytest.approx(2 / 3)
        assert out["M"]["significant_snps"] == 2  # unchanged: counts over all tested SNPs

    def test_empty_ground_truth_falls_back_to_the_significant_fraction(self):
        frame = pd.DataFrame({"P": [1e-9, 0.5]}, index=["s1", "s2"])

        out = comparison.calculate_method_power({"M": {"results": frame, "causal_snps": ["not_tested"]}})

        assert out["M"]["basis"] == "significant_fraction"
        assert out["M"]["power"] == pytest.approx(0.5)


class TestCompareAnalysisMethodsPowerPanel:
    @staticmethod
    def _methods(**extra):
        return {
            "Linear": {"results": _p_frame(4, 10), "runtime_minutes": 12.34, **extra.get("Linear", {})},
            "Mixed": {"results": _p_frame(2, 8), "runtime_minutes": 150.0, **extra.get("Mixed", {})},
        }

    def test_power_panel_draws_one_labelled_bar_per_method(self):
        pytest.importorskip("seaborn")

        fig = comparison.compare_analysis_methods(self._methods())

        panel = fig.axes[2]
        assert [bar.get_height() for bar in panel.patches] == pytest.approx([0.4, 0.25])
        assert [text.get_text() for text in panel.texts] == ["0.400", "0.250"]
        assert [label.get_text() for label in panel.get_xticklabels()] == ["Linear", "Mixed"]

    def test_panel_is_not_called_power_when_it_is_only_the_significant_fraction(self):
        pytest.importorskip("seaborn")

        panel = comparison.compare_analysis_methods(self._methods()).axes[2]

        assert "power" not in panel.get_title().lower()
        assert "power" not in panel.get_ylabel().lower()
        assert "significant" in panel.get_ylabel().lower()

    def test_panel_is_called_power_when_every_method_has_ground_truth(self):
        pytest.importorskip("seaborn")
        methods = {
            "Linear": {"results": _p_frame(4, 10, index=list("abcdefghij")), "causal_snps": list("abcd")},
            "Mixed": {"results": _p_frame(2, 8, index=list("abcdefgh")), "causal_snps": list("abcd")},
        }

        panel = comparison.compare_analysis_methods(methods).axes[2]

        assert panel.get_title() == "Method Power Comparison"
        assert panel.get_ylabel() == "Statistical Power"
        assert [bar.get_height() for bar in panel.patches] == pytest.approx([1.0, 0.5])

    def test_runtime_labels_show_the_values(self):
        pytest.importorskip("seaborn")

        panel = comparison.compare_analysis_methods(self._methods()).axes[3]

        assert [text.get_text() for text in panel.texts] == ["12.3", "150.0"]
