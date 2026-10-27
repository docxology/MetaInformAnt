"""Regression tests: ``generate_all_plots`` reads its results under every spelling and never draws a failed test.

``metainformant.gwas.visualization.general`` re-exports the plotting functions of ``_general_impl``, so the
``generate_all_plots`` imported from ``general`` (and the one ``execute_gwas_workflow`` calls) is the one in
``_general_impl``. Defects found there:

* a pandas DataFrame was silently ignored (no plots, no message), as was any input that was not a list or a path;
* a path that does not exist, or has an unsupported extension, was silently ignored too;
* chromosome, position and p-value were read under different spellings by different plots: ``manhattan_plot`` knew
  ``chrom``/``pos``/``p_value``/``pval`` only, so a record spelled ``CHR``/``BP``/``P`` or ``pvalue`` was drawn flat at
  -log10 p = 0 while the Q-Q plot read the same p-values correctly;
* a result without a p-value (a failed test: NaN, None, or "NA" in a table) was drawn at -log10 p = 50 and counted as
  genome-wide significant in the title; "NA" strings also made the Q-Q plot raise;
* when there was nothing to draw nothing said so.

Oracles: rendering is deterministic, so two inputs that mean the same thing must give byte-identical PNGs, and the
points of a figure are read straight from its axes. Hand-counted data throughout.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from metainformant.gwas.analysis.association import run_linear_model_gwas
from metainformant.gwas.analysis.summary_stats import write_summary_statistics
from metainformant.gwas.visualization import _general_impl, general
from metainformant.gwas.visualization.general import generate_all_plots, manhattan_plot, regional_plot
from metainformant.gwas.workflow.workflow_execution import execute_gwas_workflow

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
NAN = float("nan")


def rows(n: int = 24) -> List[Dict[str, Any]]:
    """Results on two chromosomes with p-values from 1e-1 down to 1e-6."""
    return [
        {"chrom": "1" if i < n // 2 else "2", "pos": 1_000 * (i + 1), "p_value": 10.0 ** -(1 + i % 6)} for i in range(n)
    ]


def with_missing(records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """The same results plus four without a usable p-value, inside the extent of chromosomes 1 and 2."""
    extra = [
        {"chrom": "1", "pos": 1_500, "p_value": NAN},
        {"chrom": "1", "pos": 2_500, "p_value": None},
        {"chrom": "2", "pos": 13_500, "p_value": "NA"},
        {"chrom": "2", "pos": 14_500},  # no p-value key at all
    ]
    return records + extra


@pytest.fixture(scope="module", autouse=True)
def settled_plot_style(tmp_path_factory: pytest.TempPathFactory) -> None:
    """The first plot of a process applies the visualization style to matplotlib globally.

    Draw one up front, so that every comparison below is between two renderings made in the same state.
    """
    generate_all_plots(rows(), tmp_path_factory.mktemp("warm_up"))


def render(data: Any, directory: Path) -> Dict[str, bytes]:
    """Plots written by ``generate_all_plots`` for ``data``, as bytes keyed by plot name."""
    written = generate_all_plots(data, directory)
    return {name: Path(path).read_bytes() for name, path in written.items()}


def assert_png(content: bytes) -> None:
    assert content.startswith(PNG_SIGNATURE) and len(content) > 5_000


def points(fig: Any) -> np.ndarray:
    """The (x, y) of every point drawn by the first scatter of the figure's first axes."""
    return np.asarray(fig.axes[0].collections[0].get_offsets())


# ---------------------------------------------------------------------------
# Which function is called
# ---------------------------------------------------------------------------


def test_the_public_generate_all_plots_is_the_implementation_in_general_impl() -> None:
    """``general.py`` re-exports ``_general_impl``; this is the function the workflow calls."""
    assert general.generate_all_plots is _general_impl.generate_all_plots
    assert general.manhattan_plot is _general_impl.manhattan_plot


# ---------------------------------------------------------------------------
# Input forms
# ---------------------------------------------------------------------------


def test_a_list_of_records_gives_a_manhattan_and_a_q_q_plot(tmp_path: Path) -> None:
    written = render(rows(), tmp_path)

    assert sorted(written) == ["manhattan", "qq"]
    for content in written.values():
        assert_png(content)
    assert sorted(p.name for p in tmp_path.glob("*.png")) == ["manhattan_plot.png", "qq_plot.png"]


def test_a_dataframe_is_read_like_the_same_list_of_records(tmp_path: Path) -> None:
    from_list = render(rows(), tmp_path / "list")

    from_frame = render(pd.DataFrame(rows()), tmp_path / "frame")

    assert sorted(from_frame) == ["manhattan", "qq"]
    assert from_frame == from_list


def test_plink_style_columns_are_read_like_the_pipelines_own_keys(tmp_path: Path) -> None:
    records = rows()
    plink = pd.DataFrame(
        {
            "CHR": [r["chrom"] for r in records],
            "BP": [r["pos"] for r in records],
            "SNP": [f"rs{i}" for i in range(len(records))],
            "P": [r["p_value"] for r in records],
        }
    )
    reference = render(records, tmp_path / "reference")

    assert render(plink, tmp_path / "frame") == reference
    assert render(plink.to_dict("records"), tmp_path / "records") == reference  # the same spelling as dicts


def test_every_spelling_of_the_chromosome_position_and_p_value_is_read_alike(tmp_path: Path) -> None:
    records = rows()
    reference = render(records, tmp_path / "reference")
    spellings = [
        ("chromosome", "position", "pvalue"),
        ("chr", "POS", "pval"),
        ("CHROM", "BP", "P"),
    ]

    for chrom_key, pos_key, p_key in spellings:
        renamed = [{chrom_key: r["chrom"], pos_key: r["pos"], p_key: r["p_value"]} for r in records]
        rendered = render(renamed, tmp_path / "-".join((chrom_key, pos_key, p_key)))
        assert rendered == reference, (chrom_key, pos_key, p_key)


def test_a_results_table_with_na_is_plotted_without_its_failed_rows(tmp_path: Path) -> None:
    """The pipeline's own ``summary_statistics.tsv`` (CHR/POS/P columns, NA for a failed test) is plottable."""
    genotypes = [[0, 1, 2, 0, 1, 2, 1, 0], [1] * 8, [2, 1, 0, 0, 1, 2, 2, 1], [0] * 8, [1, 2, 0, 1, 1, 0, 2, 1]]
    phenotype = [1.2, 2.3, 2.9, 1.8, 3.6, 2.1, 3.9, 1.1]
    info = [{"id": f"rs{i}", "chrom": "chr1" if i < 3 else "chr2", "pos": 100 * (i + 1)} for i in range(5)]
    results = run_linear_model_gwas(genotypes, phenotype, variant_info=info)
    assert [r["status"] for r in results] == ["success", "error", "success", "error", "success"]
    table = write_summary_statistics(results, info, tmp_path / "summary_statistics.tsv")

    # what the table holds (the writer keeps seven significant digits of each p-value)
    kept = [
        {"chrom": i["chrom"], "pos": i["pos"], "p_value": float(f"{r['p_value']:.6e}")}
        for r, i in zip(results, info)
        if r["status"] == "success"
    ]

    assert render(table, tmp_path / "from_table") == render(kept, tmp_path / "from_records")
    assert sorted(render(table, tmp_path / "again")) == ["manhattan", "qq"]


def test_a_json_file_of_records_is_read_too(tmp_path: Path) -> None:
    import json

    path = tmp_path / "results.json"
    path.write_text(json.dumps({"results": rows()}))

    assert render(path, tmp_path / "from_json") == render(rows(), tmp_path / "from_records")


@pytest.mark.parametrize("unsupported", [{"chrom": "1", "pos": 1, "p_value": 0.5}, 42, None])
def test_an_input_that_is_not_results_is_a_type_error_not_a_silent_no_op(tmp_path: Path, unsupported: Any) -> None:
    with pytest.raises(TypeError, match="association_results must be a path, a list of records or a pandas DataFrame"):
        generate_all_plots(unsupported, tmp_path)


def test_a_missing_file_is_reported(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="nope.tsv"):
        generate_all_plots(tmp_path / "nope.tsv", tmp_path / "plots")


def test_an_unsupported_extension_is_reported(tmp_path: Path) -> None:
    path = tmp_path / "results.xlsx"
    path.write_text("not a table")

    with pytest.raises(ValueError, match="Unsupported association results format '.xlsx'"):
        generate_all_plots(path, tmp_path / "plots")


# ---------------------------------------------------------------------------
# A result without a p-value is not drawn
# ---------------------------------------------------------------------------


def test_results_without_a_p_value_leave_the_plots_exactly_as_they_were(tmp_path: Path) -> None:
    clean = render(rows(), tmp_path / "clean")

    with_gaps = render(with_missing(rows()), tmp_path / "gaps")

    assert sorted(with_gaps) == ["manhattan", "qq"]  # the Q-Q plot used to raise on "NA" and None
    assert with_gaps == clean


def test_the_number_left_out_is_logged(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING):
        generate_all_plots(with_missing(rows()), tmp_path)

    assert "4 of 28 association results have no p-value" in caplog.text


def test_nothing_to_draw_is_said_not_silently_skipped(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING):
        nothing = generate_all_plots([], tmp_path / "empty")
        all_failed = generate_all_plots([{"chrom": "1", "pos": 5, "p_value": NAN}] * 3, tmp_path / "failed")

    assert nothing == {} and all_failed == {}
    assert "No association results were given" in caplog.text
    assert "None of the 3 association results has a p-value" in caplog.text


# ---------------------------------------------------------------------------
# Figures are closed
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("make_input", [rows, lambda: pd.DataFrame(rows()), lambda: with_missing(rows())])
def test_no_figure_is_left_open(tmp_path: Path, make_input: Any) -> None:
    before = plt.get_fignums()

    generate_all_plots(make_input(), tmp_path)

    assert plt.get_fignums() == before


def test_figures_are_closed_even_when_a_later_plot_fails(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    bad_pca = tmp_path / "pca.json"
    bad_pca.write_text('{"components": "not a matrix"}')
    before = plt.get_fignums()

    with caplog.at_level(logging.WARNING):
        written = generate_all_plots(rows(), tmp_path / "plots", pca_file=bad_pca)

    assert "manhattan" in written and "pca" not in written
    assert plt.get_fignums() == before


# ---------------------------------------------------------------------------
# manhattan_plot and regional_plot read the same records
# ---------------------------------------------------------------------------


def test_manhattan_does_not_draw_or_count_results_without_a_p_value() -> None:
    data = [{"chrom": "1", "pos": 100 * (i + 1), "p_value": 10.0 ** -(i + 1)} for i in range(5)]
    mixed = with_missing(data)  # three rows (NaN, None, "NA") and one without the key, all on chromosomes 1 and 2

    fig = manhattan_plot(mixed, significance_threshold=1e-3)
    try:
        drawn = points(fig)
        title = fig.axes[0].get_title()
    finally:
        plt.close(fig)

    assert sorted(drawn[:, 1]) == pytest.approx([1.0, 2.0, 3.0, 4.0, 5.0])  # never 50
    assert "n=5 variants, 3 GW-significant" in title  # p <= 1e-3: 1e-3, 1e-4, 1e-5; not the four without a p-value


@pytest.mark.parametrize("key", ["p_value", "pvalue", "pval", "P"])
def test_manhattan_reads_the_p_value_under_every_spelling(key: str) -> None:
    data = [{"chrom": "1", "pos": 100 * (i + 1), key: 10.0 ** -(i + 1)} for i in range(5)]

    fig = manhattan_plot(data)
    try:
        drawn = points(fig)
    finally:
        plt.close(fig)

    assert drawn[:, 1] == pytest.approx([1.0, 2.0, 3.0, 4.0, 5.0])  # used to be 0 for pvalue and P (p defaulted to 1)


def test_manhattan_accepts_a_dataframe() -> None:
    frame = pd.DataFrame({"CHR": ["1", "1", "2"], "BP": [10, 20, 30], "P": [0.1, 0.01, 0.001]})

    fig = manhattan_plot(frame)
    try:
        drawn = points(fig)
    finally:
        plt.close(fig)

    assert drawn[:, 1] == pytest.approx([1.0, 2.0, 3.0])


def test_a_zero_p_value_is_still_drawn_at_the_cap() -> None:
    fig = manhattan_plot([{"chrom": "1", "pos": 10, "p_value": 0.0}, {"chrom": "1", "pos": 20, "p_value": 0.1}])
    try:
        drawn = points(fig)
    finally:
        plt.close(fig)

    assert sorted(drawn[:, 1]) == pytest.approx([1.0, 50.0])


def test_gene_annotations_by_variant_index_still_address_the_results_as_given() -> None:
    """Indices refer to the list the caller passed, whether or not earlier results have a p-value."""
    data = [
        {"chrom": "1", "pos": 100, "p_value": NAN},
        {"chrom": "1", "pos": 200, "p_value": None},
        {"chrom": "1", "pos": 300, "p_value": 1e-2},
        {"chrom": "1", "pos": 400, "p_value": 1e-4},
        {"chrom": "1", "pos": 500, "p_value": 1e-3},
    ]

    fig = manhattan_plot(data, gene_annotations=[{"variant_index": 3, "gene_name": "G4"}])
    try:
        annotations = [t for t in fig.axes[0].texts if t.get_text() == "G4"]
        assert len(annotations) == 1
        assert annotations[0].xy == pytest.approx((400.0, 4.0))  # the fourth result: position 400, p = 1e-4
    finally:
        plt.close(fig)


def test_label_top_n_labels_the_best_result_that_has_a_p_value() -> None:
    data = [
        {"chrom": "1", "pos": 100, "p_value": 1e-3, "variant_id": "good"},
        {"chrom": "1", "pos": 200, "p_value": NAN, "variant_id": "failed"},
        {"chrom": "1", "pos": 300, "p_value": 1e-2, "variant_id": "weaker"},
    ]

    fig = manhattan_plot(data, label_top_n=1, suggestive_threshold=None)
    try:
        labels = [t.get_text() for t in fig.axes[0].texts]
    finally:
        plt.close(fig)

    assert labels == ["good"]


def test_regional_plot_leaves_results_without_a_p_value_out() -> None:
    data = [
        {"chrom": "1", "pos": 100, "p_value": 1e-2},
        {"chrom": "1", "pos": 200, "p_value": NAN},
        {"chrom": "1", "pos": 300, "pvalue": 1e-4},
        {"chrom": "1", "pos": 400, "p_value": "NA"},
    ]

    fig = regional_plot(data, "1", 0, 1_000)
    try:
        drawn = points(fig)
    finally:
        plt.close(fig)

    assert drawn[:, 0].tolist() == [100, 300]
    assert drawn[:, 1] == pytest.approx([2.0, 4.0])


def test_a_region_with_only_failed_tests_has_no_plot(caplog: pytest.LogCaptureFixture) -> None:
    with caplog.at_level(logging.WARNING):
        fig = regional_plot([{"chrom": "1", "pos": 100, "p_value": NAN}], "1", 0, 1_000)

    assert fig is None
    assert "No results with a p-value in region 1:0-1000" in caplog.text


# ---------------------------------------------------------------------------
# The workflow
# ---------------------------------------------------------------------------

# Twelve samples and seven variants; rsM and rsN are constant, so their tests fail.
GENOTYPES = {
    "rsA": [0, 1, 2, 0, 1, 2, 0, 1, 2, 1, 0, 2],
    "rsM": [0] * 12,
    "rsB": [2, 2, 1, 0, 0, 1, 2, 1, 0, 0, 1, 1],
    "rsC": [0, 0, 1, 1, 2, 2, 1, 0, 2, 1, 2, 0],
    "rsN": [1] * 12,
    "rsD": [1, 0, 2, 1, 0, 1, 2, 2, 0, 1, 1, 0],
    "rsE": [0, 2, 1, 1, 0, 2, 2, 0, 1, 0, 2, 1],
}
NOISE = [0.12, -0.31, 0.25, 0.4, -0.22, 0.05, -0.37, 0.19, -0.08, 0.33, -0.15, 0.28]
PHENOTYPE = [1.5 + 0.9 * a + 0.3 * c + e for a, c, e in zip(GENOTYPES["rsA"], GENOTYPES["rsC"], NOISE)]  # noqa: B905


def write_dataset(directory: Path, genotypes: Dict[str, List[int]]) -> Dict[str, Path]:
    calls = {0: "0/0", 1: "0/1", 2: "1/1"}
    samples = [f"S{i + 1}" for i in range(12)]
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples),
    ]
    for i, (snp, row) in enumerate(genotypes.items()):
        fields = ["chr1" if i < 4 else "chr2", str(100 * (i + 1)), snp, "A", "G", "60", "PASS", ".", "GT"]
        lines.append("\t".join(fields + [calls[g] for g in row]))
    vcf = directory / "data.vcf"
    vcf.write_text("\n".join(lines) + "\n")
    phenotypes = directory / "phenotypes.tsv"
    phenotypes.write_text("\n".join(["sample_id\tyield"] + [f"S{i + 1}\t{y}" for i, y in enumerate(PHENOTYPE)]) + "\n")
    return {"vcf": vcf, "phenotypes": phenotypes}


def workflow_config(files: Dict[str, Path], out: Path) -> Dict[str, Any]:
    return {
        "vcf_path": str(files["vcf"]),
        "phenotype_path": str(files["phenotypes"]),
        "output_dir": str(out),
        "results_dir": str(out),
        "trait": "yield",
        "model": "linear",
        "quality_control": {"min_maf": 0.0, "max_missing": 1.0, "min_hwe_p": 0.0},
        "ld_pruning": {"r2_threshold": 1.0},  # keep every variant
        "structure": {"pca_components": 1},
        "significance_threshold": 1.0,
    }


def test_execute_gwas_workflow_writes_its_manhattan_and_q_q_plots(tmp_path: Path) -> None:
    files = write_dataset(tmp_path, GENOTYPES)
    out = tmp_path / "out"
    before = plt.get_fignums()

    result = execute_gwas_workflow(workflow_config(files, out))

    assert result["status"] == "completed"
    plots = result["outputs"]["plots"]
    for name, filename in (("manhattan", "manhattan_plot.png"), ("qq", "qq_plot.png")):
        assert Path(plots[name]) == out / filename
        assert_png((out / filename).read_bytes())
    assert plt.get_fignums() == before  # every figure the workflow opened is closed


def test_the_workflow_manhattan_draws_the_tests_that_ran_and_not_the_failed_ones(tmp_path: Path) -> None:
    """Its input is the five variants with a p-value: the same plot as drawing their linregress p-values."""
    files = write_dataset(tmp_path, GENOTYPES)
    out = tmp_path / "out"
    result = execute_gwas_workflow(workflow_config(files, out))
    tested = [
        row
        for row in result["outputs"]["association_results"]
        if isinstance(row["p_value"], float) and not math.isnan(row["p_value"])
    ]
    assert len(tested) == 5 and result["outputs"]["n_tests_failed"] == 2

    redrawn = render(tested, tmp_path / "redrawn")

    assert redrawn["manhattan"] == (out / "manhattan_plot.png").read_bytes()
    assert redrawn["qq"] == (out / "qq_plot.png").read_bytes()
    fig = manhattan_plot(tested)
    try:
        assert fig.axes[0].get_title().startswith("Manhattan Plot (n=5 variants")
    finally:
        plt.close(fig)


def test_the_workflow_says_why_when_it_could_not_draw_the_core_plots(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # Identical variants: the first principal component (a covariate of every test) is collinear with the genotype,
    # so the design is singular and every test fails.
    same = GENOTYPES["rsA"]
    files = write_dataset(tmp_path, {"rsA": same, "rsB": same, "rsC": same})
    out = tmp_path / "out"

    with caplog.at_level(logging.WARNING):
        result = execute_gwas_workflow(workflow_config(files, out))

    tests = result["outputs"]["association_results"]
    assert result["status"] == "completed" and tests
    assert result["outputs"]["n_tests_failed"] == len(tests)
    assert "manhattan" not in result["outputs"]["plots"] and "qq" not in result["outputs"]["plots"]
    reason = f"none of the {len(tests)} association tests produced a p-value"
    assert f"Core plots not produced: manhattan, qq ({reason})" in caplog.text


def test_the_workflow_draws_the_maf_spectrum_from_the_mafs_of_its_results(tmp_path: Path) -> None:
    files = write_dataset(tmp_path, GENOTYPES)
    out = tmp_path / "out"

    result = execute_gwas_workflow(workflow_config(files, out))

    assert_png((out / "maf_spectrum_plot.png").read_bytes())
    assert Path(result["outputs"]["plots"]["maf_spectrum"]) == out / "maf_spectrum_plot.png"
    # the rows carry the MAF counted from the genotypes: rsA has 12 alt alleles of 24
    mafs = {i: row["maf"] for i, row in enumerate(result["outputs"]["association_results"])}
    expected = {}
    for i, genotypes in enumerate(GENOTYPES.values()):
        alt = sum(genotypes) / (2 * len(genotypes))
        expected[i] = min(alt, 1 - alt)
    assert mafs == pytest.approx(expected)
    assert stats.linregress(GENOTYPES["rsA"], PHENOTYPE).pvalue < 0.05  # the fixture still has a real signal
