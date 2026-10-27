"""Regression tests: ``scripts/eqtl/run_eqtl_real.py`` does not label synthetic data as real.

Defect: the script scans real RNA-seq expression against SYNTHETIC genotypes placed at SYNTHETIC gene positions
(``parse_gene_positions`` reads the chromosome accession from each transcript id and invents the TSS positions). It
nevertheless saved those invented positions as ``gene_positions_real.tsv``, logged "Parsing gene annotations", wrote
``"data_source": "Real RNA-seq from ENA/NCBI"`` into the summary and named every result ``*_real``. Anyone reading the
output tree could take the associations for findings on real genotypes and a real annotation.

The script is run as a separate process in a fresh directory on a small kallisto-style quantification tree that the test
writes itself (no network, no external data), and the tests read what it produced.
"""

from __future__ import annotations

import ast
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterator, List

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "eqtl" / "run_eqtl_real.py"

# Kallisto target ids as Amalgkit writes them: the chromosome accession is part of the id. The first two targets sit on
# one chromosome, as most transcripts of a real run do.
TARGETS = [
    "lcl|NC_037638.1_mrna_XM_623972.6_1",
    "lcl|NC_037638.1_mrna_XM_393171.7_2",
    "lcl|NC_037639.1_mrna_XM_624220.5_3",
    "lcl|NC_037640.1_mrna_XM_001120.4_4",
    "lcl|NC_037641.1_mrna_XM_006561.3_5",
    "lcl|NC_037642.1_mrna_XM_393300.2_6",
]
N_SAMPLES = 14
SAMPLES = [f"SRR90000{s:02d}" for s in range(N_SAMPLES)]


# Deterministic expression that varies across the samples for every target (mean well above the TPM filter)
TPM = np.random.default_rng(123).uniform(5.0, 50.0, size=(N_SAMPLES, len(TARGETS)))


def tpm(sample_index: int, target_index: int) -> float:
    return float(TPM[sample_index, target_index])


def write_quantification(root: Path) -> None:
    quant = root / "output" / "amalgkit" / "apis_mellifera_all" / "work" / "quant"
    for s, sample in enumerate(SAMPLES):
        sample_dir = quant / sample
        sample_dir.mkdir(parents=True)
        lines = ["target_id\tlength\teff_length\test_counts\ttpm"]
        for t, target in enumerate(TARGETS):
            lines.append(f"{target}\t1500\t1300.0\t{tpm(s, t) * 13:.3f}\t{tpm(s, t):.4f}")
        (sample_dir / "abundance.tsv").write_text("\n".join(lines) + "\n")


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory: pytest.TempPathFactory) -> Iterator[Path]:
    """Run the script once in a fresh directory and hand back that directory."""
    root = tmp_path_factory.mktemp("eqtl_real_run")
    write_quantification(root)
    completed = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=root,
        env={**os.environ, "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert completed.returncode == 0, completed.stderr[-3000:]
    yield root


@pytest.fixture(scope="module")
def results(run_dir: Path) -> Path:
    return run_dir / "output" / "eqtl" / "amellifera" / "results"


@pytest.fixture(scope="module")
def log_text(run_dir: Path) -> str:
    return (run_dir / "output" / "eqtl" / "amellifera" / "logs" / "eqtl_real_analysis.log").read_text()


def read_json(path: Path) -> Dict[str, object]:
    return json.loads(path.read_text())


def accession(target: str) -> str:
    match = re.search(r"NC_\d+\.\d+", target)
    assert match is not None
    return match.group(0)


# ---------------------------------------------------------------------------
# The positions file
# ---------------------------------------------------------------------------


def test_invented_gene_positions_are_saved_under_a_synthetic_name(results: Path) -> None:
    assert (results / "gene_positions_synthetic.tsv").is_file()
    assert not (results / "gene_positions_real.tsv").exists()


def test_the_positions_file_holds_accession_chromosomes_and_evenly_spaced_invented_positions(results: Path) -> None:
    positions = pd.read_csv(results / "gene_positions_synthetic.tsv", sep="\t")

    assert list(positions.columns) == ["gene_id", "chrom", "tss_position"]
    assert positions["gene_id"].tolist() == TARGETS
    assert positions["chrom"].tolist() == [accession(t) for t in TARGETS]  # the only real information in the file
    # A real annotation never puts genes at one fixed spacing: these positions depend on the row index alone.
    assert positions["tss_position"].diff().dropna().nunique() == 1


# ---------------------------------------------------------------------------
# No output claims to be real unless it is
# ---------------------------------------------------------------------------


def test_no_output_file_is_named_real_unless_it_holds_real_data(run_dir: Path, results: Path) -> None:
    outputs = [p for p in (run_dir / "output" / "eqtl").rglob("*") if p.is_file()]
    assert outputs  # guard: the run did write something

    mislabelled = sorted(p.name for p in outputs if p.stem.endswith("_real"))
    assert mislabelled == []
    # What is real keeps saying so: the expression matrix is the quantification the test wrote
    expression = pd.read_csv(results / "real_expression_matrix.tsv", sep="\t", index_col=0)
    assert expression.shape == (len(TARGETS), N_SAMPLES)
    assert expression.loc[TARGETS[2], SAMPLES[3]] == pytest.approx(tpm(3, 2), rel=1e-4)


def test_results_depending_on_synthetic_genotypes_say_so_in_their_names(results: Path, run_dir: Path) -> None:
    plots = run_dir / "output" / "eqtl" / "amellifera" / "plots"
    scan = pd.read_csv(results / "cis_eqtl_results_synthetic_genotypes.tsv", sep="\t")

    assert len(scan) > 0  # the scan really ran
    assert set(scan["gene_id"]).issubset(set(TARGETS))
    assert (results / "synthetic_genotype_matrix.tsv").is_file()
    assert (results / "summary_stats_synthetic_genotypes.json").is_file()
    assert (results / "run_summary_synthetic_genotypes.json").is_file()
    assert sorted(p.name for p in plots.iterdir()) == [
        "summary_plot_synthetic_genotypes.png",
        "top_eqtl_boxplot_synthetic_genotypes.png",
        "volcano_plot_synthetic_genotypes.png",
    ]


# ---------------------------------------------------------------------------
# Log, docstring and summaries
# ---------------------------------------------------------------------------


def test_synthetic_variant_ids_are_unique_and_every_pair_is_tested_once(results: Path) -> None:
    """Two transcripts share a chromosome here; their synthetic variants must still have distinct ids."""
    variants = pd.read_csv(results / "variant_positions.tsv", sep="\t")
    scan = pd.read_csv(results / "cis_eqtl_results_synthetic_genotypes.tsv", sep="\t")

    assert len(variants) == 3 * len(TARGETS)
    assert variants["variant_id"].is_unique
    assert not scan.duplicated(["gene_id", "variant_id"]).any()
    assert (scan["status"] == "success").all()


def test_the_log_says_the_gene_positions_and_genotypes_are_synthetic(log_text: str) -> None:
    assert "SYNTHETIC gene positions" in log_text
    assert "no GFF/GTF annotation is used" in log_text
    assert "SYNTHETIC genotypes" in log_text
    assert "Parsing gene annotations" not in log_text  # nothing is parsed from an annotation
    assert "Real A. mellifera eQTL Analysis" not in log_text  # the old banner and completion line


def test_the_log_warns_that_the_associations_are_not_findings(log_text: str) -> None:
    warnings = [line for line in log_text.splitlines() if " - WARNING - " in line]

    assert any("not biological findings" in line for line in warnings)


def test_the_module_docstring_states_that_gene_positions_are_synthetic() -> None:
    docstring = " ".join((ast.get_docstring(ast.parse(SCRIPT.read_text())) or "").split())

    assert "Gene positions are SYNTHETIC" in docstring
    assert "no GFF/GTF annotation is used" in docstring
    assert "Genotypes are SYNTHETIC" in docstring
    assert "Expression is REAL" in docstring


def test_the_summary_records_where_each_input_came_from(results: Path) -> None:
    summary = read_json(results / "summary_stats_synthetic_genotypes.json")

    assert summary["data_source"] != "Real RNA-seq from ENA/NCBI"
    assert "SYNTHETIC genotypes" in str(summary["data_source"])
    assert "SYNTHETIC gene positions" in str(summary["data_source"])
    assert str(summary["expression_data"]).startswith("real")
    assert str(summary["genotype_data"]).startswith("synthetic")
    assert str(summary["gene_position_data"]).startswith("synthetic")
    assert summary["n_samples"] == N_SAMPLES


def test_the_run_summary_does_not_call_the_analysis_real(results: Path) -> None:
    run_summary = read_json(results / "run_summary_synthetic_genotypes.json")

    assert not str(run_summary["analysis_type"]).startswith("Real")
    assert "synthetic genotypes" in str(run_summary["analysis_type"])
    assert "gene positions are synthetic" in str(run_summary["note"])
    assert str(run_summary["genotype_data"]).startswith("synthetic")
    assert str(run_summary["gene_position_data"]).startswith("synthetic")
    assert run_summary["n_genes"] == len(TARGETS)


def test_the_volcano_title_names_the_synthetic_genotypes() -> None:
    source = SCRIPT.read_text()

    assert "(real expression, synthetic genotypes)" in source
    assert "(Real Expression Data)" not in source


def test_no_script_text_calls_synthetic_positions_real() -> None:
    """Static guard on the source: the old file name and the old data-source label are gone."""
    source = SCRIPT.read_text()
    forbidden: List[str] = ["gene_positions_real", "Real RNA-seq from ENA/NCBI", "Parsing gene annotations"]

    assert [text for text in forbidden if text in source] == []
