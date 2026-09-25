"""End-to-end checks for the metainformant CLI entry point.

The CLI must stay honest: every command exercised here either performs real
work and exits 0, or fails loudly with a non-zero exit code. Fabricated
outputs (e.g. the removed ``math selection replay`` PNG stub) are pinned as
rejected.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Child processes must exercise this repository's src tree (see tests/conftest.py
# for the in-process equivalent).
CLI_ENV = {**os.environ, "PYTHONPATH": str(REPO_ROOT / "src")}


def _run_cli(
    *args: str, cwd: Path = REPO_ROOT, timeout: int = 300
) -> subprocess.CompletedProcess[str]:
    """Run ``python -m metainformant`` from the repository root."""
    return subprocess.run(
        [sys.executable, "-m", "metainformant", *args],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(cwd),
        env=CLI_ENV,
        check=False,
    )


def test_module_invocation_shows_help() -> None:
    """``python -m metainformant --help`` exits 0 and lists the real subcommands."""
    result = _run_cli("--help")

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout
    for command in (
        "protein",
        "quality",
        "rna",
        "gwas",
        "life-events",
        "simulation",
        "ontology",
        "phenotype",
        "networks",
    ):
        assert command in result.stdout


def test_no_arguments_shows_help() -> None:
    """Bare invocation prints help and exits 0."""
    result = _run_cli()

    assert result.returncode == 0
    assert "usage:" in result.stdout


def test_removed_math_subcommand_is_rejected(tmp_path: Path) -> None:
    """``math selection replay`` fabricated PNG-header files and was removed."""
    result = _run_cli("math", "selection", "replay", "--dest", str(tmp_path))

    assert result.returncode != 0
    assert "invalid choice" in (result.stderr + result.stdout).lower()


def test_bare_domain_command_is_an_error(tmp_path: Path) -> None:
    """A domain invoked without a subcommand must not exit 0 silently."""
    result = _run_cli("simulation", cwd=tmp_path)

    assert result.returncode == 1
    assert "simulation" in (result.stdout + result.stderr).lower()


def test_simulation_run_executes_real_workflow(tmp_path: Path) -> None:
    """``simulation run`` runs the simulation workflow and saves its result JSON."""
    output_dir = tmp_path / "output" / "simulation"

    result = _run_cli(
        "simulation",
        "run",
        "--model",
        "rna_expression",
        "--n",
        "3",
        "--output",
        str(output_dir),
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    result_file = output_dir / "rna_expression_result.json"
    assert result_file.exists()
    payload = json.loads(result_file.read_text())
    # --model maps to simulation_type and --n to the model's size parameter
    # (sample count for RNA expression).
    assert payload["simulation_type"] == "rna_expression"
    assert payload["config"]["n_samples"] == 3
    assert "expression_matrix_shape" in payload


def test_simulation_run_rejects_unknown_model(tmp_path: Path) -> None:
    """An unsupported --model fails loudly instead of exiting 0."""
    result = _run_cli("simulation", "run", "--model", "not_a_model", cwd=tmp_path)

    assert result.returncode == 1
    assert "Error" in result.stdout + result.stderr


def test_phenotype_run_executes_real_pipeline(tmp_path: Path) -> None:
    """``phenotype run`` runs the phenotype pipeline over a JSON dataset."""
    data = tmp_path / "specimens.json"
    data.write_text(
        json.dumps(
            [
                {"specimen": "a", "head_length": 1.2},
                {"specimen": "b", "head_length": 1.4},
            ]
        )
    )
    output_dir = tmp_path / "output" / "phenotype"

    result = _run_cli(
        "phenotype",
        "run",
        "--input",
        str(data),
        "--type",
        "morphological",
        "--output",
        str(output_dir),
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    payload = json.loads((output_dir / "pipeline_result.json").read_text())
    assert payload["success"] is True
    assert payload["steps_executed"] == ["load", "validate", "analyze", "summarize"]


def test_phenotype_run_missing_input_fails(tmp_path: Path) -> None:
    """A missing dataset fails loudly instead of exiting 0."""
    result = _run_cli(
        "phenotype", "run", "--input", str(tmp_path / "missing.json"), cwd=tmp_path
    )

    assert result.returncode == 1
    assert "Error" in result.stdout + result.stderr


def test_networks_run_builds_and_exports_network(tmp_path: Path) -> None:
    """``networks run`` builds a network from an edge list and exports results."""
    edges = tmp_path / "edges.csv"
    edges.write_text(
        "source,target,weight\na,b,0.9\na,c,0.8\nb,c,0.85\nd,e,0.9\nd,f,0.8\ne,f,0.85\n"
    )
    output_dir = tmp_path / "output" / "networks"

    result = _run_cli(
        "networks",
        "run",
        "--input",
        str(edges),
        "--output",
        str(output_dir),
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    network_summary = json.loads((output_dir / "network.json").read_text())
    assert network_summary is not None
    assert (output_dir / "metrics.json").exists()
    assert (output_dir / "communities.json").exists()
    communities = json.loads((output_dir / "communities.json").read_text())
    assert len(communities) == 6  # every node is assigned to a community


def test_networks_run_rejects_malformed_edge_list(tmp_path: Path) -> None:
    """An edge list without source/target columns fails loudly."""
    bad = tmp_path / "bad.csv"
    bad.write_text("from,to\na,b\n")

    result = _run_cli("networks", "run", "--input", str(bad), cwd=tmp_path)

    assert result.returncode == 1
    assert "source" in result.stdout + result.stderr


def test_ontology_run_skips_cleanly_without_gene_annotations(tmp_path: Path) -> None:
    """With GWAS stage inputs but no gene annotations, the ontology workflow
    writes its skipped summary and exits 0 without network access."""
    results_base = tmp_path / "results" / "ant" / "base"
    results_base.mkdir(parents=True)
    (results_base / "summary_statistics.tsv").write_text(
        "chrom\tpos\tsnp\tp_value\tbeta\tse\tmaf\n"
        "1\t100\trs1\t0.001\t0.2\t0.1\t0.3\n"
        "1\t200\trs2\t0.010\t0.1\t0.1\t0.2\n"
    )
    post_gwas = results_base / "post_gwas"
    post_gwas.mkdir()
    (post_gwas / "post_gwas_results.json").write_text(
        json.dumps({"gene_annotations": []})
    )
    config = tmp_path / "workflow.yaml"
    config.write_text("paths:\n  results_dir: results\nontology:\n  taxon_id: 7222\n")

    result = _run_cli(
        "ontology",
        "run",
        "--input",
        str(config),
        "--phenotype",
        "ant",
        "--model",
        "base",
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    summary = json.loads(
        (results_base / "ontology" / "ontology_summary.json").read_text()
    )
    assert summary["status"] == "skipped"


def test_ontology_run_requires_input_config(tmp_path: Path) -> None:
    """--input is required; argparse rejects the invocation."""
    result = _run_cli(
        "ontology", "run", "--phenotype", "ant", "--model", "base", cwd=tmp_path
    )

    assert result.returncode == 2
    assert "--input" in result.stderr


def test_ontology_run_missing_config_fails(tmp_path: Path) -> None:
    """A missing config file fails loudly instead of exiting 0."""
    result = _run_cli(
        "ontology",
        "run",
        "--input",
        str(tmp_path / "missing.yaml"),
        "--phenotype",
        "ant",
        "--model",
        "base",
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert "Error" in result.stdout + result.stderr


def test_quality_run_writes_verification_report(tmp_path: Path) -> None:
    """``quality run`` performs real docs-vs-code verification and writes a report."""
    docs = tmp_path / "docs"
    docs.mkdir()
    (docs / "readme.md").write_text("# Docs\n\nNo code examples here.\n")
    src = tmp_path / "src"
    src.mkdir()
    (src / "tiny.py").write_text("def add(a, b):\n    return a + b\n")
    report = tmp_path / "report.md"

    result = _run_cli(
        "quality",
        "run",
        "--docs-dir",
        str(docs),
        "--src-dir",
        str(src),
        "--output",
        str(report),
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    assert report.exists()
    assert "Total violations found:** 0" in report.read_text()


def test_quality_run_missing_docs_dir_fails(tmp_path: Path) -> None:
    """A missing docs directory fails loudly instead of exiting 0."""
    report = tmp_path / "report.md"

    result = _run_cli(
        "quality",
        "run",
        "--docs-dir",
        str(tmp_path / "missing"),
        "--output",
        str(report),
        cwd=tmp_path,
    )

    assert result.returncode == 1
