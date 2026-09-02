"""Handler behavior tests with real synthetic inputs (zero mocks)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from metainformant.mcp.tools import (
    core_tools,
    dna_tools,
    gwas_tools,
    math_tools,
    protein_tools,
    visualization_tools,
)


def test_dna_translate_round_trip() -> None:
    result = dna_tools._handle_translate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    assert result["protein"] == "MAIVMGR*KGAR*"
    assert result["length_aa"] == 11


def test_dna_composition_gc_content_matches_formula() -> None:
    seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    result = dna_tools._handle_composition(seq)
    gc = sum(1 for c in seq if c in "GC") / len(seq)
    assert abs(result["gc_content"] - gc) < 1e-9


def test_core_config_validate_real_yaml(tmp_path: Path) -> None:
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("a: 1\nb: [x, y]\n", encoding="utf-8")
    result = core_tools._handle_config_validate(str(cfg))
    assert result["valid"] is True
    assert result["top_level_keys"] == ["a", "b"]


def test_core_config_validate_missing_file_reports_error(tmp_path: Path) -> None:
    result = core_tools._handle_config_validate(str(tmp_path / "absent.yaml"))
    assert result["valid"] is False
    assert "not found" in result["error"]


def test_core_path_resolve_containment(tmp_path: Path) -> None:
    child = tmp_path / "sub"
    child.mkdir()
    result = core_tools._handle_path_resolve(str(child / "f.txt"), parent=str(tmp_path))
    assert result["within_parent"] is True
    assert result["exists"] is False


def test_math_popgen_summary_deterministic() -> None:
    a = math_tools._handle_popgen_summary(0.3)
    b = math_tools._handle_popgen_summary(0.3)
    assert a == b
    # genotype frequencies from p=0.3: AA=0.09, Aa=0.42, aa=0.49
    assert abs(a["freq_AA"] - 0.09) < 1e-9
    assert abs(a["freq_Aa"] - 0.42) < 1e-9
    assert abs(a["freq_aa"] - 0.49) < 1e-9


def test_protein_align_identity_identical_vs_mutated() -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLRIE"
    same = protein_tools._handle_align(seq, seq)
    mut = protein_tools._handle_align(seq, seq[:15] + "A" + seq[16:])
    assert same["identity"] == 1.0
    assert mut["identity"] < same["identity"]
    assert same["aligned_length"] == len(seq)


def test_gwas_phenotype_summary_real_table(tmp_path: Path) -> None:
    ph = tmp_path / "pheno.csv"
    pd.DataFrame({"t1": [1.0, 2.0, 3.0, 4.0], "t2": [10, 20, 30, 40]}, index=[f"s{i}" for i in range(4)]).to_csv(ph)
    out = tmp_path / "out"
    result = gwas_tools._handle_phenotype_summary(str(ph), output_dir=str(out))
    assert result["n_samples"] == 4
    assert abs(result["numeric_summary"]["t1"]["mean"] - 2.5) < 1e-9
    persisted = json.loads((out / "phenotype_summary.json").read_text())
    assert persisted == {k: v for k, v in result.items() if k != "output"}


def test_gwas_hwe_extreme_counts() -> None:
    assert gwas_tools._handle_hwe(0, 50, 50)["hwe_pvalue"] >= 0.0


def test_viz_chart_line_writes_real_png(tmp_path: Path) -> None:
    out = tmp_path / "viz"
    result = visualization_tools._handle_chart("line", [1.0, 2.0, 3.0], [2.0, 4.0, 6.0], str(out))
    png = Path(result["output"])
    assert png.exists() and png.stat().st_size > 0
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_viz_chart_rejects_unknown_kind(tmp_path: Path) -> None:
    result = visualization_tools._handle_chart("pie", [1.0], None, str(tmp_path))
    assert "error" in result
