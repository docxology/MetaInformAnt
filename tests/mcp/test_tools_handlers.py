"""Handler behavior tests with real synthetic inputs (real implementations)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from metainformant.mcp.tools import (
    core_tools,
    dna_tools,
    gwas_tools,
    math_tools,
    protein_tools,
    rna_tools,
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


def test_core_newick_parse_leaves_and_root() -> None:
    result = core_tools._handle_newick_parse("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")
    assert result["leaves"] == ["A", "B", "C", "D"]
    assert result["n_leaves"] == 4
    assert result["n_internal"] == 3
    assert result["root"] == "Internal_0"


def test_dna_transcribe_and_reverse_complement() -> None:
    fwd = dna_tools._handle_transcribe("ATGGCC")
    rev = dna_tools._handle_transcribe("ATGGCC", reverse_complement=True)
    assert fwd["rna"] == "AUGGCC"
    assert rev["rna"] == "GGCCAU"


def test_protein_composition_uniform() -> None:
    result = protein_tools._handle_composition("MKTAY")
    comp = result["amino_acid_composition"]
    assert abs(comp["M"] - 20.0) < 1e-9
    assert result["length"] == 5
    assert result["molecular_weight"] > 0


def test_rna_conservation_profiles_real_tables(tmp_path: Path) -> None:
    rng = np.random.default_rng(3)
    pdir = tmp_path / "profiles"
    pdir.mkdir()
    for sp in ("apis", "bombus"):
        pd.DataFrame(
            rng.lognormal(3, 1, size=(20, 3)),
            columns=["brain", "muscle", "gut"],
            index=[f"g{i}" for i in range(20)],
        ).to_csv(pdir / f"{sp}.csv")
    out = tmp_path / "out"
    result = rna_tools._handle_conservation_profiles(str(pdir), str(out))
    assert "error" not in result
    assert result["n_rows"] == 20
    assert (out / "profile_conservation_summary.csv").exists()


def test_rna_conservation_profiles_missing_dir(tmp_path: Path) -> None:
    result = rna_tools._handle_conservation_profiles(str(tmp_path / "nope"), str(tmp_path))
    assert "error" in result


def test_gwas_association_summary_real_table(tmp_path: Path) -> None:
    rng = np.random.default_rng(11)
    n = 100
    res = pd.DataFrame(
        {
            "snp": [f"rs{i}" for i in range(n)],
            "chrom": ["chr1"] * n,
            "pos": np.arange(n),
            "p_value": rng.uniform(0, 1, n),
            "beta": rng.normal(0, 0.3, n),
            "se": rng.uniform(0.05, 0.2, n),
            "MAF": rng.uniform(0.05, 0.5, n),
        }
    )
    p = tmp_path / "results.csv"
    res.to_csv(p, index=False)
    out = tmp_path / "out"
    result = gwas_tools._handle_association_summary(str(p), output_dir=str(out))
    assert "error" not in result
    assert result["n_variants"] == n
    assert (out / "association_summary.json").exists()


def test_gwas_association_summary_missing_columns(tmp_path: Path) -> None:
    p = tmp_path / "bad.csv"
    pd.DataFrame({"p_value": [0.1]}).to_csv(p, index=False)
    result = gwas_tools._handle_association_summary(str(p))
    assert "missing required columns" in result["error"]


def test_viz_divergence_heatmap_real_matrix(tmp_path: Path) -> None:
    vals = np.abs(np.random.default_rng(2).normal(0, 0.3, (5, 5)))
    sym = (vals + vals.T) / 2
    for i in range(5):
        sym[i, i] = 0.0
    p = tmp_path / "div.csv"
    pd.DataFrame(sym, index=list("abcde"), columns=list("abcde")).to_csv(p)
    out = tmp_path / "viz"
    result = visualization_tools._handle_divergence_heatmap(str(p), str(out))
    assert Path(result["output"]).exists()
    assert result["n_species"] == 5


def test_rna_normalize_counts_tpm_real_tables(tmp_path: Path) -> None:
    rng = np.random.default_rng(5)
    counts = pd.DataFrame(
        rng.poisson(50, size=(20, 4)), columns=[f"s{i}" for i in range(4)], index=[f"g{i}" for i in range(20)]
    )
    ct = tmp_path / "counts.csv"
    counts.to_csv(ct)
    gl = tmp_path / "lengths.csv"
    pd.DataFrame({"gene": [f"g{i}" for i in range(20)], "length": rng.integers(500, 3000, 20)}).to_csv(gl, index=False)
    out = tmp_path / "norm"
    result = rna_tools._handle_normalize_counts(str(ct), str(out), method="tpm", gene_lengths=str(gl))
    assert result["n_genes"] == 20
    norm = pd.read_csv(out / "normalized_tpm.csv", index_col=0)
    # TPM columns each sum to 1e6 (within float tolerance)
    totals = norm.sum(axis=0)
    assert (totals - 1e6).abs().max() < 1.0


def test_rna_duplication_specificity_real_bridge(tmp_path: Path) -> None:
    rng = np.random.default_rng(5)
    expr = pd.DataFrame(
        rng.lognormal(3, 1, size=(30, 3)),
        columns=["brain", "muscle", "gut"],
        index=[f"APIS_G{i:04d}" for i in range(30)],
    )
    et = tmp_path / "expr.tsv"
    expr.to_csv(et, sep="\t")
    rows = {}
    for i, og in enumerate([f"OG{i:04d}" for i in range(30)]):
        rows[og] = {
            "apis": f"APIS_G{i:04d}",
            "bombus": f"BOMB_G{i:04d},BOMB_G{i:04d}B" if i % 3 == 0 else f"BOMB_G{i:04d}",
        }
    bt = tmp_path / "bridge.tsv"
    pd.DataFrame.from_dict(rows, orient="index").to_csv(bt, sep="\t")
    out = tmp_path / "out"
    result = rna_tools._handle_duplication_specificity(str(et), str(bt), "apis", str(out))
    assert result["n_transcripts"] == 30
    assert "multicopy" in result["classes"] and "single_copy" in result["classes"]
    summary = pd.read_csv(out / "duplication_specificity_summary.csv", index_col="orthology_class")
    assert summary["count"].sum() == 30


def test_rna_duplication_specificity_unknown_species(tmp_path: Path) -> None:
    et = tmp_path / "expr.csv"
    pd.DataFrame({"a": [1.0]}, index=["g1"]).to_csv(et)
    bt = tmp_path / "bridge.csv"
    pd.DataFrame({"apis": ["g1"]}, index=["og1"]).to_csv(bt)
    result = rna_tools._handle_duplication_specificity(str(et), str(bt), "nope", str(tmp_path))
    assert "error" in result
