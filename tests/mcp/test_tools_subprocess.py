"""Subprocess round-trip tests for MCP tools (real interpreter, no mocks)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _run_handler(module: str, call: str) -> dict:
    """Execute a tool handler in a fresh subprocess and parse its JSON output."""
    code = "import json; from metainformant.mcp.tools import " f"{module} as m; print(json.dumps(m.{call}))"
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    return json.loads(proc.stdout.strip().splitlines()[-1])


def test_subprocess_round_trip_dna_translate() -> None:
    result = _run_handler("dna_tools", '_handle_translate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")')
    assert result["protein"] == "MAIVMGR*KGAR*"


def test_subprocess_round_trip_math_popgen() -> None:
    result = _run_handler("math_tools", "_handle_popgen_summary(0.3)")
    assert abs(result["freq_Aa"] - 0.42) < 1e-9


def test_subprocess_round_trip_protein_align() -> None:
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLRIE"
    result = _run_handler("protein_tools", f'_handle_align("{seq}", "{seq}")')
    assert result["identity"] == 1.0


def test_subprocess_round_trip_gwas_hwe() -> None:
    result = _run_handler("gwas_tools", "_handle_hwe(10, 20, 20)")
    assert 0.0 <= result["hwe_pvalue"] <= 1.0


def test_subprocess_round_trip_rna_tau_real_table(tmp_path: Path) -> None:
    rng = np.random.default_rng(7)
    expr = pd.DataFrame(
        rng.lognormal(3, 1, size=(30, 3)),
        columns=["a", "b", "c"],
        index=[f"g{i}" for i in range(30)],
    )
    table = tmp_path / "expr.tsv"
    expr.to_csv(table, sep="\t")
    call = f'_handle_tau("{table}", "{tmp_path / "out"}")'
    result = _run_handler("rna_tools", call)
    assert result["n_genes_input"] == 30
    assert (tmp_path / "out" / "tau_per_gene.csv").exists()
    csv = pd.read_csv(tmp_path / "out" / "tau_per_gene.csv", index_col="gene")
    assert csv["tau"].between(0, 1).dropna().all()
