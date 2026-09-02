"""RNA MCP tools: campaign status, tau tissue-specificity, atlas plots.

Read-only campaign reporting and explicit-output-dir descriptive figure
generation. All statistics are DESCRIPTIVE ONLY per the round-3 boundary.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

from metainformant.mcp.tools._spec import read_table, validate_output_dir

REPO_ROOT = Path(__file__).resolve().parents[4]
RATE_SCRIPT = REPO_ROOT / "scripts" / "rna" / "analyze_campaign_rate.py"


def _handle_campaign_status(log_paths: list[str], total: int | None = None) -> dict:
    """Safe read-only campaign rate analysis via the pinned analyzer script."""
    if not RATE_SCRIPT.exists():
        return {"error": f"analyzer script missing: {RATE_SCRIPT}"}
    cmd = [sys.executable, str(RATE_SCRIPT), "--json", *[str(p) for p in log_paths]]
    if total is not None:
        cmd[2:2] = ["--total", str(total)]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        return {"error": proc.stderr.strip()[-2000:], "returncode": proc.returncode}
    import json

    return json.loads(proc.stdout)


def _handle_tau(expression_table: str, output_dir: str, lowest_fraction: float = 0.10) -> dict:
    """Compute tau tissue-specificity from a genes x tissues TPM table (CSV/TSV).

    Writes tau_per_gene.csv into the explicit output directory. Descriptive.
    """
    from metainformant.rna.analysis.tissue_specificity import compute_tau, filter_low_expression

    frame = read_table(expression_table)
    filtered = filter_low_expression(frame, lowest_fraction=lowest_fraction)
    tau = compute_tau(filtered)
    out_dir = validate_output_dir(output_dir)
    out_path = out_dir / "tau_per_gene.csv"
    tau.rename("tau").to_csv(out_path, index_label="gene")
    valid = tau.dropna()
    return {
        "n_genes_input": int(frame.shape[0]),
        "n_genes_used": int(tau.notna().sum()),
        "tau_mean": float(valid.mean()) if len(valid) else None,
        "tau_median": float(valid.median()) if len(valid) else None,
        "output": str(out_path),
    }


def _handle_atlas_plot(
    tau_table: str,
    output_dir: str,
    plot_type: str = "heatmap",
    title: str = "Tissue specificity (tau) by species and tissue",
) -> dict:
    """Render an atlas-style tau plot (heatmap|strips) into the output dir.

    tau_table: rows = species, columns = tissues, values = mean tau in [0,1].
    """
    from metainformant.rna.analysis import atlas_plots

    frame = read_table(tau_table)
    out_dir = validate_output_dir(output_dir)
    if plot_type == "heatmap":
        path = atlas_plots.plot_tau_heatmap(frame, out_dir / "tau_heatmap.png", title=title)
    elif plot_type == "strips":
        # strips expect a per-gene table with tau + orthology_class columns
        if "tau" not in frame.columns or "orthology_class" not in frame.columns:
            return {"error": "strips plot requires columns: tau, orthology_class"}
        path = atlas_plots.plot_tau_orthology_strips(frame, out_dir / "tau_strips.png")
    else:
        return {"error": f"unknown plot_type: {plot_type} (use heatmap|strips)"}
    return {"output": str(path), "plot_type": plot_type}


TOOL_SPEC: dict[str, Any] = {
    "name": "rna_campaign_status",
    "description": "Read-only streaming-orchestrator campaign rate/ETA analysis from log files.",
    "input_schema": {
        "type": "object",
        "properties": {
            "log_paths": {"type": "array", "items": {"type": "string"}},
            "total": {"type": "integer"},
        },
        "required": ["log_paths"],
    },
    "handler": _handle_campaign_status,
    "writes": "read-only",
}

TAU_SPEC: dict[str, Any] = {
    "name": "rna_tau",
    "description": "Tau tissue-specificity (Yanai 2005; Xu & Colgan 2025 protocol) on a genes x tissues TPM table.",
    "input_schema": {
        "type": "object",
        "properties": {
            "expression_table": {"type": "string"},
            "output_dir": {"type": "string"},
            "lowest_fraction": {"type": "number"},
        },
        "required": ["expression_table", "output_dir"],
    },
    "handler": _handle_tau,
    "writes": "output-dir-only",
}

ATLAS_SPEC: dict[str, Any] = {
    "name": "rna_atlas_plot",
    "description": "Atlas-style tau heatmap or orthology strip plot (Leader+ 2024 grammar) to an explicit output dir.",
    "input_schema": {
        "type": "object",
        "properties": {
            "tau_table": {"type": "string"},
            "output_dir": {"type": "string"},
            "plot_type": {"type": "string", "enum": ["heatmap", "strips"]},
            "title": {"type": "string"},
        },
        "required": ["tau_table", "output_dir"],
    },
    "handler": _handle_atlas_plot,
    "writes": "output-dir-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC, TAU_SPEC, ATLAS_SPEC]
