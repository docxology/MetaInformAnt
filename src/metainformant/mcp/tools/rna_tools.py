"""RNA MCP tools: campaign status, tau, atlas plots, conservation profiles.

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

    heatmap: tau_table rows = species, columns = tissues, values = mean tau.
    strips: tau_table needs per-gene columns 'tau' and 'orthology_class'.
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


def _load_species_profiles(profile_dir: Path) -> dict:
    """Load per-species TPM tables (<species>.csv|tsv, genes x conditions) from a dir."""
    import pandas as pd

    profiles: dict[str, pd.DataFrame] = {}
    for path in sorted(profile_dir.iterdir()):
        if path.suffix.lower() not in {".csv", ".tsv", ".tab"}:
            continue
        sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
        profiles[path.stem] = pd.read_csv(path, sep=sep, index_col=0)
    if not profiles:
        return {"error": f"no per-species CSV/TSV tables found in {profile_dir}"}
    return profiles


def _handle_conservation_profiles(profile_dir: str, output_dir: str) -> dict:
    """Descriptive cross-species TPM profile-conservation summary to output dir.

    profile_dir must contain per-species TPM distribution CSVs as consumed by
    metainformant.rna.analysis.conservation_profiles. Descriptive only.
    """
    from pathlib import Path as _Path

    from metainformant.rna.analysis import conservation_profiles as cp

    profiles_path = _Path(profile_dir).expanduser()
    if not profiles_path.is_dir():
        return {"error": f"profile_dir not found: {profiles_path}"}
    species_profiles = _load_species_profiles(profiles_path)
    if isinstance(species_profiles, dict) and "error" in species_profiles:
        return species_profiles
    out_dir = validate_output_dir(output_dir)
    try:
        summary = cp.summarize_profile_conservation(cp.compute_profile_conservation(species_profiles))
    except Exception as exc:  # real input errors surface to the caller
        return {"error": f"{type(exc).__name__}: {exc}"}
    out_path = out_dir / "profile_conservation_summary.csv"
    summary.to_csv(out_path)
    return {
        "output": str(out_path),
        "n_rows": int(summary.shape[0]),
        "columns": list(summary.columns)[:20],
        "note": "descriptive only; no inferential tests",
    }


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

CONSERVATION_SPEC: dict[str, Any] = {
    "name": "rna_conservation_profiles",
    "description": "Descriptive cross-species TPM profile-conservation summary from per-species profile CSVs.",
    "input_schema": {
        "type": "object",
        "properties": {
            "profile_dir": {"type": "string"},
            "output_dir": {"type": "string"},
        },
        "required": ["profile_dir", "output_dir"],
    },
    "handler": _handle_conservation_profiles,
    "writes": "output-dir-only",
}


def _handle_normalize_counts(
    counts_table: str, output_dir: str, method: str = "tpm", gene_lengths: str | None = None
) -> dict:
    """Normalize a raw count matrix (genes x samples) to CPM/TPM/RPKM/quantile.

    gene_lengths: optional CSV/TSV with index=gene, single length column
    (required for tpm/rpkm).
    """
    from metainformant.rna.analysis.expression_core import normalize_counts

    frame = read_table(counts_table)
    lengths = None
    if gene_lengths is not None:
        lengths = read_table(gene_lengths, index_col=None)
        if lengths.shape[1] == 1:
            lengths = lengths.set_index(lengths.columns[0]).iloc[:, 0]
        else:
            lengths = lengths.set_index(lengths.columns[0])
        lengths = lengths.iloc[:, 0] if lengths.ndim > 1 else lengths
    result = normalize_counts(frame, method=method, gene_lengths=lengths)
    out_dir = validate_output_dir(output_dir)
    out_path = out_dir / f"normalized_{method}.csv"
    result.to_csv(out_path)
    return {
        "output": str(out_path),
        "method": method,
        "n_genes": int(result.shape[0]),
        "n_samples": int(result.shape[1]),
    }


def _handle_duplication_specificity(expression_table: str, bridge_table: str, species: str, output_dir: str) -> dict:
    """Descriptive tau-by-orthology-class summary (Xu & Colgan 2025 methods).

    Joins the species expression table with the ortholog bridge, computes tau
    per transcript, and summarizes by 1:1 vs multicopy orthology class.
    """
    from metainformant.rna.analysis.tissue_specificity import (
        compute_tau,
        duplication_specificity_summary,
        join_expression_with_orthology,
    )

    expression = read_table(expression_table)
    bridge = read_table(bridge_table)
    try:
        joined = join_expression_with_orthology(expression, bridge, species)
        tissue_cols = [c for c in joined.columns if c not in ("orthogroup", "orthology_class")]
        tau = compute_tau(joined[tissue_cols]).reindex(joined.index)
        summary = duplication_specificity_summary(tau, joined["orthology_class"])
    except (ValueError, KeyError) as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}
    out_dir = validate_output_dir(output_dir)
    out_path = out_dir / "duplication_specificity_summary.csv"
    summary.to_csv(out_path)
    return {
        "output": str(out_path),
        "n_transcripts": int(len(joined)),
        "classes": sorted(str(c) for c in summary.index),
        "note": "descriptive only; inferential comparison gated post-freeze",
    }


NORMALIZE_SPEC: dict[str, Any] = {
    "name": "rna_normalize_counts",
    "description": "Normalize a raw gene x sample count matrix (cpm|tpm|rpkm|quantile).",
    "input_schema": {
        "type": "object",
        "properties": {
            "counts_table": {"type": "string"},
            "output_dir": {"type": "string"},
            "method": {"type": "string", "enum": ["cpm", "tpm", "rpkm", "quantile"]},
            "gene_lengths": {"type": "string"},
        },
        "required": ["counts_table", "output_dir"],
    },
    "handler": _handle_normalize_counts,
    "writes": "output-dir-only",
}

DUP_SPEC: dict[str, Any] = {
    "name": "rna_duplication_specificity",
    "description": "Descriptive tau summary across 1:1 vs multicopy orthology classes for one species.",
    "input_schema": {
        "type": "object",
        "properties": {
            "expression_table": {"type": "string"},
            "bridge_table": {"type": "string"},
            "species": {"type": "string"},
            "output_dir": {"type": "string"},
        },
        "required": ["expression_table", "bridge_table", "species", "output_dir"],
    },
    "handler": _handle_duplication_specificity,
    "writes": "output-dir-only",
}

ALL_SPECS: list[dict[str, Any]] = [
    TOOL_SPEC,
    TAU_SPEC,
    ATLAS_SPEC,
    CONSERVATION_SPEC,
    NORMALIZE_SPEC,
    DUP_SPEC,
]
