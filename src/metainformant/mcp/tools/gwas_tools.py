"""GWAS MCP tools: phenotype summary and HWE/titv descriptive stats (read-only)."""

from __future__ import annotations

from typing import Any

from metainformant.mcp.tools._spec import read_table, validate_output_dir


def _handle_phenotype_summary(phenotype_table: str, output_dir: str | None = None) -> dict:
    """Descriptive summary of a phenotype table (index/sample column + trait columns)."""
    frame = read_table(phenotype_table)
    numeric = frame.select_dtypes(include=["number"])
    summary: dict[str, Any] = {
        "n_samples": int(frame.shape[0]),
        "n_traits": int(frame.shape[1]),
        "trait_columns": list(frame.columns),
    }
    if not numeric.empty:
        desc = numeric.describe()
        summary["numeric_summary"] = {
            col: {stat: float(desc.loc[stat, col]) for stat in ("mean", "std", "min", "max")} for col in desc.columns
        }
    if output_dir is not None:
        out_dir = validate_output_dir(output_dir)
        out_path = out_dir / "phenotype_summary.json"
        from metainformant.mcp.tools._spec import dump_json

        dump_json(summary, out_path)
        summary["output"] = str(out_path)
    return summary


def _handle_hwe(obs_hets: int, obs_hom1: int, obs_hom2: int) -> dict:
    """Descriptive Hardy-Weinberg p-value (chi-square approximation) from counts."""
    from metainformant.gwas.analysis.summary_stats import calculate_hwe_pvalue

    return {"hwe_pvalue": calculate_hwe_pvalue(obs_hets, obs_hom1, obs_hom2)}


TOOL_SPEC: dict[str, Any] = {
    "name": "gwas_phenotype_summary",
    "description": "Descriptive per-trait summary statistics for a phenotype table.",
    "input_schema": {
        "type": "object",
        "properties": {
            "phenotype_table": {"type": "string"},
            "output_dir": {"type": "string"},
        },
        "required": ["phenotype_table"],
    },
    "handler": _handle_phenotype_summary,
    "writes": "output-dir-only",
}

HWE_SPEC: dict[str, Any] = {
    "name": "gwas_hwe",
    "description": "Hardy-Weinberg equilibrium p-value from genotype counts.",
    "input_schema": {
        "type": "object",
        "properties": {
            "obs_hets": {"type": "integer"},
            "obs_hom1": {"type": "integer"},
            "obs_hom2": {"type": "integer"},
        },
        "required": ["obs_hets", "obs_hom1", "obs_hom2"],
    },
    "handler": _handle_hwe,
    "writes": "read-only",
}

ALL_SPECS: list[dict[str, Any]] = [TOOL_SPEC, HWE_SPEC]
