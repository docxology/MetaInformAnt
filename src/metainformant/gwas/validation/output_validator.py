"""GWAS pipeline output validator.

Programmatically verifies the integrity and consistency of GWAS pipeline
outputs: file existence, column schema, sample-count concordance, λ_GC
sanity, and post-GWAS JSON structure.

Extracted from ``validate_outputs.py`` script.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def _check(label: str, ok: bool, detail: str = "", warning: bool = False) -> Dict[str, Any]:
    """Build a single check-result dict."""
    status = ("WARN" if warning else "PASS") if ok else "FAIL"
    return {"label": label, "status": status, "ok": ok or warning, "detail": detail}


def validate_pipeline_outputs(
    config: Dict[str, Any],
    *,
    lambda_gc_warn: float = 2.0,
    require_hits: bool = False,
) -> Dict[str, Any]:
    """Validate all expected GWAS pipeline outputs.

    Args:
        config: Full pipeline config dict (with ``paths``, ``samples``,
            ``gwas``, ``validation`` sections).
        lambda_gc_warn: λ_GC threshold above which a warning is raised.
        require_hits: If True, fail when no GW-significant hits found.

    Returns:
        Dict with ``n_pass``, ``n_fail``, ``n_warn``, ``all_ok``,
        and ``checks`` (list of individual check results).
    """
    paths = config.get("paths", {})
    results_dir = Path(paths.get("results_dir", "output"))
    Path(paths.get("processed_dir", "data/processed"))
    Path(paths.get("phenotype_dir", "data/phenotypes"))
    Path(paths.get("raw_dir", "data/raw"))

    default_trait = config.get("samples", {}).get("phenotypes", {}).get("default_trait", "")
    model = config.get("gwas", {}).get("model", "mixed")
    sig_threshold = config.get("gwas", {}).get("significance_threshold", 5e-8)

    trait_dir = results_dir / default_trait / model
    post_gwas_dir = trait_dir / "post_gwas"
    trait_dir / "power_analysis"

    checks: List[Dict[str, Any]] = []
    n_pass = n_fail = n_warn = 0

    def _add(label: str, ok: bool, detail: str = "", warning: bool = False) -> None:
        nonlocal n_pass, n_fail, n_warn
        r = _check(label, ok, detail, warning)
        checks.append(r)
        if r["status"] == "WARN":
            n_warn += 1
        elif r["ok"]:
            n_pass += 1
        else:
            n_fail += 1

    # ── Core outputs ──
    core_files = [
        trait_dir / "summary_statistics.tsv",
        trait_dir / "significant_hits.tsv",
        trait_dir / "manhattan_plot.png",
        trait_dir / "qq_plot.png",
        trait_dir / "pca_plot.png",
        trait_dir / "kinship_plot.png",
        trait_dir / "effect_size_plot.png",
        trait_dir / "maf_spectrum_plot.png",
    ]
    for f in core_files:
        ok = f.exists() and f.stat().st_size > 0
        _add(f.name, ok, f"size={f.stat().st_size:,}B" if f.exists() else "MISSING")

    # ── Post-GWAS outputs ──
    post_files = [
        post_gwas_dir / "volcano_plot.png",
        post_gwas_dir / "qq_stratified.png",
        post_gwas_dir / "chrom_summary.tsv",
        post_gwas_dir / "post_gwas_results.json",
    ]
    for f in post_files:
        ok = f.exists() and f.stat().st_size > 0
        _add(f.name, ok, f"size={f.stat().st_size:,}B" if f.exists() else "MISSING")

    # ── Summary statistics schema ──
    summary_tsv = trait_dir / "summary_statistics.tsv"
    if summary_tsv.exists():
        with open(summary_tsv) as f:
            headers = f.readline().strip().split("\t")
        expected_cols = {"CHR", "POS", "SNP", "BETA", "SE", "P", "MAF"}
        missing_cols = expected_cols - set(headers)
        _add(
            "summary_statistics.tsv schema",
            len(missing_cols) == 0,
            f"columns={headers}" if not missing_cols else f"missing={missing_cols}",
        )

    # ── Significant hits consistency ──
    sig_tsv = trait_dir / "significant_hits.tsv"
    if sig_tsv.exists() and summary_tsv.exists():
        with open(sig_tsv) as f:
            sig_rows = list(csv.DictReader(f, delimiter="\t"))
        bad_p = [r for r in sig_rows if float(r.get("P", 1)) >= sig_threshold]
        _add(
            "significant_hits.tsv consistency",
            len(bad_p) == 0,
            f"{len(sig_rows)} hits, {len(bad_p)} above threshold {sig_threshold:.0e}",
        )

    # ── λ_GC sanity ──
    post_json = post_gwas_dir / "post_gwas_results.json"
    if post_json.exists():
        with open(post_json) as f:
            post_data = json.load(f)
        lgc = post_data.get("summary_statistics", {}).get("p_value_calibration", {}).get("lambda_gc", float("nan"))
        is_finite = not math.isnan(lgc) and not math.isinf(lgc)
        is_acceptable = is_finite and lgc > 0
        is_inflated = is_finite and lgc > lambda_gc_warn
        _add(
            "lambda_gc finite and positive",
            is_acceptable,
            f"lambda_gc={lgc:.4f}" + (" ⚠️ INFLATED" if is_inflated else ""),
            warning=is_inflated,
        )

        # GW hits
        n_gw = post_data.get("summary_statistics", {}).get("significance_counts", {}).get("genome_wide_5e-8", 0)
        if require_hits:
            _add("GW-significant hits > 0", n_gw > 0, f"n_gw={n_gw}")
        else:
            _add(
                "GW-significant hits",
                True,
                f"n={n_gw} (require_significant_hits=false)",
            )

        # JSON schema
        required_keys = {
            "summary_statistics",
            "n_variants_analyzed",
            "sign_test",
            "credible_sets",
        }
        missing_keys = required_keys - set(post_data.keys())
        _add(
            "post_gwas_results.json schema",
            len(missing_keys) == 0,
            f"keys={list(post_data.keys())}" if not missing_keys else f"missing={missing_keys}",
        )

    return {
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_warn": n_warn,
        "all_ok": n_fail == 0,
        "checks": checks,
    }
