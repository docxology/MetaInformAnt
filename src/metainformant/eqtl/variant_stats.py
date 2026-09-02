"""bcftools stats parsing and per-sample/population variant summaries.

Example:
    >>> from metainformant.eqtl.variant_stats import parse_bcftools_stats
    >>> stats = parse_bcftools_stats("SN\\t0\\tnumber of records:\\t12\\n")
    >>> stats["n_records"]
    12
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def parse_bcftools_stats(stats_text: str) -> dict[str, Any]:
    """Parse record/SNP/indel counts and Ti/Tv ratio from bcftools stats text.

    Args:
        stats_text: Raw stdout of ``bcftools stats``.

    Returns:
        Dict with keys ``n_records``, ``n_snps``, ``n_indels``,
        ``ts_tv_ratio`` (0.0 when absent).
    """
    stats: dict[str, Any] = {
        "n_records": 0,
        "n_snps": 0,
        "n_indels": 0,
        "ts_tv_ratio": 0.0,
    }

    for line in stats_text.split("\n"):
        if not line.startswith("SN"):
            continue
        if "number of records:" in line:
            stats["n_records"] = int(line.split("\t")[-1])
        elif "number of SNPs:" in line:
            stats["n_snps"] = int(line.split("\t")[-1])
        elif "number of indels:" in line:
            stats["n_indels"] = int(line.split("\t")[-1])

    # Ti/Tv from bcftools stats
    for line in stats_text.split("\n"):
        if line.startswith("TSTV"):
            parts = line.split("\t")
            if len(parts) >= 5:
                try:
                    stats["ts_tv_ratio"] = float(parts[4])
                except (ValueError, IndexError):
                    pass

    return stats


def _run_bcftools_stats(vcf_path: Path) -> str:
    """Run ``bcftools stats`` on a VCF and return stdout.

    Raises:
        RuntimeError: If bcftools is missing or fails.
    """
    try:
        result = subprocess.run(
            ["bcftools", "stats", str(vcf_path)],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("bcftools not found on PATH; install bcftools to compute variant stats") from exc
    if result.returncode != 0:
        raise RuntimeError(f"bcftools stats failed for {vcf_path}: {result.stderr[:300]}")
    return result.stdout


def compute_sample_stats(vcf_path: Path, output_json: Path) -> dict[str, Any]:
    """Compute per-sample variant statistics from a VCF and write JSON.

    Args:
        vcf_path: Input compressed VCF.
        output_json: Destination JSON file (parent dirs created).

    Returns:
        Stats dict with ``vcf_file`` plus parsed bcftools fields.
    """
    stdout = _run_bcftools_stats(vcf_path)
    stats = parse_bcftools_stats(stdout)
    stats["vcf_file"] = str(vcf_path)

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(stats, f, indent=2)

    return stats


def compute_popgen_summary(
    merged_vcf: Path,
    sample_stats: list[dict[str, Any]],
    output_json: Path,
) -> dict[str, Any]:
    """Compute population-level genetics summary from a merged VCF.

    Args:
        merged_vcf: Merged population VCF.
        sample_stats: Per-sample stats dicts (embedded in the summary).
        output_json: Destination JSON file (parent dirs created).

    Returns:
        Summary dict with totals, Ti/Tv ratio, and per-sample stats.
    """
    stdout = _run_bcftools_stats(merged_vcf)
    parsed = parse_bcftools_stats(stdout)

    summary: dict[str, Any] = {
        "n_samples": len(sample_stats),
        "total_variants": parsed["n_records"],
        "total_snps": parsed["n_snps"],
        "total_indels": parsed["n_indels"],
        "ts_tv_ratio": parsed["ts_tv_ratio"],
        "per_sample_stats": sample_stats,
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(summary, f, indent=2)

    return summary


def compute_allele_frequencies(merged_vcf: Path, output_tsv: Path) -> int:
    """Extract per-site allele frequencies from a merged VCF via bcftools query.

    Args:
        merged_vcf: Merged population VCF.
        output_tsv: Destination TSV with header ``chrom pos ref alt af``.

    Returns:
        Number of variant rows written.

    Raises:
        RuntimeError: If bcftools query fails.
    """
    cmd = [
        "bcftools",
        "query",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n",
        str(merged_vcf),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"bcftools query failed: {result.stderr[:300]}")

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as f:
        f.write("chrom\tpos\tref\talt\taf\n")
        f.write(result.stdout)

    n_lines = result.stdout.count("\n")
    logger.info(f"Wrote allele frequencies for {n_lines} variants -> {output_tsv}")
    return n_lines


def load_real_expression_data_hint() -> str:
    """Return pointer to the real-expression loader (kept for API symmetry).

    The real-expression loader lives in
    :func:`metainformant.eqtl.synthetic.load_real_expression_data`.
    """
    return "metainformant.eqtl.synthetic.load_real_expression_data"
