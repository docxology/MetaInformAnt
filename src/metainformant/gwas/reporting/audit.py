"""GWAS pipeline audit and provenance utilities.

Provides functions for SHA-256 checksumming, library version detection,
and structured audit log generation.  These were previously embedded in
the ``09_audit.py`` orchestrator script.
"""

from __future__ import annotations

import datetime
import hashlib
import platform
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def sha256_file(path: str | Path, *, chunk_size: int = 1 << 20) -> str:
    """Compute the SHA-256 hex-digest of a file.

    Args:
        path: File to hash.
        chunk_size: Read chunk size in bytes (default 1 MB).

    Returns:
        Lowercase hex-digest string.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            data = f.read(chunk_size)
            if not data:
                break
            h.update(data)
    return h.hexdigest()


def checksum_directory(
    root: str | Path,
    *,
    extensions: Sequence[str] = (
        ".png",
        ".tsv",
        ".json",
        ".html",
        ".yaml",
        ".vcf",
        ".csv",
        ".fna",
    ),
) -> Dict[str, str]:
    """Compute SHA-256 checksums for all matching files under *root*.

    Args:
        root: Root directory to walk recursively.
        extensions: File extensions to include (case-insensitive match).

    Returns:
        Dict mapping ``relative_path → sha256_hex``.
    """
    root = Path(root)
    ext_set = set(e.lower() for e in extensions)
    checksums: Dict[str, str] = {}
    for p in sorted(root.rglob("*")):
        if p.is_file() and p.suffix.lower() in ext_set:
            try:
                checksums[str(p.relative_to(root))] = sha256_file(p)
            except Exception as exc:
                checksums[str(p.relative_to(root))] = f"ERROR: {exc}"
    return checksums


def library_versions(
    libs: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    """Detect installed versions of key scientific libraries.

    Args:
        libs: Library names to check.  Defaults to numpy, matplotlib,
              yaml, scipy, metainformant.

    Returns:
        Dict mapping ``library_name → version_string``.
    """
    if libs is None:
        libs = ("numpy", "matplotlib", "yaml", "scipy", "metainformant")

    versions: Dict[str, str] = {}
    for lib in libs:
        try:
            m = __import__(lib)
            versions[lib] = getattr(m, "__version__", "unknown")
        except ImportError:
            versions[lib] = "not installed"
    return versions


def python_environment() -> Dict[str, str]:
    """Capture the Python runtime environment.

    Returns:
        Dict with ``python``, ``platform``, and ``libraries`` keys.
    """
    return {
        "python": f"Python {platform.python_version()} on {platform.system()} {platform.release()}",
        "platform": platform.platform(),
        "libraries": library_versions(),
    }


def generate_audit_log(
    *,
    config: Dict[str, Any],
    results_dir: str | Path,
    processed_dir: str | Path,
    raw_dirs: Optional[List[str | Path]] = None,
    variant_count: int = 0,
    sig_count: int = 0,
    lambda_gc: Optional[float] = None,
) -> Dict[str, Any]:
    """Build a structured audit log dictionary.

    The returned dict is suitable for serialising to YAML.

    Args:
        config: Full pipeline config dict.
        results_dir: Path to results directory.
        processed_dir: Path to processed data directory.
        raw_dirs: Optional list of raw/phenotype dirs to checksum.
        variant_count: Number of variants tested.
        sig_count: Number of GW-significant variants.
        lambda_gc: Genomic inflation factor.

    Returns:
        Structured audit dictionary.
    """
    results_dir = Path(results_dir)
    processed_dir = Path(processed_dir)

    # Checksums
    result_checksums = checksum_directory(results_dir) if results_dir.exists() else {}
    proc_checksums = checksum_directory(processed_dir) if processed_dir.exists() else {}
    raw_checksums: Dict[str, str] = {}
    for d in raw_dirs or []:
        d = Path(d)
        if d.exists():
            for p in sorted(d.glob("*")):
                if p.is_file() and p.suffix in (".yaml", ".csv", ".tsv"):
                    raw_checksums[str(p)] = sha256_file(p)

    # Parameters snapshot
    gwas_cfg = config.get("gwas", {})
    qc_cfg = config.get("qc", {})
    pg_cfg = config.get("post_gwas", {})
    params = {
        "trait": config.get("samples", {}).get("phenotypes", {}).get("default_trait", "unknown"),
        "model": gwas_cfg.get("model", "mixed"),
        "significance_threshold": gwas_cfg.get("significance_threshold", 5e-8),
        "maf_threshold": qc_cfg.get("min_maf", 0.01),
        "max_missing": qc_cfg.get("max_missing_rate", 0.1),
        "hwe_threshold": qc_cfg.get("hwe_p_threshold", 1e-6),
        "ld_r2_prune": config.get("ld_pruning", {}).get("r2_threshold", 0.2),
        "n_pca_components": gwas_cfg.get("n_pca_components", 10),
        "top_hits_to_label": pg_cfg.get("top_hits_to_label", 5),
        "gene_annotation_radius_bp": pg_cfg.get("gene_annotation_radius_bp", 100000),
    }

    return {
        "run_timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "pipeline_version": "metainformant-gwas",
        "environment": python_environment(),
        "parameters": params,
        "results_summary": {
            "n_variants_tested": variant_count,
            "n_gw_significant": sig_count,
            "lambda_gc": lambda_gc,
        },
        "file_counts": {
            "result_files": len(result_checksums),
            "processed_files": len(proc_checksums),
            "raw_metadata_files": len(raw_checksums),
            "total_audited": len(result_checksums) + len(proc_checksums) + len(raw_checksums),
        },
    }
