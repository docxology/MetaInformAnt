"""VCF file utilities for GWAS pipeline operations.

Provides reusable functions for VCF discovery, compression, indexing,
merging, subsampling, and variant counting. These are the shared
building blocks used by all orchestrator scripts that operate on VCF files.

All functions delegate to ``bcftools``, ``bgzip``, and ``tabix`` CLI tools,
which must be available on ``$PATH``.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, Optional

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def discover_sample_vcfs(
    processed_dir: str | Path,
    *,
    exclude_patterns: tuple[str, ...] = ("merged", "subset", "subsample", "partial"),
) -> List[Path]:
    """Auto-discover per-sample VCF files in a directory.

    Scans for ``*.vcf`` files and excludes any whose filename contains
    a pattern from *exclude_patterns* (e.g. merged artifacts).

    Args:
        processed_dir: Directory to scan.
        exclude_patterns: Substrings to skip (case-sensitive match on stem).

    Returns:
        Sorted list of per-sample VCF paths.
    """
    processed_dir = Path(processed_dir)
    vcfs = sorted(processed_dir.glob("*.vcf"))
    result = [v for v in vcfs if not any(pat in v.name for pat in exclude_patterns)]
    logger.info("Discovered %d per-sample VCFs in %s", len(result), processed_dir)
    return result


def count_variants(vcf_path: str | Path) -> int:
    """Count data lines (variants) in a VCF, excluding the header.

    Args:
        vcf_path: Path to VCF or VCF.gz file.

    Returns:
        Number of variant records.
    """
    result = subprocess.run(
        ["bcftools", "view", "-H", str(vcf_path)],
        capture_output=True,
        text=True,
    )
    return sum(1 for _ in result.stdout.splitlines())


def bgzip_and_index(vcf_path: str | Path) -> Optional[Path]:
    """Compress a VCF with bgzip and create a tabix index.

    If the ``.vcf.gz`` file already exists and is non-trivial (>1 KB),
    the compression step is skipped.  Indexing is always attempted if
    ``*.tbi`` is absent.

    Args:
        vcf_path: Path to an uncompressed ``.vcf`` file.

    Returns:
        Path to the ``.vcf.gz`` file, or ``None`` on failure.
    """
    vcf_path = Path(vcf_path)
    gz = vcf_path.with_suffix(".vcf.gz") if vcf_path.suffix == ".vcf" else vcf_path
    try:
        if vcf_path.suffix == ".vcf":
            if not gz.exists() or gz.stat().st_size < 1000:
                tmp = gz.with_suffix(".tmp")
                subprocess.run(
                    f"bgzip -c {vcf_path} > {tmp} && mv {tmp} {gz}",
                    shell=True,
                    check=True,
                )
        tbi = Path(str(gz) + ".tbi")
        if not tbi.exists():
            subprocess.run(["tabix", "-p", "vcf", str(gz)], check=True)
        return gz
    except Exception as exc:
        logger.warning("Failed to compress/index %s: %s", vcf_path, exc)
        return None


def merge_vcfs(
    vcf_gz_paths: List[Path],
    output_path: str | Path,
    *,
    threads: int = 1,
) -> Path:
    """Multi-sample VCF merge via ``bcftools merge``.

    Args:
        vcf_gz_paths: Bgzipped + indexed VCF files to merge.
        output_path: Output path (should end in ``.vcf.gz``).
        threads: bcftools thread count.

    Returns:
        Path to the merged VCF.
    """
    output_path = Path(output_path)
    vcf_list = " ".join(str(v) for v in vcf_gz_paths)
    logger.info("Merging %d VCFs → %s", len(vcf_gz_paths), output_path)
    subprocess.run(
        f"bcftools merge -0 --threads {threads} -O z -o {output_path} {vcf_list}",
        shell=True,
        check=True,
    )
    subprocess.run(["tabix", "-p", "vcf", str(output_path)], check=True)
    return output_path


def subsample_vcf(
    input_vcf: str | Path,
    output_vcf: str | Path,
    *,
    fraction: Optional[float] = None,
    n_variants: Optional[int] = None,
    seed: int = 42,
) -> int:
    """Random-fraction or fixed-count VCF subsampling.

    Uses ``bcftools view`` + ``awk`` for speed.  Exactly one of
    *fraction* or *n_variants* must be provided.

    Args:
        input_vcf: Input VCF (plain or bgzipped).
        output_vcf: Output path (will be bgzipped + indexed).
        fraction: Fraction of variants to keep, in (0, 1].
        n_variants: Exact number of variants to target.
        seed: Random seed for reproducibility.

    Returns:
        Number of variants in the output.

    Raises:
        ValueError: If neither or both of fraction/n_variants are set.
    """
    if (fraction is None) == (n_variants is None):
        raise ValueError("Specify exactly one of fraction or n_variants")

    if fraction is not None:
        if not 0 < fraction <= 1:
            raise ValueError(f"fraction must be in (0, 1], got {fraction}")
        effective_fraction = fraction
    else:
        total = count_variants(str(input_vcf))
        effective_fraction = min(1.0, n_variants / max(total, 1))  # type: ignore[operator]
        logger.info(
            "Total variants: %d, requesting %d → fraction=%.4f",
            total,
            n_variants,
            effective_fraction,
        )

    output_vcf = Path(output_vcf)
    tmp = Path(str(output_vcf) + ".tmp")
    logger.info(
        "Subsampling %s at %.1f%% (seed=%d)",
        input_vcf,
        effective_fraction * 100,
        seed,
    )
    cmd = (
        f"bcftools view -h {input_vcf} > {tmp} && "
        f"bcftools view -H {input_vcf} | "
        f"awk 'BEGIN{{srand({seed})}} rand() < {effective_fraction}' >> {tmp} && "
        f"bgzip -c {tmp} > {output_vcf} && "
        f"tabix -p vcf {output_vcf} && "
        f"rm {tmp}"
    )
    subprocess.run(cmd, shell=True, check=True)
    n = count_variants(str(output_vcf))
    logger.info("Subsampled VCF contains %d variants", n)
    return n


def extract_sample_ids(vcf_path: str | Path) -> List[str]:
    """Extract sample IDs from a VCF header using ``bcftools query -l``.

    Args:
        vcf_path: Path to VCF or VCF.gz.

    Returns:
        List of sample ID strings.
    """
    result = subprocess.run(
        ["bcftools", "query", "-l", str(vcf_path)],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip().splitlines()
