"""Synthetic GWAS dataset generator with statistical calibration.

Produces a realistic synthetic VCF + phenotype CSV for pipeline validation
and benchmarking.  Supports NULL, LOW, and HIGH calibration modes that
control the expected number of genome-wide significant variants.

Extracted from ``99_create_synthetic_data.py`` to make the generator
importable by tests and notebooks.
"""

from __future__ import annotations

import csv
import math
import random
from pathlib import Path
from typing import Any, Dict, List, Optional

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


# ── Pure-Python statistical helpers (no scipy) ───────────────────────────


def normal_cdf(x: float) -> float:
    """Standard normal CDF via Abramowitz & Stegun."""
    sign = 1 if x >= 0 else -1
    x_abs = abs(x) / math.sqrt(2)
    t = 1.0 / (1.0 + 0.3275911 * x_abs)
    y = 1.0 - (
        (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t
        + 0.254829592
    ) * t * math.exp(-x_abs * x_abs)
    return 0.5 * (1.0 + sign * y)


def inv_normal_cdf(p: float) -> float:
    """Inverse normal CDF (probit) via rational approximation."""
    if p <= 0 or p >= 1:
        raise ValueError(f"p must be in (0,1), got {p}")
    if p < 0.5:
        return -inv_normal_cdf(1 - p)
    t = math.sqrt(-2 * math.log(1 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t)


def calculate_power(
    beta: float,
    maf: float,
    n: int,
    noise_sd: float = 1.0,
    alpha: float = 5e-8,
) -> float:
    """Analytical statistical power for a linear association test.

    Power = Φ(|β| × √(2 × n × MAF × (1−MAF)) / σ − z_{1−α/2})

    Args:
        beta: Effect size.
        maf: Minor allele frequency.
        n: Sample size.
        noise_sd: Phenotype noise standard deviation.
        alpha: Significance threshold.

    Returns:
        Power in [0, 1].
    """
    try:
        z_alpha = inv_normal_cdf(1 - alpha / 2)
        var_g = 2 * maf * (1 - maf)
        ncp = abs(beta) * math.sqrt(n * var_g) / noise_sd
        return 1.0 - normal_cdf(z_alpha - ncp)
    except (ValueError, OverflowError):
        return 0.0


# ── Apis mellifera karyotype constants ───────────────────────────────────

APIS_CHROMOSOMES = [
    ("NC_037638.1", 27_754_200, "CM009931.2"),
    ("NC_037639.1", 16_089_512, "CM009932.2"),
    ("NC_037640.1", 13_619_445, "CM009933.2"),
    ("NC_037641.1", 13_404_451, "CM009934.2"),
    ("NC_037642.1", 13_896_941, "CM009935.2"),
    ("NC_037643.1", 17_789_102, "CM009936.2"),
    ("NC_037644.1", 14_198_698, "CM009937.2"),
    ("NC_037645.1", 12_717_210, "CM009938.2"),
    ("NC_037646.1", 12_354_651, "CM009939.2"),
    ("NC_037647.1", 12_360_052, "CM009940.2"),
    ("NC_037648.1", 16_352_600, "CM009941.2"),
    ("NC_037649.1", 11_514_234, "CM009942.2"),
    ("NC_037650.1", 11_279_722, "CM009943.2"),
    ("NC_037651.1", 10_670_842, "CM009944.2"),
    ("NC_037652.1", 9_534_514, "CM009945.2"),
    ("NC_037653.1", 7_238_532, "CM009946.2"),
]


# ── Known Apis mellifera gene coordinates for realistic causal sites ────

APIS_GENES = [
    ("NC_037638.1", 1_200_000, "LOC551115"),
    ("NC_037638.1", 3_500_000, "LOC726601"),
    ("NC_037639.1", 2_800_000, "LOC552286"),
    ("NC_037640.1", 5_100_000, "LOC411577"),
    ("NC_037641.1", 7_200_000, "LOC726244"),
    ("NC_037643.1", 4_300_000, "LOC412025"),
    ("NC_037644.1", 8_600_000, "LOC725056"),
    ("NC_037648.1", 6_100_000, "LOC551840"),
]


# ── Calibration modes ───────────────────────────────────────────────────

CALIBRATION_MODES = {
    "NULL": {
        "target_gw_sig_pct": 0.0,
        "h2": 0.0,
        "label": "Pure null (no true effects)",
    },
    "LOW": {
        "target_gw_sig_pct": 0.1,
        "h2": 0.05,
        "label": "Low-heritability polygenic",
    },
    "HIGH": {
        "target_gw_sig_pct": 1.0,
        "h2": 0.30,
        "label": "High-heritability with major QTL",
    },
}


def create_synthetic_data(
    output_dir: str | Path,
    *,
    sample_ids: Optional[List[str]] = None,
    n_samples: int = 33,
    n_variants: int = 2000,
    calibration_mode: str = "LOW",
    noise_sd: float = 1.0,
    seed: int = 42,
) -> Dict[str, Any]:
    """Generate a complete synthetic GWAS dataset.

    Produces:
      - ``synthetic.vcf`` with realistic Apis karyotype placement
      - ``phenotype.csv`` with sample_id and trait values
      - ``synthetic_params.json`` documenting all generation parameters

    Args:
        output_dir: Directory to write output files.
        sample_ids: Explicit sample IDs (default: C/I/M/R pattern).
        n_samples: Number of samples (only used if sample_ids is None).
        n_variants: Number of variants to generate.
        calibration_mode: One of NULL, LOW, HIGH.
        noise_sd: Phenotype noise standard deviation.
        seed: Random seed for reproducibility.

    Returns:
        Dict with paths to generated files and calibration metadata.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    random.seed(seed)

    mode = CALIBRATION_MODES.get(calibration_mode, CALIBRATION_MODES["LOW"])
    h2 = mode["h2"]

    # Generate sample IDs if not provided
    if sample_ids is None:
        strains = ["C", "I", "M", "R"]
        sample_ids = []
        per_strain = n_samples // len(strains)
        for s in strains:
            for i in range(1, per_strain + 1):
                for caste in ["ITQ", "ITW", "WORK", "G"]:
                    if len(sample_ids) < n_samples:
                        sample_ids.append(f"{s}{i}{caste}")
        while len(sample_ids) < n_samples:
            sample_ids.append(f"C{len(sample_ids)}WORK")
    else:
        n_samples = len(sample_ids)

    # Assign causal variant positions near known genes
    n_causal = max(1, int(n_variants * mode["target_gw_sig_pct"] / 100))
    causal_indices = set(random.sample(range(n_variants), min(n_causal, n_variants)))

    # Generate variant coordinates across karyotype
    total_genome = sum(length for _, length, _ in APIS_CHROMOSOMES)
    variants = []
    for _ in range(n_variants):
        pos_in_genome = random.randint(0, total_genome - 1)
        cumulative = 0
        for chrom, length, _ in APIS_CHROMOSOMES:
            cumulative += length
            if pos_in_genome < cumulative:
                chrom_pos = pos_in_genome - (cumulative - length)
                variants.append((chrom, chrom_pos + 1))
                break
    variants.sort()

    # Generate genotypes and effect sizes
    genotype_matrix = []  # n_variants × n_samples
    effect_sizes = []
    mafs = []

    for vi in range(n_variants):
        maf = random.uniform(0.05, 0.49)
        mafs.append(maf)

        # Generate genotypes under HWE
        genos = []
        for _ in range(n_samples):
            a1 = 1 if random.random() < maf else 0
            a2 = 1 if random.random() < maf else 0
            genos.append(a1 + a2)
        genotype_matrix.append(genos)

        if vi in causal_indices and h2 > 0:
            # Calibrate effect size to achieve detectable power
            var_g = 2 * maf * (1 - maf)
            if var_g > 0:
                beta = math.sqrt(h2 / (n_causal * var_g)) * noise_sd
                if random.random() < 0.5:
                    beta = -beta
            else:
                beta = 0.0
        else:
            beta = 0.0
        effect_sizes.append(beta)

    # Generate phenotypes: y = Σ(β_i × g_i) + noise
    phenotypes = []
    for si in range(n_samples):
        genetic_value = sum(
            effect_sizes[vi] * genotype_matrix[vi][si] for vi in range(n_variants)
        )
        noise = random.gauss(0, noise_sd)
        phenotypes.append(genetic_value + noise)

    # Write VCF
    vcf_path = output_dir / "synthetic.vcf"
    _write_vcf(vcf_path, variants, genotype_matrix, sample_ids)

    # Write phenotype CSV
    phenotype_path = output_dir / "phenotype.csv"
    with open(phenotype_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample_id", "trait_value"])
        for sid, pheno in zip(sample_ids, phenotypes):
            writer.writerow([sid, f"{pheno:.6f}"])

    # Compute power for causal variants
    causal_powers = []
    for vi in causal_indices:
        pwr = calculate_power(effect_sizes[vi], mafs[vi], n_samples, noise_sd)
        causal_powers.append(pwr)
    mean_power = sum(causal_powers) / len(causal_powers) if causal_powers else 0.0

    import json

    params = {
        "calibration_mode": calibration_mode,
        "mode_label": mode["label"],
        "h2": h2,
        "n_samples": n_samples,
        "n_variants": n_variants,
        "n_causal": n_causal,
        "noise_sd": noise_sd,
        "seed": seed,
        "mean_causal_power": round(mean_power, 4),
        "vcf_path": str(vcf_path),
        "phenotype_path": str(phenotype_path),
    }
    with open(output_dir / "synthetic_params.json", "w") as f:
        json.dump(params, f, indent=2)

    logger.info(
        "Generated synthetic dataset: %d samples × %d variants, "
        "%d causal (mode=%s, h²=%.2f, mean_power=%.3f)",
        n_samples,
        n_variants,
        n_causal,
        calibration_mode,
        h2,
        mean_power,
    )
    return params


def _write_vcf(
    vcf_path: Path,
    variants: List[tuple],
    genotype_matrix: List[List[int]],
    sample_ids: List[str],
) -> None:
    """Write a minimal VCF 4.2 file."""
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sid in sample_ids:
            f.write(f"\t{sid}")
        f.write("\n")

        for vi, (chrom, pos) in enumerate(variants):
            snp_id = f"syn_{vi:06d}"
            genos = genotype_matrix[vi]
            gt_strs = []
            for g in genos:
                if g == 0:
                    gt_strs.append("0/0")
                elif g == 1:
                    gt_strs.append("0/1")
                else:
                    gt_strs.append("1/1")
            f.write(
                f"{chrom}\t{pos}\t{snp_id}\tA\tT\t999\tPASS\t.\tGT\t"
                + "\t".join(gt_strs)
                + "\n"
            )
