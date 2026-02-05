#!/usr/bin/env python3
"""GWAS simulation script.

This script generates synthetic GWAS data including variant VCF files,
genotype matrices, phenotypes, and population structure data.

Usage:
    python3 scripts/simulation/simulate_gwas.py --type variants --n-variants 1000 --n-samples 100
    python3 scripts/simulation/simulate_gwas.py --type genotypes --n-samples 200 --n-sites 500
    python3 scripts/simulation/simulate_gwas.py --type phenotype --n-samples 150
    python3 scripts/simulation/simulate_gwas.py --type population --n-samples 100 --n-populations 2
"""

import argparse
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation
from metainformant.simulation.popgen import generate_genotype_matrix

logger = logging.get_logger(__name__)


def simulate_variants(
    output_dir: Path,
    n_variants: int,
    n_samples: int,
    chromosome: str,
    seed: int,
) -> dict:
    """Simulate variant VCF file.

    Args:
        output_dir: Output directory for results
        n_variants: Number of variants to generate
        n_samples: Number of samples
        chromosome: Chromosome name
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_variants, min_val=1, name="n_variants")
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_type(chromosome, str, "chromosome")

    logger.info(f"Generating variants: {n_variants} variants, {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate VCF header
    vcf_lines = [
        "##fileformat=VCFv4.3",
        f"##contig=<ID={chromosome},length=1000000>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join([f"sample_{i:03d}" for i in range(n_samples)]),
    ]

    # Generate variants
    positions = sorted(rng.sample(range(1, 1000000), n_variants))
    bases = ["A", "C", "G", "T"]

    for pos_idx, pos in enumerate(positions):
        ref = rng.choice(bases)
        alt = rng.choice([b for b in bases if b != ref])

        # Generate genotypes for each sample
        genotypes = []
        for _ in range(n_samples):
            # Random genotype: 0/0, 0/1, or 1/1
            gt = rng.choice(["0/0", "0/1", "1/1"], p=[0.5, 0.3, 0.2])
            dp = rng.randint(10, 50)
            genotypes.append(f"{gt}:{dp}")

        vcf_line = f"{chromosome}\t{pos}\tvariant_{pos_idx:05d}\t{ref}\t{alt}\t.\tPASS\t.\tGT:DP\t" + "\t".join(
            genotypes
        )
        vcf_lines.append(vcf_line)

    # Save VCF
    vcf_file = output_dir / "variants.vcf"
    with open(vcf_file, "w") as f:
        f.write("\n".join(vcf_lines))

    logger.info(f"VCF file saved to {vcf_file}")

    return {
        "type": "variants",
        "n_variants": n_variants,
        "n_samples": n_samples,
        "chromosome": chromosome,
        "output_file": str(vcf_file),
    }


def simulate_genotypes(
    output_dir: Path,
    n_samples: int,
    n_sites: int,
    seed: int,
) -> dict:
    """Simulate genotype matrix.

    Args:
        output_dir: Output directory for results
        n_samples: Number of samples
        n_sites: Number of variant sites
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_range(n_sites, min_val=1, name="n_sites")

    logger.info(f"Generating genotype matrix: {n_samples} samples, {n_sites} sites")
    rng = random.Random(seed)

    # Generate genotype matrix
    genotypes = generate_genotype_matrix(n_samples, n_sites, rng=rng)

    # Create DataFrame
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    site_ids = [f"site_{i:05d}" for i in range(n_sites)]
    df = pd.DataFrame(genotypes, columns=site_ids, index=sample_ids)

    csv_file = output_dir / "genotype_matrix.csv"
    df.to_csv(csv_file)

    logger.info(f"Genotype matrix saved to {csv_file}")

    return {
        "type": "genotypes",
        "n_samples": n_samples,
        "n_sites": n_sites,
        "output_file": str(csv_file),
    }


def simulate_phenotype(
    output_dir: Path,
    n_samples: int,
    n_causal_variants: int,
    heritability: float,
    seed: int,
) -> dict:
    """Simulate phenotype with genetic effects.

    Args:
        output_dir: Output directory for results
        n_samples: Number of samples
        n_causal_variants: Number of causal variants
        heritability: Heritability (0.0-1.0)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_range(n_causal_variants, min_val=1, name="n_causal_variants")
    validation.validate_range(heritability, min_val=0.0, max_val=1.0, name="heritability")

    logger.info(f"Generating phenotype: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Generate genotypes for causal variants
    genotypes = generate_genotype_matrix(n_samples, n_causal_variants, rng=rng)

    # Generate effect sizes
    effect_sizes = np.random.normal(0, 1, n_causal_variants)

    # Calculate genetic component
    genetic_component = np.array(genotypes) @ effect_sizes

    # Scale to achieve heritability
    genetic_var = np.var(genetic_component)
    if genetic_var > 0:
        genetic_component = genetic_component * np.sqrt(heritability / genetic_var)

    # Add environmental noise
    environmental_component = np.random.normal(0, np.sqrt(1 - heritability), n_samples)

    # Phenotype
    phenotype = genetic_component + environmental_component

    # Create DataFrame
    df = pd.DataFrame(
        {
            "sample_id": [f"sample_{i:03d}" for i in range(n_samples)],
            "phenotype": phenotype,
        }
    )

    csv_file = output_dir / "phenotype.csv"
    df.to_csv(csv_file, index=False)

    # Save effect sizes
    effects_df = pd.DataFrame(
        {
            "variant_id": [f"variant_{i:05d}" for i in range(n_causal_variants)],
            "effect_size": effect_sizes,
        }
    )
    effects_file = output_dir / "causal_effects.csv"
    effects_df.to_csv(effects_file, index=False)

    logger.info(f"Phenotype data saved to {csv_file}")

    return {
        "type": "phenotype",
        "n_samples": n_samples,
        "n_causal_variants": n_causal_variants,
        "heritability": heritability,
        "output_file": str(csv_file),
        "effects_file": str(effects_file),
    }


def simulate_population(
    output_dir: Path,
    n_samples: int,
    n_populations: int,
    fst: float,
    seed: int,
) -> dict:
    """Simulate population structure.

    Args:
        output_dir: Output directory for results
        n_samples: Number of samples
        n_populations: Number of populations
        fst: Fst value (0.0-1.0)
        seed: Random seed for reproducibility

    Returns:
        Dictionary with simulation results and metadata

    Raises:
        ValidationError: If parameters are invalid
    """
    # Validate parameters
    validation.validate_range(n_samples, min_val=1, name="n_samples")
    validation.validate_range(n_populations, min_val=1, max_val=n_samples, name="n_populations")
    validation.validate_range(fst, min_val=0.0, max_val=1.0, name="fst")

    logger.info(f"Generating population structure: {n_samples} samples, {n_populations} populations")
    rng = random.Random(seed)
    np.random.seed(seed)

    # Assign samples to populations
    samples_per_pop = n_samples // n_populations
    assignments = []

    for pop_idx in range(n_populations):
        for _ in range(samples_per_pop):
            assignments.append(pop_idx)

    # Add remaining samples
    for _ in range(n_samples - len(assignments)):
        assignments.append(rng.randint(0, n_populations - 1))

    # Generate PCA coordinates (simplified)
    pca_coords = []
    for pop_idx in assignments:
        # Each population has distinct coordinates
        x = pop_idx * 2.0 + np.random.normal(0, 0.5)
        y = np.random.normal(0, 0.5)
        pca_coords.append([x, y])

    # Create DataFrame
    df = pd.DataFrame(
        {
            "sample_id": [f"sample_{i:03d}" for i in range(n_samples)],
            "population": [f"pop_{a:02d}" for a in assignments],
            "pca1": [c[0] for c in pca_coords],
            "pca2": [c[1] for c in pca_coords],
        }
    )

    csv_file = output_dir / "population_structure.csv"
    df.to_csv(csv_file, index=False)

    logger.info(f"Population structure saved to {csv_file}")

    return {
        "type": "population",
        "n_samples": n_samples,
        "n_populations": n_populations,
        "fst": fst,
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="GWAS simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate variants (VCF)
  %(prog)s --type variants --n-variants 1000 --n-samples 100

  # Simulate genotype matrix
  %(prog)s --type genotypes --n-samples 200 --n-sites 500

  # Simulate phenotype
  %(prog)s --type phenotype --n-samples 150 --n-causal 10 --heritability 0.5

  # Simulate population structure
  %(prog)s --type population --n-samples 100 --n-populations 2
        """,
    )

    parser.add_argument(
        "--type",
        required=True,
        choices=["variants", "genotypes", "phenotype", "population"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/gwas"),
        help="Output directory (default: output/simulation/gwas)",
    )
    parser.add_argument("--n-variants", type=int, default=1000, help="Number of variants (variants type)")
    parser.add_argument("--n-samples", type=int, default=100, help="Number of samples")
    parser.add_argument("--n-sites", type=int, default=500, help="Number of sites (genotypes type)")
    parser.add_argument("--chromosome", type=str, default="chr1", help="Chromosome name (variants type)")
    parser.add_argument("--n-causal", type=int, default=10, help="Number of causal variants (phenotype type)")
    parser.add_argument("--heritability", type=float, default=0.5, help="Heritability (phenotype type, 0-1)")
    parser.add_argument("--n-populations", type=int, default=2, help="Number of populations (population type)")
    parser.add_argument("--fst", type=float, default=0.1, help="Fst value (population type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")

    args = parser.parse_args()

    if args.verbose:
        import logging as std_logging

        logger.setLevel(std_logging.DEBUG)

    # Validate output directory
    output_dir = paths.ensure_directory(args.output)

    # Validate common parameters
    validation.validate_range(args.n_samples, min_val=1, name="n_samples")
    if hasattr(args, "n_variants"):
        validation.validate_range(args.n_variants, min_val=1, name="n_variants")
    if hasattr(args, "n_sites"):
        validation.validate_range(args.n_sites, min_val=1, name="n_sites")
    if hasattr(args, "n_causal"):
        validation.validate_range(args.n_causal, min_val=1, name="n_causal")
    if hasattr(args, "heritability"):
        validation.validate_range(args.heritability, min_val=0.0, max_val=1.0, name="heritability")
    if hasattr(args, "n_populations"):
        validation.validate_range(args.n_populations, min_val=1, max_val=args.n_samples, name="n_populations")
    if hasattr(args, "fst"):
        validation.validate_range(args.fst, min_val=0.0, max_val=1.0, name="fst")

    try:
        if args.type == "variants":
            results = simulate_variants(output_dir, args.n_variants, args.n_samples, args.chromosome, args.seed)
        elif args.type == "genotypes":
            results = simulate_genotypes(output_dir, args.n_samples, args.n_sites, args.seed)
        elif args.type == "phenotype":
            results = simulate_phenotype(output_dir, args.n_samples, args.n_causal, args.heritability, args.seed)
        elif args.type == "population":
            results = simulate_population(output_dir, args.n_samples, args.n_populations, args.fst, args.seed)

        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")

        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
