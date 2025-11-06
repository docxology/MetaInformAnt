#!/usr/bin/env python3
"""Multi-omics simulation script.

This script generates synthetic multi-omics data including cross-platform
datasets and integrated omics data.

Usage:
    python3 scripts/simulation/simulate_multiomics.py --type cross-platform --n-samples 20
    python3 scripts/simulation/simulate_multiomics.py --type integrated --n-samples 30
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
from metainformant.simulation.rna import simulate_counts_negative_binomial
from metainformant.simulation.popgen import generate_genotype_matrix
from metainformant.simulation.sequences import generate_random_protein

logger = logging.get_logger(__name__)


def simulate_crossplatform(
    output_dir: Path,
    n_samples: int,
    n_genes: int,
    n_variants: int,
    n_proteins: int,
    seed: int,
) -> dict:
    """Simulate cross-platform multi-omics data."""
    logger.info(f"Generating cross-platform data: {n_samples} samples")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    
    # Generate genomics (genotypes)
    logger.info("Generating genomics data...")
    genotypes = generate_genotype_matrix(n_samples, n_variants, rng=rng)
    genomics_df = pd.DataFrame(genotypes, columns=[f"variant_{i:05d}" for i in range(n_variants)])
    genomics_df.insert(0, "sample_id", sample_ids)
    genomics_file = output_dir / "genomics.csv"
    genomics_df.to_csv(genomics_file, index=False)
    
    # Generate transcriptomics (expression)
    logger.info("Generating transcriptomics data...")
    expression = simulate_counts_negative_binomial(n_genes, n_samples, mean_expression=100.0, dispersion=0.1, rng=rng)
    transcriptomics_df = pd.DataFrame(expression, columns=sample_ids)
    transcriptomics_df.index = [f"gene_{i:05d}" for i in range(n_genes)]
    transcriptomics_file = output_dir / "transcriptomics.csv"
    transcriptomics_df.to_csv(transcriptomics_file)
    
    # Generate proteomics (protein abundance)
    logger.info("Generating proteomics data...")
    protein_abundance = np.random.lognormal(mean=5.0, sigma=1.0, size=(n_proteins, n_samples))
    proteomics_df = pd.DataFrame(protein_abundance, columns=sample_ids)
    proteomics_df.index = [f"protein_{i:05d}" for i in range(n_proteins)]
    proteomics_file = output_dir / "proteomics.csv"
    proteomics_df.to_csv(proteomics_file)
    
    # Create sample metadata
    metadata = pd.DataFrame({
        "sample_id": sample_ids,
        "batch": [f"batch_{i % 3}" for i in range(n_samples)],
        "condition": [f"condition_{i % 2}" for i in range(n_samples)],
    })
    metadata_file = output_dir / "sample_metadata.csv"
    metadata.to_csv(metadata_file, index=False)
    
    logger.info(f"Cross-platform data saved to {output_dir}")
    
    return {
        "type": "cross-platform",
        "n_samples": n_samples,
        "genomics_file": str(genomics_file),
        "transcriptomics_file": str(transcriptomics_file),
        "proteomics_file": str(proteomics_file),
        "metadata_file": str(metadata_file),
    }


def simulate_integrated(
    output_dir: Path,
    n_samples: int,
    n_features_per_omic: int,
    n_omics: int,
    seed: int,
) -> dict:
    """Simulate integrated multi-omics dataset."""
    logger.info(f"Generating integrated data: {n_samples} samples, {n_omics} omics types")
    rng = random.Random(seed)
    np.random.seed(seed)
    
    sample_ids = [f"sample_{i:03d}" for i in range(n_samples)]
    omics_types = ["genomics", "transcriptomics", "proteomics", "metabolomics", "epigenomics"]
    
    integrated_data = {}
    
    for omic_idx in range(min(n_omics, len(omics_types))):
        omic_type = omics_types[omic_idx]
        logger.info(f"Generating {omic_type} data...")
        
        # Generate omic-specific data
        if omic_type == "genomics":
            data = generate_genotype_matrix(n_samples, n_features_per_omic, rng=rng)
        elif omic_type == "transcriptomics":
            data = simulate_counts_negative_binomial(
                n_features_per_omic, n_samples, mean_expression=100.0, dispersion=0.1, rng=rng
            )
        else:
            # Other omics: log-normal distribution
            data = np.random.lognormal(mean=5.0, sigma=1.0, size=(n_features_per_omic, n_samples))
            data = data.tolist()
        
        # Create DataFrame
        feature_names = [f"{omic_type}_feature_{i:05d}" for i in range(n_features_per_omic)]
        if omic_type == "genomics":
            df = pd.DataFrame(data, columns=feature_names)
            df.insert(0, "sample_id", sample_ids)
        else:
            df = pd.DataFrame(data, columns=sample_ids)
            df.index = feature_names
        
        integrated_data[omic_type] = df
        
        # Save individual omic
        omic_file = output_dir / f"{omic_type}.csv"
        df.to_csv(omic_file)
    
    # Create integrated matrix (samples x all features)
    all_features = []
    for omic_type, df in integrated_data.items():
        if "sample_id" in df.columns:
            # Transpose genomics
            df_t = df.set_index("sample_id").T
            all_features.append(df_t)
        else:
            # Transcriptomics, etc. are already features x samples
            all_features.append(df.T)
    
    integrated_df = pd.concat(all_features, axis=1)
    integrated_file = output_dir / "integrated_data.csv"
    integrated_df.to_csv(integrated_file)
    
    logger.info(f"Integrated data saved to {integrated_file}")
    
    return {
        "type": "integrated",
        "n_samples": n_samples,
        "n_omics": n_omics,
        "n_features_per_omic": n_features_per_omic,
        "output_file": str(integrated_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Multi-omics simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate cross-platform data
  %(prog)s --type cross-platform --n-samples 20 --n-genes 1000 --n-variants 500

  # Simulate integrated data
  %(prog)s --type integrated --n-samples 30 --n-features 200 --n-omics 3
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["cross-platform", "integrated"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/multiomics"),
        help="Output directory (default: output/simulation/multiomics)",
    )
    parser.add_argument("--n-samples", type=int, default=20, help="Number of samples")
    parser.add_argument("--n-genes", type=int, default=1000, help="Number of genes (cross-platform type)")
    parser.add_argument("--n-variants", type=int, default=500, help="Number of variants (cross-platform type)")
    parser.add_argument("--n-proteins", type=int, default=500, help="Number of proteins (cross-platform type)")
    parser.add_argument("--n-features", type=int, default=200, help="Features per omic (integrated type)")
    parser.add_argument("--n-omics", type=int, default=3, help="Number of omics types (integrated type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "cross-platform":
            results = simulate_crossplatform(
                output_dir,
                args.n_samples,
                args.n_genes,
                args.n_variants,
                args.n_proteins,
                args.seed,
            )
        elif args.type == "integrated":
            results = simulate_integrated(
                output_dir, args.n_samples, args.n_features, args.n_omics, args.seed
            )
        
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

