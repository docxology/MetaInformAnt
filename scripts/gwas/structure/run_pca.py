#!/usr/bin/env python3
"""Run PCA for population structure analysis.

Calculates Principal Components from genomic variants.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core import logging
from metainformant.gwas import compute_pca, parse_vcf_full

# logging.setup_logging()
logger = logging.get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run PCA")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--n-components", type=int, default=10, help="Number of PC components")
    parser.add_argument("--output", required=True, help="Output JSON results file")
    parser.add_argument("--plot", help="Output plot file (optional)")

    args = parser.parse_args()

    logger.info(f"Parsing VCF: {args.vcf}")
    try:
        vcf_data = parse_vcf_full(str(args.vcf))
        genotypes = vcf_data.get("genotypes")

        if not genotypes:
            logger.error("No genotypes found in VCF")
            sys.exit(1)

        logger.info(f"Computing PCA with {args.n_components} components on {len(genotypes)} samples...")
        result = compute_pca(genotypes, n_components=args.n_components)

        if result["status"] == "success":
            logger.info("PCA completed successfully")

            # Convert numpy arrays to lists for JSON
            # (Result usually has lists, but verifying)
            if "pcs" in result and isinstance(result["pcs"], np.ndarray):
                result["pcs"] = result["pcs"].tolist()
            if "explained_variance_ratio" in result and isinstance(result["explained_variance_ratio"], np.ndarray):
                result["explained_variance_ratio"] = result["explained_variance_ratio"].tolist()

            with open(args.output, "w") as f:
                json.dump(result, f, indent=2)
            logger.info(f"Results saved to: {args.output}")

            if args.plot:
                # Optional: Import plotting logic
                # For now, just logging that plotting is separate or use visualization module
                logger.info("Plotting not fully integrated in CLI yet, use visualization scripts.")
        else:
            logger.error(f"PCA failed: {result.get('error')}")
            sys.exit(1)

    except Exception as e:
        logger.error(f"Error executing PCA: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
