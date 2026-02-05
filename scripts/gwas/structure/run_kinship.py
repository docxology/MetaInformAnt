#!/usr/bin/env python3
"""Run Kinship Matrix computation.

Calculates the kinship matrix (relatedness) between samples.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core import logging
from metainformant.gwas import compute_kinship_matrix, parse_vcf_full

# logging.setup_logging()
logger = logging.get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run Kinship Matrix Computation")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--method", default="vanraden", choices=["vanraden", "ibs"], help="Kinship method")
    parser.add_argument("--max-variants", type=int, default=5000, help="Max variants to use (subsampling)")
    parser.add_argument("--output", required=True, help="Output JSON results file")

    args = parser.parse_args()

    logger.info(f"Parsing VCF: {args.vcf}")
    try:
        vcf_data = parse_vcf_full(str(args.vcf))
        genotypes = vcf_data.get("genotypes")

        if not genotypes:
            logger.error("No genotypes found in VCF")
            sys.exit(1)

        n_variants = len(genotypes[0])
        logger.info(f"Total variants available: {n_variants}")

        # Subsampling optimization
        if n_variants > args.max_variants:
            step = n_variants // args.max_variants
            indices = list(range(0, n_variants, step))[: args.max_variants]
            logger.info(f"Subsampling to {len(indices)} variants for efficiency...")
            genotypes = [[row[i] for i in indices] for row in genotypes]

        logger.info(f"Computing Kinship Matrix using {args.method}...")
        result = compute_kinship_matrix(genotypes, method=args.method)

        if result["status"] == "success":
            logger.info("Kinship matrix computed successfully")

            if "kinship_matrix" in result and isinstance(result["kinship_matrix"], np.ndarray):
                result["kinship_matrix"] = result["kinship_matrix"].tolist()

            with open(args.output, "w") as f:
                json.dump(result, f, indent=2)
            logger.info(f"Results saved to: {args.output}")
        else:
            logger.error(f"Kinship computation failed: {result.get('error')}")
            sys.exit(1)

    except Exception as e:
        logger.error(f"Error executing Kinship: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
