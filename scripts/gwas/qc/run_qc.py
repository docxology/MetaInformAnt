#!/usr/bin/env python3
"""Run Quality Control on VCF files.

This script applies standard GWAS quality control filters:
- Minor Allele Frequency (MAF)
- Missingness rate
- Hardy-Weinberg Equilibrium (HWE)
- Quality score
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core import logging
from metainformant.gwas import apply_qc_filters

# logging.setup_logging()
logger = logging.get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run GWAS Quality Control")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--output", help="Output filtered VCF file (optional)")
    parser.add_argument("--min-maf", type=float, default=0.05, help="Minimum Minor Allele Frequency")
    parser.add_argument("--max-missing", type=float, default=0.10, help="Maximum missingness rate")
    parser.add_argument("--hwe-pval", type=float, default=1e-6, help="HWE p-value threshold")
    parser.add_argument("--min-qual", type=float, default=30.0, help="Minimum quality score")
    parser.add_argument("--report", help="Output QC report JSON file")

    args = parser.parse_args()

    qc_config = {
        "min_maf": args.min_maf,
        "max_missing": args.max_missing,
        "hwe_pval": args.hwe_pval,
        "min_qual": args.min_qual,
    }

    logger.info(f"Running QC on {args.vcf}")
    logger.info(f"Parameters: {qc_config}")

    result = apply_qc_filters(vcf_path=args.vcf, config=qc_config, output_vcf=args.output)

    if result["status"] == "success":
        logger.info("QC completed successfully")
        logger.info(f"Variants passing: {result.get('num_variants_after', 0)}")
        if args.output:
            logger.info(f"Filtered VCF saved to: {args.output}")

        if args.report:
            with open(args.report, "w") as f:
                json.dump(result, f, indent=2)
            logger.info(f"Report saved to: {args.report}")
    else:
        logger.error(f"QC failed: {result.get('error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()
