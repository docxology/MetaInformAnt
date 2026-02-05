#!/usr/bin/env python3
"""Run GWAS Association Tests.

Performs linear or logistic regression for phenotype-genotype association.
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

from metainformant.core import logging
from metainformant.gwas import parse_vcf_full
from metainformant.gwas.association import association_test_linear, association_test_logistic
from metainformant.gwas.correction import bonferroni_correction, fdr_correction

# logging.setup_logging()
logger = logging.get_logger(__name__)


def main():
    parser = argparse.ArgumentParser(description="Run GWAS Association Tests")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--phenotypes", required=True, help="Phenotype TSV file")
    parser.add_argument("--trait", required=True, help="Trait column name to test")
    parser.add_argument("--test-type", choices=["linear", "logistic", "auto"], default="auto", help="Test type")
    parser.add_argument("--output", required=True, help="Output results TSV file")
    parser.add_argument("--max-variants", type=int, help="Limit number of variants (for testing)")

    args = parser.parse_args()

    # Load Phenotypes
    logger.info(f"Loading phenotypes from {args.phenotypes}")
    pheno_df = pd.read_csv(args.phenotypes, sep="\t")

    if args.trait not in pheno_df.columns:
        logger.error(f"Trait '{args.trait}' not found in phenotype file")
        sys.exit(1)

    phenotype_values = pheno_df[args.trait].tolist()

    # Determine test method
    test_type = args.test_type
    if test_type == "auto":
        unique_vals = len(set(phenotype_values))
        if unique_vals <= 2:
            test_type = "logistic"
        else:
            test_type = "linear"
    logger.info(f"Using test type: {test_type}")

    # Parse VCF
    logger.info(f"Parsing VCF: {args.vcf}")
    vcf_data = parse_vcf_full(str(args.vcf))
    genotypes = vcf_data.get("genotypes")
    variants = vcf_data.get("variants")

    if not genotypes:
        logger.error("No genotypes found")
        sys.exit(1)

    # Limit variants if requested
    if args.max_variants and args.max_variants < len(variants):
        logger.info(f"Limiting to first {args.max_variants} variants")
        variants = variants[: args.max_variants]
        # Transpose constraint implicitly assumed by loop order?
        # Genotypes usually [sample][variant].
        # So we just iterate up to max_variants index.

    results = []

    logger.info(f"Testing {len(variants)} variants...")
    for idx, variant in enumerate(variants):
        # Extract variant genotype vector
        variant_geno = [row[idx] for row in genotypes]

        try:
            if test_type == "linear":
                res = association_test_linear(variant_geno, phenotype_values)
            else:
                # Logistic requires int/float 0/1 usually
                res = association_test_logistic(variant_geno, [float(x) for x in phenotype_values])

            if res["status"] == "success":
                results.append(
                    {
                        "variant_id": variant.get("ID", f"var_{idx}"),
                        "chrom": variant.get("CHROM"),
                        "pos": variant.get("POS"),
                        "beta": res.get("beta"),
                        "se": res.get("se"),
                        "p_value": res.get("p_value"),
                    }
                )
        except Exception as e:
            # logger.debug(f"Test failed for variant {idx}: {e}")
            pass

        if (idx + 1) % 1000 == 0:
            logger.info(f"Processed {idx+1}/{len(variants)}")

    if not results:
        logger.warning("No successful association tests.")
        sys.exit(0)

    # Correction
    p_values = [r["p_value"] for r in results]
    bonf = bonferroni_correction(p_values)
    fdr = fdr_correction(p_values)

    # Annotate results
    bonf_sig = set(bonf.get("significant_indices", []))
    fdr_sig = set(fdr.get("significant_indices", []))

    for i, r in enumerate(results):
        r["bonferroni_sig"] = i in bonf_sig
        r["fdr_sig"] = i in fdr_sig

    # Save
    out_df = pd.DataFrame(results)
    out_df.to_csv(args.output, sep="\t", index=False)
    logger.info(f"Saved {len(results)} associations to {args.output}")


if __name__ == "__main__":
    main()
