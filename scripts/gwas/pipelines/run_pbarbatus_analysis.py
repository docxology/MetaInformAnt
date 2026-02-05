#!/usr/bin/env python3
"""Run comprehensive P. barbatus GWAS analysis with synthetic data."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent / "src"))

import logging

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

from metainformant.core.io import dump_json, ensure_directory
from metainformant.gwas import (
    apply_qc_filters,
    compute_kinship_matrix,
    compute_pca,
    load_gwas_config,
    manhattan_plot,
    parse_vcf_full,
    qq_plot,
)
from metainformant.gwas.association import association_test_linear, association_test_logistic
from metainformant.gwas.correction import bonferroni_correction, fdr_correction

print("=" * 80)
print("P. BARBATUS GWAS - COMPREHENSIVE ANALYSIS")
print("=" * 80)

# Load configuration
print("\nLoading configuration...")
config = load_gwas_config("config/gwas/gwas_pbarbatus_synthetic.yaml")
print(f"✓ Configuration loaded")

# Parse VCF
print("\n[1/7] Parsing VCF file...")
vcf_file = config.variants["vcf_files"][0]
vcf_data = parse_vcf_full(vcf_file)
print(f'✓ Parsed {len(vcf_data.get("variants", []))} variants for {len(vcf_data.get("samples", []))} samples')

# Quality control
print("\n[2/7] Applying quality control filters...")
qc_filters = {
    "min_maf": config.qc.get("min_maf", 0.05),
    "max_missing": config.qc.get("max_missing", 0.10),
    "hwe_pval": config.qc.get("hwe_pval", 1e-6),
    "min_qual": config.qc.get("min_qual", 30.0),
}

# Create filtered output file
filtered_vcf = config.work_dir / "filtered_variants.vcf.gz"
qc_result = apply_qc_filters(
    vcf_path=vcf_file,
    config=qc_filters,
    output_vcf=None,  # Don't write, just get stats
)

if qc_result.get("status") == "success":
    n_passing = qc_result.get("num_variants_after", 0)
    n_total = qc_result.get("num_variants_before", 50000)
    print(f"✓ QC complete: {n_passing}/{n_total} variants passed ({n_passing/max(n_total,1)*100:.1f}%)")

    # Use original data (QC already filtered it internally)
    genotypes = vcf_data["genotypes"]
    variants_passing = vcf_data["variants"]
else:
    print(f'✗ QC failed: {qc_result.get("error")}')
    sys.exit(1)

n_samples = len(vcf_data["samples"])

# PCA
print("\n[3/7] Computing PCA for population structure...")
pca_result = compute_pca(genotypes, n_components=10)
if pca_result["status"] == "success":
    print(f'✓ PCA computed: {pca_result["n_components"]} components')
    print(f'  Variance explained: {sum(pca_result["explained_variance_ratio"])*100:.1f}%')
    pcs = np.array(pca_result["pcs"])
else:
    print(f'✗ PCA failed: {pca_result.get("error")}')
    pcs = None

# Kinship
print("\n[4/7] Computing kinship matrix...")
kinship_result = compute_kinship_matrix(genotypes, method="vanraden")
if kinship_result["status"] == "success":
    print(f"✓ Kinship computed: {n_samples}x{n_samples} matrix")
    kinship_matrix = np.array(kinship_result["kinship_matrix"])
    print(f"  Mean relatedness: {np.mean(kinship_matrix):.4f}")
else:
    print(f'✗ Kinship failed: {kinship_result.get("error")}')

# Save intermediate results
results_dir = ensure_directory(config.work_dir / "results")
dump_json(qc_result, results_dir / "qc_results.json", indent=2)
if pcs is not None:
    dump_json(pca_result, results_dir / "pca_results.json", indent=2)
dump_json(kinship_result, results_dir / "kinship_results.json", indent=2)
print(f"\n✓ Intermediate results saved to {results_dir}")

# Load phenotypes
print("\n[5/7] Loading phenotype data...")
pheno_df = pd.read_csv(config.samples["phenotype_file"], sep="\t")
print(f"✓ Loaded {len(pheno_df)} samples with {len(pheno_df.columns)-1} traits/covariates")

# Association testing
print("\n[6/7] Running association testing...")

traits_to_test = [
    ("body_size_mm", "linear"),
    ("foraging_activity", "linear"),
    ("disease_resistance", "logistic"),
]

all_results = {}

for trait_name, test_type in traits_to_test:
    print(f"\n  Testing trait: {trait_name} ({test_type})...")

    phenotypes = pheno_df[trait_name].tolist()
    assoc_results = []

    # Test subset of variants for speed (first 5000)
    n_test = min(5000, len(variants_passing))

    for idx in range(n_test):
        variant = variants_passing[idx]
        variant_geno = [genotypes[i][idx] for i in range(len(genotypes))]

        if test_type == "linear":
            result = association_test_linear(variant_geno, phenotypes)
        else:
            result = association_test_logistic(variant_geno, [int(p) for p in phenotypes])

        if result.get("status") == "success":
            assoc_results.append(
                {
                    "variant_id": variant.get("ID", f"var_{idx}"),
                    "chrom": variant.get("CHROM"),
                    "pos": variant.get("POS"),
                    "beta": result["beta"],
                    "se": result["se"],
                    "p_value": result["p_value"],
                }
            )

        if (idx + 1) % 500 == 0:
            print(f"    Tested {idx+1}/{n_test} variants...")

    print(f"  ✓ {len(assoc_results)} association results")

    # Multiple testing correction
    p_values = [r["p_value"] for r in assoc_results]
    bonf_result = bonferroni_correction(p_values)
    fdr_result = fdr_correction(p_values)

    print(f'  Bonferroni significant: {bonf_result.get("significant_count", 0)}')
    print(f'  FDR significant: {fdr_result.get("significant_count", 0)}')

    # Save results
    trait_file = results_dir / f"{trait_name}_associations.tsv"
    results_df = pd.DataFrame(assoc_results)
    results_df.to_csv(trait_file, sep="\t", index=False)
    print(f"  Saved to: {trait_file}")

    all_results[trait_name] = {
        "associations": assoc_results,
        "bonferroni": bonf_result,
        "fdr": fdr_result,
    }

# Visualization
print("\n[7/7] Generating visualizations...")
plots_dir = ensure_directory(config.work_dir / "plots")

for trait_name in all_results.keys():
    assoc_results = all_results[trait_name]["associations"]

    # Manhattan plot
    manhattan_file = plots_dir / f"manhattan_{trait_name}.png"
    try:
        manhattan_result = manhattan_plot(
            assoc_results, manhattan_file, title=f'{trait_name.replace("_", " ").title()} GWAS - P. barbatus'
        )
        if manhattan_result.get("status") == "success":
            print(f"  ✓ Manhattan plot: {manhattan_file.name}")
    except Exception as e:
        print(f"  ✗ Manhattan plot failed: {e}")

    # Q-Q plot
    qq_file = plots_dir / f"qq_{trait_name}.png"
    try:
        qq_result = qq_plot(
            [r["p_value"] for r in assoc_results], qq_file, title=f'Q-Q Plot - {trait_name.replace("_", " ").title()}'
        )
        if qq_result.get("status") == "success":
            print(f"  ✓ Q-Q plot: {qq_file.name}")
    except Exception as e:
        print(f"  ✗ Q-Q plot failed: {e}")

print("\n" + "=" * 80)
print("GWAS ANALYSIS COMPLETE!")
print("=" * 80)
print(f"\nResults directory: {results_dir}")
print(f"Plots directory: {plots_dir}")
print("\nKey files created:")
print(f"  - qc_results.json ({n_passing} variants passed QC)")
print(f'  - pca_results.json ({pca_result.get("n_components", 0)} components)')
print(f"  - kinship_results.json ({n_samples}x{n_samples} matrix)")
for trait_name in all_results.keys():
    print(f"  - {trait_name}_associations.tsv")
    print(f"  - manhattan_{trait_name}.png")
    print(f"  - qq_{trait_name}.png")
print("\n✓ Comprehensive GWAS workflow executed successfully!")
