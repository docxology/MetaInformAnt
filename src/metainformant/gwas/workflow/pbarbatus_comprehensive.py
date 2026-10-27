"""Run the comprehensive P. barbatus GWAS analysis with synthetic data.

The CLI wrapper lives in ``scripts/gwas/pipelines/run_pbarbatus_analysis.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from metainformant.core.io import dump_json, ensure_directory
from metainformant.gwas.analysis.association import association_test_linear, association_test_logistic
from metainformant.gwas.analysis.correction import bonferroni_correction, fdr_correction
from metainformant.gwas.analysis.quality import apply_qc_filters, parse_vcf_full
from metainformant.gwas.analysis.structure import compute_kinship_matrix, compute_pca
from metainformant.gwas.data.config import load_gwas_config
from metainformant.gwas.visualization._general_impl import manhattan_plot, qq_plot
from metainformant.gwas.workflow.workflow_config import _phenotypes_in_sample_order


def run(config_path: str = "config/gwas/gwas_pbarbatus_synthetic.yaml") -> None:
    """Run the comprehensive P. barbatus GWAS analysis."""

    print("=" * 80)
    print("P. BARBATUS GWAS - COMPREHENSIVE ANALYSIS")
    print("=" * 80)

    # Load configuration
    print("\nLoading configuration...")
    config = load_gwas_config(config_path)
    print("✓ Configuration loaded")

    # Parse VCF
    print("\n[1/7] Parsing VCF file...")
    vcf_file = config.variants["vcf_files"][0]
    vcf_data = parse_vcf_full(vcf_file)
    print(f'✓ Parsed {len(vcf_data.get("variants", []))} variants for {len(vcf_data.get("samples", []))} samples')

    # Quality control
    print("\n[2/7] Applying quality control filters...")
    # float(): YAML 1.1 loads "1e-6" (no decimal point) as a string
    qc_config = {
        "min_maf": float(config.qc.get("min_maf", 0.05)),
        "max_missing": float(config.qc.get("max_missing", 0.10)),
        "min_hwe_p": float(config.qc.get("hwe_pval", 1e-6)),
        "min_qual": float(config.qc.get("min_qual", 30.0)),
    }

    qc_result = apply_qc_filters(vcf_data, qc_config)

    if qc_result.get("status") == "success":
        n_passing = qc_result.get("num_variants_after", 0)
        n_total = qc_result.get("num_variants_before", 50000)
        print(f"✓ QC complete: {n_passing}/{n_total} variants passed ({n_passing/max(n_total, 1)*100:.1f}%)")
        if n_passing == 0:
            print("✗ No variants passed QC; nothing to analyze")
            sys.exit(1)

        # Everything downstream runs on the QC-filtered variants and genotypes
        filtered_data = qc_result["filtered_data"]
        genotypes = filtered_data["genotypes"]
        variants_passing = filtered_data["variants"]
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
    results_dir = ensure_directory(Path(config.work_dir) / "results")
    qc_summary = {key: value for key, value in qc_result.items() if key != "filtered_data"}
    dump_json(qc_summary, results_dir / "qc_results.json", indent=2)
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

    all_results: Dict[str, Dict[str, Any]] = {}

    for trait_name, test_type in traits_to_test:
        print(f"\n  Testing trait: {trait_name} ({test_type})...")

        phenotypes = _phenotypes_in_sample_order(pheno_df, trait_name, vcf_data["samples"])
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
                        "variant_id": variant.get("id", f"var_{idx}"),
                        "chrom": variant.get("chrom"),
                        "pos": variant.get("pos"),
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
    plots_dir = ensure_directory(Path(config.work_dir) / "plots")

    for trait_name in all_results.keys():
        assoc_results = all_results[trait_name]["associations"]

        # Manhattan plot
        manhattan_file = plots_dir / f"manhattan_{trait_name}.png"
        try:
            manhattan_fig = manhattan_plot(
                assoc_results, manhattan_file, title=f'{trait_name.replace("_", " ").title()} GWAS - P. barbatus'
            )
            if manhattan_fig is not None:
                plt.close(manhattan_fig)
                print(f"  ✓ Manhattan plot: {manhattan_file.name}")
        except Exception as e:
            print(f"  ✗ Manhattan plot failed: {e}")

        # Q-Q plot
        qq_file = plots_dir / f"qq_{trait_name}.png"
        try:
            qq_fig = qq_plot(
                [r["p_value"] for r in assoc_results],
                qq_file,
                title=f'Q-Q Plot - {trait_name.replace("_", " ").title()}',
            )
            if qq_fig is not None:
                plt.close(qq_fig)
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
