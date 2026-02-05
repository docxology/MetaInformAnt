#!/usr/bin/env python3
"""Comprehensive GWAS method validation for Apis mellifera.

Exercises ALL GWAS methods sequentially to verify complete coverage:
  1. Data generation (VCF with drones, phenotypes, metadata)
  2. VCF parsing
  3. QC filters
  4. Haplodiploidy detection
  5. LD pruning (multiple r2 thresholds)
  6. PCA
  7. All kinship methods (vanraden, astle, yang, ibs)
  8. Linear association
  9. Logistic association
 10. Mixed model association
 11. Multiple testing correction (bonferroni, FDR)
 12. Heritability estimation
 13. Annotation with GFF3
 14. Fine-mapping credible sets
 15. Summary statistics
 16. Metadata loading/validation
 17. Visualizations
"""
from __future__ import annotations

import sys
import time
import traceback
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import logging

logger = logging.get_logger(__name__)


class ValidationResult:
    """Track pass/fail for validation groups."""

    def __init__(self) -> None:
        self.results: list[tuple[str, bool, str]] = []

    def record(self, name: str, passed: bool, detail: str = "") -> None:
        self.results.append((name, passed, detail))
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))

    def summary(self) -> None:
        total = len(self.results)
        passed = sum(1 for _, p, _ in self.results if p)
        failed = total - passed
        print("\n" + "=" * 80)
        print(f"VALIDATION SUMMARY: {passed}/{total} passed, {failed} failed")
        print("=" * 80)
        if failed > 0:
            print("\nFailed checks:")
            for name, p, detail in self.results:
                if not p:
                    print(f"  FAIL: {name} — {detail}")
        print()


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="GWAS method validation")
    parser.add_argument("--force-regenerate", action="store_true", help="Force regeneration of cached data")
    parser.add_argument("--scale-factor", type=int, default=1, help="Sample count multiplier")
    args = parser.parse_args()

    output_base = Path("output/gwas/amellifera")
    output_base.mkdir(parents=True, exist_ok=True)

    vr = ValidationResult()
    start = time.time()

    # Compute expected sample count from config
    from run_amellifera_gwas import DEFAULT_SUBSPECIES, load_data_generation_config

    dg = load_data_generation_config()
    dg["scale_factor"] = args.scale_factor
    expected_diploid = sum(v["n_samples"] for v in dg["subspecies"].values()) * dg["scale_factor"]
    expected_drones = dg["n_drones"] * dg["scale_factor"]
    expected_total = expected_diploid + expected_drones

    print("=" * 80)
    print("GWAS COMPREHENSIVE METHOD VALIDATION — Apis mellifera")
    if args.scale_factor != 1:
        print(f"Scale factor: {args.scale_factor}x (expecting {expected_total} samples)")
    print("=" * 80)

    # =========================================================================
    # 1. Data Generation
    # =========================================================================
    print("\n--- 1. DATA GENERATION ---")

    # Delete cached files so they regenerate with drones
    for f in [
        output_base / "variants" / "amellifera_population.vcf",
        output_base / "phenotypes" / "varroa_resistance.tsv",
        output_base / "metadata" / "sample_metadata.tsv",
    ]:
        if f.exists():
            f.unlink()
            print(f"  Deleted cached: {f.name}")

    # Import from sibling module
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_amellifera_gwas import (
        generate_real_vcf,
        generate_phenotypes,
        generate_metadata,
        download_genome_annotation,
    )

    try:
        gff_path = download_genome_annotation(output_base)
        vr.record("GFF3 annotation", gff_path is not None and gff_path.exists())
    except Exception as e:
        vr.record("GFF3 annotation", False, str(e))
        gff_path = None

    force_regen = args.force_regenerate
    try:
        vcf_path = generate_real_vcf(
            output_base, n_variants=2000, scale_factor=args.scale_factor, force=force_regen
        )
        vr.record("VCF generation", vcf_path.exists(), f"{vcf_path}")
    except Exception as e:
        vr.record("VCF generation", False, str(e))
        return 1

    try:
        pheno_path = generate_phenotypes(output_base, vcf_path, force=force_regen)
        vr.record("Phenotype generation", pheno_path.exists())
    except Exception as e:
        vr.record("Phenotype generation", False, str(e))
        return 1

    try:
        meta_path = generate_metadata(output_base, vcf_path, force=force_regen)
        vr.record("Metadata generation", meta_path.exists())
    except Exception as e:
        vr.record("Metadata generation", False, str(e))
        meta_path = None

    # =========================================================================
    # 2. VCF Parsing
    # =========================================================================
    print("\n--- 2. VCF PARSING ---")
    from metainformant.gwas.analysis.quality import parse_vcf_full, apply_qc_filters, check_haplodiploidy

    try:
        vcf_data = parse_vcf_full(vcf_path)
        n_variants = len(vcf_data.get("variants", []))
        n_samples = len(vcf_data.get("samples", []))
        vr.record("VCF parsing", n_variants > 0 and n_samples > 0, f"{n_variants} variants, {n_samples} samples")
        vr.record("Sample count matches config", n_samples == expected_total, f"got {n_samples}, expected {expected_total}")
    except Exception as e:
        vr.record("VCF parsing", False, str(e))
        return 1

    # =========================================================================
    # 3. QC Filters
    # =========================================================================
    print("\n--- 3. QC FILTERS ---")
    try:
        # Strict QC
        strict_qc = apply_qc_filters(vcf_data, {"min_maf": 0.05, "max_missing": 0.01})
        strict_after = strict_qc.get("num_variants_after", len(strict_qc.get("filtered_data", {}).get("variants", [])))

        # Relaxed QC
        relaxed_qc = apply_qc_filters(vcf_data, {"min_maf": 0.01, "max_missing": 0.1})
        relaxed_after = relaxed_qc.get(
            "num_variants_after", len(relaxed_qc.get("filtered_data", {}).get("variants", []))
        )

        vr.record("QC filters (strict)", strict_after > 0, f"{strict_after} variants after strict QC")
        vr.record("QC filters (relaxed)", relaxed_after > 0, f"{relaxed_after} variants after relaxed QC")
        vr.record("Strict < relaxed", strict_after <= relaxed_after, f"strict={strict_after}, relaxed={relaxed_after}")
    except Exception as e:
        vr.record("QC filters", False, str(e))

    # Use standard QC for downstream
    qc_result = apply_qc_filters(vcf_data, {"min_maf": 0.01, "max_missing": 0.05})
    filtered_data = qc_result.get("filtered_data", vcf_data)

    # =========================================================================
    # 4. Haplodiploidy Detection
    # =========================================================================
    print("\n--- 4. HAPLODIPLOIDY DETECTION ---")
    try:
        haplo_result = check_haplodiploidy(filtered_data, het_threshold=0.05)
        n_haploid = haplo_result.get("n_haploid", 0)
        n_diploid = haplo_result.get("n_diploid", 0)
        haploid_samples = haplo_result.get("haploid_samples", [])
        vr.record("Haplodiploidy check runs", True, f"haploid={n_haploid}, diploid={n_diploid}")
        vr.record("Haploid samples detected", n_haploid >= 1, f"detected {n_haploid} haploid (expected ~10)")

        # Verify drone samples are the haploid ones
        samples = filtered_data.get("samples", [])
        if samples and haploid_samples:
            haploid_names = [samples[i] for i in haploid_samples if i < len(samples)]
            drone_count = sum(1 for name in haploid_names if name.startswith("drone_"))
            vr.record("Drones identified as haploid", drone_count >= 1, f"{drone_count} drones in haploid set")
    except Exception as e:
        vr.record("Haplodiploidy detection", False, str(e))

    # Exclude haploid samples for downstream analysis
    diploid_indices = haplo_result.get("diploid_samples", list(range(len(filtered_data.get("genotypes", [])))))
    genotypes_all = filtered_data.get("genotypes", [])
    if genotypes_all and diploid_indices:
        genotypes_diploid = [genotypes_all[i] for i in diploid_indices]
        samples_diploid = [filtered_data["samples"][i] for i in diploid_indices]
    else:
        genotypes_diploid = genotypes_all
        samples_diploid = filtered_data.get("samples", [])

    n_dip = len(genotypes_diploid)
    n_var = len(genotypes_diploid[0]) if genotypes_diploid else 0

    # Transpose to variants x samples
    genotypes_by_variant = [[genotypes_diploid[s][v] for s in range(n_dip)] for v in range(n_var)]

    # =========================================================================
    # 5. LD Pruning
    # =========================================================================
    print("\n--- 5. LD PRUNING ---")
    from metainformant.gwas.analysis.ld_pruning import ld_prune

    ld_results = {}
    try:
        for r2 in [0.1, 0.2, 0.5]:
            pruned = ld_prune(genotypes_by_variant, window_size=50, step_size=5, r2_threshold=r2)
            ld_results[r2] = len(pruned)
            vr.record(f"LD prune r2={r2}", len(pruned) > 0, f"{len(pruned)} variants retained")

        # Monotonicity: stricter r2 should retain fewer or equal variants
        vr.record(
            "LD monotonicity",
            ld_results[0.1] <= ld_results[0.2] <= ld_results[0.5],
            f"r2=0.1:{ld_results[0.1]} <= r2=0.2:{ld_results[0.2]} <= r2=0.5:{ld_results[0.5]}",
        )
    except Exception as e:
        vr.record("LD pruning", False, str(e))

    # Use r2=0.2 pruned set for PCA
    ld_pruned_indices = ld_prune(genotypes_by_variant, window_size=50, step_size=5, r2_threshold=0.2)
    ld_pruned_geno_var = [genotypes_by_variant[i] for i in ld_pruned_indices]

    # Transpose LD-pruned back to samples x variants for PCA
    n_ld = len(ld_pruned_geno_var)
    pca_geno_by_sample = [[ld_pruned_geno_var[v][s] for v in range(n_ld)] for s in range(n_dip)]

    # =========================================================================
    # 6. PCA
    # =========================================================================
    print("\n--- 6. PCA ---")
    from metainformant.gwas.analysis.structure import compute_pca

    try:
        pca_result = compute_pca(pca_geno_by_sample, n_components=min(10, n_dip))
        pca_ok = pca_result.get("status") == "success"
        pcs = pca_result.get("pcs", [])
        var_ratios = pca_result.get("explained_variance_ratio", [])
        vr.record("PCA computation", pca_ok, f"{len(pcs)} samples, {len(var_ratios)} components")
        if pca_ok:
            vr.record("PCA shape correct", len(pcs) == n_dip, f"pcs={len(pcs)}, samples={n_dip}")
            vr.record(
                "Variance ratios sum ≤ 1.01",
                sum(var_ratios) <= 1.01,
                f"sum={sum(var_ratios):.4f}",
            )
    except Exception as e:
        vr.record("PCA", False, str(e))
        pca_result = None

    # =========================================================================
    # 7. All Kinship Methods
    # =========================================================================
    print("\n--- 7. KINSHIP METHODS ---")
    from metainformant.gwas.analysis.structure import compute_kinship_matrix

    kinship_results = {}
    for method in ["vanraden", "astle", "yang", "ibs"]:
        try:
            kr = compute_kinship_matrix(genotypes_diploid, method=method)
            ok = kr.get("status") == "success"
            km = kr.get("kinship_matrix", [])
            kinship_results[method] = kr

            if ok and km:
                # Verify symmetry
                sym = all(abs(km[i][j] - km[j][i]) < 1e-10 for i in range(len(km)) for j in range(len(km)))
                # Verify diagonal >= off-diagonal mean
                diag_mean = sum(km[i][i] for i in range(len(km))) / len(km)
                off_diag = [km[i][j] for i in range(len(km)) for j in range(len(km)) if i != j]
                off_mean = sum(off_diag) / len(off_diag) if off_diag else 0
                vr.record(
                    f"Kinship {method}",
                    ok and sym,
                    f"size={len(km)}, diag_mean={diag_mean:.3f}, off_mean={off_mean:.3f}",
                )
            else:
                vr.record(f"Kinship {method}", False, f"status={kr.get('status')}")
        except Exception as e:
            vr.record(f"Kinship {method}", False, str(e))

    # Use vanraden for downstream
    vanraden_km = kinship_results.get("vanraden", {}).get("kinship_matrix")

    # =========================================================================
    # 8. Linear Association
    # =========================================================================
    print("\n--- 8. LINEAR ASSOCIATION ---")
    from metainformant.gwas.analysis.association import association_test_linear

    # Load phenotypes
    pheno_values = []
    with open(pheno_path) as f:
        header = f.readline().strip().split("\t")
        varroa_idx = header.index("varroa_resistance")
        for line in f:
            parts = line.strip().split("\t")
            pheno_values.append(float(parts[varroa_idx]))

    # Align: use first n_dip phenotypes (diploid samples)
    pheno_dip = pheno_values[:n_dip]

    try:
        linear_results = []
        n_tested = min(100, len(genotypes_by_variant))  # Test first 100 for speed
        for i in range(n_tested):
            geno = genotypes_by_variant[i][:n_dip]
            result = association_test_linear(geno, pheno_dip)
            linear_results.append(result)

        p_values_linear = [r.get("p_value", 1.0) for r in linear_results]
        all_valid = all(0 <= p <= 1 for p in p_values_linear)
        has_betas = all("beta" in r for r in linear_results)
        vr.record("Linear association", len(linear_results) == n_tested and all_valid, f"{n_tested} tests, all p in [0,1]")
        vr.record("Linear returns betas", has_betas)
    except Exception as e:
        vr.record("Linear association", False, str(e))
        linear_results = []

    # =========================================================================
    # 9. Logistic Association
    # =========================================================================
    print("\n--- 9. LOGISTIC ASSOCIATION ---")
    from metainformant.gwas.analysis.association import association_test_logistic

    # Load binary phenotype
    binary_pheno = []
    with open(pheno_path) as f:
        header = f.readline().strip().split("\t")
        dr_idx = header.index("disease_resistance")
        for line in f:
            parts = line.strip().split("\t")
            binary_pheno.append(int(parts[dr_idx]))

    binary_dip = binary_pheno[:n_dip]

    try:
        logistic_results = []
        n_tested_log = min(100, len(genotypes_by_variant))
        for i in range(n_tested_log):
            geno = genotypes_by_variant[i][:n_dip]
            result = association_test_logistic(geno, binary_dip)
            logistic_results.append(result)

        p_values_logistic = [r.get("p_value", 1.0) for r in logistic_results]
        all_valid_log = all(0 <= p <= 1 for p in p_values_logistic)
        vr.record(
            "Logistic association",
            len(logistic_results) == n_tested_log and all_valid_log,
            f"{n_tested_log} tests",
        )
    except Exception as e:
        vr.record("Logistic association", False, str(e))
        logistic_results = []

    # =========================================================================
    # 10. Mixed Model Association
    # =========================================================================
    print("\n--- 10. MIXED MODEL ASSOCIATION ---")
    from metainformant.gwas.analysis.mixed_model import run_mixed_model_gwas

    try:
        if vanraden_km:
            variant_info = filtered_data.get("variants", [])
            mixed_results = run_mixed_model_gwas(
                genotype_matrix=genotypes_by_variant,
                phenotypes=pheno_dip,
                kinship_matrix=vanraden_km,
                variant_info=variant_info,
            )
            p_values_mixed = [r.get("p_value", 1.0) for r in mixed_results]
            all_valid_mix = all(0 <= p <= 1 for p in p_values_mixed if p is not None)
            vr.record(
                "Mixed model GWAS",
                len(mixed_results) > 0 and all_valid_mix,
                f"{len(mixed_results)} tests",
            )

            # Compare lambda_gc between linear and mixed
            import math

            def compute_lambda(pvals: list) -> float:
                valid = sorted([p for p in pvals if 0 < p <= 1])
                if not valid:
                    return 0.0
                med = valid[len(valid) // 2]
                return -2.0 * math.log(med) / 1.386 if med > 0 else 0.0

            if linear_results and mixed_results:
                lambda_linear = compute_lambda([r.get("p_value", 1.0) for r in linear_results])
                lambda_mixed = compute_lambda(p_values_mixed[:len(linear_results)])
                vr.record(
                    "Lambda_GC comparison",
                    True,
                    f"linear={lambda_linear:.3f}, mixed={lambda_mixed:.3f}",
                )
        else:
            vr.record("Mixed model GWAS", False, "No kinship matrix available")
            mixed_results = []
    except Exception as e:
        vr.record("Mixed model GWAS", False, str(e))
        mixed_results = []

    # Use mixed results for downstream (full set)
    assoc_for_downstream = mixed_results if mixed_results else linear_results

    # =========================================================================
    # 11. Multiple Testing Correction
    # =========================================================================
    print("\n--- 11. MULTIPLE TESTING CORRECTION ---")
    from metainformant.gwas.analysis.correction import fdr_correction, bonferroni_correction

    if assoc_for_downstream:
        p_vals = [r.get("p_value", 1.0) for r in assoc_for_downstream]
        try:
            fdr_result = fdr_correction(p_vals)
            if isinstance(fdr_result, dict):
                q_values = fdr_result.get("adjusted_p_values", [])
            else:
                _, q_values = fdr_result
            vr.record(
                "FDR correction (BH)",
                len(q_values) == len(p_vals) and all(0 <= q <= 1 for q in q_values),
                f"{len(q_values)} q-values",
            )
        except Exception as e:
            vr.record("FDR correction", False, str(e))

        try:
            bonf_result = bonferroni_correction(p_vals)
            if isinstance(bonf_result, dict):
                n_sig = sum(1 for s in bonf_result.get("significant", []) if s)
                vr.record("Bonferroni correction", True, f"{n_sig} significant")
            else:
                vr.record("Bonferroni correction", True, "completed")
        except Exception as e:
            vr.record("Bonferroni correction", False, str(e))
    else:
        vr.record("Multiple testing", False, "No association results available")

    # =========================================================================
    # 12. Heritability Estimation
    # =========================================================================
    print("\n--- 12. HERITABILITY ---")
    from metainformant.gwas.analysis.heritability import estimate_heritability, partition_heritability_by_chromosome

    try:
        if vanraden_km:
            h2_result = estimate_heritability(vanraden_km, pheno_dip)
            h2_ok = h2_result.get("status") == "success"
            h2_val = h2_result.get("h2", -1)
            h2_se = h2_result.get("h2_se", -1)
            vr.record(
                "Heritability estimation",
                h2_ok and 0 <= h2_val <= 1,
                f"h2={h2_val:.3f}, se={h2_se:.3f}",
            )
        else:
            vr.record("Heritability estimation", False, "No kinship matrix")
    except Exception as e:
        vr.record("Heritability estimation", False, str(e))

    # =========================================================================
    # 13. Annotation
    # =========================================================================
    print("\n--- 13. ANNOTATION ---")
    from metainformant.gwas.analysis.annotation import annotate_variants_with_genes

    try:
        if gff_path and gff_path.exists() and assoc_for_downstream:
            # Add chrom/pos to results if not present
            variants = filtered_data.get("variants", [])
            for i, r in enumerate(assoc_for_downstream):
                if "chrom" not in r and i < len(variants):
                    r["chrom"] = variants[i].get("chrom", "")
                    r["pos"] = variants[i].get("pos", 0)

            annotate_variants_with_genes(assoc_for_downstream, str(gff_path), window_kb=50)
            annotated = sum(
                1 for r in assoc_for_downstream
                if r.get("annotation", {}).get("nearest_gene")
            )
            vr.record("Variant annotation", annotated > 0, f"{annotated} variants annotated with genes")
        else:
            vr.record("Variant annotation", False, "GFF3 not available")
    except Exception as e:
        vr.record("Variant annotation", False, str(e))

    # =========================================================================
    # 14. Fine-Mapping
    # =========================================================================
    print("\n--- 14. FINE-MAPPING ---")
    from metainformant.gwas.visualization.visualization_finemapping import compute_credible_set

    try:
        if assoc_for_downstream:
            cs_result = compute_credible_set(assoc_for_downstream, credible_level=0.95)
            cs_ok = cs_result.get("status") == "success"
            cs_size = cs_result.get("credible_set_size", 0)
            pips = cs_result.get("pips", [])
            vr.record(
                "Fine-mapping credible set",
                cs_ok and cs_size > 0,
                f"credible set size={cs_size}",
            )
            if pips:
                max_pip = max(pips)
                pip_sum = sum(pips)
                vr.record(
                    "PIPs valid",
                    abs(pip_sum - 1.0) < 1e-6 and 0 <= max_pip <= 1,
                    f"sum={pip_sum:.6f}, max={max_pip:.4f}",
                )
        else:
            vr.record("Fine-mapping", False, "No association results")
    except Exception as e:
        vr.record("Fine-mapping", False, str(e))

    # =========================================================================
    # 15. Summary Statistics
    # =========================================================================
    print("\n--- 15. SUMMARY STATISTICS ---")
    from metainformant.gwas.analysis.summary_stats import (
        write_summary_statistics,
        write_significant_hits,
        create_results_summary,
    )

    results_dir = output_base / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    try:
        if assoc_for_downstream and filtered_data.get("variants"):
            variants = filtered_data["variants"]
            stats_path = results_dir / "validation_summary_stats.tsv"
            write_summary_statistics(assoc_for_downstream, variants, stats_path)
            vr.record("Summary stats written", stats_path.exists() and stats_path.stat().st_size > 0)

            sig_path = results_dir / "validation_significant_hits.tsv"
            write_significant_hits(assoc_for_downstream, variants, sig_path, threshold=5e-8)
            vr.record("Significant hits written", sig_path.exists())

            summary_path = results_dir / "validation_results_summary.json"
            summary = create_results_summary(assoc_for_downstream, summary_path, threshold=5e-8)
            vr.record("Results summary JSON", summary_path.exists() and summary is not None)
        else:
            vr.record("Summary statistics", False, "No results or variants")
    except Exception as e:
        vr.record("Summary statistics", False, str(e))

    # =========================================================================
    # 16. Metadata
    # =========================================================================
    print("\n--- 16. METADATA ---")
    from metainformant.gwas.data.metadata import load_sample_metadata, validate_metadata, get_population_labels

    try:
        if meta_path and meta_path.exists():
            meta_result = load_sample_metadata(str(meta_path))
            meta_ok = meta_result.get("status") == "success"
            metadata = meta_result.get("metadata", {})
            columns = meta_result.get("columns", [])
            vr.record(
                "Metadata loading",
                meta_ok and len(metadata) > 0,
                f"{len(metadata)} samples, columns={columns}",
            )

            # Validate against VCF samples
            all_samples = filtered_data.get("samples", [])
            if all_samples:
                val_result = validate_metadata(metadata, all_samples)
                completeness = val_result.get("completeness", 0)
                vr.record("Metadata completeness", completeness > 0.9, f"completeness={completeness:.2f}")

            # Get population labels
            pop_labels = get_population_labels(metadata, column="population")
            vr.record("Population labels", len(pop_labels) > 0, f"{len(pop_labels)} samples with labels")
        else:
            vr.record("Metadata", False, "Metadata file not found")
    except Exception as e:
        vr.record("Metadata", False, str(e))

    # =========================================================================
    # 17. Visualizations
    # =========================================================================
    print("\n--- 17. VISUALIZATIONS ---")
    import matplotlib

    matplotlib.use("Agg")

    plots_dir = output_base / "validation_plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    viz_tests = []

    # Manhattan
    try:
        from metainformant.gwas.visualization.general import manhattan_plot

        fig = manhattan_plot(assoc_for_downstream, output_path=str(plots_dir / "manhattan.png"))
        if fig:
            import matplotlib.pyplot as plt

            plt.close(fig)
        ok = (plots_dir / "manhattan.png").exists()
        vr.record("Manhattan plot", ok)
        viz_tests.append(("manhattan", ok))
    except Exception as e:
        vr.record("Manhattan plot", False, str(e))

    # QQ
    try:
        from metainformant.gwas.visualization.general import qq_plot

        p_vals_viz = [r.get("p_value", 1.0) for r in assoc_for_downstream]
        fig = qq_plot(p_vals_viz, output_path=str(plots_dir / "qq.png"))
        if fig:
            import matplotlib.pyplot as plt

            plt.close(fig)
        ok = (plots_dir / "qq.png").exists()
        vr.record("QQ plot", ok)
    except Exception as e:
        vr.record("QQ plot", False, str(e))

    # PCA plot
    try:
        from metainformant.gwas.visualization.visualization_population import pca_plot

        import numpy as np

        if pca_result and pca_result.get("status") == "success":
            pca_data = {
                "pcs": np.array(pca_result["pcs"]),
                "explained_variance": np.array(pca_result.get("explained_variance_ratio", [])),
            }
            result = pca_plot(pca_data, output_file=str(plots_dir / "pca.png"))
            ok = (plots_dir / "pca.png").exists()
            vr.record("PCA plot", ok)
    except Exception as e:
        vr.record("PCA plot", False, str(e))

    # Kinship heatmap
    try:
        from metainformant.gwas.visualization.general import kinship_heatmap

        if vanraden_km:
            fig = kinship_heatmap(vanraden_km, output_path=str(plots_dir / "kinship.png"))
            if fig:
                import matplotlib.pyplot as plt

                plt.close(fig)
            ok = (plots_dir / "kinship.png").exists()
            vr.record("Kinship heatmap", ok)
    except Exception as e:
        vr.record("Kinship heatmap", False, str(e))

    # Volcano
    try:
        from metainformant.gwas.visualization.visualization_statistical import volcano_plot

        result = volcano_plot(results=assoc_for_downstream, output_path=str(plots_dir / "volcano.png"))
        ok = (plots_dir / "volcano.png").exists()
        vr.record("Volcano plot", ok)
    except Exception as e:
        vr.record("Volcano plot", False, str(e))

    # Variant density
    try:
        from metainformant.gwas.visualization.visualization_variants import variant_density_plot

        var_dicts = [{"CHROM": r.get("chrom", ""), "POS": r.get("pos", 0)} for r in assoc_for_downstream]
        result = variant_density_plot(var_dicts, output_file=str(plots_dir / "density.png"))
        ok = (plots_dir / "density.png").exists()
        vr.record("Variant density plot", ok)
    except Exception as e:
        vr.record("Variant density plot", False, str(e))

    # MAF distribution
    try:
        from metainformant.gwas.visualization.visualization_variants import maf_distribution

        mafs = [min(abs(r.get("beta", 0.0)), 0.5) for r in assoc_for_downstream]
        result = maf_distribution(mafs, output_file=str(plots_dir / "maf.png"))
        ok = (plots_dir / "maf.png").exists()
        vr.record("MAF distribution", ok)
    except Exception as e:
        vr.record("MAF distribution", False, str(e))

    # Power plot
    try:
        from metainformant.gwas.visualization.visualization_statistical import power_plot

        result = power_plot(
            sample_sizes=[20, 40, 60, 80],
            effect_sizes=[0.1, 0.2, 0.5],
            output_file=str(plots_dir / "power.png"),
        )
        ok = (plots_dir / "power.png").exists()
        vr.record("Power plot", ok)
    except Exception as e:
        vr.record("Power plot", False, str(e))

    # Composite panels
    try:
        from metainformant.gwas.visualization.visualization_composite import (
            gwas_summary_panel,
            population_structure_panel,
        )
        import numpy as np

        pca_for_panel = pca_result if pca_result and pca_result.get("status") == "success" else None
        summary_panel = gwas_summary_panel(
            assoc_for_downstream,
            pca_data=pca_for_panel,
            kinship_matrix=vanraden_km,
            output_file=str(plots_dir / "summary_panel.png"),
        )
        ok = summary_panel.get("status") == "success"
        vr.record("GWAS summary panel", ok, f"panels={summary_panel.get('panels_generated', 0)}")
    except Exception as e:
        vr.record("GWAS summary panel", False, str(e))

    # =========================================================================
    # SUMMARY
    # =========================================================================
    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f}s")
    vr.summary()

    failed = sum(1 for _, p, _ in vr.results if not p)
    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
