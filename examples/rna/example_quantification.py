#!/usr/bin/env python3
"""RNA expression quantification example.

This example demonstrates RNA-seq expression quantification and analysis using METAINFORMANT's RNA analysis toolkit.

Usage:
    python examples/rna/example_quantification.py

Output:
    output/examples/rna/expression_analysis.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate RNA expression quantification."""
    # Setup output directory
    output_dir = Path("output/examples/rna")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT RNA Quantification Example ===")

    # 1. Create simulated RNA-seq count data
    print("\n1. Creating simulated RNA-seq count data...")

    # Simulated gene expression data
    genes = [f"Gene_{i:03d}" for i in range(1, 21)]  # 20 genes

    # Simulated count data for 3 samples
    raw_counts = {
        "sample_1": [150, 300, 450, 120, 800, 90, 650, 220, 1100, 330,
                     780, 95, 420, 180, 560, 290, 1340, 410, 670, 185],
        "sample_2": [180, 280, 420, 140, 750, 110, 680, 250, 950, 310,
                     820, 85, 380, 200, 590, 270, 1280, 390, 690, 205],
        "sample_3": [120, 350, 480, 100, 900, 75, 620, 190, 1200, 360,
                     740, 105, 460, 160, 530, 310, 1400, 430, 650, 170]
    }

    # Gene lengths for normalization
    gene_lengths = [1200, 1800, 900, 2100, 1500, 800, 1900, 1100, 2500, 1400,
                    1700, 950, 1300, 1000, 1600, 1200, 2800, 1350, 1750, 1150]

    print(f"✓ Created expression data for {len(genes)} genes across {len(raw_counts)} samples")

    # 2. Basic count statistics
    print("\n2. Calculating basic count statistics...")

    count_stats = {}

    for sample, counts in raw_counts.items():
        total_reads = sum(counts)
        mean_expression = total_reads / len(counts)
        max_expression = max(counts)
        expressed_genes = sum(1 for c in counts if c > 0)

        count_stats[sample] = {
            "total_reads": total_reads,
            "mean_expression": mean_expression,
            "max_expression": max_expression,
            "expressed_genes": expressed_genes,
            "expression_range": f"{min(counts)} - {max(counts)}"
        }

        print(f"  {sample}: {total_reads} total reads, {expressed_genes} expressed genes")

    # 3. TPM (Transcripts Per Million) normalization
    print("\n3. Calculating TPM normalization...")

    def calculate_tpm(counts, lengths):
        """Calculate TPM (Transcripts Per Million) values."""
        # Step 1: Calculate reads per kilobase (RPK)
        rpk = [count / (length / 1000) for count, length in zip(counts, lengths)]

        # Step 2: Calculate total RPK
        total_rpk = sum(rpk)

        # Step 3: Calculate TPM
        if total_rpk > 0:
            tpm = [(rpk_value / total_rpk) * 1_000_000 for rpk_value in rpk]
        else:
            tpm = [0] * len(rpk)

        return tpm

    tpm_normalized = {}

    for sample, counts in raw_counts.items():
        tpm_values = calculate_tpm(counts, gene_lengths)
        tpm_normalized[sample] = tpm_values

        # Calculate TPM statistics
        total_tpm = sum(tpm_values)
        mean_tpm = total_tpm / len(tpm_values)

        print(f"  {sample}: Total TPM = {total_tpm:.1f}, Mean TPM = {mean_tpm:.2f}")

    # 4. FPKM (Fragments Per Kilobase Million) normalization
    print("\n4. Calculating FPKM normalization...")

    def calculate_fpkm(counts, lengths, total_mapped_reads=None):
        """Calculate FPKM (Fragments Per Kilobase Million) values."""
        # For demonstration, assume total mapped reads if not provided
        if total_mapped_reads is None:
            total_mapped_reads = sum(counts)  # Simplified

        fpkm_values = []
        for count, length in zip(counts, lengths):
            # FPKM = (count * 1e9) / (length * total_mapped_reads)
            fpkm = (count * 1_000_000_000) / (length * total_mapped_reads)
            fpkm_values.append(fpkm)

        return fpkm_values

    fpkm_normalized = {}

    for sample, counts in raw_counts.items():
        fpkm_values = calculate_fpkm(counts, gene_lengths)
        fpkm_normalized[sample] = fpkm_values

        mean_fpkm = sum(fpkm_values) / len(fpkm_values)
        print(f"  {sample}: Mean FPKM = {mean_fpkm:.3f}")

    # 5. Differential expression analysis (simplified)
    print("\n5. Performing simplified differential expression analysis...")

    # Compare sample 1 vs sample 3 (assuming different conditions)
    sample_1_tpm = tpm_normalized["sample_1"]
    sample_3_tpm = tpm_normalized["sample_3"]

    fold_changes = []
    diff_expressed = []

    for i, (tpm1, tpm3) in enumerate(zip(sample_1_tpm, sample_3_tpm)):
        if tpm1 > 0 and tpm3 > 0:
            fc = tpm3 / tpm1
            log2_fc = log2(fc) if fc > 0 else 0
            fold_changes.append((genes[i], fc, log2_fc))

            # Simple threshold for "differential expression"
            if abs(log2_fc) > 1:  # 2-fold change
                diff_expressed.append((genes[i], log2_fc))

    print(f"  Analyzed {len(fold_changes)} genes")
    print(f"  Found {len(diff_expressed)} differentially expressed genes (|log2FC| > 1)")

    if diff_expressed:
        # Show top 3 by absolute fold change
        top_diff = sorted(diff_expressed, key=lambda x: abs(x[1]), reverse=True)[:3]
        for gene, log2fc in top_diff:
            direction = "up" if log2fc > 0 else "down"
            print(f"    {gene}: {direction}-regulated (log2FC = {log2fc:.2f})")

    # 6. Quality control metrics
    print("\n6. Calculating quality control metrics...")

    qc_metrics = {}

    for sample, counts in raw_counts.items():
        # Percentage of reads from top-expressed genes
        sorted_counts = sorted(counts, reverse=True)
        top_10_percent = sum(sorted_counts[:2])  # Top 10% of 20 genes = 2 genes
        total_reads = sum(counts)
        pct_top_expressed = (top_10_percent / total_reads) * 100 if total_reads > 0 else 0

        # Number of undetected genes (count = 0)
        undetected = sum(1 for c in counts if c == 0)

        # Expression distribution
        low_expression = sum(1 for c in counts if c < 100)
        medium_expression = sum(1 for c in counts if 100 <= c < 500)
        high_expression = sum(1 for c in counts if c >= 500)

        qc_metrics[sample] = {
            "percent_reads_top_expressed": pct_top_expressed,
            "undetected_genes": undetected,
            "expression_distribution": {
                "low (< 100)": low_expression,
                "medium (100-499)": medium_expression,
                "high (≥ 500)": high_expression
            },
            "qc_passed": pct_top_expressed < 80  # Arbitrary threshold
        }

        print(f"  {sample}: {undetected} undetected genes, {pct_top_expressed:.1f}% reads from top expressed")

    # 7. Create comprehensive expression analysis results
    print("\n7. Creating comprehensive expression analysis results...")

    expression_results = {
        "rna_expression_quantification_demo": {
            "timestamp": "2024-12-26T10:00:00Z",
            "input_data": {
                "genes": genes,
                "samples": list(raw_counts.keys()),
                "gene_lengths": gene_lengths
            },
            "raw_counts": raw_counts,
            "normalization_methods": {
                "tpm": {
                    "description": "Transcripts Per Million - normalizes for gene length and library size",
                    "formula": "TPM = (RPK / total_RPK) × 1,000,000",
                    "use_case": "Comparing expression within a sample",
                    "results": tpm_normalized
                },
                "fpkm": {
                    "description": "Fragments Per Kilobase Million - normalizes for gene length and library size",
                    "formula": "FPKM = (count × 1e9) / (length × total_mapped_reads)",
                    "use_case": "Comparing expression across samples",
                    "results": fpkm_normalized
                }
            },
            "differential_expression": {
                "comparison": "sample_1 vs sample_3",
                "method": "Simplified fold change analysis",
                "threshold": "|log2FC| > 1",
                "differentially_expressed_genes": len(diff_expressed),
                "top_changes": diff_expressed[:5] if len(diff_expressed) > 5 else diff_expressed
            },
            "quality_control": qc_metrics,
            "summary_statistics": {
                "total_genes": len(genes),
                "total_samples": len(raw_counts),
                "average_reads_per_sample": sum(sum(counts) for counts in raw_counts.values()) / len(raw_counts),
                "genes_with_expression": sum(1 for gene_idx in range(len(genes))
                                           if any(counts[gene_idx] > 0 for counts in raw_counts.values())),
                "normalization_methods_applied": ["TPM", "FPKM"],
                "qc_checks_performed": ["Expression distribution", "Top gene contribution", "Undetected genes"]
            },
            "key_concepts_demonstrated": [
                "Raw count processing and statistics",
                "TPM and FPKM normalization methods",
                "Differential expression analysis",
                "Quality control metrics",
                "Expression distribution analysis"
            ]
        }
    }

    results_file = output_dir / "expression_analysis.json"
    io.dump_json(expression_results, results_file, indent=2)

    print(f"✓ Comprehensive expression quantification analysis saved to: {results_file}")

    print("\n=== RNA Quantification Example Complete ===")
    print("This example demonstrated METAINFORMANT's RNA expression analysis capabilities:")
    print("- Raw count processing and quality assessment")
    print("- TPM and FPKM normalization techniques")
    print("- Differential expression analysis")
    print("- Quality control metrics calculation")

    print(f"\nAll outputs saved to: {output_dir}")

def log2(x):
    """Simple log2 calculation."""
    import math
    return math.log2(x) if x > 0 else 0

if __name__ == "__main__":
    main()
