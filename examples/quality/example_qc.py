#!/usr/bin/env python3
"""Quality control analysis example.

This example demonstrates data quality assessment using METAINFORMANT's quality toolkit.

Usage:
    python examples/quality/example_qc.py

Output:
    output/examples/quality/qc_analysis.json
"""

from __future__ import annotations

from pathlib import Path
from metainformant.core import io

def main():
    """Demonstrate quality control analysis."""
    # Setup output directory
    output_dir = Path("output/examples/quality")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Quality Control Example ===")

    # Simulated quality metrics
    qc_metrics = {
        "sequencing_quality": {
            "mean_phred_score": 35.2,
            "bases_above_q30": 0.92,
            "adapter_contamination": 0.02,
            "status": "good"
        },
        "data_completeness": {
            "total_reads": 1000000,
            "mapped_reads": 950000,
            "mapping_rate": 0.95,
            "duplication_rate": 0.15
        },
        "biological_qc": {
            "genes_detected": 15000,
            "housekeeping_gene_coverage": 0.98,
            "gender_consistency": True,
            "contamination_estimate": 0.01
        }
    }

    # Overall quality assessment
    overall_quality = "pass" if (
        qc_metrics["sequencing_quality"]["bases_above_q30"] > 0.8 and
        qc_metrics["data_completeness"]["mapping_rate"] > 0.9 and
        qc_metrics["biological_qc"]["contamination_estimate"] < 0.05
    ) else "fail"

    results = {
        "qc_metrics": qc_metrics,
        "overall_assessment": overall_quality,
        "recommendations": [
            "Data quality is acceptable for downstream analysis" if overall_quality == "pass"
            else "Data quality issues detected - review before proceeding"
        ]
    }

    print(f"✓ Quality assessment: {overall_quality.upper()}")
    print(f"Mapping rate: {qc_metrics['data_completeness']['mapping_rate']:.1%}")
    print(f"Bases > Q30: {qc_metrics['sequencing_quality']['bases_above_q30']:.1%}")

    # Save results
    results_file = output_dir / "qc_analysis.json"
    io.dump_json({
        "quality_control": results
    }, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Quality Control Example Complete ===")

if __name__ == "__main__":
    main()
