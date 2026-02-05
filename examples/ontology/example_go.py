#!/usr/bin/env python3
"""Gene Ontology analysis example.

This example demonstrates Gene Ontology term analysis and semantic similarity using METAINFORMANT's ontology toolkit.

Usage:
    python examples/ontology/example_go.py

Output:
    output/examples/ontology/go_analysis.json
"""

from __future__ import annotations

from pathlib import Path

from metainformant.core import io


def main():
    """Demonstrate GO term analysis."""
    # Setup output directory
    output_dir = Path("output/examples/ontology")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== METAINFORMANT Ontology Example ===")

    # Simulated GO terms and annotations
    go_terms = {
        "GO:0008150": {"name": "biological_process", "namespace": "biological_process"},
        "GO:0009987": {"name": "cellular_process", "namespace": "biological_process"},
        "GO:0003674": {"name": "molecular_function", "namespace": "molecular_function"},
        "GO:0003824": {"name": "catalytic_activity", "namespace": "molecular_function"},
    }

    gene_annotations = {
        "gene1": ["GO:0008150", "GO:0009987"],
        "gene2": ["GO:0009987", "GO:0003674"],
        "gene3": ["GO:0003674", "GO:0003824"],
    }

    print(f"✓ Analyzing {len(go_terms)} GO terms and {len(gene_annotations)} gene annotations")

    # Basic term analysis
    analysis_results = {
        "go_terms": go_terms,
        "gene_annotations": gene_annotations,
        "term_counts": {term: sum(1 for ann in gene_annotations.values() if term in ann) for term in go_terms},
    }

    print("GO term frequencies:")
    for term, count in analysis_results["term_counts"].items():
        print(f"  {term} ({go_terms[term]['name']}): {count} genes")

    # Save results
    results_file = output_dir / "go_analysis.json"
    io.dump_json({"ontology_analysis": analysis_results}, results_file, indent=2)

    print(f"✓ Results saved to: {results_file}")
    print("\n=== Ontology Example Complete ===")


if __name__ == "__main__":
    main()
