#!/usr/bin/env python3
"""Ontology simulation script.

This script generates synthetic GO term annotations and enrichment test data.

Usage:
    python3 scripts/simulation/simulate_ontology.py --type annotations --n-genes 1000 --n-terms 100
    python3 scripts/simulation/simulate_ontology.py --type enrichment --n-genes 500 --n-terms 50
"""

import argparse
import random
import sys
from pathlib import Path

import pandas as pd

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core import io, logging, paths, validation

logger = logging.get_logger(__name__)


def simulate_annotations(
    output_dir: Path,
    n_genes: int,
    n_terms: int,
    annotations_per_gene: int,
    seed: int,
) -> dict:
    """Simulate gene-to-GO-term annotations."""
    logger.info(f"Generating annotations: {n_genes} genes, {n_terms} GO terms")
    rng = random.Random(seed)
    
    # Generate GO term IDs
    go_terms = [f"GO:{rng.randint(1000000, 9999999):07d}" for _ in range(n_terms)]
    
    # Generate gene IDs
    genes = [f"gene_{i:05d}" for i in range(n_genes)]
    
    # Generate annotations
    annotations = []
    for gene in genes:
        n_annotations = rng.randint(1, annotations_per_gene)
        selected_terms = rng.sample(go_terms, min(n_annotations, len(go_terms)))
        
        for term in selected_terms:
            # Randomly assign namespace
            namespace = rng.choice(["biological_process", "molecular_function", "cellular_component"])
            
            annotations.append({
                "gene_id": gene,
                "go_term": term,
                "namespace": namespace,
                "evidence": rng.choice(["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "ISA", "ISM", "ISO"]),
            })
    
    # Save as CSV
    df = pd.DataFrame(annotations)
    csv_file = output_dir / "gene_annotations.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Annotations saved to {csv_file}")
    
    return {
        "type": "annotations",
        "n_genes": n_genes,
        "n_terms": n_terms,
        "n_annotations": len(annotations),
        "output_file": str(csv_file),
    }


def simulate_enrichment(
    output_dir: Path,
    n_genes: int,
    n_terms: int,
    n_significant: int,
    seed: int,
) -> dict:
    """Simulate enrichment test results."""
    logger.info(f"Generating enrichment results: {n_terms} terms")
    rng = random.Random(seed)
    
    # Generate GO terms
    go_terms = [f"GO:{rng.randint(1000000, 9999999):07d}" for _ in range(n_terms)]
    
    # Generate enrichment results
    results = []
    for term in go_terms:
        is_significant = len(results) < n_significant
        
        if is_significant:
            pvalue = rng.uniform(0.0001, 0.05)
            odds_ratio = rng.uniform(1.5, 10.0)
        else:
            pvalue = rng.uniform(0.05, 1.0)
            odds_ratio = rng.uniform(0.1, 1.5)
        
        results.append({
            "go_term": term,
            "term_name": f"simulated_term_{term.split(':')[1]}",
            "namespace": rng.choice(["biological_process", "molecular_function", "cellular_component"]),
            "n_genes_in_term": rng.randint(5, 100),
            "n_genes_in_list": rng.randint(10, n_genes),
            "n_genes_annotated": rng.randint(50, n_genes),
            "pvalue": pvalue,
            "adjusted_pvalue": pvalue * rng.uniform(1.0, 5.0),
            "odds_ratio": odds_ratio,
        })
    
    # Save as CSV
    df = pd.DataFrame(results)
    csv_file = output_dir / "enrichment_results.csv"
    df.to_csv(csv_file, index=False)
    
    logger.info(f"Enrichment results saved to {csv_file}")
    
    return {
        "type": "enrichment",
        "n_terms": n_terms,
        "n_significant": n_significant,
        "output_file": str(csv_file),
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Ontology simulation for testing and validation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Simulate gene annotations
  %(prog)s --type annotations --n-genes 1000 --n-terms 100

  # Simulate enrichment results
  %(prog)s --type enrichment --n-genes 500 --n-terms 50 --n-significant 10
        """,
    )
    
    parser.add_argument(
        "--type",
        required=True,
        choices=["annotations", "enrichment"],
        help="Type of simulation to run",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/simulation/ontology"),
        help="Output directory (default: output/simulation/ontology)",
    )
    parser.add_argument("--n-genes", type=int, default=1000, help="Number of genes")
    parser.add_argument("--n-terms", type=int, default=100, help="Number of GO terms")
    parser.add_argument("--annotations-per-gene", type=int, default=3, help="Average annotations per gene (annotations type)")
    parser.add_argument("--n-significant", type=int, default=10, help="Number of significant terms (enrichment type)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    args = parser.parse_args()
    
    if args.verbose:
        import logging as std_logging
        logger.setLevel(std_logging.DEBUG)
    
    output_dir = paths.ensure_directory(args.output)
    
    try:
        if args.type == "annotations":
            results = simulate_annotations(
                output_dir, args.n_genes, args.n_terms, args.annotations_per_gene, args.seed
            )
        elif args.type == "enrichment":
            results = simulate_enrichment(
                output_dir, args.n_genes, args.n_terms, args.n_significant, args.seed
            )
        
        # Save summary
        summary_file = output_dir / "simulation_summary.json"
        io.dump_json(results, summary_file, indent=2)
        logger.info(f"Simulation complete. Summary saved to {summary_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Simulation failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())

