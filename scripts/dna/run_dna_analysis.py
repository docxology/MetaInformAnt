#!/usr/bin/env python3
"""DNA analysis workflow orchestrator.

This script provides a comprehensive orchestration layer for DNA analysis workflows,
including sequence processing, composition analysis, quality control, and reporting.

Usage:
    python3 scripts/dna/run_dna_analysis.py --input data/dna/sequences.fasta --output output/dna/results
    python3 scripts/dna/run_dna_analysis.py --input sequences.fasta --analyze-composition --analyze-population
    python3 scripts/dna/run_dna_analysis.py --help
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json
from metainformant.core.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="DNA analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic sequence analysis
  %(prog)s --input sequences.fasta --output output/dna/basic

  # Full analysis with all modules
  %(prog)s --input sequences.fasta --analyze-composition --analyze-population --analyze-phylogeny

  # Quality control and composition only
  %(prog)s --input sequences.fasta --analyze-composition --min-length 100
        """
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input DNA sequences (FASTA format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/dna"),
        help="Output directory (default: output/dna)",
    )
    parser.add_argument(
        "--analyze-composition",
        action="store_true",
        help="Perform sequence composition analysis (GC content, dinucleotides, etc.)",
    )
    parser.add_argument(
        "--analyze-population",
        action="store_true",
        help="Perform population genetics analysis (diversity, FST, etc.)",
    )
    parser.add_argument(
        "--analyze-phylogeny",
        action="store_true",
        help="Perform phylogenetic analysis (distance matrices, trees)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="Minimum sequence length to include (default: 0)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without executing",
    )
    return parser.parse_args()


def analyze_composition(sequences: dict[str, str], output_dir: Path) -> dict[str, Any]:
    """Analyze sequence composition.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with composition analysis results
    """
    logger.info("Analyzing sequence composition...")
    
    from metainformant.dna import sequences as dna_sequences, composition
    
    results = {
        "gc_content": {},
        "sequence_lengths": {},
        "composition_stats": {},
    }
    
    gc_values = []
    lengths = []
    
    for seq_id, seq in sequences.items():
        # GC content
        gc = dna_sequences.gc_content(seq)
        results["gc_content"][seq_id] = gc
        gc_values.append(gc)
        
        # Sequence length
        length = len(seq)
        results["sequence_lengths"][seq_id] = length
        lengths.append(length)
    
    # Overall statistics
    if gc_values:
        results["composition_stats"] = {
            "mean_gc_content": sum(gc_values) / len(gc_values),
            "min_gc_content": min(gc_values),
            "max_gc_content": max(gc_values),
            "mean_length": sum(lengths) / len(lengths),
            "min_length": min(lengths),
            "max_length": max(lengths),
            "total_sequences": len(sequences),
        }
    
    # Save results
    output_file = output_dir / "composition_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Composition analysis saved to {output_file}")
    
    return results


def analyze_population(sequences: dict[str, str], output_dir: Path) -> dict[str, Any]:
    """Analyze population genetics metrics.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with population genetics results
    """
    logger.info("Analyzing population genetics...")
    
    from metainformant.dna import population
    
    # Convert to list for population functions
    seq_list = list(sequences.values())
    
    results = {}
    
    if len(seq_list) >= 2:
        # Nucleotide diversity
        try:
            pi = population.nucleotide_diversity(seq_list)
            results["nucleotide_diversity"] = pi
            logger.info(f"Nucleotide diversity (π): {pi:.6f}")
        except Exception as e:
            logger.warning(f"Could not calculate nucleotide diversity: {e}")
        
        # Segregating sites
        try:
            segregating = population.segregating_sites(seq_list)
            results["segregating_sites"] = segregating
            logger.info(f"Segregating sites: {segregating}")
        except Exception as e:
            logger.warning(f"Could not calculate segregating sites: {e}")
    else:
        logger.warning("Need at least 2 sequences for population genetics analysis")
        results["note"] = "Insufficient sequences for population analysis"
    
    # Save results
    output_file = output_dir / "population_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Population genetics analysis saved to {output_file}")
    
    return results


def analyze_phylogeny(sequences: dict[str, str], output_dir: Path) -> dict[str, Any]:
    """Perform phylogenetic analysis.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with phylogenetic results
    """
    logger.info("Analyzing phylogeny...")
    
    from metainformant.dna import distances, phylogeny
    
    # Convert to list for distance calculations
    seq_list = list(sequences.values())
    seq_ids = list(sequences.keys())
    
    results = {}
    
    if len(seq_list) >= 3:
        try:
            # Calculate distance matrix
            dist_matrix = distances.distance_matrix(seq_list)
            results["distance_matrix"] = {
                "ids": seq_ids,
                "distances": dist_matrix.tolist() if hasattr(dist_matrix, 'tolist') else dist_matrix,
            }
            logger.info(f"Calculated distance matrix for {len(seq_list)} sequences")
            
            # Build phylogenetic tree (if possible)
            try:
                tree = phylogeny.build_tree_neighbor_joining(dist_matrix, seq_ids)
                results["tree_newick"] = tree
                logger.info("Built neighbor-joining tree")
                
                # Save tree to file
                tree_file = output_dir / "phylogenetic_tree.newick"
                with open(tree_file, 'w') as f:
                    f.write(tree)
                logger.info(f"Tree saved to {tree_file}")
            except Exception as e:
                logger.warning(f"Could not build tree: {e}")
                results["tree_error"] = str(e)
        except Exception as e:
            logger.warning(f"Could not calculate distance matrix: {e}")
            results["error"] = str(e)
    else:
        logger.warning("Need at least 3 sequences for phylogenetic analysis")
        results["note"] = "Insufficient sequences for phylogenetic analysis"
    
    # Save results
    output_file = output_dir / "phylogeny_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Phylogenetic analysis saved to {output_file}")
    
    return results


def run_workflow(args):
    """Execute DNA analysis workflow.
    
    Args:
        args: Parsed command-line arguments
    """
    logger.info("Starting DNA analysis workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")
    
    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        logger.info(f"Would analyze: {args.input}")
        logger.info(f"Would write to: {args.output}")
        if args.analyze_composition:
            logger.info("Would perform composition analysis")
        if args.analyze_population:
            logger.info("Would perform population genetics analysis")
        if args.analyze_phylogeny:
            logger.info("Would perform phylogenetic analysis")
        return
    
    # Ensure output directory exists
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load sequences
    logger.info(f"Loading sequences from {args.input}")
    from metainformant.dna import sequences as dna_sequences
    
    try:
        all_sequences = dna_sequences.read_fasta(str(args.input))
        logger.info(f"Loaded {len(all_sequences)} sequences")
    except Exception as e:
        logger.error(f"Failed to load sequences: {e}")
        raise
    
    # Filter by minimum length
    if args.min_length > 0:
        filtered = {sid: seq for sid, seq in all_sequences.items() if len(seq) >= args.min_length}
        logger.info(f"Filtered to {len(filtered)} sequences (min length: {args.min_length})")
        sequences = filtered
    else:
        sequences = all_sequences
    
    if not sequences:
        logger.warning("No sequences remaining after filtering")
        return
    
    # Run analyses
    workflow_results = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "total_sequences": len(sequences),
        "analyses": {},
    }
    
    # Composition analysis (always run if no specific analysis requested, or if explicitly requested)
    if args.analyze_composition or not any([args.analyze_population, args.analyze_phylogeny]):
        try:
            comp_results = analyze_composition(sequences, output_dir)
            workflow_results["analyses"]["composition"] = comp_results
        except Exception as e:
            logger.error(f"Composition analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["composition"] = {"error": str(e)}
    
    # Population genetics analysis
    if args.analyze_population:
        try:
            pop_results = analyze_population(sequences, output_dir)
            workflow_results["analyses"]["population"] = pop_results
        except Exception as e:
            logger.error(f"Population genetics analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["population"] = {"error": str(e)}
    
    # Phylogenetic analysis
    if args.analyze_phylogeny:
        try:
            phylo_results = analyze_phylogeny(sequences, output_dir)
            workflow_results["analyses"]["phylogeny"] = phylo_results
        except Exception as e:
            logger.error(f"Phylogenetic analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["phylogeny"] = {"error": str(e)}
    
    # Save workflow summary
    summary_file = output_dir / "workflow_summary.json"
    dump_json(workflow_results, summary_file, indent=2)
    logger.info(f"Workflow summary saved to {summary_file}")
    
    logger.info("Workflow complete")


def main():
    """Main entry point."""
    args = parse_args()
    
    # Setup logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    try:
        run_workflow(args)
        return 0
    except Exception as e:
        logger.error(f"Workflow failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())




