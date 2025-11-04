#!/usr/bin/env python3
"""Protein analysis workflow orchestrator.

This script provides comprehensive orchestration for protein analysis workflows,
including sequence validation, composition analysis, secondary structure prediction,
and functional annotation.

Usage:
    python3 scripts/protein/run_protein_analysis.py --input sequences.fasta --output output/protein/results
    python3 scripts/protein/run_protein_analysis.py --input sequences.fasta --analyze-composition --analyze-secondary
    python3 scripts/protein/run_protein_analysis.py --help
"""

import argparse
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
        description="Protein analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic sequence analysis
  %(prog)s --input sequences.fasta --output output/protein/basic

  # Full analysis with all modules
  %(prog)s --input sequences.fasta --analyze-composition --analyze-secondary --analyze-alignment

  # Composition and validation only
  %(prog)s --input sequences.fasta --analyze-composition --min-length 50
        """
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input protein sequences (FASTA format)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/protein"),
        help="Output directory (default: output/protein)",
    )
    parser.add_argument(
        "--analyze-composition",
        action="store_true",
        help="Perform amino acid composition analysis",
    )
    parser.add_argument(
        "--analyze-secondary",
        action="store_true",
        help="Analyze secondary structure propensity",
    )
    parser.add_argument(
        "--analyze-alignment",
        action="store_true",
        help="Perform pairwise alignment analysis",
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
    """Analyze amino acid composition.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with composition analysis results
    """
    logger.info("Analyzing amino acid composition...")
    
    from metainformant.protein import calculate_aa_composition, is_valid_protein_sequence
    
    results = {
        "aa_composition": {},
        "sequence_lengths": {},
        "validation": {},
        "composition_stats": {},
    }
    
    all_compositions = []
    lengths = []
    valid_count = 0
    
    for seq_id, seq in sequences.items():
        # Validate sequence
        is_valid = is_valid_protein_sequence(seq)
        results["validation"][seq_id] = is_valid
        if is_valid:
            valid_count += 1
        
        # Amino acid composition
        comp = calculate_aa_composition(seq)
        results["aa_composition"][seq_id] = comp
        all_compositions.append(comp)
        
        # Sequence length
        length = len(seq)
        results["sequence_lengths"][seq_id] = length
        lengths.append(length)
    
    # Overall statistics
    if all_compositions:
        # Calculate mean composition across all sequences
        mean_comp = {}
        aa_set = set()
        for comp in all_compositions:
            aa_set.update(comp.keys())
        
        for aa in aa_set:
            mean_comp[aa] = sum(c.get(aa, 0) for c in all_compositions) / len(all_compositions)
        
        results["composition_stats"] = {
            "mean_aa_composition": mean_comp,
            "mean_length": sum(lengths) / len(lengths) if lengths else 0,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "total_sequences": len(sequences),
            "valid_sequences": valid_count,
        }
    
    # Save results
    output_file = output_dir / "composition_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Composition analysis saved to {output_file}")
    
    return results


def analyze_secondary(sequences: dict[str, str], output_dir: Path) -> dict[str, Any]:
    """Analyze secondary structure propensity.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with secondary structure analysis results
    """
    logger.info("Analyzing secondary structure propensity...")
    
    from metainformant.protein import simple_helix_coil_propensity, is_valid_protein_sequence
    
    results = {
        "propensity_profiles": {},
        "mean_propensity": {},
    }
    
    all_profiles = []
    valid_sequences = {}
    
    for seq_id, seq in sequences.items():
        if not is_valid_protein_sequence(seq):
            logger.warning(f"Skipping invalid sequence: {seq_id}")
            continue
        
        try:
            propensity = simple_helix_coil_propensity(seq)
            results["propensity_profiles"][seq_id] = propensity
            all_profiles.append(propensity)
            valid_sequences[seq_id] = seq
        except Exception as e:
            logger.warning(f"Could not calculate propensity for {seq_id}: {e}")
    
    # Calculate mean propensity across all sequences
    if all_profiles:
        # Find maximum length
        max_len = max(len(p) for p in all_profiles)
        
        # Calculate mean at each position
        mean_propensity = []
        for i in range(max_len):
            values = [p[i] for p in all_profiles if i < len(p)]
            if values:
                mean_propensity.append(sum(values) / len(values))
        
        results["mean_propensity"] = mean_propensity
        results["num_analyzed"] = len(all_profiles)
    
    # Save results
    output_file = output_dir / "secondary_structure_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Secondary structure analysis saved to {output_file}")
    
    return results


def analyze_alignment(sequences: dict[str, str], output_dir: Path) -> dict[str, Any]:
    """Perform pairwise alignment analysis.
    
    Args:
        sequences: Dictionary of sequence ID to sequence string
        output_dir: Output directory for results
        
    Returns:
        Dictionary with alignment results
    """
    logger.info("Analyzing pairwise alignments...")
    
    from metainformant.protein import alignment, is_valid_protein_sequence
    
    # Filter to valid sequences
    valid_seqs = {sid: seq for sid, seq in sequences.items() if is_valid_protein_sequence(seq)}
    
    results = {
        "pairwise_identity": {},
        "alignments": {},
    }
    
    seq_ids = list(valid_seqs.keys())
    
    if len(seq_ids) >= 2:
        # Calculate pairwise identity for all pairs
        for i, id1 in enumerate(seq_ids):
            for id2 in seq_ids[i+1:]:
                seq1 = valid_seqs[id1]
                seq2 = valid_seqs[id2]
                
                try:
                    identity = alignment.pairwise_identity(seq1, seq2)
                    pair_key = f"{id1}_vs_{id2}"
                    results["pairwise_identity"][pair_key] = identity
                    
                    # Perform full alignment for first few pairs
                    if len(results["alignments"]) < 10:
                        score, align1, align2 = alignment.needleman_wunsch(seq1, seq2)
                        results["alignments"][pair_key] = {
                            "score": score,
                            "alignment1": align1,
                            "alignment2": align2,
                        }
                except Exception as e:
                    logger.warning(f"Could not align {id1} vs {id2}: {e}")
        
        logger.info(f"Calculated {len(results['pairwise_identity'])} pairwise identities")
    else:
        logger.warning("Need at least 2 valid sequences for alignment analysis")
        results["note"] = "Insufficient sequences for alignment analysis"
    
    # Save results
    output_file = output_dir / "alignment_analysis.json"
    dump_json(results, output_file)
    logger.info(f"Alignment analysis saved to {output_file}")
    
    return results


def run_workflow(args):
    """Execute protein analysis workflow.
    
    Args:
        args: Parsed command-line arguments
    """
    logger.info("Starting protein analysis workflow")
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
        if args.analyze_secondary:
            logger.info("Would perform secondary structure analysis")
        if args.analyze_alignment:
            logger.info("Would perform alignment analysis")
        return
    
    # Ensure output directory exists
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    # Load sequences
    logger.info(f"Loading sequences from {args.input}")
    from metainformant.protein import parse_fasta
    
    try:
        all_sequences = parse_fasta(args.input)
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
    if args.analyze_composition or not any([args.analyze_secondary, args.analyze_alignment]):
        try:
            comp_results = analyze_composition(sequences, output_dir)
            workflow_results["analyses"]["composition"] = comp_results
        except Exception as e:
            logger.error(f"Composition analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["composition"] = {"error": str(e)}
    
    # Secondary structure analysis
    if args.analyze_secondary:
        try:
            sec_results = analyze_secondary(sequences, output_dir)
            workflow_results["analyses"]["secondary"] = sec_results
        except Exception as e:
            logger.error(f"Secondary structure analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["secondary"] = {"error": str(e)}
    
    # Alignment analysis
    if args.analyze_alignment:
        try:
            align_results = analyze_alignment(sequences, output_dir)
            workflow_results["analyses"]["alignment"] = align_results
        except Exception as e:
            logger.error(f"Alignment analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["alignment"] = {"error": str(e)}
    
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




