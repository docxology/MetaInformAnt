#!/usr/bin/env python3
"""Network analysis workflow orchestrator.

This script provides comprehensive orchestration for network analysis workflows,
including network construction, metrics calculation, community detection, and analysis.

Usage:
    python3 scripts/networks/run_network_analysis.py --input interactions.tsv --output output/networks/results
    python3 scripts/networks/run_network_analysis.py --input interactions.tsv --analyze-metrics --detect-communities
    python3 scripts/networks/run_network_analysis.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import dump_json, ensure_directory, read_csv
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Network analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic network construction and metrics
  %(prog)s --input interactions.tsv --output output/networks/basic

  # Full analysis with all modules
  %(prog)s --input interactions.tsv --analyze-metrics --detect-communities --analyze-centrality

  # Metrics and centrality only
  %(prog)s --input interactions.tsv --analyze-metrics --analyze-centrality
        """,
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input interaction data (TSV/CSV format with columns: node1, node2, [weight])",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/networks"),
        help="Output directory (default: output/networks)",
    )
    parser.add_argument(
        "--analyze-metrics",
        action="store_true",
        help="Calculate network metrics (density, clustering, etc.)",
    )
    parser.add_argument(
        "--detect-communities",
        action="store_true",
        help="Detect network communities using Louvain algorithm",
    )
    parser.add_argument(
        "--analyze-centrality",
        action="store_true",
        help="Calculate centrality measures (degree, betweenness, etc.)",
    )
    parser.add_argument(
        "--directed",
        action="store_true",
        help="Treat network as directed (default: undirected)",
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


def load_interactions(input_path: Path) -> list[tuple[str, str, float]]:
    """Load interaction data from file.

    Args:
        input_path: Path to input file

    Returns:
        List of (node1, node2, weight) tuples
    """
    logger.info(f"Loading interactions from {input_path}")

    try:
        df = read_csv(input_path, sep="\t" if input_path.suffix == ".tsv" else ",")

        # Find column names
        if "node1" in df.columns and "node2" in df.columns:
            node1_col = "node1"
            node2_col = "node2"
        elif len(df.columns) >= 2:
            node1_col = df.columns[0]
            node2_col = df.columns[1]
        else:
            raise ValueError("Need at least 2 columns for node1 and node2")

        # Get weight column if present
        weight_col = None
        if "weight" in df.columns:
            weight_col = "weight"
        elif len(df.columns) >= 3:
            weight_col = df.columns[2]

        interactions = []
        for _, row in df.iterrows():
            node1 = str(row[node1_col])
            node2 = str(row[node2_col])
            weight = float(row[weight_col]) if weight_col and weight_col in row else 1.0
            interactions.append((node1, node2, weight))

        logger.info(f"Loaded {len(interactions)} interactions")
        return interactions
    except Exception as e:
        logger.error(f"Failed to load interactions: {e}")
        raise


def analyze_metrics(network, output_dir: Path) -> dict[str, Any]:
    """Calculate network metrics.

    Args:
        network: BiologicalNetwork object
        output_dir: Output directory for results

    Returns:
        Dictionary with network metrics
    """
    logger.info("Calculating network metrics...")

    from metainformant.networks import network_metrics

    try:
        metrics = network_metrics(network)

        # Save results
        output_file = output_dir / "network_metrics.json"
        dump_json(metrics, output_file)
        logger.info(f"Network metrics saved to {output_file}")

        return metrics
    except Exception as e:
        logger.error(f"Failed to calculate metrics: {e}")
        return {"error": str(e)}


def analyze_centrality(network, output_dir: Path) -> dict[str, Any]:
    """Calculate centrality measures.

    Args:
        network: BiologicalNetwork object
        output_dir: Output directory for results

    Returns:
        Dictionary with centrality measures
    """
    logger.info("Calculating centrality measures...")

    from metainformant.networks import centrality_measures

    try:
        centralities = centrality_measures(network)

        # Save results
        output_file = output_dir / "centrality_measures.json"
        dump_json(centralities, output_file)
        logger.info(f"Centrality measures saved to {output_file}")

        return centralities
    except Exception as e:
        logger.error(f"Failed to calculate centrality: {e}")
        return {"error": str(e)}


def detect_communities_analysis(network, output_dir: Path) -> dict[str, Any]:
    """Detect network communities.

    Args:
        network: BiologicalNetwork object
        output_dir: Output directory for results

    Returns:
        Dictionary with community detection results
    """
    logger.info("Detecting network communities...")

    from metainformant.networks import community_metrics, detect_communities, modularity

    try:
        # Detect communities
        communities = detect_communities(network, method="louvain")

        # Calculate modularity
        mod = modularity(network, communities)

        # Calculate community metrics
        comm_metrics = community_metrics(network, communities)

        results = {
            "communities": communities,
            "modularity": mod,
            "community_metrics": comm_metrics,
            "num_communities": len(set(communities.values())),
        }

        # Save results
        output_file = output_dir / "community_detection.json"
        dump_json(results, output_file)
        logger.info(f"Community detection saved to {output_file}")
        logger.info(f"Detected {results['num_communities']} communities (modularity: {mod:.3f})")

        return results
    except Exception as e:
        logger.error(f"Failed to detect communities: {e}")
        return {"error": str(e)}


def run_workflow(args):
    """Execute network analysis workflow.

    Args:
        args: Parsed command-line arguments
    """
    logger.info("Starting network analysis workflow")
    logger.info(f"Input: {args.input}")
    logger.info(f"Output: {args.output}")

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        logger.info(f"Would analyze: {args.input}")
        logger.info(f"Would write to: {args.output}")
        if args.analyze_metrics:
            logger.info("Would calculate network metrics")
        if args.analyze_centrality:
            logger.info("Would calculate centrality measures")
        if args.detect_communities:
            logger.info("Would detect communities")
        return

    # Ensure output directory exists
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")

    # Load interactions
    try:
        interactions = load_interactions(args.input)
    except Exception as e:
        logger.error(f"Failed to load interactions: {e}")
        raise

    if not interactions:
        logger.warning("No interactions loaded")
        return

    # Construct network
    logger.info("Constructing network...")
    from metainformant.networks import add_edges_from_interactions, create_network

    network = create_network(directed=args.directed)
    add_edges_from_interactions(network, interactions)

    logger.info(f"Network constructed: {network.num_nodes()} nodes, {network.num_edges()} edges")

    # Run analyses
    workflow_results = {
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "network_size": {
            "nodes": network.num_nodes(),
            "edges": network.num_edges(),
            "directed": args.directed,
        },
        "analyses": {},
    }

    # Metrics analysis (always run if no specific analysis requested, or if explicitly requested)
    if args.analyze_metrics or not any([args.analyze_centrality, args.detect_communities]):
        try:
            metrics_results = analyze_metrics(network, output_dir)
            workflow_results["analyses"]["metrics"] = metrics_results
        except Exception as e:
            logger.error(f"Metrics analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["metrics"] = {"error": str(e)}

    # Centrality analysis
    if args.analyze_centrality:
        try:
            centrality_results = analyze_centrality(network, output_dir)
            workflow_results["analyses"]["centrality"] = centrality_results
        except Exception as e:
            logger.error(f"Centrality analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["centrality"] = {"error": str(e)}

    # Community detection
    if args.detect_communities:
        try:
            community_results = detect_communities_analysis(network, output_dir)
            workflow_results["analyses"]["communities"] = community_results
        except Exception as e:
            logger.error(f"Community detection failed: {e}", exc_info=True)
            workflow_results["analyses"]["communities"] = {"error": str(e)}

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
