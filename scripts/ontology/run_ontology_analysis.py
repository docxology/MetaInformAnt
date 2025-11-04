#!/usr/bin/env python3
"""Ontology analysis workflow orchestrator.

This script provides comprehensive orchestration for ontology analysis workflows,
including GO loading, term queries, and ontology summaries.

Usage:
    python3 scripts/ontology/run_ontology_analysis.py --go go.obo --output output/ontology/results
    python3 scripts/ontology/run_ontology_analysis.py --go go.obo --query-term GO:0008150 --ancestors --descendants
    python3 scripts/ontology/run_ontology_analysis.py --help
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
        description="Ontology analysis workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Load GO and generate summary
  %(prog)s --go go.obo --write-summary --output output/ontology/go_summary

  # Query specific term
  %(prog)s --go go.obo --query-term GO:0008150 --ancestors --descendants

  # Extract subgraph
  %(prog)s --go go.obo --subgraph GO:0008150,GO:0003674 --output output/ontology/subgraph
        """
    )
    parser.add_argument(
        "--go",
        type=Path,
        help="Input Gene Ontology OBO file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/ontology"),
        help="Output directory (default: output/ontology)",
    )
    parser.add_argument(
        "--write-summary",
        action="store_true",
        help="Write GO summary to file",
    )
    parser.add_argument(
        "--query-term",
        type=str,
        help="Query specific GO term ID",
    )
    parser.add_argument(
        "--ancestors",
        action="store_true",
        help="Get ancestor terms for query term",
    )
    parser.add_argument(
        "--descendants",
        action="store_true",
        help="Get descendant terms for query term",
    )
    parser.add_argument(
        "--subgraph",
        type=str,
        help="Extract subgraph from root terms (comma-separated)",
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


def run_workflow(args):
    """Execute ontology analysis workflow."""
    logger.info("Starting ontology analysis workflow")
    logger.info(f"Output: {args.output}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        if args.go:
            logger.info(f"Would load GO: {args.go}")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    workflow_results = {
        "output_dir": str(output_dir),
        "analyses": {},
    }
    
    # Load GO ontology
    onto = None
    if args.go:
        if not args.go.exists():
            raise FileNotFoundError(f"GO file not found: {args.go}")
        
        try:
            logger.info(f"Loading Gene Ontology from {args.go}")
            from metainformant.ontology import load_go_obo
            
            onto = load_go_obo(args.go)
            logger.info(f"Loaded {onto.num_terms()} GO terms")
            workflow_results["go_file"] = str(args.go)
            workflow_results["n_terms"] = onto.num_terms()
        except Exception as e:
            logger.error(f"Failed to load GO: {e}", exc_info=True)
            workflow_results["analyses"]["load_go"] = {"error": str(e)}
    
    if onto is None:
        logger.warning("No ontology loaded - cannot perform analysis")
        return
    
    # Write summary
    if args.write_summary:
        try:
            logger.info("Writing GO summary...")
            from metainformant.ontology import write_go_summary
            
            summary_file = write_go_summary(onto, dest=output_dir / "go_summary.txt")
            logger.info(f"GO summary written to {summary_file}")
            workflow_results["analyses"]["summary"] = {"file": str(summary_file)}
        except Exception as e:
            logger.error(f"Summary writing failed: {e}", exc_info=True)
            workflow_results["analyses"]["summary"] = {"error": str(e)}
    
    # Query term
    if args.query_term:
        try:
            logger.info(f"Querying term: {args.query_term}")
            from metainformant.ontology.query import ancestors, descendants
            
            query_results = {}
            
            if onto.has_term(args.query_term):
                term = onto.terms[args.query_term]
                query_results["term_info"] = {
                    "id": term.term_id,
                    "name": term.name,
                    "namespace": term.namespace,
                    "definition": term.definition[:200] if term.definition else None,  # Truncate
                }
                
                if args.ancestors:
                    anc = ancestors(onto, args.query_term)
                    query_results["ancestors"] = {
                        "count": len(anc),
                        "terms": list(anc)[:50],  # Limit output
                    }
                    logger.info(f"Found {len(anc)} ancestor terms")
                
                if args.descendants:
                    desc = descendants(onto, args.query_term)
                    query_results["descendants"] = {
                        "count": len(desc),
                        "terms": list(desc)[:50],  # Limit output
                    }
                    logger.info(f"Found {len(desc)} descendant terms")
            else:
                query_results["error"] = f"Term {args.query_term} not found in ontology"
                logger.warning(f"Term {args.query_term} not found")
            
            output_file = output_dir / "term_query.json"
            dump_json(query_results, output_file)
            workflow_results["analyses"]["term_query"] = query_results
            logger.info(f"Term query results saved to {output_file}")
        except Exception as e:
            logger.error(f"Term query failed: {e}", exc_info=True)
            workflow_results["analyses"]["term_query"] = {"error": str(e)}
    
    # Subgraph extraction
    if args.subgraph:
        try:
            logger.info("Extracting subgraph...")
            from metainformant.ontology.query import subgraph
            
            root_terms = [t.strip() for t in args.subgraph.split(",")]
            sub_onto = subgraph(onto, root_terms)
            
            subgraph_results = {
                "root_terms": root_terms,
                "n_terms_in_subgraph": sub_onto.num_terms(),
            }
            
            output_file = output_dir / "subgraph_analysis.json"
            dump_json(subgraph_results, output_file)
            workflow_results["analyses"]["subgraph"] = subgraph_results
            logger.info(f"Subgraph extracted: {sub_onto.num_terms()} terms")
        except Exception as e:
            logger.error(f"Subgraph extraction failed: {e}", exc_info=True)
            workflow_results["analyses"]["subgraph"] = {"error": str(e)}
    
    # Save summary
    summary_file = output_dir / "workflow_summary.json"
    dump_json(workflow_results, summary_file, indent=2)
    logger.info(f"Workflow summary saved to {summary_file}")
    
    logger.info("Workflow complete")


def main():
    """Main entry point."""
    args = parse_args()
    
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
