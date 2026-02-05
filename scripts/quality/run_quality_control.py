#!/usr/bin/env python3
"""Quality control workflow orchestrator.

This script provides comprehensive orchestration for quality control workflows,
including FASTQ quality analysis, contamination detection, and quality reporting.

Usage:
    python3 scripts/quality/run_quality_control.py --fastq reads.fastq --output output/quality/report
    python3 scripts/quality/run_quality_control.py --fastq reads.fastq --analyze-fastq --detect-contamination
    python3 scripts/quality/run_quality_control.py --help
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Any

# Add project to path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.core.io import ensure_directory, dump_json
from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Quality control workflow orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # FASTQ quality analysis
  %(prog)s --fastq reads.fastq --analyze-fastq --output output/quality/fastq

  # Full QC with contamination detection
  %(prog)s --fastq reads.fastq --analyze-fastq --detect-contamination --generate-report

  # Subsampled analysis for large files
  %(prog)s --fastq large_file.fastq --analyze-fastq --subsample 10000
        """
    )
    parser.add_argument(
        "--fastq",
        type=Path,
        help="Input FASTQ file for quality analysis",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/quality"),
        help="Output directory (default: output/quality)",
    )
    parser.add_argument(
        "--analyze-fastq",
        action="store_true",
        help="Perform comprehensive FASTQ quality analysis",
    )
    parser.add_argument(
        "--detect-contamination",
        action="store_true",
        help="Detect sequence contamination",
    )
    parser.add_argument(
        "--generate-report",
        action="store_true",
        help="Generate comprehensive quality report",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        help="Subsample reads for analysis (for large files)",
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
    """Execute quality control workflow."""
    logger.info("Starting quality control workflow")
    logger.info(f"Output: {args.output}")
    
    if args.dry_run:
        logger.info("DRY RUN - no changes will be made")
        if args.fastq:
            logger.info(f"Would analyze: {args.fastq}")
        return
    
    output_dir = ensure_directory(args.output)
    logger.info(f"Output directory: {output_dir}")
    
    workflow_results = {
        "output_dir": str(output_dir),
        "analyses": {},
    }
    
    # FASTQ quality analysis
    if args.analyze_fastq and args.fastq:
        if not args.fastq.exists():
            raise FileNotFoundError(f"FASTQ file not found: {args.fastq}")
        
        try:
            logger.info(f"Analyzing FASTQ quality from {args.fastq}...")
            from metainformant.quality import analyze_fastq_quality
            
            quality_results = analyze_fastq_quality(
                args.fastq,
                subsample=args.subsample,
                seed=42
            )
            
            output_file = output_dir / "fastq_quality_analysis.json"
            dump_json(quality_results, output_file)
            workflow_results["analyses"]["fastq_quality"] = quality_results
            workflow_results["input_file"] = str(args.fastq)
            logger.info(f"FASTQ quality analysis saved to {output_file}")
        except Exception as e:
            logger.error(f"FASTQ analysis failed: {e}", exc_info=True)
            workflow_results["analyses"]["fastq_quality"] = {"error": str(e)}
    
    # Contamination detection
    if args.detect_contamination and args.fastq:
        try:
            logger.info("Detecting contamination...")
            from metainformant.quality import detect_adapter_contamination
            from metainformant.dna.fastq import iter_fastq
            
            # Read sequences
            sequences = []
            for read_id, seq, qual in iter_fastq(args.fastq):
                sequences.append(seq)
                if len(sequences) >= 10000:  # Limit for performance
                    break
            
            # Detect adapter contamination
            adapter_results = detect_adapter_contamination(sequences)
            
            output_file = output_dir / "contamination_detection.json"
            dump_json(adapter_results, output_file)
            workflow_results["analyses"]["contamination"] = adapter_results
            logger.info(f"Contamination detection saved to {output_file}")
        except Exception as e:
            logger.error(f"Contamination detection failed: {e}", exc_info=True)
            workflow_results["analyses"]["contamination"] = {"error": str(e)}
    
    # Generate quality report
    if args.generate_report:
        try:
            logger.info("Generating quality report...")
            from metainformant.quality import generate_quality_report
            
            # Collect quality data
            quality_data = {}
            if "fastq_quality" in workflow_results["analyses"]:
                fastq_data = workflow_results["analyses"]["fastq_quality"]
                if "quality_metrics" in fastq_data:
                    quality_data["quality"] = fastq_data["quality_metrics"]
                if "gc_metrics" in fastq_data:
                    quality_data["gc"] = fastq_data["gc_metrics"]
                if "length_metrics" in fastq_data:
                    quality_data["length"] = fastq_data["length_metrics"]
            
            if quality_data:
                report_text = generate_quality_report(quality_data, sample_name=args.fastq.stem if args.fastq else "Unknown")
                
                report_file = output_dir / "quality_report.txt"
                with open(report_file, 'w') as f:
                    f.write(report_text)
                logger.info(f"Quality report saved to {report_file}")
                workflow_results["analyses"]["report"] = {"file": str(report_file)}
            else:
                logger.warning("No quality data available for report generation")
        except Exception as e:
            logger.error(f"Report generation failed: {e}", exc_info=True)
            workflow_results["analyses"]["report"] = {"error": str(e)}
    
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
