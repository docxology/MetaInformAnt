#!/usr/bin/env python3
"""
Unified Species Processing Script.

Consolidates logic for processing RNA-seq samples using the core
metainformant.rna.engine.orchestrator.StreamingPipeline.

Usage:
    python3 scripts/rna/process_species.py \
        --species apis_mellifera \
        --strategy ena \
        --workers 20

Strategies:
    - ena: Direct download from ENA FTP/HTTP (fastest, no extraction).
    - sra: (Not yet implemented in new engine) Use local SRA Toolkit.
"""

import argparse
import sys
import logging
from pathlib import Path

# Add src to path for imports
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "src"))

from metainformant.core.utils.logging import get_logger, setup_logger
from metainformant.rna.engine.orchestrator import StreamingPipeline

logger = get_logger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="Process species RNA-seq samples.")
    parser.add_argument("--species", required=True, help="Species name (e.g., apis_mellifera)")
    parser.add_argument("--strategy", choices=["ena", "sra"], default="ena", help="Download strategy")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel workers")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without executing")
    
    parser.add_argument("--watchdog", action="store_true", help="Enable stalled process watchdog")
    parser.add_argument("--watchdog-timeout", type=int, default=3600, help="Watchdog stall timeout in seconds")
    
    # Optional overrides
    parser.add_argument("--work-dir", help="Override working directory")
    parser.add_argument("--index-file", help="Override index file path")
    parser.add_argument("--fastq-dir", help="Override download directory")
    
    return parser.parse_args()

def resolve_paths(species: str, args) -> dict:
    """Resolve paths based on species or overrides."""
    # Defaults for internal cluster
    base_dir = Path("output/amalgkit")
    
    # 1. Work Dir
    if args.work_dir:
        work_dir = Path(args.work_dir)
    else:
        work_dir = base_dir / f"{species}_all" / "work"
        
    # 2. Index File
    if args.index_file:
        index_file = Path(args.index_file)
    else:
        # Heuristic: work_dir/index/Species_name_transcripts.idx
        # Need capitalized species name? "Apis_mellifera"
        species_cap = species.capitalize() # "Apis_mellifera" if input is "apis_mellifera"?
        # Actually usually it's "Apis_mellifera" from "apis_mellifera"
        parts = species.split("_")
        species_cap = "_".join([p.capitalize() for p in parts])
        
        index_file = work_dir / "index" / f"{species_cap}_transcripts.idx"

    # 3. FASTQ Dir
    if args.fastq_dir:
        fastq_dir = Path(args.fastq_dir)
    else:
        # Default to external drive output/amalgkit/{species}
        fastq_dir = Path("output/amalgkit") / species

    return {
        "work_dir": work_dir,
        "index_file": index_file,
        "fastq_dir": fastq_dir
    }

def load_metadata_samples(work_dir: Path) -> list[str]:
    """Load sample IDs from Amalgkit metadata."""
    import csv
    metadata_file = work_dir / "metadata" / "metadata_selected.tsv"
    
    if not metadata_file.exists():
        logger.warning(f"Metadata not found at {metadata_file}. Using dummy list for dry-run if applicable.")
        return []
    
    samples = []
    try:
        with open(metadata_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if "run" not in (reader.fieldnames or []):
                return []
            for row in reader:
                if row.get("run"):
                    samples.append(row["run"])
    except Exception as e:
        logger.error(f"Failed to read metadata: {e}")
        
    return samples

def main():
    args = parse_args()
    setup_logger(__name__, level="INFO")
    
    logger.info(f"Species: {args.species}")
    logger.info(f"Strategy: {args.strategy}")
    
    paths = resolve_paths(args.species, args)
    logger.info(f"Work Dir: {paths['work_dir']}")
    logger.info(f"Index: {paths['index_file']}")
    logger.info(f"FASTQ Dir: {paths['fastq_dir']}")

    if args.strategy == "sra":
        logger.error("SRA strategy not yet implemented in StreamingPipeline.")
        sys.exit(1)

    # Initialize Engine
    pipeline = StreamingPipeline(
        species=args.species,
        index_file=paths["index_file"],
        work_dir=paths["work_dir"],
        fastq_dir=paths["fastq_dir"],
        workers=args.workers,
        dry_run=args.dry_run,
        use_watchdog=args.watchdog,
        watchdog_timeout=args.watchdog_timeout
    )

    # Load Samples
    # In a dry run without metadata, we might mock samples
    samples = load_metadata_samples(paths["work_dir"])
    
    if not samples:
        if args.dry_run:
            logger.info("[DRY RUN] No metadata found. Using mock sample list.")
            samples = ["SRR11092056", "SRR_MOCK_2"]
        else:
             logger.error("No samples found in metadata. Cannot proceed.")
             sys.exit(1)
             
    logger.info(f"Found {len(samples)} samples to process.")
    
    # Run
    pipeline.run(samples)

if __name__ == "__main__":
    main()
