#!/usr/bin/env python3
"""Prepare transcriptome FASTA files for species with downloaded genomes.

This script:
1. Scans all amalgkit config files
2. Identifies species with downloaded genomes but missing prepared transcriptomes
3. Extracts RNA FASTA from genome packages
4. Prepares them in the expected location for kallisto
5. Logs progress and results

Usage:
    python3 scripts/rna/prepare_transcriptomes.py
    python3 scripts/rna/prepare_transcriptomes.py --species camponotus_floridanus
    python3 scripts/rna/prepare_transcriptomes.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

# Add src directory to path for imports
script_dir = Path(__file__).parent
repo_root = script_dir.parent.parent
src_dir = repo_root / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from metainformant.core.config import load_mapping_from_file
from metainformant.core.io import dump_json, ensure_directory
from metainformant.rna.genome_prep import (
    find_rna_fasta_in_genome_dir,
    prepare_transcriptome_for_kallisto,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format."""
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.2f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.2f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB"


def check_transcriptome_prepared(work_dir: Path, species_name: str) -> bool:
    """Check if transcriptome FASTA is already prepared.
    
    Args:
        work_dir: Work directory for species
        species_name: Species name (with underscores)
        
    Returns:
        True if transcriptome is prepared, False otherwise
    """
    fasta_dir = work_dir / "fasta"
    expected_name = species_name.replace(" ", "_") + "_rna.fasta"
    fasta_path = fasta_dir / expected_name
    
    return fasta_path.exists() and fasta_path.stat().st_size > 0


def prepare_transcriptome_for_species(
    config_path: Path,
    repo_root: Path,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Prepare transcriptome for a single species.
    
    Args:
        config_path: Path to species config file
        repo_root: Repository root directory
        dry_run: If True, only report what would be done
        
    Returns:
        Dictionary with preparation results
    """
    try:
        config = load_mapping_from_file(config_path)
    except Exception as e:
        return {
            "config_file": str(config_path),
            "success": False,
            "error": f"Failed to load config: {e}",
        }
    
    species_name = config.get("species_list", ["unknown"])[0]
    genome = config.get("genome", {})
    
    result: dict[str, Any] = {
        "config_file": str(config_path),
        "species_name": species_name,
        "success": False,
        "prepared": False,
        "dry_run": dry_run,
    }
    
    if not genome:
        result["error"] = "No genome configuration found"
        return result
    
    accession = genome.get("accession")
    if not accession:
        result["error"] = "No genome accession in config"
        return result
    
    result["accession"] = accession
    
    # Determine directories
    dest_dir_str = genome.get("dest_dir", "")
    if not dest_dir_str:
        work_dir_str = config.get("work_dir", "")
        if work_dir_str:
            work_dir_path = Path(work_dir_str).expanduser()
            if not work_dir_path.is_absolute():
                work_dir_path = repo_root / work_dir_path
            dest_dir = work_dir_path.parent / "genome"
            work_dir = work_dir_path
        else:
            result["error"] = "Cannot determine directories"
            return result
    else:
        dest_dir = Path(dest_dir_str).expanduser()
        if not dest_dir.is_absolute():
            dest_dir = repo_root / dest_dir
        
        work_dir_str = config.get("work_dir", "")
        if work_dir_str:
            work_dir = Path(work_dir_str).expanduser()
            if not work_dir.is_absolute():
                work_dir = repo_root / work_dir
        else:
            result["error"] = "Cannot determine work directory"
            return result
    
    result["genome_dir"] = str(dest_dir)
    result["work_dir"] = str(work_dir)
    
    # Check if transcriptome already prepared
    if check_transcriptome_prepared(work_dir, species_name):
        logger.info(f"  ✓ {species_name}: Transcriptome already prepared")
        result["prepared"] = True
        result["already_prepared"] = True
        result["success"] = True
        return result
    
    # Check if genome is downloaded (will check for RNA FASTA and CDS fallback in prepare_transcriptome_for_kallisto)
    # Just verify genome directory exists
    if not dest_dir.exists():
        result["error"] = "Genome directory not found"
        logger.warning(f"  ⚠ {species_name}: {result['error']}")
        return result
    
    if dry_run:
        expected_fasta = work_dir / "fasta" / f"{species_name.replace(' ', '_')}_rna.fasta"
        logger.info(f"  [DRY RUN] Would prepare transcriptome for {species_name}")
        logger.info(f"    Will search for RNA FASTA, fallback to CDS if not found")
        logger.info(f"    Destination: {expected_fasta}")
        result["success"] = True
        result["would_prepare"] = True
        return result
    
    # Prepare transcriptome (will try RNA FASTA first, then CDS fallback)
    logger.info(f"  Preparing transcriptome for {species_name}...")
    # Log resolved work_dir path for clarity
    logger.debug(f"    Work directory: {work_dir}")
    expected_output = work_dir / "fasta" / f"{species_name.replace(' ', '_')}_rna.fasta"
    logger.debug(f"    Expected output: {expected_output}")
    
    try:
        import time
        start_time = time.time()
        
        fasta_path = prepare_transcriptome_for_kallisto(
            dest_dir,
            species_name,
            work_dir,
            accession=accession,
            use_cds_fallback=True,  # Enable CDS fallback
        )
        
        elapsed_time = time.time() - start_time
        
        if fasta_path:
            result["success"] = True
            result["prepared"] = True
            result["newly_prepared"] = True
            result["fasta_path"] = str(fasta_path)
            result["elapsed_seconds"] = round(elapsed_time, 2)
            
            # Get source file info (for logging)
            source_file = fasta_path  # Will be updated from prepare_transcriptome_for_kallisto if needed
            result["source_rna_fasta"] = str(source_file)
            
            # Get output file size
            output_size = fasta_path.stat().st_size
            output_size_mb = round(output_size / (1024 * 1024), 2)
            result["fasta_size_bytes"] = output_size
            result["fasta_size_mb"] = output_size_mb
            
            logger.info(f"    Destination: {fasta_path.name}")
            logger.debug(f"    Destination path: {fasta_path}")
            logger.info(f"    Output size: {format_file_size(output_size)}")
            logger.info(f"  ✓ {species_name}: Transcriptome prepared ({elapsed_time:.1f}s)")
        else:
            result["error"] = "Failed to prepare transcriptome (RNA FASTA and CDS not found)"
            logger.error(f"  ✗ {species_name}: {result['error']}")
        
        return result
        
    except Exception as e:
        result["error"] = str(e)
        logger.error(f"  ✗ {species_name}: Exception during preparation - {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return result


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Prepare transcriptome FASTA files for species with downloaded genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config-dir",
        type=str,
        default="config/amalgkit",
        help="Directory containing amalgkit config files",
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Process only this species (config filename without prefix)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/transcriptome_preparation_results.json",
        help="Output JSON file path",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=".",
        help="Repository root directory",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only report what would be done, don't actually prepare",
    )
    
    args = parser.parse_args()
    
    # Determine repo root from script location if not provided
    if args.repo_root == ".":
        script_dir = Path(__file__).parent
        repo_root = script_dir.parent.parent.resolve()
    else:
        repo_root = Path(args.repo_root).resolve()
    
    config_dir = repo_root / args.config_dir
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = repo_root / output_path
    
    if not config_dir.exists():
        logger.error(f"Config directory not found: {config_dir}")
        sys.exit(1)
    
    # Find config files
    if args.species:
        config_files = [config_dir / f"amalgkit_{args.species}.yaml"]
        config_files = [f for f in config_files if f.exists()]
    else:
        config_files = sorted(config_dir.glob("amalgkit_*.yaml"))
        config_files = [
            f for f in config_files
            if f.name not in {"amalgkit_template.yaml", "amalgkit_test.yaml"}
        ]
    
    if not config_files:
        logger.error("No config files found")
        sys.exit(1)
    
    logger.info(f"Found {len(config_files)} config file(s)")
    if args.dry_run:
        logger.info("DRY RUN MODE - No transcriptomes will be prepared")
    logger.info("")
    
    results: list[dict[str, Any]] = []
    prepared = 0
    already_prepared = 0
    failed = 0
    skipped = 0
    total_size_mb = 0.0
    
    for idx, config_file in enumerate(config_files, 1):
        logger.info(f"[{idx}/{len(config_files)}] Processing {config_file.name}...")
        result = prepare_transcriptome_for_species(
            config_file,
            repo_root,
            dry_run=args.dry_run,
        )
        results.append(result)
        
        # Track size
        if result.get("fasta_size_mb"):
            total_size_mb += result["fasta_size_mb"]
        
        # Count results correctly
        if result.get("would_prepare"):
            prepared += 1  # Dry run - would prepare
        elif result.get("already_prepared"):
            already_prepared += 1
        elif result.get("newly_prepared"):
            prepared += 1  # Actually prepared in this run
        elif result.get("success"):
            # Success but not clearly categorized
            if result.get("elapsed_seconds"):
                prepared += 1  # Likely newly prepared
            else:
                already_prepared += 1  # Likely already existed
        elif result.get("error") and "not downloaded" in result.get("error", "").lower():
            skipped += 1
        else:
            failed += 1
    
    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("PREPARATION SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total species: {len(results)}")
    logger.info(f"Already prepared: {already_prepared}")
    if args.dry_run:
        logger.info(f"Would prepare: {prepared}")
    else:
        logger.info(f"Prepared: {prepared}")
        if total_size_mb > 0:
            logger.info(f"Total size: {total_size_mb:.2f} MB")
    logger.info(f"Skipped (no genome): {skipped}")
    logger.info(f"Failed: {failed}")
    logger.info("=" * 80)
    
    # Write output
    ensure_directory(output_path.parent)
    output_data = {
        "summary": {
            "total": len(results),
            "already_prepared": already_prepared,
            "prepared": prepared if not args.dry_run else 0,
            "would_prepare": prepared if args.dry_run else 0,
            "skipped": skipped,
            "failed": failed,
        },
        "results": results,
    }
    dump_json(output_data, output_path, indent=2)
    logger.info(f"\nResults written to: {output_path}")
    
    if failed > 0:
        sys.exit(1)


if __name__ == "__main__":
    main()

