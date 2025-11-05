#!/usr/bin/env python3
"""Build kallisto indexes for species with prepared transcriptomes.

This script:
1. Scans all amalgkit config files
2. Identifies species with prepared transcriptomes but missing kallisto indexes
3. Builds kallisto indexes from transcriptome FASTA files
4. Logs progress and results

Usage:
    python3 scripts/rna/build_kallisto_indexes.py
    python3 scripts/rna/build_kallisto_indexes.py --species camponotus_floridanus
    python3 scripts/rna/build_kallisto_indexes.py --kmer-size 31 --dry-run
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
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
    build_kallisto_index,
    get_expected_index_path,
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


def find_transcriptome_fasta(work_dir: Path, species_name: str) -> Path | None:
    """Find transcriptome FASTA file for species.
    
    Args:
        work_dir: Work directory for species
        species_name: Species name (with underscores)
        
    Returns:
        Path to FASTA file if found, None otherwise
    """
    fasta_dir = work_dir / "fasta"
    if not fasta_dir.exists():
        return None
    
    # Try expected name
    expected_name = species_name.replace(" ", "_") + "_rna.fasta"
    fasta_path = fasta_dir / expected_name
    if fasta_path.exists():
        return fasta_path
    
    # Search for any matching FASTA files
    pattern = species_name.replace(" ", "_").replace("_", "*")
    matches = list(fasta_dir.glob(f"{pattern}*rna*.fasta"))
    if matches:
        return matches[0]
    
    # Search for any .fasta files
    all_fastas = list(fasta_dir.glob("*.fasta"))
    if len(all_fastas) == 1:
        return all_fastas[0]
    
    return None


def build_index_for_species(
    config_path: Path,
    repo_root: Path,
    *,
    kmer_size: int = 31,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Build kallisto index for a single species.
    
    Args:
        config_path: Path to species config file
        repo_root: Repository root directory
        kmer_size: k-mer size for kallisto index
        dry_run: If True, only report what would be done
        
    Returns:
        Dictionary with build results
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
    
    result: dict[str, Any] = {
        "config_file": str(config_path),
        "species_name": species_name,
        "success": False,
        "index_built": False,
        "dry_run": dry_run,
        "kmer_size": kmer_size,
    }
    
    # Check if kallisto is available
    kallisto_exe = shutil.which("kallisto")
    if not kallisto_exe:
        result["error"] = "kallisto not found on PATH"
        logger.error(f"  ✗ {species_name}: {result['error']}")
        return result
    
    result["kallisto_path"] = kallisto_exe
    
    # Determine work directory
    work_dir_str = config.get("work_dir", "")
    if not work_dir_str:
        result["error"] = "Cannot determine work directory"
        return result
    
    work_dir = Path(work_dir_str).expanduser()
    if not work_dir.is_absolute():
        work_dir = repo_root / work_dir
    
    result["work_dir"] = str(work_dir)
    
    # Check if index already exists
    index_path = get_expected_index_path(work_dir, species_name)
    if index_path.exists():
        logger.info(f"  ✓ {species_name}: Index already exists at {index_path}")
        result["index_built"] = True
        result["success"] = True
        result["index_path"] = str(index_path)
        return result
    
    result["index_path"] = str(index_path)
    
    # Find transcriptome FASTA
    fasta_path = find_transcriptome_fasta(work_dir, species_name)
    if not fasta_path:
        result["error"] = "Transcriptome FASTA not found"
        logger.warning(f"  ⚠ {species_name}: {result['error']}")
        return result
    
    result["fasta_path"] = str(fasta_path)
    
    if dry_run:
        logger.info(f"  [DRY RUN] Would build kallisto index for {species_name}")
        logger.info(f"    FASTA: {fasta_path}")
        logger.info(f"    Index: {index_path}")
        logger.info(f"    k-mer size: {kmer_size}")
        result["success"] = True
        result["would_build"] = True
        return result
    
    # Build index
    logger.info(f"  Building kallisto index for {species_name} (k={kmer_size})...")
    logger.info(f"    FASTA: {fasta_path.name}")
    
    # Get FASTA file size
    fasta_size = fasta_path.stat().st_size
    fasta_size_mb = round(fasta_size / (1024 * 1024), 2)
    logger.info(f"    FASTA size: {format_file_size(fasta_size)}")
    logger.info(f"    k-mer size: {kmer_size}")
    
    try:
        import time
        start_time = time.time()
        
        success = build_kallisto_index(
            fasta_path,
            index_path,
            kmer_size=kmer_size,
            check_existing=False,  # Already checked above
        )
        
        elapsed_time = time.time() - start_time
        
        if success and index_path.exists():
            result["success"] = True
            result["index_built"] = True
            result["elapsed_seconds"] = round(elapsed_time, 2)
            
            # Get index file size
            index_size = index_path.stat().st_size
            index_size_mb = round(index_size / (1024 * 1024), 2)
            result["index_size_bytes"] = index_size
            result["index_size_mb"] = index_size_mb
            
            logger.info(f"    Index: {index_path.name}")
            logger.info(f"    Index size: {format_file_size(index_size)}")
            logger.info(f"  ✓ {species_name}: Index built successfully ({elapsed_time:.1f}s)")
        else:
            result["error"] = "Failed to build index" if not success else "Index file not created"
            logger.error(f"  ✗ {species_name}: {result['error']}")
        
        return result
        
    except Exception as e:
        result["error"] = str(e)
        logger.error(f"  ✗ {species_name}: Exception during index building - {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return result


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Build kallisto indexes for species with prepared transcriptomes",
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
        default="output/kallisto_index_build_results.json",
        help="Output JSON file path",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=".",
        help="Repository root directory",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=31,
        help="k-mer size for kallisto index (default: 31, use 23 for short reads)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only report what would be done, don't actually build indexes",
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
    
    # Check kallisto availability
    kallisto_exe = shutil.which("kallisto")
    if not kallisto_exe:
        logger.error("kallisto not found on PATH. Please install kallisto first.")
        logger.error("  Install: conda install -c bioconda kallisto")
        sys.exit(1)
    
    logger.info(f"Using kallisto: {kallisto_exe}")
    
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
    logger.info(f"k-mer size: {args.kmer_size}")
    if args.dry_run:
        logger.info("DRY RUN MODE - No indexes will be built")
    logger.info("")
    
    results: list[dict[str, Any]] = []
    built = 0
    already_built = 0
    failed = 0
    skipped = 0
    total_size_mb = 0.0
    
    for idx, config_file in enumerate(config_files, 1):
        logger.info(f"[{idx}/{len(config_files)}] Processing {config_file.name}...")
        result = build_index_for_species(
            config_file,
            repo_root,
            kmer_size=args.kmer_size,
            dry_run=args.dry_run,
        )
        results.append(result)
        
        # Track size
        if result.get("index_size_mb"):
            total_size_mb += result["index_size_mb"]
        
        logger.info("")  # Blank line between species
        
        if result.get("index_built"):
            if result.get("would_build"):
                built += 1
            elif result.get("success"):
                already_built += 1
        elif result.get("success"):
            built += 1
        elif result.get("error") and "not found" in result.get("error", "").lower():
            skipped += 1
        else:
            failed += 1
    
    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("INDEX BUILD SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total species: {len(results)}")
    logger.info(f"Already built: {already_built}")
    if args.dry_run:
        logger.info(f"Would build: {built}")
    else:
        logger.info(f"Built: {built}")
        if total_size_mb > 0:
            logger.info(f"Total index size: {total_size_mb:.2f} MB")
    logger.info(f"Skipped (no transcriptome): {skipped}")
    logger.info(f"Failed: {failed}")
    logger.info("=" * 80)
    
    # Write output
    ensure_directory(output_path.parent)
    output_data = {
        "summary": {
            "total": len(results),
            "already_built": already_built,
            "built": built if not args.dry_run else 0,
            "would_build": built if args.dry_run else 0,
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

