#!/usr/bin/env python3
"""Download genomes for species that are missing them.

This script:
1. Scans all amalgkit config files
2. Identifies species with missing genome downloads
3. Downloads genomes using the workflow's download mechanism
4. Logs progress and results

Usage:
    python3 scripts/rna/download_missing_genomes.py
    python3 scripts/rna/download_missing_genomes.py --species camponotus_floridanus
    python3 scripts/rna/download_missing_genomes.py --threads 4 --dry-run
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
from metainformant.dna.ncbi import download_genome_package_best_effort

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def format_file_size(size_bytes: int) -> str:
    """Format file size in human-readable format.
    
    Args:
        size_bytes: Size in bytes
        
    Returns:
        Formatted string (e.g., "1.5 MB", "2.3 GB")
    """
    if size_bytes < 1024:
        return f"{size_bytes} B"
    elif size_bytes < 1024 * 1024:
        return f"{size_bytes / 1024:.2f} KB"
    elif size_bytes < 1024 * 1024 * 1024:
        return f"{size_bytes / (1024 * 1024):.2f} MB"
    else:
        return f"{size_bytes / (1024 * 1024 * 1024):.2f} GB"


def check_genome_downloaded(genome_dir: Path, accession: str) -> bool:
    """Check if genome is already downloaded.
    
    Args:
        genome_dir: Genome directory path
        accession: NCBI assembly accession
        
    Returns:
        True if genome is downloaded, False otherwise
    """
    # Check download record
    download_record = genome_dir / "download_record.json"
    if download_record.exists():
        try:
            with open(download_record, "r") as f:
                record = json.load(f)
                if record.get("return_code") == 0:
                    return True
        except Exception:
            pass
    
    # Check if extracted directory exists
    extracted_dirs = [
        genome_dir / "ncbi_dataset_api_extracted",
        genome_dir / "ncbi_dataset_extracted",
    ]
    return any(d.exists() for d in extracted_dirs)


def download_genome_for_species(
    config_path: Path,
    repo_root: Path,
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Download genome for a single species.
    
    Args:
        config_path: Path to species config file
        repo_root: Repository root directory
        dry_run: If True, only report what would be done
        
    Returns:
        Dictionary with download results
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
        "downloaded": False,
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
    
    # Determine genome directory
    dest_dir_str = genome.get("dest_dir", "")
    if not dest_dir_str:
        work_dir_str = config.get("work_dir", "")
        if work_dir_str:
            work_dir_path = Path(work_dir_str).expanduser()
            if not work_dir_path.is_absolute():
                work_dir_path = repo_root / work_dir_path
            dest_dir = work_dir_path.parent / "genome"
        else:
            result["error"] = "Cannot determine genome directory"
            return result
    else:
        dest_dir = Path(dest_dir_str).expanduser()
        if not dest_dir.is_absolute():
            dest_dir = repo_root / dest_dir
    
    result["genome_dir"] = str(dest_dir)
    
    # Check if already downloaded
    if check_genome_downloaded(dest_dir, accession):
        logger.info(f"  ✓ {species_name}: Genome already downloaded")
        result["downloaded"] = True
        result["success"] = True
        return result
    
    if dry_run:
        logger.info(f"  [DRY RUN] Would download {species_name} ({accession}) to {dest_dir}")
        result["success"] = True
        result["would_download"] = True
        return result
    
    # Download genome
    logger.info(f"  Downloading {species_name} ({accession})...")
    logger.info(f"    Method: Best effort (CLI → API → FTP)")
    logger.info(f"    Include: {', '.join(genome.get('include', []))}")
    
    include = genome.get("include") or ["gff3", "rna", "cds", "protein", "genome", "seq-report"]
    ftp_url = genome.get("ftp_url")
    
    ensure_directory(dest_dir)
    
    try:
        import time
        start_time = time.time()
        
        dl_rec = download_genome_package_best_effort(
            accession,
            dest_dir,
            include=include,
            ftp_url=ftp_url,
        )
        
        elapsed_time = time.time() - start_time
        
        return_code = dl_rec.get("return_code", 1)
        result["return_code"] = return_code
        result["method"] = dl_rec.get("method", "unknown")
        result["extracted_dir"] = dl_rec.get("extracted_dir", "")
        result["elapsed_seconds"] = round(elapsed_time, 2)
        
        # Get file information
        zip_path_str = dl_rec.get("zip_path", "")
        if zip_path_str:
            zip_path = Path(zip_path_str)
            if zip_path.exists():
                zip_size = zip_path.stat().st_size
                result["zip_path"] = str(zip_path)
                result["zip_size_bytes"] = zip_size
                result["zip_size_mb"] = round(zip_size / (1024 * 1024), 2)
                logger.info(f"    Downloaded: {zip_path.name} ({format_file_size(zip_size)})")
        
        # Get extracted directory information
        extracted_dir_str = dl_rec.get("extracted_dir", "")
        if extracted_dir_str:
            extracted_dir = Path(extracted_dir_str)
            if extracted_dir.exists():
                # Calculate total size of extracted directory
                total_size = sum(f.stat().st_size for f in extracted_dir.rglob("*") if f.is_file())
                result["extracted_size_bytes"] = total_size
                result["extracted_size_mb"] = round(total_size / (1024 * 1024), 2)
                logger.info(f"    Extracted: {extracted_dir.name} ({format_file_size(total_size)})")
                
                # Find RNA FASTA file if present
                rna_files = list(extracted_dir.rglob("rna.fna*"))
                if rna_files:
                    rna_file = rna_files[0]
                    rna_size = rna_file.stat().st_size
                    result["rna_fasta_found"] = True
                    result["rna_fasta_path"] = str(rna_file)
                    result["rna_fasta_size_mb"] = round(rna_size / (1024 * 1024), 2)
                    logger.info(f"    RNA FASTA: {rna_file.name} ({format_file_size(rna_size)})")
        
        if return_code == 0:
            result["success"] = True
            result["downloaded"] = True
            logger.info(f"  ✓ {species_name}: Download successful ({result['method']}, {elapsed_time:.1f}s)")
        else:
            result["error"] = dl_rec.get("error", "Download failed")
            error_msg = result["error"]
            
            # Provide more helpful error messages
            if "not a zip file" in error_msg.lower():
                logger.error(f"  ✗ {species_name}: Download failed - Invalid file format")
                logger.error(f"    Error: {error_msg}")
                logger.error(f"    This may indicate the accession is invalid or the API returned an error page")
                logger.error(f"    Check: {dl_rec.get('url', 'N/A')}")
                logger.error(f"    Try: Verify accession {accession} is correct")
            else:
                logger.error(f"  ✗ {species_name}: Download failed - {error_msg}")
            
            # Log additional error details if available
            if dl_rec.get("stderr"):
                logger.debug(f"    stderr: {dl_rec['stderr'][:200]}")
        
        return result
        
    except Exception as e:
        result["error"] = str(e)
        logger.error(f"  ✗ {species_name}: Exception during download - {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return result


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Download genomes for species missing them",
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
        help="Process only this species (config filename without prefix, e.g., 'camponotus_floridanus')",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/genome_download_results.json",
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
        help="Only report what would be done, don't actually download",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of parallel downloads (not yet implemented, downloads are sequential)",
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
        logger.info("DRY RUN MODE - No downloads will be performed")
    logger.info("")
    
    results: list[dict[str, Any]] = []
    downloaded = 0
    already_downloaded = 0
    failed = 0
    total_size_mb = 0.0
    
    for idx, config_file in enumerate(config_files, 1):
        logger.info(f"[{idx}/{len(config_files)}] Processing {config_file.name}...")
        result = download_genome_for_species(
            config_file,
            repo_root,
            dry_run=args.dry_run,
        )
        results.append(result)
        
        # Track size
        if result.get("zip_size_mb"):
            total_size_mb += result["zip_size_mb"]
        
        if result.get("downloaded"):
            if result.get("would_download"):
                pass  # Dry run
            elif result.get("success"):
                already_downloaded += 1
            else:
                downloaded += 1
        elif result.get("success"):
            downloaded += 1
        else:
            failed += 1
        
        logger.info("")  # Blank line between species
    
    # Summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("DOWNLOAD SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total species: {len(results)}")
    logger.info(f"Already downloaded: {already_downloaded}")
    if args.dry_run:
        logger.info(f"Would download: {downloaded}")
    else:
        logger.info(f"Downloaded: {downloaded}")
        logger.info(f"Failed: {failed}")
        if total_size_mb > 0:
            total_bytes = int(total_size_mb * 1024 * 1024)
            logger.info(f"Total size downloaded: {format_file_size(total_bytes)}")
    logger.info("")
    
    # Show failed downloads
    if failed > 0:
        logger.info("Failed downloads:")
        for result in results:
            if not result.get("success") and not result.get("downloaded"):
                species = result.get("species_name", "unknown")
                accession = result.get("accession", "N/A")
                error = result.get("error", "Unknown error")
                logger.info(f"  ✗ {species} ({accession}): {error}")
        logger.info("")
    
    logger.info("=" * 80)
    
    # Write output
    ensure_directory(output_path.parent)
    output_data = {
        "summary": {
            "total": len(results),
            "already_downloaded": already_downloaded,
            "downloaded": downloaded if not args.dry_run else 0,
            "would_download": downloaded if args.dry_run else 0,
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

