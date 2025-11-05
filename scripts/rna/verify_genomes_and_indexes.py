#!/usr/bin/env python3
"""Verify genome downloads and kallisto indexes for all amalgkit species configs.

This script scans all amalgkit configuration files and checks:
1. Genome download status
2. RNA FASTA file availability
3. Kallisto index existence

Usage:
    python3 scripts/rna/verify_genomes_and_indexes.py
    python3 scripts/rna/verify_genomes_and_indexes.py --output output/genome_status.json
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

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def find_rna_fasta_in_genome_dir(genome_dir: Path, accession: str) -> Path | None:
    """Find RNA FASTA file in extracted genome directory.
    
    NCBI datasets extracts to:
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/rna.fna.gz
    - ncbi_dataset_api_extracted/ncbi_dataset/data/{accession}/transcript_sequences.fna
    
    Also checks for alternative paths and naming conventions.
    
    Args:
        genome_dir: Base genome directory
        accession: NCBI assembly accession
        
    Returns:
        Path to RNA FASTA file if found, None otherwise
    """
    # Common extraction patterns
    patterns = [
        # API extraction pattern
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna",
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna.gz",
        genome_dir / "ncbi_dataset_api_extracted" / "ncbi_dataset" / "data" / accession / "transcript_sequences.fna",
        # CLI extraction pattern
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna",
        genome_dir / "ncbi_dataset_extracted" / "ncbi_dataset" / "data" / accession / "rna.fna.gz",
        # Direct file patterns
        genome_dir / "rna.fna",
        genome_dir / "rna.fna.gz",
        # Check for any rna*.fna files
    ]
    
    for pattern in patterns:
        if pattern.exists() and pattern.is_file():
            return pattern
    
    # Search recursively for RNA FASTA files
    for pattern in ["rna.fna", "rna.fna.gz", "transcript_sequences.fna", "*_rna_from_genomic.fna.gz"]:
        matches = list(genome_dir.rglob(pattern))
        if matches:
            # Prefer unzipped files
            unzipped = [m for m in matches if not m.name.endswith(".gz")]
            if unzipped:
                return unzipped[0]
            return matches[0]
    
    return None


def get_species_name_from_config(config: dict[str, Any]) -> str:
    """Extract species name from config.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Species name (with underscores)
    """
    species_list = config.get("species_list", [])
    if species_list:
        return species_list[0]
    
    # Fallback: try to infer from genome accession or other fields
    genome = config.get("genome", {})
    if genome:
        accession = genome.get("accession", "")
        if accession:
            # Could parse from accession, but better to use species_list
            pass
    
    return "unknown"


def check_kallisto_index(work_dir: Path, species_name: str) -> Path | None:
    """Check if kallisto index exists for species.
    
    Expected location: work_dir/index/{Species_Name}_transcripts.idx
    
    Args:
        work_dir: Work directory for species
        species_name: Species name (with underscores)
        
    Returns:
        Path to index file if found, None otherwise
    """
    index_dir = work_dir / "index"
    if not index_dir.exists():
        return None
    
    # Expected naming: Species_Name_transcripts.idx
    index_file = index_dir / f"{species_name}_transcripts.idx"
    if index_file.exists():
        return index_file
    
    # Check for any .idx files matching species name pattern
    pattern = species_name.replace("_", "*")
    matches = list(index_dir.glob(f"{pattern}*.idx"))
    if matches:
        return matches[0]
    
    # Check for any .idx files (fallback)
    all_indices = list(index_dir.glob("*.idx"))
    if len(all_indices) == 1:
        # Single index file - likely the one we want
        return all_indices[0]
    
    return None


def verify_species_config(config_path: Path, repo_root: Path) -> dict[str, Any]:
    """Verify genome and index status for a single species config.
    
    Args:
        config_path: Path to amalgkit config file
        repo_root: Repository root directory
        
    Returns:
        Dictionary with verification results
    """
    try:
        config = load_mapping_from_file(config_path)
    except Exception as e:
        return {
            "config_file": str(config_path),
            "error": f"Failed to load config: {e}",
            "genome_downloaded": False,
            "rna_fasta_found": False,
            "kallisto_index_found": False,
        }
    
    species_name = get_species_name_from_config(config)
    genome = config.get("genome", {})
    
    result: dict[str, Any] = {
        "config_file": str(config_path),
        "species_name": species_name,
        "genome_downloaded": False,
        "rna_fasta_found": False,
        "rna_fasta_path": None,
        "kallisto_index_found": False,
        "kallisto_index_path": None,
        "genome_accession": None,
        "work_dir": None,
    }
    
    if not genome:
        result["error"] = "No genome configuration found"
        return result
    
    accession = genome.get("accession")
    result["genome_accession"] = accession
    
    if not accession:
        result["error"] = "No genome accession in config"
        return result
    
    # Check genome download
    dest_dir_str = genome.get("dest_dir", "")
    if not dest_dir_str:
        # Infer from work_dir
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
    
    # Check download record
    download_record = dest_dir / "download_record.json"
    if download_record.exists():
        try:
            with open(download_record, "r") as f:
                record = json.load(f)
                if record.get("return_code") == 0:
                    result["genome_downloaded"] = True
        except Exception:
            pass
    
    # Check if extracted directory exists (alternative check)
    extracted_dirs = [
        dest_dir / "ncbi_dataset_api_extracted",
        dest_dir / "ncbi_dataset_extracted",
    ]
    if any(d.exists() for d in extracted_dirs):
        result["genome_downloaded"] = True
    
    # Find RNA FASTA
    rna_fasta = find_rna_fasta_in_genome_dir(dest_dir, accession)
    if rna_fasta:
        result["rna_fasta_found"] = True
        result["rna_fasta_path"] = str(rna_fasta)
    
    # Check kallisto index
    work_dir_str = config.get("work_dir", "")
    if work_dir_str:
        work_dir = Path(work_dir_str).expanduser()
        if not work_dir.is_absolute():
            work_dir = repo_root / work_dir
        result["work_dir"] = str(work_dir)
        
        index_path = check_kallisto_index(work_dir, species_name)
        if index_path:
            result["kallisto_index_found"] = True
            result["kallisto_index_path"] = str(index_path)
    
    return result


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(description="Verify genome downloads and kallisto indexes")
    parser.add_argument(
        "--config-dir",
        type=str,
        default="config/amalgkit",
        help="Directory containing amalgkit config files",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="output/genome_index_status.json",
        help="Output JSON file path",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=".",
        help="Repository root directory",
    )
    
    args = parser.parse_args()
    
    # Determine repo root from script location if not provided
    if args.repo_root == ".":
        script_dir = Path(__file__).parent
        repo_root = script_dir.parent.parent.resolve()
    else:
        repo_root = Path(args.repo_root).resolve()
    
    config_dir = repo_root / args.config_dir
    # Handle relative output paths
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = repo_root / output_path
    
    if not config_dir.exists():
        logger.error(f"Config directory not found: {config_dir}")
        return
    
    # Find all config files
    config_files = sorted(config_dir.glob("amalgkit_*.yaml"))
    config_files = [f for f in config_files if f.name != "amalgkit_template.yaml" and f.name != "amalgkit_test.yaml"]
    
    logger.info(f"Found {len(config_files)} config file(s)")
    logger.info("")
    
    results: list[dict[str, Any]] = []
    total = len(config_files)
    
    for idx, config_file in enumerate(config_files, 1):
        logger.info(f"[{idx}/{total}] Checking {config_file.name}...")
        result = verify_species_config(config_file, repo_root)
        results.append(result)
        
        # Print detailed summary
        species_name = result.get("species_name", "unknown")
        status_parts = []
        
        if result.get("genome_downloaded"):
            status_parts.append("✓ Genome")
            if result.get("genome_dir"):
                logger.debug(f"    Genome dir: {result['genome_dir']}")
        else:
            status_parts.append("✗ Genome")
            if result.get("genome_accession"):
                logger.debug(f"    Missing: {result['genome_accession']}")
        
        if result.get("rna_fasta_found"):
            status_parts.append("✓ RNA FASTA")
            if result.get("rna_fasta_path"):
                logger.debug(f"    RNA FASTA: {result['rna_fasta_path']}")
        else:
            status_parts.append("✗ RNA FASTA")
        
        if result.get("kallisto_index_found"):
            status_parts.append("✓ Index")
            if result.get("kallisto_index_path"):
                logger.debug(f"    Index: {result['kallisto_index_path']}")
        else:
            status_parts.append("✗ Index")
        
        logger.info(f"  {species_name}: {' | '.join(status_parts)}")
        
        if result.get("error"):
            logger.warning(f"    Warning: {result['error']}")
        
        logger.info("")  # Blank line between species
    
    # Generate detailed summary
    total = len(results)
    genomes_downloaded = sum(1 for r in results if r.get("genome_downloaded"))
    rna_fasta_found = sum(1 for r in results if r.get("rna_fasta_found"))
    kallisto_indexes_found = sum(1 for r in results if r.get("kallisto_index_found"))
    
    # Calculate missing counts
    genomes_missing = total - genomes_downloaded
    rna_fasta_missing = genomes_downloaded - rna_fasta_found
    indexes_missing = rna_fasta_found - kallisto_indexes_found
    
    # Calculate percentages
    genome_pct = (genomes_downloaded / total * 100) if total > 0 else 0
    rna_pct = (rna_fasta_found / genomes_downloaded * 100) if genomes_downloaded > 0 else 0
    index_pct = (kallisto_indexes_found / rna_fasta_found * 100) if rna_fasta_found > 0 else 0
    
    summary = {
        "total_species": total,
        "genomes_downloaded": genomes_downloaded,
        "genomes_missing": genomes_missing,
        "rna_fasta_found": rna_fasta_found,
        "rna_fasta_missing": rna_fasta_missing,
        "kallisto_indexes_found": kallisto_indexes_found,
        "indexes_missing": indexes_missing,
        "percentages": {
            "genomes": round(genome_pct, 1),
            "rna_fasta": round(rna_pct, 1),
            "indexes": round(index_pct, 1),
        },
    }
    
    logger.info("=" * 80)
    logger.info("VERIFICATION SUMMARY")
    logger.info("=" * 80)
    logger.info(f"Total species: {total}")
    logger.info("")
    logger.info("Genome Downloads:")
    logger.info(f"  ✓ Downloaded: {genomes_downloaded}/{total} ({genome_pct:.1f}%)")
    logger.info(f"  ✗ Missing: {genomes_missing}/{total}")
    logger.info("")
    logger.info("RNA FASTA Files (of downloaded genomes):")
    logger.info(f"  ✓ Found: {rna_fasta_found}/{genomes_downloaded} ({rna_pct:.1f}%)")
    logger.info(f"  ✗ Missing: {rna_fasta_missing}/{genomes_downloaded}")
    logger.info("")
    logger.info("Kallisto Indexes (of prepared transcriptomes):")
    logger.info(f"  ✓ Built: {kallisto_indexes_found}/{rna_fasta_found} ({index_pct:.1f}%)")
    logger.info(f"  ✗ Missing: {indexes_missing}/{rna_fasta_found}")
    logger.info("")
    logger.info("=" * 80)
    logger.info("")
    logger.info("Next Steps:")
    if genomes_missing > 0:
        logger.info(f"  → Download {genomes_missing} missing genomes:")
        logger.info(f"     python3 scripts/rna/download_missing_genomes.py")
    if rna_fasta_missing > 0:
        logger.info(f"  → Prepare {rna_fasta_missing} transcriptomes:")
        logger.info(f"     python3 scripts/rna/prepare_transcriptomes.py")
    if indexes_missing > 0:
        logger.info(f"  → Build {indexes_missing} kallisto indexes:")
        logger.info(f"     python3 scripts/rna/build_kallisto_indexes.py")
    if genomes_missing == 0 and rna_fasta_missing == 0 and indexes_missing == 0:
        logger.info("  ✓ All species are fully prepared!")
    logger.info("")
    
    # Write output
    ensure_directory(output_path.parent)
    output_data = {
        "summary": summary,
        "results": results,
    }
    dump_json(output_data, output_path, indent=2)
    logger.info(f"\nResults written to: {output_path}")


if __name__ == "__main__":
    main()

