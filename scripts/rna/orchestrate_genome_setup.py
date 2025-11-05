#!/usr/bin/env python3
"""Master orchestrator for complete genome setup pipeline.

This script orchestrates the complete genome setup process:
1. Verify current status
2. Download missing genomes
3. Prepare transcriptomes
4. Build kallisto indexes
5. Verify final status

All steps are logged and progress is tracked.

Usage:
    python3 scripts/rna/orchestrate_genome_setup.py
    python3 scripts/rna/orchestrate_genome_setup.py --species camponotus_floridanus
    python3 scripts/rna/orchestrate_genome_setup.py --skip-download --skip-prepare
    python3 scripts/rna/orchestrate_genome_setup.py --dry-run
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def run_script(script_name: str, args: list[str], repo_root: Path) -> int:
    """Run a script and return its exit code.
    
    Args:
        script_name: Name of script to run (without .py extension)
        args: Additional arguments to pass to script
        repo_root: Repository root directory
        
    Returns:
        Exit code from script
    """
    script_path = repo_root / "scripts" / "rna" / f"{script_name}.py"
    
    if not script_path.exists():
        logger.error(f"Script not found: {script_path}")
        return 1
    
    # Ensure script path is absolute and exists
    script_path = script_path.resolve()
    if not script_path.exists():
        logger.error(f"Script not found: {script_path}")
        return 1
    
    cmd = [sys.executable, str(script_path)] + args
    logger.info(f"Running: {' '.join(cmd)}")
    
    try:
        # Run from repo root to ensure relative paths work correctly
        result = subprocess.run(
            cmd,
            cwd=str(repo_root),
            check=False,
            capture_output=False,  # Let output stream to console
        )
        return result.returncode
    except Exception as e:
        logger.error(f"Error running script: {e}")
        import traceback
        logger.debug(traceback.format_exc())
        return 1


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Orchestrate complete genome setup pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--species",
        type=str,
        help="Process only this species (config filename without prefix)",
    )
    parser.add_argument(
        "--kmer-size",
        type=int,
        default=31,
        help="k-mer size for kallisto index (default: 31)",
    )
    parser.add_argument(
        "--skip-verify-initial",
        action="store_true",
        help="Skip initial verification step",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip genome download step",
    )
    parser.add_argument(
        "--skip-prepare",
        action="store_true",
        help="Skip transcriptome preparation step",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip kallisto index building step",
    )
    parser.add_argument(
        "--skip-verify-final",
        action="store_true",
        help="Skip final verification step",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run mode (don't actually perform operations)",
    )
    parser.add_argument(
        "--repo-root",
        type=str,
        default=".",
        help="Repository root directory",
    )
    
    args = parser.parse_args()
    
    # Determine repo root
    if args.repo_root == ".":
        script_dir = Path(__file__).parent
        repo_root = script_dir.parent.parent.resolve()
    else:
        repo_root = Path(args.repo_root).resolve()
    
    logger.info("=" * 80)
    logger.info("Genome Setup Orchestrator")
    logger.info("=" * 80)
    logger.info(f"Repository root: {repo_root}")
    if args.species:
        logger.info(f"Species filter: {args.species}")
    if args.dry_run:
        logger.info("DRY RUN MODE")
    logger.info("")
    
    # Build common arguments
    common_args = []
    if args.species:
        common_args.extend(["--species", args.species])
    if args.dry_run:
        common_args.append("--dry-run")
    
    exit_code = 0
    
    # Step 1: Initial verification
    if not args.skip_verify_initial:
        logger.info("=" * 80)
        logger.info("Step 1: Initial Verification")
        logger.info("=" * 80)
        verify_args = common_args + ["--output", "output/genome_index_status_initial.json"]
        code = run_script("verify_genomes_and_indexes", verify_args, repo_root)
        if code != 0:
            logger.warning(f"Verification script exited with code {code}")
        logger.info("")
    
    # Step 2: Download missing genomes
    if not args.skip_download:
        logger.info("=" * 80)
        logger.info("Step 2: Download Missing Genomes")
        logger.info("=" * 80)
        download_args = common_args + ["--output", "output/genome_download_results.json"]
        code = run_script("download_missing_genomes", download_args, repo_root)
        if code != 0:
            logger.error(f"Download script exited with code {code}")
            exit_code = max(exit_code, code)
        logger.info("")
    else:
        logger.info("Skipping genome download step")
        logger.info("")
    
    # Step 3: Prepare transcriptomes
    if not args.skip_prepare:
        logger.info("=" * 80)
        logger.info("Step 3: Prepare Transcriptomes")
        logger.info("=" * 80)
        prepare_args = common_args + ["--output", "output/transcriptome_preparation_results.json"]
        code = run_script("prepare_transcriptomes", prepare_args, repo_root)
        if code != 0:
            logger.error(f"Prepare script exited with code {code}")
            exit_code = max(exit_code, code)
        logger.info("")
    else:
        logger.info("Skipping transcriptome preparation step")
        logger.info("")
    
    # Step 4: Build kallisto indexes
    if not args.skip_build:
        logger.info("=" * 80)
        logger.info("Step 4: Build Kallisto Indexes")
        logger.info("=" * 80)
        build_args = common_args + [
            "--kmer-size", str(args.kmer_size),
            "--output", "output/kallisto_index_build_results.json",
        ]
        code = run_script("build_kallisto_indexes", build_args, repo_root)
        if code != 0:
            logger.error(f"Build script exited with code {code}")
            exit_code = max(exit_code, code)
        logger.info("")
    else:
        logger.info("Skipping kallisto index building step")
        logger.info("")
    
    # Step 5: Final verification
    if not args.skip_verify_final:
        logger.info("=" * 80)
        logger.info("Step 5: Final Verification")
        logger.info("=" * 80)
        verify_args = common_args + ["--output", "output/genome_index_status_final.json"]
        code = run_script("verify_genomes_and_indexes", verify_args, repo_root)
        if code != 0:
            logger.warning(f"Verification script exited with code {code}")
        logger.info("")
    else:
        logger.info("Skipping final verification step")
        logger.info("")
    
    # Summary
    logger.info("=" * 80)
    logger.info("Pipeline Complete")
    logger.info("=" * 80)
    logger.info(f"Exit code: {exit_code}")
    logger.info("")
    logger.info("Results files:")
    logger.info("  - output/genome_index_status_initial.json (if initial verify ran)")
    logger.info("  - output/genome_download_results.json (if download ran)")
    logger.info("  - output/transcriptome_preparation_results.json (if prepare ran)")
    logger.info("  - output/kallisto_index_build_results.json (if build ran)")
    logger.info("  - output/genome_index_status_final.json (if final verify ran)")
    logger.info("")
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()

