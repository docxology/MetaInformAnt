#!/usr/bin/env python3
"""Process samples sequentially: download → quantify → delete per sample.

This script runs the unified download-quant workflow in sequential mode,
processing one sample at a time for maximum disk efficiency.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

# Import setup utilities (must be before other imports)
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Auto-setup and activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

from metainformant.rna.workflow import load_workflow_config
from metainformant.rna import steps as _steps_mod
from metainformant.core.io import read_delimited
from metainformant.core.utils.logging import get_logger

logger = get_logger("process_samples_sequential")


def main() -> int:
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Process samples sequentially: download → quantify → delete",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to species workflow config file",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Maximum number of samples to process (default: all)",
    )
    
    args = parser.parse_args()
    
    config_path = args.config.resolve()
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1
    
    # Load config
    cfg = load_workflow_config(config_path)
    
    # Get parameters
    getfastq_params = cfg.per_step.get("getfastq", {}).copy()
    quant_params = cfg.per_step.get("quant", {}).copy()
    
    # Force sequential mode (num_workers=1)
    getfastq_params["num_download_workers"] = 1
    
    # Find metadata file
    metadata_paths = [
        cfg.work_dir / "metadata" / "metadata.filtered.clean.tsv",
        cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv",
        cfg.work_dir / "metadata" / "metadata.tsv",
    ]
    metadata_file = None
    for mp in metadata_paths:
        if mp.exists():
            try:
                rows = list(read_delimited(mp, delimiter="\t"))
                if rows and "run" in rows[0]:
                    metadata_file = mp
                    logger.info(f"Using metadata file: {metadata_file}")
                    break
            except Exception:
                continue
    
    if not metadata_file:
        logger.error("No metadata file with 'run' column found")
        return 1
    
    # Run unified workflow in sequential mode
    logger.info("🚀 Starting sequential processing (one sample at a time)")
    logger.info("   Each sample: download → quantify → delete FASTQs → next sample")
    logger.info("   Maximum disk efficiency: only one sample's FASTQs exist at a time")
    
    try:
        stats = _steps_mod.run_download_quant_workflow(
            metadata_path=metadata_file,
            getfastq_params=getfastq_params,
            quant_params=quant_params,
            work_dir=cfg.work_dir,
            log_dir=(cfg.log_dir or (cfg.work_dir / "logs")),
            num_workers=1,  # Sequential mode: one sample at a time
            max_samples=args.max_samples,
            skip_completed=True,  # Skip already quantified samples
        )
        
        logger.info("=" * 80)
        logger.info("Processing complete!")
        logger.info(f"  Total samples: {stats['total_samples']}")
        logger.info(f"  Processed: {stats['processed']}")
        logger.info(f"  Skipped (already done): {stats['skipped']}")
        logger.info(f"  Failed: {stats['failed']}")
        if stats.get("failed_runs"):
            logger.warning(f"  Failed runs: {stats['failed_runs'][:10]}")
        logger.info("=" * 80)
        
        return 0 if stats.get("failed", 0) == 0 else 1
        
    except Exception as e:
        logger.error(f"Sequential processing failed: {e}", exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())


