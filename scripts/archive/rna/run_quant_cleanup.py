#!/usr/bin/env python3
"""Quantify completed samples and delete their FASTQ files.

This script finds samples with completed FASTQ files, quantifies them,
and deletes the FASTQ files to free disk space.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Import setup utilities (must be before other imports) - same as run_workflow.py
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import ensure_venv_activated, check_environment_or_exit

# Ensure virtual environment is activated and check environment
ensure_venv_activated()
check_environment_or_exit()

from metainformant.rna.orchestration import cleanup_unquantified_samples
from metainformant.rna.workflow import load_workflow_config


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/rna/run_quant_cleanup.py <config_file>")
        print("Example: python3 scripts/rna/run_quant_cleanup.py config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
        sys.exit(1)
    
    config_path = Path(sys.argv[1])
    if not config_path.exists():
        print(f"❌ Config file not found: {config_path}")
        sys.exit(1)
    
    print("=" * 80)
    print("🔬 Quantify and Cleanup Completed Samples")
    print("=" * 80)
    print()
    
    # Load config to get paths
    config = load_workflow_config(config_path)
    fastq_dir = Path(config.per_step.get("getfastq", {}).get("out_dir", config.work_dir / "fastq"))
    quant_dir = Path(config.per_step.get("quant", {}).get("out_dir", config.work_dir / "quant"))
    
    print(f"📁 FASTQ directory: {fastq_dir}")
    print(f"📁 Quant directory: {quant_dir}")
    print()
    
    # Use cleanup function directly - it will find unquantified samples
    from metainformant.rna.monitoring import find_unquantified_samples
    
    unquantified = find_unquantified_samples(config_path)
    
    if not unquantified:
        print("✅ No unquantified samples found")
        print("   (All downloaded samples are already quantified)")
        return
    
    print(f"📋 Found {len(unquantified)} downloaded but unquantified samples")
    print()
    print("🚀 Starting quantification and cleanup...")
    print("   (This will quantify samples and delete FASTQ files after quantification)")
    print()
    
    quantified, failed = cleanup_unquantified_samples(
        config_path,
        log_dir=config.log_dir or (config.work_dir.parent / "logs"),
    )
    
    print()
    print("=" * 80)
    print("✅ Summary")
    print("=" * 80)
    print(f"   Quantified: {quantified}")
    print(f"   Failed: {failed}")
    print("=" * 80)


if __name__ == "__main__":
    main()

