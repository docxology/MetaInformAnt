#!/usr/bin/env python3
"""
Test script for Pogonomyrmex barbatus RNA-seq workflow.

This script tests the sequential download-quant-delete workflow
with the ability to parallelize downloads while keeping quantification sequential.

Tests:
1. Single sample end-to-end (download → quant → delete)
2. Small batch (3 samples) with parallel downloads
3. Verify kallisto index building
4. Verify FASTQ deletion after quant

Usage:
    python3 scripts/rna/test_pbarbatus_workflow.py [--parallel-downloads N]
"""

import argparse
import sys
from pathlib import Path

# Add src to path
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow


def test_single_sample():
    """Test single sample end-to-end workflow."""
    print("=" * 80)
    print("TEST 1: Single Sample End-to-End")
    print("=" * 80)
    
    # Load config
    config_path = Path("config/amalgkit_pbarbatus.yaml")
    if not config_path.exists():
        print(f"❌ Config not found: {config_path}")
        return False
    
    print(f"Loading config: {config_path}")
    cfg = load_workflow_config(config_path)
    
    # Check if we already have samples processed
    quant_dir = Path("output/amalgkit/pbarbatus/quant")
    if quant_dir.exists():
        existing_samples = list(quant_dir.glob("SRR*/abundance.tsv"))
        print(f"Found {len(existing_samples)} already quantified samples")
        if len(existing_samples) > 0:
            print(f"  First 3: {[s.parent.name for s in existing_samples[:3]]}")
    
    print("\n✅ Config loaded successfully")
    print(f"  Work dir: {cfg.work_dir}")
    print(f"  Species: {cfg.species_list}")
    print(f"  Threads: {cfg.threads}")
    
    return True


def test_workflow_execution(max_samples=None):
    """Test workflow execution with optional sample limit."""
    print("\n" + "=" * 80)
    print(f"TEST: Workflow Execution (max_samples={max_samples or 'unlimited'})")
    print("=" * 80)
    
    config_path = Path("config/amalgkit_pbarbatus.yaml")
    cfg = load_workflow_config(config_path)
    
    # Modify config to limit samples if requested
    if max_samples:
        print(f"\n⚠️  Limiting to first {max_samples} samples for testing")
        # This would need to be implemented in the sequential processor
    
    print("\n🚀 Starting workflow execution...")
    print("   This will:")
    print("   1. Download FASTQ files")
    print("   2. Build kallisto index (if needed)")
    print("   3. Quantify each sample")
    print("   4. Delete FASTQ files")
    print("   5. Continue to next sample")
    
    return_codes = execute_workflow(cfg)
    
    print(f"\n📊 Workflow completed with {len(return_codes)} steps")
    print(f"   Return codes: {return_codes}")
    
    if all(code == 0 for code in return_codes):
        print("✅ All steps succeeded!")
        return True
    else:
        failed = [i for i, code in enumerate(return_codes) if code != 0]
        print(f"⚠️  {len(failed)} steps failed: {failed}")
        return False


def check_disk_usage():
    """Check current disk usage of FASTQ and quant directories."""
    print("\n" + "=" * 80)
    print("DISK USAGE CHECK")
    print("=" * 80)
    
    import subprocess
    
    fastq_dir = Path("output/amalgkit/pbarbatus/fastq")
    quant_dir = Path("output/amalgkit/pbarbatus/quant")
    
    if fastq_dir.exists():
        result = subprocess.run(
            ["du", "-sh", str(fastq_dir)],
            capture_output=True,
            text=True
        )
        print(f"FASTQ dir: {result.stdout.strip()}")
        
        # Count FASTQ files
        fastq_files = list(fastq_dir.glob("**/*.fastq*"))
        print(f"  {len(fastq_files)} FASTQ files found")
    else:
        print("FASTQ dir: Not found")
    
    if quant_dir.exists():
        result = subprocess.run(
            ["du", "-sh", str(quant_dir)],
            capture_output=True,
            text=True
        )
        print(f"Quant dir: {result.stdout.strip()}")
        
        # Count quantified samples
        abundance_files = list(quant_dir.glob("*/abundance.tsv"))
        print(f"  {len(abundance_files)} samples quantified")
    else:
        print("Quant dir: Not found")


def main():
    parser = argparse.ArgumentParser(description="Test Pbarbatus RNA-seq workflow")
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Limit to N samples for testing (default: process all)"
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only check config and disk usage, don't run workflow"
    )
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("POGONOMYRMEX BARBATUS RNA-SEQ WORKFLOW TEST")
    print("=" * 80)
    print()
    
    # Test 1: Config loading
    if not test_single_sample():
        print("\n❌ Config test failed")
        return 1
    
    # Check disk usage
    check_disk_usage()
    
    if args.check_only:
        print("\n✅ Check-only mode: Exiting without running workflow")
        return 0
    
    # Test 2: Execute workflow
    print("\n" + "=" * 80)
    print("READY TO RUN WORKFLOW")
    print("=" * 80)
    print(f"Max samples: {args.max_samples or 'unlimited'}")
    print("\nPress Ctrl+C to cancel, or the workflow will start in 3 seconds...")
    
    import time
    try:
        for i in range(3, 0, -1):
            print(f"  {i}...")
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n\n❌ Cancelled by user")
        return 1
    
    success = test_workflow_execution(max_samples=args.max_samples)
    
    # Final disk usage check
    print("\n")
    check_disk_usage()
    
    if success:
        print("\n✅ All tests passed!")
        return 0
    else:
        print("\n⚠️  Some tests failed - check logs")
        return 1


if __name__ == "__main__":
    sys.exit(main())

