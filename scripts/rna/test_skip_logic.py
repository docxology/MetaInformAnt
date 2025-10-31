#!/usr/bin/env python3
"""Test script to verify skip logic works correctly for already-quantified samples."""

import sys
from pathlib import Path

# Add src to path
repo_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_root / "src"))

from metainformant.rna.steps.batched_process import run_batched_download_quant, _sample_already_quantified
from metainformant.rna.workflow import load_workflow_config

def test_pbarbatus_skip_logic():
    """Test that Pbarbatus skip logic works."""
    print("=" * 80)
    print("TESTING SKIP LOGIC FOR PBARBATUS")
    print("=" * 80)
    print()
    
    # Load Pbarbatus config
    config_path = repo_root / "config" / "amalgkit_pbarbatus.yaml"
    print(f"Loading config: {config_path}")
    config = load_workflow_config(config_path)
    print(f"  Work dir: {config.work_dir}")
    
    # Check quant directory
    quant_params = config.per_step.get("quant", {})
    quant_dir = Path(quant_params.get("out_dir", config.work_dir.parent / "quant")).absolute()
    print(f"  Quant dir: {quant_dir}")
    
    # Count quantified samples
    quantified = sum(1 for _ in quant_dir.glob("SRR*/abundance.tsv") if _sample_already_quantified(_.parent.name, quant_dir))
    print(f"  Quantified samples: {quantified}")
    print()
    
    # Find metadata file
    metadata_paths = [
        config.work_dir / "metadata" / "metadata.filtered.tissue.tsv",
        config.work_dir / "metadata" / "metadata.filtered.clean.tsv",
        config.work_dir / "metadata" / "metadata.tsv",
    ]
    
    metadata_file = None
    for mp in metadata_paths:
        if mp.exists():
            metadata_file = mp
            break
    
    if not metadata_file:
        print("❌ No metadata file found. Running workflow first to generate metadata...")
        print("   Expected paths:")
        for mp in metadata_paths:
            print(f"     - {mp}")
        return False
    
    print(f"Using metadata: {metadata_file}")
    print()
    
    # Get params
    getfastq_params = config.per_step.get("getfastq", {})
    quant_params = config.per_step.get("quant", {})
    
    print("Running batched processing (should skip all quantified samples)...")
    print()
    
    try:
        stats = run_batched_download_quant(
            metadata_path=metadata_file,
            getfastq_params=getfastq_params,
            quant_params=quant_params,
            work_dir=config.work_dir,
            log_dir=config.work_dir / "logs",
            batch_size=8,
        )
        
        print()
        print("=" * 80)
        print("RESULTS")
        print("=" * 80)
        print(f"Total samples: {stats['total_samples']}")
        print(f"Skipped: {stats['skipped']}")
        print(f"Processed: {stats['processed']}")
        print(f"Failed: {stats['failed']}")
        print(f"Batches: {stats['batches']}")
        print()
        
        # Verify results
        if stats['skipped'] == quantified and stats['processed'] == 0 and stats['batches'] == 0:
            print("✅ SUCCESS: All samples properly skipped!")
            print("   - No downloads initiated")
            print("   - No batches created")
            print("   - Skip logic working correctly")
            return True
        else:
            print("⚠️  UNEXPECTED RESULTS:")
            if stats['skipped'] != quantified:
                print(f"   - Expected {quantified} skipped, got {stats['skipped']}")
            if stats['processed'] > 0:
                print(f"   - Expected 0 processed, got {stats['processed']}")
            if stats['batches'] > 0:
                print(f"   - Expected 0 batches, got {stats['batches']}")
            return False
            
    except Exception as e:
        print(f"❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_pbarbatus_skip_logic()
    sys.exit(0 if success else 1)

