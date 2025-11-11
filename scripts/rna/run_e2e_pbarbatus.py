#!/usr/bin/env python3
"""Run end-to-end workflow for Pogonomyrmex barbatus with single sample (TESTING ONLY).

This script runs the complete amalgkit pipeline from metadata to sanity check
for a single sample to verify the entire system works correctly.

**NOTE**: This script is for testing purposes only (limits to 1 sample).
For production end-to-end workflows with all samples, use:
    python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml
"""

from __future__ import annotations

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, plan_workflow, execute_workflow
from metainformant.core.logging import get_logger

logger = get_logger("e2e_pbarbatus")


def main() -> int:
    """Run end-to-end workflow for P. barbatus."""
    # Load config
    config_path = Path("config/amalgkit/amalgkit_pogonomyrmex_barbatus.yaml")
    if not config_path.exists():
        logger.error(f"Config file not found: {config_path}")
        return 1
    
    logger.info("=" * 80)
    logger.info("END-TO-END WORKFLOW: Pogonomyrmex barbatus")
    logger.info("=" * 80)
    logger.info(f"Loading config: {config_path}")
    
    cfg = load_workflow_config(config_path)
    
    # Limit to 1 sample for end-to-end testing
    logger.info("Limiting to 1 sample for comprehensive end-to-end test")
    if "metadata" not in cfg.per_step:
        cfg.per_step["metadata"] = {}
    cfg.per_step["metadata"]["max_samples"] = 1
    logger.info("Set max_samples=1 in metadata step")
    
    # Create filtered single-sample metadata file BEFORE workflow execution
    # This ensures only 1 sample is processed through the entire pipeline
    from metainformant.core.io import read_delimited, write_delimited
    
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if metadata_file.exists():
        rows = list(read_delimited(metadata_file, delimiter="\t"))
        if len(rows) > 1:  # Has header + at least 1 sample
            # Create single-sample metadata
            filtered_rows = [rows[0]] + [rows[1]]  # Header + first sample
            single_metadata = metadata_file.parent / "metadata.single_sample.tsv"
            write_delimited(filtered_rows, single_metadata, delimiter="\t")
            logger.info(f"✅ Created single-sample metadata: {single_metadata}")
            logger.info(f"   Sample to process: {rows[1].get('run', 'unknown')}")
            logger.info(f"   Total samples in full metadata: {len(rows)-1}")
            
            # Update config to use filtered metadata for ALL downstream steps
            # This ensures the workflow processes only 1 sample
            if "select" not in cfg.per_step:
                cfg.per_step["select"] = {}
            cfg.per_step["select"]["metadata"] = str(single_metadata)
            
            if "getfastq" not in cfg.per_step:
                cfg.per_step["getfastq"] = {}
            cfg.per_step["getfastq"]["metadata"] = str(single_metadata)
            
            logger.info("✅ Updated config to use single-sample metadata for select and getfastq steps")
        else:
            logger.warning("Not enough rows in metadata file to create single-sample version")
    else:
        logger.warning("Metadata file not found - workflow will create it during execution")
    
    # Plan workflow
    logger.info("\nPlanning workflow...")
    steps = plan_workflow(cfg)
    step_names = [name for name, _ in steps]
    logger.info(f"Planned {len(steps)} steps:")
    for i, (name, _) in enumerate(steps, 1):
        logger.info(f"  {i}. {name}")
    
    # Execute workflow
    logger.info("\n" + "=" * 80)
    logger.info("EXECUTING WORKFLOW")
    logger.info("=" * 80)
    
    return_codes = execute_workflow(cfg, check=False)
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("WORKFLOW EXECUTION SUMMARY")
    logger.info("=" * 80)
    
    for i, (step_name, return_code) in enumerate(zip(step_names, return_codes), 1):
        status = "✅ PASSED" if return_code == 0 else f"❌ FAILED (code {return_code})"
        logger.info(f"{i:2d}. {step_name:15s} {status}")
    
    # Overall status
    all_passed = all(rc == 0 for rc in return_codes)
    if all_passed:
        logger.info("\n✅ All steps completed successfully!")
        return 0
    else:
        failed_count = sum(1 for rc in return_codes if rc != 0)
        logger.warning(f"\n⚠️  {failed_count} step(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

