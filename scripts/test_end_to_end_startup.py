#!/usr/bin/env python3
"""Test end-to-end workflow startup and continuation for all species.

This script performs a dry-run test of the workflow startup capability
without actually executing the full workflow.
"""

import os
import sys
import yaml
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

def test_all_species_startup():
    """Test that all species configs can be loaded and workflows planned."""
    from metainformant.rna.workflow import load_workflow_config, plan_workflow
    
    config_dir = REPO_ROOT / "config" / "amalgkit"
    configs = sorted(config_dir.glob("amalgkit_*.yaml"))
    configs = [c for c in configs if "template" not in c.stem.lower() and "test" not in c.stem.lower()]
    
    print("=" * 80)
    print("TESTING END-TO-END WORKFLOW STARTUP FOR ALL SPECIES")
    print("=" * 80)
    print()
    
    results = []
    
    for config_path in configs:
        species_name = config_path.stem.replace("amalgkit_", "")
        
        try:
            # Load config
            cfg = load_workflow_config(config_path)
            
            # Plan workflow
            steps = plan_workflow(cfg)
            
            # Test with environment variable override
            os.environ["AK_THREADS"] = "12"
            cfg_override = load_workflow_config(config_path)
            os.environ.pop("AK_THREADS", None)
            
            results.append({
                "species": species_name,
                "config_path": str(config_path),
                "threads": cfg.threads,
                "threads_override": cfg_override.threads,
                "steps": len(steps),
                "work_dir": str(cfg.work_dir),
                "status": "OK"
            })
            
            print(f"✅ {species_name:30s} - {len(steps)} steps, threads: {cfg.threads} (override: {cfg_override.threads})")
        
        except Exception as e:
            results.append({
                "species": species_name,
                "config_path": str(config_path),
                "status": f"ERROR: {e}"
            })
            print(f"❌ {species_name:30s} - ERROR: {e}")
    
    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total species: {len(configs)}")
    print(f"Successful: {sum(1 for r in results if r['status'] == 'OK')}")
    print(f"Failed: {sum(1 for r in results if r['status'] != 'OK')}")
    
    if all(r['status'] == 'OK' for r in results):
        print()
        print("✅ ALL SPECIES CAN START END-TO-END WORKFLOWS")
        print()
        print("To start all species with configurable threads:")
        print("  export AK_THREADS=12")
        print("  python3 scripts/rna/run_multi_species.py")
        return True
    else:
        print()
        print("❌ SOME SPECIES HAVE ISSUES - Review errors above")
        return False

if __name__ == "__main__":
    success = test_all_species_startup()
    sys.exit(0 if success else 1)

