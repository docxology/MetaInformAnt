#!/usr/bin/env python3
"""
Test workflow for a single species to verify fixes.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.rna.workflow import load_workflow_config, execute_workflow
from metainformant.core.logging import get_logger

def main():
    logger = get_logger("test_single_species")
    
    # Test with P. barbatus (main test species)
    config_path = Path("config/amalgkit/amalgkit_pbarbatus.yaml")
    
    # Fall back to old location if new structure not yet adopted
    if not config_path.exists():
        config_path = Path("config/amalgkit_pbarbatus.yaml")
    
    logger.info("="*80)
    logger.info(f"Testing workflow with: {config_path}")
    logger.info("="*80)
    
    config = load_workflow_config(config_path)
    logger.info(f"Species: {config.species_list}")
    logger.info(f"Work dir: {config.work_dir}")
    logger.info(f"Filters: {config.filters}")
    logger.info("")
    
    # Run workflow
    return_codes = execute_workflow(config, check=False)
    
    logger.info("")
    logger.info("="*80)
    logger.info("RESULTS")
    logger.info("="*80)
    
    steps = ["genome", "metadata", "integrate", "config", "select", 
             "getfastq", "quant", "merge", "cstmm", "curate", "csca", "sanity"]
    
    for i, (step, code) in enumerate(zip(steps, return_codes)):
        status = "✅ PASS" if code == 0 else f"❌ FAIL (code {code})"
        logger.info(f"{step:15s}: {status}")
    
    logger.info("="*80)

if __name__ == "__main__":
    main()

