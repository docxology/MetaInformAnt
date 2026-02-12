#!/usr/bin/env python3
"""
RNA-seq Multi-Species Orchestrator

This script provides robust orchestration for running the amalgkit pipeline across multiple species.
It handles:
- Sequential execution of pipeline steps (metadata -> getfastq -> ... -> curate)
- Error handling and retries
- Skipping species that fail specific steps
- Centralized logging
- Resuming from last successful step

Usage:
    uv run python scripts/rna/orchestrate_species.py --config config/amalgkit/cross_species.yaml
"""

import argparse
import sys
from pathlib import Path

# Add src to path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from metainformant.rna.engine.orchestration_multi_species import PipelineOrchestrator


def main():
    parser = argparse.ArgumentParser(description="RNA-seq Orchestrator")
    parser.add_argument("--config", type=Path, required=True, help="Path to YAML config file")
    args = parser.parse_args()
    
    orchestrator = PipelineOrchestrator(args.config)
    orchestrator.run()

if __name__ == "__main__":
    main()
