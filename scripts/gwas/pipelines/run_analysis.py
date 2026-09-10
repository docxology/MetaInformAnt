#!/usr/bin/env python3
"""Thin CLI wrapper for the full-scale P. barbatus GWAS pipeline (all variants, all traits).

Business logic lives in ``metainformant.gwas.workflow.pbarbatus_full_scale``;
this script only bootstraps the src tree and delegates.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from metainformant.gwas.workflow.pbarbatus_full_scale import run


def main() -> int:
    parser = argparse.ArgumentParser(description="Full-scale P. barbatus GWAS analysis (all variants, all traits)")
    parser.add_argument(
        "--config",
        type=str,
        default="config/gwas/gwas_pbarbatus_synthetic.yaml",
        help="GWAS YAML config path (default: config/gwas/gwas_pbarbatus_synthetic.yaml)",
    )
    args = parser.parse_args()
    run(args.config)
    return 0


if __name__ == "__main__":
    sys.exit(main())
