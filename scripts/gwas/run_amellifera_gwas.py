#!/usr/bin/env python3
"""Thin CLI wrapper for the Apis mellifera end-to-end GWAS pipeline.

Business logic lives in ``metainformant.gwas.workflow.amellifera_pipeline``;
this script only bootstraps the src tree and delegates.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.gwas.workflow.amellifera_pipeline import main

if __name__ == "__main__":
    main()
