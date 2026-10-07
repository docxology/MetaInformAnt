#!/usr/bin/env python3
"""Thin entry point for the shared local/AWS Amalgkit manifest worker."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"src"))
from metainformant.rna.engine.acquisition_worker_cli import main  # noqa: E402
if __name__ == "__main__":
    raise SystemExit(main())
