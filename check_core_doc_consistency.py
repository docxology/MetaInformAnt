#!/usr/bin/env python3
"""Backward-compatible entry point for the core documentation checker.

The implementation lives under ``scripts/`` with the other repository
validation tools. Keeping this small wrapper prevents the historical root
command from drifting away from the source-derived API contract.
"""

from scripts.core_docs_cross_check import main


if __name__ == "__main__":
    raise SystemExit(main())
