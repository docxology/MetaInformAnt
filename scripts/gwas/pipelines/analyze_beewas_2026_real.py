"""Compatibility wrapper for the BeeWAS real-cohort reporter.

The BeeWAS project implementation now lives in
``projects/apis_gwas/beewas/real/analysis.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module, run_project_main

reexport_project_module("analyze_beewas_2026_real", globals())


def main() -> int:
    """Run the canonical BeeWAS real-cohort analysis entry point."""
    return run_project_main("analyze_beewas_2026_real")


if __name__ == "__main__":
    raise SystemExit(main())
