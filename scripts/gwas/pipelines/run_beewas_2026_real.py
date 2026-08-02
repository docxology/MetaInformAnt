"""Compatibility wrapper for the BeeWAS real-read runner.

The BeeWAS project implementation now lives in
``projects/apis_gwas/beewas/real/processing.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module, run_project_main

reexport_project_module("run_beewas_2026_real", globals())


def main() -> int:
    """Run the canonical BeeWAS real-read processing entry point."""
    return run_project_main("run_beewas_2026_real")


if __name__ == "__main__":
    raise SystemExit(main())
