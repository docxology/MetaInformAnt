"""Compatibility wrapper for the BeeWAS real-cohort reporter.

The BeeWAS project implementation now lives in
``projects/apis_gwas/scripts/beewas/analyze_beewas_2026_real.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module


reexport_project_module("analyze_beewas_2026_real", globals())


if __name__ == "__main__":
    raise SystemExit(main())
