"""Compatibility wrapper for the BeeWAS genomic-estimator validator.

The BeeWAS project implementation now lives in
``projects/apis_gwas/beewas/validation/genomic_estimators.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module, run_project_main

reexport_project_module("validate_beewas_genomic_estimators", globals())


def main() -> int:
    """Run the canonical BeeWAS genomic-estimator validation entry point."""
    return run_project_main("validate_beewas_genomic_estimators")


if __name__ == "__main__":
    raise SystemExit(main())
