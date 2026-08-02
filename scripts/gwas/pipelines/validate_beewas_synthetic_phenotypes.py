"""Compatibility wrapper for BeeWAS synthetic-phenotype validation.

The BeeWAS project implementation now lives in
``projects/apis_gwas/beewas/validation/synthetic_phenotypes.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module, run_project_main

reexport_project_module("validate_beewas_synthetic_phenotypes", globals())


def main() -> int:
    """Run the canonical BeeWAS synthetic-phenotype validation entry point."""
    return run_project_main("validate_beewas_synthetic_phenotypes")


if __name__ == "__main__":
    raise SystemExit(main())
