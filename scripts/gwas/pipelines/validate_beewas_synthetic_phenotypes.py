"""Compatibility wrapper for BeeWAS synthetic-phenotype validation.

The BeeWAS project implementation now lives in
``projects/apis_gwas/scripts/beewas/validate_beewas_synthetic_phenotypes.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module


reexport_project_module("validate_beewas_synthetic_phenotypes", globals())


if __name__ == "__main__":
    raise SystemExit(main())
