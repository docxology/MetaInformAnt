"""Compatibility wrapper for BeeWAS phenotype curation.

The BeeWAS project implementation now lives in
``projects/apis_gwas/scripts/beewas/beewas_phenotypes.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module


reexport_project_module("beewas_phenotypes", globals())


if __name__ == "__main__":
    raise SystemExit(main())
