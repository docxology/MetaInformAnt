"""Compatibility wrapper for BeeWAS phenotype curation.

The BeeWAS project implementation now lives in
``projects/apis_gwas/beewas/phenotypes.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module, run_project_main

reexport_project_module("beewas_phenotypes", globals())


def main() -> int:
    """Run the canonical BeeWAS phenotype-curation entry point."""
    return run_project_main("beewas_phenotypes")


if __name__ == "__main__":
    raise SystemExit(main())
