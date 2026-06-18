"""Compatibility wrapper for BeeWAS reporting helpers.

The BeeWAS project implementation now lives in
``projects/apis_gwas/scripts/beewas/beewas_reporting.py``.
"""

from __future__ import annotations

from _beewas_project_loader import reexport_project_module


reexport_project_module("beewas_reporting", globals())
