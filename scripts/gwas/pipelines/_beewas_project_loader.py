"""Compatibility loader for BeeWAS modules moved into ``projects/apis_gwas``.

The ``apis_gwas`` submodule keeps its canonical implementation as a Python
package under ``beewas/``.  The parent repository exposes historical
top-level script names, so this loader bridges those names without copying or
duplicating the project implementation.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType

PROJECT_ROOT = Path(__file__).resolve().parents[3] / "projects" / "apis_gwas"
SOURCE_ROOT = PROJECT_ROOT.parent.parent / "src"
MODULE_ALIASES = {
    "beewas_phenotypes": "beewas.phenotypes",
    "beewas_reporting": "beewas.reporting",
    "analyze_beewas_2026_real": "beewas.real.analysis",
    "run_beewas_2026_real": "beewas.real.processing",
    "validate_beewas_synthetic_phenotypes": "beewas.validation.synthetic_phenotypes",
    "validate_beewas_genomic_estimators": "beewas.validation.genomic_estimators",
}


def load_project_module(module_name: str) -> ModuleType:
    """Import a BeeWAS project module without requiring a separate install."""
    import_name = MODULE_ALIASES.get(module_name, module_name)
    if not PROJECT_ROOT.is_dir():
        raise ImportError(f"BeeWAS project checkout not found: {PROJECT_ROOT}")
    if SOURCE_ROOT.is_dir() and str(SOURCE_ROOT) not in sys.path:
        sys.path.insert(0, str(SOURCE_ROOT))
    if str(PROJECT_ROOT) not in sys.path:
        sys.path.insert(0, str(PROJECT_ROOT))
    try:
        return importlib.import_module(import_name)
    except ModuleNotFoundError as exc:
        missing_name = exc.name or ""
        if missing_name == import_name or missing_name.startswith(f"{import_name}."):
            raise ImportError(
                f"BeeWAS project module not found: " f"{PROJECT_ROOT / import_name.replace('.', '/')}.py"
            ) from exc
        raise


def reexport_project_module(module_name: str, namespace: dict[str, object]) -> None:
    """Re-export public names from the migrated BeeWAS project module."""
    module = load_project_module(module_name)
    namespace["__doc__"] = module.__doc__
    for name in dir(module):
        if name.startswith("__"):
            continue
        namespace[name] = getattr(module, name)


def run_project_main(module_name: str) -> int:
    """Run a migrated project's CLI entry point and normalize its exit code."""
    module = load_project_module(module_name)
    entrypoint = getattr(module, "main", None)
    if entrypoint is None:
        raise RuntimeError(f"BeeWAS project module {module_name!r} has no CLI entrypoint")
    result = entrypoint()
    return 0 if result is None else int(result)
