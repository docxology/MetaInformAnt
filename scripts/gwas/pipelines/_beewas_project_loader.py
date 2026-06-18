"""Compatibility loader for BeeWAS scripts moved into projects/apis_gwas."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType


PROJECT_SCRIPT_DIR = Path(__file__).resolve().parents[3] / "projects" / "apis_gwas" / "scripts" / "beewas"


def load_project_module(module_name: str) -> ModuleType:
    """Load a BeeWAS project module by filename without requiring package install."""
    module_path = PROJECT_SCRIPT_DIR / f"{module_name}.py"
    if not module_path.exists():
        raise ImportError(f"BeeWAS project module not found: {module_path}")
    if str(PROJECT_SCRIPT_DIR) not in sys.path:
        sys.path.insert(0, str(PROJECT_SCRIPT_DIR))
    alias = f"_apis_gwas_project_{module_name}"
    cached = sys.modules.get(alias)
    if cached is not None:
        return cached
    spec = importlib.util.spec_from_file_location(alias, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load BeeWAS project module: {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


def reexport_project_module(module_name: str, namespace: dict[str, object]) -> None:
    """Re-export public names from the migrated BeeWAS project module."""
    module = load_project_module(module_name)
    namespace["__doc__"] = module.__doc__
    for name in dir(module):
        if name.startswith("__"):
            continue
        namespace[name] = getattr(module, name)
