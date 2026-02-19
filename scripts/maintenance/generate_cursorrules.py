#!/usr/bin/env python3
"""
Generate module-specific cursorrules files.

This script automates the creation of .cursorrules files for each domain module
in the MetaInformAnt repository, ensuring consistent application of:
1. Zero Mock Policy
2. UV Package Management
3. Output Isolation
4. Domain-specific configuration
"""

import os
from pathlib import Path
from typing import Dict

# Configuration
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CURSORRULES_DIR = REPO_ROOT / "cursorrules"
SRC_DIR = REPO_ROOT / "src" / "metainformant"

# map module name to env var prefix
MODULE_PREFIXES: Dict[str, str] = {
    "core": "CORE",
    "dna": "DNA",
    "rna": "AK",
    "protein": "PROT",
    "epigenome": "EPI",
    "ontology": "ONT",
    "phenotype": "PHEN",
    "ecology": "ECO",
    "math": "MATH",
    "gwas": "GWAS",
    "information": "INFO",
    "life_events": "LE",
    "visualization": "VIZ",
    "simulation": "SIM",
    "singlecell": "SC",
    "quality": "QC",
    "networks": "NET",
    "ml": "ML",
    "multiomics": "MULTI",
    "longread": "LR",
    "metagenomics": "META",
    "structural_variants": "SV",
    "spatial": "SPATIAL",
    "pharmacogenomics": "PHARMA",
    "menu": "MENU",
    "metabolomics": "METAB",
}

TEMPLATE = """# Cursor Rules for {module_name} Module ({prefix}_)

## 🧠 Context & Intent
- **Domain**: {module_name}
- **Path**: `src/metainformant/{module_name}/`
- **Config Prefix**: `{prefix}_` (e.g., `{prefix}_THREADS`)

## 🧱 Directory Structure
- **Source**: `src/metainformant/{module_name}/`
- **Tests**: `tests/test_{module_name}_*.py`
- **Docs**: `docs/{module_name}/`
- **Output**: `output/{module_name}/` (Default for all artifacts)

## ⚙️ Configuration
- Use `metainformant.core.config` to load config.
- Support env checks with `{prefix}_` prefix.
- Example:
  ```python
  config = load_domain_config("config/{module_name}/default.yaml", prefix="{prefix}")
  ```

## 🧪 Testing (ZERO MOCK POLICY)
- **STRICTLY NO MOCKS**: Use real implementations only.
- **External Tools**: Skip if not found (don't mock).
- **Network**: Use real API calls with `@pytest.mark.network`.
- **Output**: Write all test artifacts to `output/` via `tmp_path`.

## 📦 Dependencies
- Manage via `uv`: `uv add`, `uv remove`.
- Import pattern:
  ```python
  from metainformant.core import io, logging
  from metainformant.{module_name} import utils
  ```
"""

def main():
    """Generate cursorrules files."""
    CURSORRULES_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating cursorrules in {CURSORRULES_DIR}...")
    
    modules = [d.name for d in SRC_DIR.iterdir() if d.is_dir() and not d.name.startswith("__")]
    
    for module in modules:
        if module not in MODULE_PREFIXES:
            print(f"Warning: No prefix definition for {module}, skipping.")
            continue
            
        prefix = MODULE_PREFIXES[module]
        content = TEMPLATE.format(module_name=module, prefix=prefix)
        
        output_file = CURSORRULES_DIR / f"{module}.cursorrules"
        with open(output_file, "w") as f:
            f.write(content)
            
        print(f"Created {output_file.name}")

    print("Done.")

if __name__ == "__main__":
    main()
