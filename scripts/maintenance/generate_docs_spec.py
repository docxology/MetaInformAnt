#!/usr/bin/env python3
"""
Generate missing SPEC.md files in docs/ subdirectories.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

# List of modules to check and generate for
MODULES = [
    "core", "dna", "rna", "protein", "epigenome", "ontology", "phenotype",
    "ecology", "math", "gwas", "information", "life_events", "visualization",
    "simulation", "singlecell", "quality", "networks", "ml", "multiomics",
    "longread", "metagenomics", "structural_variants", "spatial",
    "pharmacogenomics", "menu"
]

TEMPLATE = """# Specification: {module_name}

## 🎯 Scope
Documentation for the {module_name} domain in MetaInformAnt.

## 🧱 Architecture
- **Dependency Level**: Documentation
- **Component Type**: Guide
- **Location**: `docs/{module_name}/`

## 💾 Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## 🔌 Integration
- **Source**: `src/metainformant/{module_name}/`
- **Tests**: `tests/test_{module_name}_*.py`
"""

def main():
    print(f"Checking docs directories in {DOCS_DIR}...")
    
    for module in MODULES:
        module_dir = DOCS_DIR / module
        if not module_dir.exists():
            print(f"Skipping {module} (directory not found)")
            continue
            
        spec_file = module_dir / "SPEC.md"
        if not spec_file.exists():
            print(f"Generating SPEC.md for {module}...")
            content = TEMPLATE.format(module_name=module)
            with open(spec_file, "w") as f:
                f.write(content)
        else:
            print(f"SPEC.md exists for {module}")

    print("Done.")

if __name__ == "__main__":
    main()
