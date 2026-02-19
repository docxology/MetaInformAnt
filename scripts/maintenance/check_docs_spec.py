#!/usr/bin/env python3
"""
Check for missing SPEC.md files in docs/ subdirectories.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

# List of expected modules that should have documentation
MODULES = [
    "core", "dna", "rna", "protein", "epigenome", "ontology", "phenotype",
    "ecology", "math", "gwas", "information", "life_events", "visualization",
    "simulation", "singlecell", "quality", "networks", "ml", "multiomics",
    "longread", "metagenomics", "structural_variants", "spatial",
    "pharmacogenomics", "menu"
]

def main():
    missing = []
    for module in MODULES:
        spec_file = DOCS_DIR / module / "SPEC.md"
        if not spec_file.exists():
            # Check if directory exists first
            if (DOCS_DIR / module).exists():
                missing.append(str(spec_file))
            else:
                print(f"Warning: docs/{module} directory does not exist")
    
    if missing:
        print("Missing SPEC.md files:")
        for m in missing:
            print(m)
    else:
        print("All docs subdirectories have SPEC.md")

if __name__ == "__main__":
    main()
