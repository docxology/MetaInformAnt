import os
from pathlib import Path

MODULES = {
    "core": "Shared utilities and infrastructure used across all domains",
    "dna": "Comprehensive DNA sequence analysis and manipulation",
    "ecology": "Ecological metadata and community analysis",
    "epigenome": "Epigenetic modification analysis",
    "gwas": "Genome-Wide Association Studies workflow",
    "information": "Information-theoretic analysis for biological data",
    "life_events": "Life course and event sequence analysis",
    "menu": "Interactive CLI menu system",
    "ml": "Statistical and machine learning methods",
    "multiomics": "Integrated multi-omic data analysis",
    "networks": "Biological network analysis",
    "phenotype": "Phenotypic trait analysis and curation",
    "protein": "Proteomic analysis and sequence manipulation",
    "quality": "Data quality assessment and control",
    "simulation": "Synthetic data generation and agent-based modeling",
    "singlecell": "Single-cell transcriptomic analysis",
    "visualization": "Unified plotting and animation framework",
}

TEMPLATE = """# {title} Technical Specification

## Architecture
The `{module}` module is a component of the `metainformant` package, designed to encapsulate {description}.

## Components
*   **`__init__.py`**: Exposes public API.
*   **Submodules**: Specialized components for domain logic.

## Dependencies
*   **Internal**: `metainformant.core` and related domain modules.
*   **External**: Standard scientific stack (numpy, pandas, scipy) and domain-specific tools.

## Standards
*   **Code Style**: PEP 8
*   **Docstrings**: Google Style
*   **Type Hints**: Full coverage required
*   **Testing**: Pytest with 100% coverage target
"""


def generate_specs():
    base_dir = Path("src/metainformant")
    created = []

    for module, desc in MODULES.items():
        spec_path = base_dir / module / "SPEC.md"
        if not spec_path.exists():
            print(f"Creating SPEC.md for {module}...")
            title = f"{module.replace('_', ' ').title()} Module"
            if module == "gwas":
                title = "GWAS Module"
            if module == "dna":
                title = "DNA Module"
            if module == "rna":
                title = "RNA Module"
            if module == "ml":
                title = "Machine Learning Module"

            content = TEMPLATE.format(title=title, module=module, description=desc.lower())
            spec_path.write_text(content)
            created.append(module)
        else:
            print(f"SPEC.md already exists for {module}")

    if created:
        print(f"\nCreated SPEC.md for: {', '.join(created)}")
    else:
        print("\nNo new SPEC.md files created.")


if __name__ == "__main__":
    generate_specs()
