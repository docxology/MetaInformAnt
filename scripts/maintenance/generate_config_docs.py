#!/usr/bin/env python3
"""
Generate documentation for config/ subdirectories.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
CONFIG_DIR = REPO_ROOT / "config"

README_TEMPLATE = """# Configuration: {module_name}

## 📝 Overview
Configuration files for the {module_name} module.

## 📂 Structure
- `*.yaml`: Configuration templates
- Environment overrides supported via `{prefix}_` prefix.

## 🔗 Links
- [Module Docs](../../docs/{module_name}/README.md)
"""

AGENTS_TEMPLATE = """# AI Agents - {module_name} Config

## 🤖 Agents
- **Config Validator**: Validates YAML structure against schema.
- **Migration Assistant**: Helps port legacy configs.

## 🛠️ Workflows
- **Validation**: Schema validation on load.
"""

SPEC_TEMPLATE = """# Specification: {module_name} Config

## 🎯 Scope
Configuration templates and schemas for {module_name}.

## 🧱 Architecture
- **Component**: Configuration
- **Location**: `config/{module_name}/`

## 💾 Data Structures
- **Format**: YAML
- **Schema**: `metainformant.{module_name}.config`
"""

# Map module to prefix
PREFIXES = {
    "amalgkit": "AK",
    "config_base": "CORE",
    "eqtl": "GWAS",
    "gwas": "GWAS",
    "life_events": "LE",
    "longread": "LR",
    "multiomics": "MULTI",
    "ncbi": "DNA",
    "networks": "NET",
    "phenotype": "PHEN",
    "singlecell": "SC"
}

def main():
    print(f"Checking config directories in {CONFIG_DIR}...")
    
    for module_dir in CONFIG_DIR.iterdir():
        if not module_dir.is_dir() or module_dir.name.startswith("__"):
            continue
            
        module = module_dir.name
        print(f"Processing {module}...")
        
        prefix = PREFIXES.get(module, module.upper())
        
        # README.md
        readme = module_dir / "README.md"
        if not readme.exists():
            with open(readme, "w") as f:
                f.write(README_TEMPLATE.format(module_name=module, prefix=prefix))
            print(f"Created {readme.name}")
            
        # AGENTS.md
        agents = module_dir / "AGENTS.md"
        if not agents.exists():
            with open(agents, "w") as f:
                f.write(AGENTS_TEMPLATE.format(module_name=module))
            print(f"Created {agents.name}")
            
        # SPEC.md
        spec = module_dir / "SPEC.md"
        if not spec.exists():
            with open(spec, "w") as f:
                f.write(SPEC_TEMPLATE.format(module_name=module))
            print(f"Created {spec.name}")

    print("Done.")

if __name__ == "__main__":
    main()
