
import os
import ast
from pathlib import Path
import textwrap

# Exclusion list
EXCLUDE_DIRS = {
    "__pycache__", "_build", ".git", ".github", ".vscode", ".idea", 
    ".pytest_cache", ".uv-cache", ".venv", ".tmp", "metainformant.egg-info",
    "html", "doctrees" # Sphinx build artifacts
}

REQUIRED_FILES = ["README.md", "AGENTS.md", "SPEC.md", "PAI.md"]

def get_dir_context(path: Path) -> dict:
    """
    Analyze the directory to extract context for documentation.
    Scans __init__.py or the first available .py file for docstrings.
    """
    context = {
        "description": f"Functionality for {path.name}.",
        "modules": [],
        "classes": [],
        "functions": []
    }
    
    # Try to read __init__.py or other py files
    py_files = sorted(list(path.glob("*.py")))
    
    for py_file in py_files:
        if py_file.name == "__init__.py":
            try:
                with open(py_file, "r") as f:
                    tree = ast.parse(f.read())
                    docstring = ast.get_docstring(tree)
                    if docstring:
                        context["description"] = docstring.strip().split('\n')[0]
            except Exception:
                pass
        
        # Collect basic stats
        context["modules"].append(py_file.name)
        
    return context

TEMPLATES = {
    "PAI.md": """# Personal AI Infrastructure (PAI) - {name}

## 🧠 Context & Intent
- **Path**: `{path}`
- **Purpose**: {description}
- **Domain**: {parent}

## 🏗️ Virtual Hierarchy
- **Type**: {type}
- **Parent**: `{parent}`

## 📝 Maintenance Notes
- **System**: Part of the METAINFORMANT {dep_level} layer.
- **Style**: Strict type hinting, no mocks in tests.
- **Stability**: API boundaries should be respected.

## 🔄 AI Workflows
- **Modification**: Run functional tests in `{test_dir}` before committing.
- **Documentation**: Update `SPEC.md` if architectural patterns change.
""",

    "SPEC.md": """# Specification: {name}

## 🎯 Scope
{description}

## 🧱 Architecture
- **Dependency Level**: {dep_level}
- **Component Type**: {type}

## 💾 Data Structures
- **Modules**: {module_count} Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
{module_list}
""",

    "AGENTS.md": """# Agent Directives: {name}

## 🤖 Role
Specialized agent context for the `{name}` component.

## 🛠️ Tools & Capabilities
- **Context**: {description}
- **Pattern**: {type} Pattern

## ⚠️ Rules & Constraints
- **Imports**: Prefer absolute imports from `metainformant`.
- **I/O**: Use `metainformant.core.io` for all file operations.
- **Logging**: Use `metainformant.core.logging`.
""",

    "README.md": """# {name_capitalized}

## Overview
{description}

## 📦 Contents
{contents_list}

## 📊 Structure

```mermaid
graph TD
    {name}[{name}]
    style {name} fill:#f9f,stroke:#333,stroke-width:2px
```

## Usage
Import module:
```python
from metainformant.{import_path} import ...
```
"""
}

def get_dir_type(path: Path) -> str:
    parts = path.parts
    if "src" in parts: return "Source Code"
    if "scripts" in parts: return "Orchestration Script"
    if "tests" in parts: return "Test Suite"
    if "docs" in parts: return "Documentation"
    if "config" in parts: return "Configuration"
    if "examples" in parts: return "Example/Demo"
    return "Directory"

def get_import_path(path: Path) -> str:
    # Try to construct python import path from src/metainformant/...
    try:
        parts = path.parts
        if "src" in parts:
            idx = parts.index("src")
            return ".".join(parts[idx+1:])
        return path.name
    except:
        return path.name

def get_contents_list(path: Path) -> str:
    items = []
    try:
        for item in path.iterdir():
            if item.name.startswith(".") or item.name == "__pycache__":
                continue
            if item.is_dir():
                items.append(f"- **[{item.name}/]({item.name}/)**")
            elif item.name.endswith(".py"):
                items.append(f"- `[{item.name}]({item.name})`")
    except Exception:
        pass
    return "\n".join(sorted(items))

def generate_for_directory(path: Path):
    if not path.is_dir():
        return

    # Check exclusions
    if any(ex in path.parts for ex in EXCLUDE_DIRS):
        return
    
    # Context
    name = path.name
    parent = path.parent.name
    rel_path = path.resolve()
    
    # Heuristics
    dep_level = "Core" if "core" in path.parts else "Domain"
    dir_type = get_dir_type(path)
    contents = get_contents_list(path)
    ctx = get_dir_context(path)
    import_path = get_import_path(path)
    
    # Locate likely test dir
    test_dir = "tests/"
    
    for filename, template in TEMPLATES.items():
        target_file = path / filename
        
        # For readability, we replace newlines in list
        module_list_str = "\n".join([f"- `{m}`" for m in ctx["modules"][:10]])
        if len(ctx["modules"]) > 10:
            module_list_str += "\n- ..."

        content = template.format(
            name=name,
            name_capitalized=name.upper(),
            parent=parent,
            path=str(rel_path),
            type=dir_type,
            dep_level=dep_level,
            contents_list=contents,
            description=ctx["description"],
            module_count=len(ctx["modules"]),
            module_list=module_list_str,
            import_path=import_path,
            test_dir=test_dir
        )
        
        # Always overwrite to ensure accuracy as per user request
        with open(target_file, "w") as f:
            f.write(content)
        # print(f"Updated {target_file}")

def walk_and_generate(root_path: Path):
    # Only walk relevant dirs
    relevant_roots = ["src", "scripts", "docs", "config", "tests", "examples"]
    
    for r in relevant_roots:
        p = root_path / r
        if not p.exists():
            continue
            
        for root, dirs, files in os.walk(p):
            # Modify dirs in-place to skip excluded
            dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
            
            current_path = Path(root)
            
            # Additional filter: if it's a python package or specific folders
            has_init = (current_path / "__init__.py").exists()
            is_valid_dir = any(x in current_path.parts for x in relevant_roots)
            
            if is_valid_dir:
                generate_for_directory(current_path)

if __name__ == "__main__":
    walk_and_generate(Path("."))
