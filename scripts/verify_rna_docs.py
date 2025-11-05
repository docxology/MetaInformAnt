#!/usr/bin/env python3
"""Comprehensive verification script for RNA documentation and code.

Verifies that all RNA documentation, methods, and scripts are accurate,
complete, and functional.
"""

import ast
import importlib
import inspect
import subprocess
import sys
from pathlib import Path
from typing import Any

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

REPO_ROOT = Path(__file__).parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "rna"
SRC_DIR = REPO_ROOT / "src" / "metainformant" / "rna"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "rna"

issues = []
warnings = []


def check_module_docstrings(module_path: Path, module_name: str) -> list[str]:
    """Check that all public functions/classes have docstrings."""
    issues_found = []
    
    try:
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            return [f"Cannot load module {module_name}"]
        
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
        # Get all public members
        for name, obj in inspect.getmembers(module):
            if name.startswith("_"):
                continue
            
            # Check functions
            if inspect.isfunction(obj):
                if obj.__doc__ is None or not obj.__doc__.strip():
                    issues_found.append(f"{module_name}.{name}: Missing docstring")
                elif len(obj.__doc__.strip()) < 20:
                    issues_found.append(f"{module_name}.{name}: Docstring too short")
            
            # Check classes
            elif inspect.isclass(obj):
                if obj.__doc__ is None or not obj.__doc__.strip():
                    issues_found.append(f"{module_name}.{name}: Missing class docstring")
                
                # Check methods
                for method_name, method in inspect.getmembers(obj, inspect.isfunction):
                    if method_name.startswith("_"):
                        continue
                    if method.__doc__ is None or not method.__doc__.strip():
                        issues_found.append(f"{module_name}.{name}.{method_name}: Missing docstring")
    
    except Exception as e:
        issues_found.append(f"{module_name}: Error checking - {e}")
    
    return issues_found


def check_imports_exist(module_name: str) -> bool:
    """Check if a module can be imported."""
    try:
        importlib.import_module(module_name)
        return True
    except Exception:
        return False


def check_doc_examples(file_path: Path) -> list[str]:
    """Check code examples in documentation files."""
    issues_found = []
    
    if not file_path.exists():
        return [f"Documentation file missing: {file_path}"]
    
    content = file_path.read_text()
    
    # Find code blocks
    in_code_block = False
    code_lines = []
    code_start_line = 0
    
    for i, line in enumerate(content.split("\n"), 1):
        if line.strip().startswith("```python"):
            in_code_block = True
            code_lines = []
            code_start_line = i
            continue
        elif line.strip().startswith("```") and in_code_block:
            # Check syntax of code block
            try:
                ast.parse("\n".join(code_lines))
            except SyntaxError as e:
                issues_found.append(
                    f"{file_path}:{code_start_line}: Syntax error in code example: {e.msg}"
                )
            in_code_block = False
            code_lines = []
        elif in_code_block:
            code_lines.append(line)
    
    return issues_found


def check_script_syntax(script_path: Path) -> list[str]:
    """Check that a script has valid Python syntax."""
    issues_found = []
    
    if not script_path.exists():
        return [f"Script missing: {script_path}"]
    
    try:
        with open(script_path, "r") as f:
            ast.parse(f.read(), filename=str(script_path))
    except SyntaxError as e:
        issues_found.append(f"{script_path}: Syntax error: {e.msg} at line {e.lineno}")
    
    return issues_found


def check_file_exists(file_path: Path) -> bool:
    """Check if a file exists."""
    return file_path.exists()


def check_doc_links(doc_file: Path) -> list[str]:
    """Check internal markdown links in documentation."""
    issues_found = []
    
    if not doc_file.exists():
        return [f"Documentation file missing: {doc_file}"]
    
    content = doc_file.read_text()
    doc_dir = doc_file.parent
    
    # Find markdown links [text](path)
    import re
    link_pattern = re.compile(r'\[([^\]]+)\]\(([^)]+)\)')
    
    for match in link_pattern.finditer(content):
        link_text, link_path = match.groups()
        
        # Skip external links
        if link_path.startswith("http"):
            continue
        
        # Resolve relative paths
        if link_path.startswith("/"):
            target = REPO_ROOT / link_path.lstrip("/")
        else:
            target = (doc_dir / link_path).resolve()
        
        # Check if file exists
        if not target.exists():
            issues_found.append(f"{doc_file}: Broken link to {link_path} ({link_text})")
    
    return issues_found


def main():
    """Run comprehensive verification."""
    print("=" * 80)
    print("RNA DOCUMENTATION AND CODE VERIFICATION")
    print("=" * 80)
    print()
    
    # Phase 1: Method Documentation Verification
    print("Phase 1: Method Documentation Verification")
    print("-" * 80)
    
    rna_modules = [
        ("metainformant.rna.amalgkit", SRC_DIR / "amalgkit.py"),
        ("metainformant.rna.workflow", SRC_DIR / "workflow.py"),
        ("metainformant.rna.configs", SRC_DIR / "configs.py"),
        ("metainformant.rna.monitoring", SRC_DIR / "monitoring.py"),
        ("metainformant.rna.environment", SRC_DIR / "environment.py"),
    ]
    
    for module_name, module_path in rna_modules:
        if module_path.exists():
            module_issues = check_module_docstrings(module_path, module_name)
            issues.extend(module_issues)
            if module_issues:
                print(f"  ❌ {module_name}: {len(module_issues)} issues")
                for issue in module_issues[:3]:
                    print(f"     - {issue}")
            else:
                print(f"  ✅ {module_name}: All public methods documented")
        else:
            warnings.append(f"Module file not found: {module_path}")
    
    # Check step modules
    steps_dir = SRC_DIR / "steps"
    if steps_dir.exists():
        for step_file in steps_dir.glob("*.py"):
            if step_file.name.startswith("_") or step_file.name == "__init__.py":
                continue
            step_name = f"metainformant.rna.steps.{step_file.stem}"
            step_issues = check_module_docstrings(step_file, step_name)
            issues.extend(step_issues)
            if step_issues:
                print(f"  ❌ {step_name}: {len(step_issues)} issues")
            else:
                print(f"  ✅ {step_name}: Documented")
    
    print()
    
    # Phase 2: Documentation Accuracy Verification
    print("Phase 2: Documentation Accuracy Verification")
    print("-" * 80)
    
    key_docs = [
        DOCS_DIR / "README.md",
        DOCS_DIR / "WORKFLOW.md",
        DOCS_DIR / "STEPS.md",
        DOCS_DIR / "CONFIGURATION.md",
        DOCS_DIR / "ORCHESTRATION" / "ENA_WORKFLOW.md",
        DOCS_DIR / "ORCHESTRATION" / "MULTI_SPECIES.md",
        DOCS_DIR / "ORCHESTRATION" / "BATCH_DOWNLOAD.md",
    ]
    
    for doc_file in key_docs:
        if doc_file.exists():
            # Check syntax in examples
            example_issues = check_doc_examples(doc_file)
            issues.extend(example_issues)
            
            # Check links
            link_issues = check_doc_links(doc_file)
            issues.extend(link_issues)
            
            if example_issues or link_issues:
                print(f"  ❌ {doc_file.name}: {len(example_issues) + len(link_issues)} issues")
            else:
                print(f"  ✅ {doc_file.name}: Examples and links valid")
        else:
            warnings.append(f"Documentation file missing: {doc_file}")
    
    print()
    
    # Phase 3: Script Verification
    print("Phase 3: Script Verification")
    print("-" * 80)
    
    orchestrator_scripts = [
        SCRIPTS_DIR / "workflow_ena_integrated.py",
        SCRIPTS_DIR / "run_multi_species.py",
        SCRIPTS_DIR / "batch_download_species.py",
    ]
    
    for script in orchestrator_scripts:
        if script.exists():
            syntax_issues = check_script_syntax(script)
            issues.extend(syntax_issues)
            
            # Check for config block
            content = script.read_text()
            if "# ============================================================================" in content and "CONFIGURATION" in content:
                config_ok = True
            else:
                config_ok = False
                warnings.append(f"{script.name}: Missing config block")
            
            if syntax_issues:
                print(f"  ❌ {script.name}: Syntax errors")
            elif not config_ok:
                print(f"  ⚠️  {script.name}: Missing config block")
            else:
                print(f"  ✅ {script.name}: Valid syntax and config block")
        else:
            warnings.append(f"Script missing: {script}")
    
    print()
    
    # Phase 4: Import Verification
    print("Phase 4: Import Verification")
    print("-" * 80)
    
    # Check __all__ exports
    try:
        from metainformant.rna import __all__ as rna_all
        
        for export_name in rna_all:
            try:
                exec(f"from metainformant.rna import {export_name}")
                print(f"  ✅ {export_name}: Importable")
            except Exception as e:
                issues.append(f"Cannot import {export_name}: {e}")
                print(f"  ❌ {export_name}: Import failed")
    except Exception as e:
        issues.append(f"Cannot check __all__ exports: {e}")
        print(f"  ❌ Cannot check exports: {e}")
    
    print()
    
    # Phase 5: Documentation Completeness
    print("Phase 5: Documentation Completeness")
    print("-" * 80)
    
    # Check module README exists
    module_readme = SRC_DIR / "README.md"
    if module_readme.exists():
        print(f"  ✅ Module README exists: {module_readme}")
    else:
        warnings.append(f"Module README missing: {module_readme}")
        print(f"  ⚠️  Module README missing")
    
    # Check orchestrator docs exist
    orchestrator_docs = [
        DOCS_DIR / "ORCHESTRATION" / "ENA_WORKFLOW.md",
        DOCS_DIR / "ORCHESTRATION" / "MULTI_SPECIES.md",
        DOCS_DIR / "ORCHESTRATION" / "BATCH_DOWNLOAD.md",
    ]
    
    for doc in orchestrator_docs:
        if doc.exists():
            print(f"  ✅ {doc.name}: Exists")
        else:
            warnings.append(f"Orchestrator doc missing: {doc}")
            print(f"  ❌ {doc.name}: Missing")
    
    print()
    print("=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Issues found: {len(issues)}")
    print(f"Warnings: {len(warnings)}")
    
    if issues:
        print("\nIssues:")
        for issue in issues[:20]:
            print(f"  - {issue}")
        if len(issues) > 20:
            print(f"  ... and {len(issues) - 20} more")
    
    if warnings:
        print("\nWarnings:")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more")
    
    if not issues and not warnings:
        print("\n✅ All verifications passed!")
        return 0
    elif issues:
        print("\n❌ Verification failed with issues")
        return 1
    else:
        print("\n⚠️  Verification passed with warnings")
        return 0


if __name__ == "__main__":
    import importlib.util
    sys.exit(main())

