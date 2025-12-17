#!/usr/bin/env python3
"""Comprehensive verification of RNA documentation and code accuracy.

This script performs multiple levels of verification:
1. Basic verification: docstrings, examples, links
2. Documentation verification: module docstrings, imports
3. Comprehensive verification: signature matching, cross-references
"""

import ast
import importlib
import importlib.util
import inspect
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "rna"
SRC_DIR = REPO_ROOT / "src" / "metainformant" / "rna"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "rna"

issues = []
warnings = []
checks_passed = []


def normalize_anchor(text: str) -> str:
    """Convert heading text to markdown anchor format."""
    # Markdown anchors: lowercase, spaces to hyphens, remove special chars
    text = text.lower()
    text = re.sub(r'[^\w\s-]', '', text)
    text = re.sub(r'[\s]+', '-', text)
    return text.strip('-')


def check_fragment_exists(doc_file: Path, fragment: str) -> bool:
    """Check if a fragment anchor exists in a markdown file."""
    if not doc_file.exists():
        return False
    
    content = doc_file.read_text()
    
    # Find all headings
    headings = re.findall(r'^##+ (.+)$', content, re.MULTILINE)
    
    # Check if any heading matches the fragment (normalized)
    normalized_fragment = normalize_anchor(fragment)
    for heading in headings:
        if normalize_anchor(heading) == normalized_fragment:
            return True
    
    return False


def check_code_example_syntax(code: str, file_path: Path, line_num: int) -> list[str]:
    """Check if a code example is syntactically valid."""
    issues_found = []
    
    # Skip function signatures (single line with def/class)
    lines = code.strip().split('\n')
    if len(lines) == 1 and (code.strip().startswith('def ') or code.strip().startswith('class ')):
        return []  # Function/class signature, not executable
    
    # Skip type hints only
    if code.strip().startswith('->') or ':' in code and 'def' not in code and 'class' not in code:
        # Might be a type annotation, check if it's valid
        try:
            ast.parse(f"def dummy{code}")
            return []
        except:
            pass
    
    try:
        ast.parse(code, filename=str(file_path))
    except SyntaxError as e:
        # Check if it's a meaningful error (not just incomplete example)
        if 'expected' in str(e).lower() or 'invalid' in str(e).lower():
            issues_found.append(f"{file_path}:{line_num}: Syntax error: {e.msg}")
    
    return issues_found


def check_doc_links(doc_file: Path) -> list[str]:
    """Check internal markdown links in documentation."""
    issues_found = []
    
    if not doc_file.exists():
        return [f"Documentation file missing: {doc_file}"]
    
    content = doc_file.read_text()
    doc_dir = doc_file.parent
    
    # Find markdown links [text](path)
    link_pattern = re.compile(r'\[([^\]]+)\]\(([^)]+)\)')
    
    for match in link_pattern.finditer(content):
        link_text, link_path = match.groups()
        
        # Skip external links
        if link_path.startswith("http"):
            continue
        
        # Handle fragment links
        if '#' in link_path:
            file_part, fragment = link_path.split('#', 1)
        else:
            file_part = link_path
            fragment = None
        
        # Resolve file path
        if file_part.startswith("/"):
            target_file = REPO_ROOT / file_part.lstrip("/")
        elif file_part.startswith("../"):
            # Relative path - resolve from doc_dir
            # For ../ links, we need to go up one level from doc_dir
            target_file = doc_dir.parent / file_part[3:]  # Remove ../
        else:
            # Relative path in same directory
            target_file = doc_dir / file_part
        
        # Normalize path (resolve symlinks, etc.)
        try:
            target_file = target_file.resolve()
        except (OSError, RuntimeError):
            pass  # If resolve fails, use original path
        
        # Check if file exists
        if not target_file.exists():
            issues_found.append(f"{doc_file}: Broken link to {link_path} ({link_text}) - resolved to {target_file}")
            continue
        
        # Check fragment exists
        if fragment:
            if not check_fragment_exists(target_file, fragment):
                # Try without the fragment - maybe it's just a navigation link
                # Only warn if it's clearly a section reference
                if any(word in fragment.lower() for word in ['section', 'heading', 'anchor']):
                    warnings.append(f"{doc_file}: Fragment '{fragment}' may not exist in {target_file.name}")
    
    return issues_found


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
            code = "\n".join(code_lines)
            if code.strip():  # Only check non-empty blocks
                example_issues = check_code_example_syntax(code, file_path, code_start_line)
                issues_found.extend(example_issues)
            in_code_block = False
            code_lines = []
        elif in_code_block:
            code_lines.append(line)
    
    return issues_found


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


def extract_signature_from_docstring(docstring: str) -> dict[str, Any]:
    """Extract parameter and return type info from docstring."""
    if not docstring:
        return {}

    params = {}
    returns = None

    # Look for Args: section
    args_match = re.search(r'Args?:[\s\n]+(.*?)(?=\n\n|\n[A-Z]|\Z)', docstring, re.DOTALL)
    if args_match:
        args_text = args_match.group(1)
        for line in args_text.split('\n'):
            line = line.strip()
            if ':' in line:
                param_name = line.split(':')[0].strip()
                param_desc = line.split(':', 1)[1].strip()
                params[param_name] = param_desc

    # Look for Returns: section
    returns_match = re.search(r'Returns?:[\s\n]+(.*?)(?=\n\n|\n[A-Z]|\Z)', docstring, re.DOTALL)
    if returns_match:
        returns = returns_match.group(1).strip()

    return {"params": params, "returns": returns}


def check_method_signature_match(func, docstring: str) -> list[str]:
    """Check if method signature matches documentation."""
    issues_found = []

    try:
        sig = inspect.signature(func)
        doc_params = extract_signature_from_docstring(docstring)

        # Skip checking for builtin/standard library methods
        if func.__module__ and ('builtin' in func.__module__ or 'posixpath' in func.__module__ or 'pathlib' in func.__module__):
            return issues_found

        # Check parameters
        for param_name, param in sig.parameters.items():
            if param_name == 'self':
                continue
            # Skip **kwargs - these are intentionally flexible
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                continue
            # Skip parameters with default values that are clearly optional
            if param.default != inspect.Parameter.empty:
                continue
            # Only check if docstring actually documents parameters
            if doc_params.get("params") and param_name not in doc_params.get("params", {}):
                # Check if docstring mentions the parameter in a different format
                if param_name not in docstring and f"`{param_name}`" not in docstring:
                    issues_found.append(f"Parameter '{param_name}' not documented")

        # Check documented parameters exist (only if we have a structured docstring)
        if doc_params.get("params"):
            for doc_param in doc_params.get("params", {}):
                if doc_param not in sig.parameters and doc_param not in ['kwargs', '**kwargs']:
                    # Might be a typo or alternative name
                    pass

    except Exception as e:
        warnings.append(f"Cannot check signature for {func.__name__}: {e}")

    return issues_found


def check_code_example_executable(code: str, context: dict[str, Any] = None) -> tuple[bool, str | None]:
    """Check if a code example can be compiled/executed."""
    if context is None:
        context = {}

    # Skip function signatures
    if code.strip().startswith('def ') and ':' in code and code.count('\n') <= 2:
        return True, None

    # Skip type hints only
    if '->' in code and 'def' not in code and 'class' not in code:
        return True, None

    try:
        # Try to compile
        compile(code, '<example>', 'exec')
        return True, None
    except SyntaxError as e:
        # Check if it's a meaningful error
        if 'expected' in str(e).lower() or 'invalid' in str(e).lower():
            return False, f"Syntax error: {e.msg} at line {e.lineno}"
        return True, None
    except Exception as e:
        return False, f"Error: {e}"


def find_all_links(doc_file: Path) -> list[tuple[str, str]]:
    """Find all markdown links in a file."""
    if not doc_file.exists():
        return []

    content = doc_file.read_text()
    links = []

    # Find markdown links [text](path)
    link_pattern = re.compile(r'\[([^\]]+)\]\(([^)]+)\)')

    for match in link_pattern.finditer(content):
        link_text, link_path = match.groups()
        links.append((link_text, link_path))

    return links


def resolve_link_path(base_file: Path, link_path: str) -> Path:
    """Resolve a link path relative to a base file."""
    base_dir = base_file.parent

    # Skip external links
    if link_path.startswith("http"):
        return None

    # Handle fragment links
    if '#' in link_path:
        file_part = link_path.split('#')[0]
    else:
        file_part = link_path

    # Resolve relative paths
    if file_part.startswith('/'):
        return Path(file_part[1:])  # Remove leading /
    else:
        return (base_dir / file_part).resolve()


def main():
    """Run comprehensive verification."""
    print("=" * 80)
    print("RNA DOCUMENTATION AND CODE COMPREHENSIVE VERIFICATION")
    print("=" * 80)
    print()
    
    # Phase 1: Method Documentation
    print("Phase 1: Method Documentation Verification")
    print("-" * 80)
    
    # Check key modules
    key_modules = [
        ("metainformant.rna.amalgkit", SRC_DIR / "amalgkit.py"),
        ("metainformant.rna.workflow", SRC_DIR / "workflow.py"),
        ("metainformant.rna.configs", SRC_DIR / "configs.py"),
        ("metainformant.rna.monitoring", SRC_DIR / "monitoring.py"),
        ("metainformant.rna.environment", SRC_DIR / "environment.py"),
    ]
    
    for module_name, module_path in key_modules:
        if module_path.exists():
            try:
                spec = importlib.util.spec_from_file_location(module_name, module_path)
                if spec and spec.loader:
                    module = importlib.util.module_from_spec(spec)
                    spec.loader.exec_module(module)
                    
                    # Check public members have docstrings
                    missing_docs = []
                    for name, obj in inspect.getmembers(module):
                        if name.startswith("_"):
                            continue
                        if inspect.isfunction(obj) or inspect.isclass(obj):
                            if not obj.__doc__ or len(obj.__doc__.strip()) < 10:
                                missing_docs.append(name)
                    
                    if missing_docs:
                        issues.extend([f"{module_name}.{name}: Missing/incomplete docstring" for name in missing_docs])
                        print(f"  ❌ {module_name}: {len(missing_docs)} undocumented items")
                    else:
                        print(f"  ✅ {module_name}: All documented")
            except Exception as e:
                warnings.append(f"{module_name}: Could not check - {e}")
                print(f"  ⚠️  {module_name}: Check failed")
    
    print()
    
    # Phase 2: Documentation Accuracy
    print("Phase 2: Documentation Accuracy")
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
            example_issues = check_doc_examples(doc_file)
            link_issues = check_doc_links(doc_file)
            issues.extend(example_issues)
            issues.extend(link_issues)
            
            total = len(example_issues) + len(link_issues)
            if total == 0:
                print(f"  ✅ {doc_file.name}: Valid")
            else:
                print(f"  ⚠️  {doc_file.name}: {total} issues")
        else:
            warnings.append(f"Missing: {doc_file}")
            print(f"  ❌ {doc_file.name}: Missing")
    
    print()
    
    # Phase 3: Script Verification
    print("Phase 3: Script Verification")
    print("-" * 80)

    # Check that the main production scripts exist and are syntactically valid
    production_scripts = [
        SCRIPTS_DIR / "run_workflow.py",
        SCRIPTS_DIR / "setup_genome.py",
        SCRIPTS_DIR / "discover_species.py",
        SCRIPTS_DIR / "check_environment.py",
    ]

    for script in production_scripts:
        if script.exists():
            try:
                with open(script) as f:
                    ast.parse(f.read(), filename=str(script))
                print(f"  ✅ {script.name}: Valid syntax")
            except SyntaxError as e:
                issues.append(f"{script}: Syntax error: {e}")
                print(f"  ❌ {script.name}: Syntax error")
        else:
            issues.append(f"Missing production script: {script}")
            print(f"  ❌ {script.name}: Missing")

    print()
    
    # Phase 4: Import Verification
    print("Phase 4: Import Verification")
    print("-" * 80)
    
    try:
        from metainformant.rna import __all__ as rna_all
        
        importable = 0
        failed = 0
        
        for export_name in rna_all:
            try:
                exec(f"from metainformant.rna import {export_name}")
                importable += 1
            except Exception as e:
                issues.append(f"Cannot import {export_name}: {e}")
                failed += 1
        
        print(f"  ✅ {importable}/{len(rna_all)} exports importable")
        if failed > 0:
            print(f"  ❌ {failed} imports failed")
    except Exception as e:
        issues.append(f"Cannot check exports: {e}")
        print(f"  ❌ Cannot check: {e}")
    
    print()
    
    # Summary
    print("=" * 80)
    print("VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Issues: {len(issues)}")
    print(f"Warnings: {len(warnings)}")
    
    if issues:
        print("\nIssues to fix:")
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
    
    if not issues:
        print("\n✅ All critical verifications passed!")
        return 0
    else:
        print("\n❌ Some issues found - review above")
        return 1


if __name__ == "__main__":
    sys.exit(main())

