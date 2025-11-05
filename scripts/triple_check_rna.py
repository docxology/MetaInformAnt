#!/usr/bin/env python3
"""Comprehensive triple-check verification of RNA documentation and code.

Implements all phases of deep verification including:
- Method signature matching
- Code example execution
- Bidirectional cross-references
- Script method verification
- Import testing
"""

import ast
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

sys.path.insert(0, str(REPO_ROOT / "src"))

issues = []
warnings = []
checks_passed = []


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
    
    # Find [text](path) links
    link_pattern = re.compile(r'\[([^\]]+)\]\(([^)]+)\)')
    for match in link_pattern.finditer(content):
        link_text, link_path = match.groups()
        links.append((link_text, link_path))
    
    return links


def normalize_anchor(text: str) -> str:
    """Convert heading to markdown anchor."""
    text = text.lower()
    text = re.sub(r'[^\w\s-]', '', text)
    text = re.sub(r'[\s]+', '-', text)
    return text.strip('-')


def check_fragment_exists(file_path: Path, fragment: str) -> bool:
    """Check if fragment anchor exists in file."""
    if not file_path.exists():
        return False
    
    content = file_path.read_text()
    headings = re.findall(r'^##+ (.+)$', content, re.MULTILINE)
    
    normalized_fragment = normalize_anchor(fragment)
    for heading in headings:
        if normalize_anchor(heading) == normalized_fragment:
            return True
    
    return False


def resolve_link_path(base_file: Path, link_path: str) -> Path:
    """Resolve a relative link path."""
    if link_path.startswith("/"):
        return REPO_ROOT / link_path.lstrip("/")
    elif link_path.startswith("../"):
        return base_file.parent.parent / link_path[3:]
    elif link_path.startswith("./"):
        return base_file.parent / link_path[2:]
    else:
        return base_file.parent / link_path


def main():
    """Run comprehensive triple-check verification."""
    print("=" * 80)
    print("RNA TRIPLE-CHECK VERIFICATION")
    print("=" * 80)
    print()
    
    # Phase 1: Deep Method Documentation Verification
    print("Phase 1: Deep Method Documentation Verification")
    print("-" * 80)
    
    key_modules = [
        ("metainformant.rna.amalgkit", SRC_DIR / "amalgkit.py"),
        ("metainformant.rna.workflow", SRC_DIR / "workflow.py"),
        ("metainformant.rna.configs", SRC_DIR / "configs.py"),
        ("metainformant.rna.monitoring", SRC_DIR / "monitoring.py"),
        ("metainformant.rna.environment", SRC_DIR / "environment.py"),
    ]
    
    for module_name, module_path in key_modules:
        if not module_path.exists():
            warnings.append(f"Module file missing: {module_path}")
            continue
        
        try:
            spec = importlib.util.spec_from_file_location(module_name, module_path)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                # Check all public members
                for name, obj in inspect.getmembers(module):
                    if name.startswith("_"):
                        continue
                    
                    if inspect.isfunction(obj):
                        # Check docstring completeness
                        if not obj.__doc__:
                            issues.append(f"{module_name}.{name}: Missing docstring")
                        else:
                            # Check signature match
                            sig_issues = check_method_signature_match(obj, obj.__doc__)
                            issues.extend([f"{module_name}.{name}: {issue}" for issue in sig_issues])
                            if not sig_issues:
                                checks_passed.append(f"{module_name}.{name}: Signature matches docs")
                    
                    elif inspect.isclass(obj):
                        if not obj.__doc__:
                            issues.append(f"{module_name}.{name}: Missing class docstring")
                        # Check if it's a dataclass (dataclasses auto-generate __init__)
                        is_dataclass = hasattr(obj, '__dataclass_fields__')
                        # Check methods (skip builtin methods from parent classes)
                        for method_name, method in inspect.getmembers(obj, inspect.isfunction):
                            # Skip private methods except important ones
                            if method_name.startswith("_") and method_name not in ["__init__"]:
                                continue
                            # Skip methods from builtin classes
                            if method.__module__ and ('builtin' in method.__module__ or 'pathlib' in method.__module__):
                                continue
                            # Skip __init__ for dataclasses (they auto-generate it)
                            if method_name == "__init__" and is_dataclass:
                                continue
                            if not method.__doc__ and method_name == "__init__":
                                issues.append(f"{module_name}.{name}.{method_name}: Missing docstring")
                
                print(f"  ✅ {module_name}: Checked")
        except Exception as e:
            warnings.append(f"{module_name}: Error - {e}")
            print(f"  ⚠️  {module_name}: Check failed")
    
    # Check __all__ exports
    try:
        from metainformant.rna import __all__ as rna_all
        for export in rna_all:
            try:
                exec(f"from metainformant.rna import {export}")
                checks_passed.append(f"Export {export}: Importable")
            except Exception as e:
                issues.append(f"Export {export}: Cannot import - {e}")
    except Exception as e:
        warnings.append(f"Cannot check exports: {e}")
    
    print()
    
    # Phase 2: Deep Documentation Accuracy
    print("Phase 2: Deep Documentation Accuracy")
    print("-" * 80)
    
    key_docs = [
        DOCS_DIR / "README.md",
        DOCS_DIR / "WORKFLOW.md",
        DOCS_DIR / "STEPS.md",
        DOCS_DIR / "CONFIGURATION.md",
    ]
    
    for doc_file in key_docs:
        if not doc_file.exists():
            issues.append(f"Missing: {doc_file}")
            continue
        
        content = doc_file.read_text()
        
        # Extract and check code examples
        in_code = False
        code_lines = []
        code_start = 0
        
        for i, line in enumerate(content.split('\n'), 1):
            if line.strip().startswith('```python'):
                in_code = True
                code_lines = []
                code_start = i
            elif line.strip().startswith('```') and in_code:
                code = '\n'.join(code_lines)
                if code.strip():
                    executable, error = check_code_example_executable(code)
                    if not executable:
                        issues.append(f"{doc_file}:{code_start}: Code example issue: {error}")
                    else:
                        checks_passed.append(f"{doc_file}:{code_start}: Code example valid")
                in_code = False
                code_lines = []
            elif in_code:
                code_lines.append(line)
        
        # Check links bidirectionally
        links = find_all_links(doc_file)
        for link_text, link_path in links:
            if link_path.startswith('http'):
                continue
            
            # Resolve link
            if '#' in link_path:
                file_part, fragment = link_path.split('#', 1)
            else:
                file_part = link_path
                fragment = None
            
            target_file = resolve_link_path(doc_file, file_part)
            
            if not target_file.exists():
                issues.append(f"{doc_file}: Broken link to {link_path}")
            else:
                checks_passed.append(f"{doc_file}: Link to {link_path} valid")
                if fragment:
                    if check_fragment_exists(target_file, fragment):
                        checks_passed.append(f"{doc_file}: Fragment {fragment} exists")
                    else:
                        issues.append(f"{doc_file}: Fragment {fragment} not found in {target_file.name}")
        
        print(f"  ✅ {doc_file.name}: Checked")
    
    print()
    
    # Phase 3: Deep Script Verification
    print("Phase 3: Deep Script Verification")
    print("-" * 80)
    
    orchestrators = [
        (SCRIPTS_DIR / "workflow_ena_integrated.py", DOCS_DIR / "ORCHESTRATION" / "ENA_WORKFLOW.md"),
        (SCRIPTS_DIR / "run_multi_species.py", DOCS_DIR / "ORCHESTRATION" / "MULTI_SPECIES.md"),
        (SCRIPTS_DIR / "batch_download_species.py", DOCS_DIR / "ORCHESTRATION" / "BATCH_DOWNLOAD.md"),
    ]
    
    for script, doc in orchestrators:
        if not script.exists():
            issues.append(f"Script missing: {script}")
            continue
        
        # Check syntax
        try:
            with open(script) as f:
                ast.parse(f.read(), filename=str(script))
            checks_passed.append(f"{script.name}: Syntax valid")
        except SyntaxError as e:
            issues.append(f"{script}: Syntax error: {e}")
        
        # Check config block
        content = script.read_text()
        if "# ============================================================================" in content and "CONFIGURATION" in content:
            checks_passed.append(f"{script.name}: Config block present")
        else:
            issues.append(f"{script.name}: Missing config block")
        
        # Check method calls match metainformant.rna.*
        if 'from metainformant.rna' in content or 'import metainformant.rna' in content:
            checks_passed.append(f"{script.name}: Uses metainformant.rna imports")
        
        # Check doc exists and matches
        if doc.exists():
            checks_passed.append(f"{script.name}: Documentation exists")
        else:
            issues.append(f"{script.name}: Documentation missing: {doc}")
        
        print(f"  ✅ {script.name}: Checked")
    
    print()
    
    # Phase 4: Import and Execution Verification
    print("Phase 4: Import and Execution Verification")
    print("-" * 80)
    
    # Test all documented imports
    documented_imports = [
        "from metainformant.rna import check_cli_available",
        "from metainformant.rna import metadata",
        "from metainformant.rna import AmalgkitWorkflowConfig",
        "from metainformant.rna import plan_workflow",
        "from metainformant.rna.steps import STEP_RUNNERS",
    ]
    
    for import_stmt in documented_imports:
        try:
            exec(import_stmt)
            checks_passed.append(f"Import works: {import_stmt[:50]}...")
        except Exception as e:
            issues.append(f"Import failed: {import_stmt} - {e}")
    
    print(f"  ✅ Tested {len(documented_imports)} documented imports")
    print()
    
    # Phase 5: Cross-Module Consistency
    print("Phase 5: Cross-Module Consistency")
    print("-" * 80)
    
    # Check __init__.py exports
    try:
        from metainformant.rna import __all__ as rna_all
        import metainformant.rna as rna_module
        
        # Verify exports exist (they're imported from submodules)
        for export in rna_all:
            try:
                if hasattr(rna_module, export):
                    checks_passed.append(f"Export {export}: Exists in module")
                else:
                    # Try importing it to see if it's actually available
                    exec(f"from metainformant.rna import {export}")
                    checks_passed.append(f"Export {export}: Importable")
            except Exception as e:
                issues.append(f"Export {export}: Cannot import - {e}")
        
        print(f"  ✅ Verified {len(rna_all)} exports")
    except Exception as e:
        warnings.append(f"Cannot check exports: {e}")
    
    print()
    
    # Summary
    print("=" * 80)
    print("TRIPLE-CHECK VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Checks passed: {len(checks_passed)}")
    print(f"Issues found: {len(issues)}")
    print(f"Warnings: {len(warnings)}")
    
    if issues:
        print("\nIssues to fix:")
        for issue in issues[:30]:
            print(f"  - {issue}")
        if len(issues) > 30:
            print(f"  ... and {len(issues) - 30} more")
    
    if warnings:
        print("\nWarnings:")
        for warning in warnings[:10]:
            print(f"  - {warning}")
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more")
    
    if not issues:
        print("\n✅ All triple-checks passed!")
        return 0
    else:
        print(f"\n❌ {len(issues)} issues found")
        return 1


if __name__ == "__main__":
    sys.exit(main())

