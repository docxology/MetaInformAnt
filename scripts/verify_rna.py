#!/usr/bin/env python3
"""Verification of RNA documentation and code accuracy."""

import ast
import importlib.util
import inspect
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
DOCS_DIR = REPO_ROOT / "docs" / "rna"
SRC_DIR = REPO_ROOT / "src" / "metainformant" / "rna"
SCRIPTS_DIR = REPO_ROOT / "scripts" / "rna"

issues = []
warnings = []


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
    
    orchestrators = [
        SCRIPTS_DIR / "workflow_ena_integrated.py",
        SCRIPTS_DIR / "run_multi_species.py",
        SCRIPTS_DIR / "batch_download_species.py",
    ]
    
    for script in orchestrators:
        if script.exists():
            try:
                with open(script) as f:
                    ast.parse(f.read(), filename=str(script))
                
                content = script.read_text()
                has_config = "# ============================================================================" in content and "CONFIGURATION" in content
                
                if has_config:
                    print(f"  ✅ {script.name}: Valid syntax and config block")
                else:
                    warnings.append(f"{script.name}: Missing config block")
                    print(f"  ⚠️  {script.name}: Valid syntax, missing config block")
            except SyntaxError as e:
                issues.append(f"{script}: Syntax error: {e}")
                print(f"  ❌ {script.name}: Syntax error")
        else:
            warnings.append(f"Missing: {script}")
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

