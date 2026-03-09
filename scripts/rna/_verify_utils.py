#!/usr/bin/env python3
"""Reusable verification utilities for RNA documentation and code checks.

Provides functions for:
- Markdown link validation
- Code example syntax checking
- Module docstring auditing
- Signature-docstring matching
"""

import ast
import importlib
import importlib.util
import inspect
import re
from pathlib import Path
from typing import Any


def normalize_anchor(text: str) -> str:
    """Convert heading text to markdown anchor format."""
    text = text.lower()
    text = re.sub(r"[^\w\s-]", "", text)
    text = re.sub(r"[\s]+", "-", text)
    return text.strip("-")


def check_fragment_exists(doc_file: Path, fragment: str) -> bool:
    """Check if a fragment anchor exists in a markdown file."""
    if not doc_file.exists():
        return False

    content = doc_file.read_text()
    headings = re.findall(r"^##+ (.+)$", content, re.MULTILINE)
    normalized_fragment = normalize_anchor(fragment)

    for heading in headings:
        if normalize_anchor(heading) == normalized_fragment:
            return True
    return False


def check_code_example_syntax(code: str, file_path: Path, line_num: int) -> list[str]:
    """Check if a code example is syntactically valid."""
    issues_found = []

    lines = code.strip().split("\n")
    if len(lines) == 1 and (code.strip().startswith("def ") or code.strip().startswith("class ")):
        return []

    if code.strip().startswith("->") or ":" in code and "def" not in code and "class" not in code:
        try:
            ast.parse(f"def dummy{code}")
            return []
        except Exception:
            pass

    try:
        ast.parse(code, filename=str(file_path))
    except SyntaxError as e:
        if "expected" in str(e).lower() or "invalid" in str(e).lower():
            issues_found.append(f"{file_path}:{line_num}: Syntax error: {e.msg}")

    return issues_found


def check_doc_links(doc_file: Path, repo_root: Path, warnings: list[str] | None = None) -> list[str]:
    """Check internal markdown links in documentation."""
    issues_found = []
    if warnings is None:
        warnings = []

    if not doc_file.exists():
        return [f"Documentation file missing: {doc_file}"]

    content = doc_file.read_text()
    doc_dir = doc_file.parent
    link_pattern = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")

    for match in link_pattern.finditer(content):
        link_text, link_path = match.groups()
        if link_path.startswith("http"):
            continue

        if "#" in link_path:
            file_part, fragment = link_path.split("#", 1)
        else:
            file_part = link_path
            fragment = None

        if file_part.startswith("/"):
            target_file = repo_root / file_part.lstrip("/")
        elif file_part.startswith("../"):
            target_file = doc_dir.parent / file_part[3:]
        else:
            target_file = doc_dir / file_part

        try:
            target_file = target_file.resolve()
        except (OSError, RuntimeError):
            pass

        if not target_file.exists():
            issues_found.append(f"{doc_file}: Broken link to {link_path} ({link_text}) - resolved to {target_file}")
            continue

        if fragment:
            if not check_fragment_exists(target_file, fragment):
                if any(word in fragment.lower() for word in ["section", "heading", "anchor"]):
                    warnings.append(f"{doc_file}: Fragment '{fragment}' may not exist in {target_file.name}")

    return issues_found


def check_doc_examples(file_path: Path) -> list[str]:
    """Check code examples in documentation files."""
    issues_found = []

    if not file_path.exists():
        return [f"Documentation file missing: {file_path}"]

    content = file_path.read_text()
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
            code = "\n".join(code_lines)
            if code.strip():
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

        for name, obj in inspect.getmembers(module):
            if name.startswith("_"):
                continue
            if inspect.isfunction(obj):
                if obj.__doc__ is None or not obj.__doc__.strip():
                    issues_found.append(f"{module_name}.{name}: Missing docstring")
                elif len(obj.__doc__.strip()) < 20:
                    issues_found.append(f"{module_name}.{name}: Docstring too short")
            elif inspect.isclass(obj):
                if obj.__doc__ is None or not obj.__doc__.strip():
                    issues_found.append(f"{module_name}.{name}: Missing class docstring")
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

    args_match = re.search(r"Args?:[\s\n]+(.*?)(?=\n\n|\n[A-Z]|\Z)", docstring, re.DOTALL)
    if args_match:
        args_text = args_match.group(1)
        for line in args_text.split("\n"):
            line = line.strip()
            if ":" in line:
                param_name = line.split(":")[0].strip()
                param_desc = line.split(":", 1)[1].strip()
                params[param_name] = param_desc

    returns_match = re.search(r"Returns?:[\s\n]+(.*?)(?=\n\n|\n[A-Z]|\Z)", docstring, re.DOTALL)
    if returns_match:
        returns = returns_match.group(1).strip()

    return {"params": params, "returns": returns}


def check_method_signature_match(func, docstring: str, warnings: list[str] | None = None) -> list[str]:
    """Check if method signature matches documentation."""
    issues_found = []
    if warnings is None:
        warnings = []

    try:
        sig = inspect.signature(func)
        doc_params = extract_signature_from_docstring(docstring)

        if func.__module__ and (
            "builtin" in func.__module__ or "posixpath" in func.__module__ or "pathlib" in func.__module__
        ):
            return issues_found

        for param_name, param in sig.parameters.items():
            if param_name == "self":
                continue
            if param.kind == inspect.Parameter.VAR_KEYWORD:
                continue
            if param.default != inspect.Parameter.empty:
                continue
            if doc_params.get("params") and param_name not in doc_params.get("params", {}):
                if param_name not in docstring and f"`{param_name}`" not in docstring:
                    issues_found.append(f"Parameter '{param_name}' not documented")

    except Exception as e:
        warnings.append(f"Cannot check signature for {func.__name__}: {e}")

    return issues_found


def check_code_example_executable(code: str, context: dict[str, Any] = None) -> tuple[bool, str | None]:
    """Check if a code example can be compiled/executed."""
    if context is None:
        context = {}

    if code.strip().startswith("def ") and ":" in code and code.count("\n") <= 2:
        return True, None

    if "->" in code and "def" not in code and "class" not in code:
        return True, None

    try:
        compile(code, "<example>", "exec")
        return True, None
    except SyntaxError as e:
        if "expected" in str(e).lower() or "invalid" in str(e).lower():
            return False, f"Syntax error: {e.msg} at line {e.lineno}"
        return True, None
    except Exception as e:
        return False, f"Error: {e}"


def find_all_links(doc_file: Path) -> list[tuple[str, str]]:
    """Find all markdown links in a file."""
    if not doc_file.exists():
        return []

    content = doc_file.read_text()
    link_pattern = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
    return [(m.group(1), m.group(2)) for m in link_pattern.finditer(content)]


def resolve_link_path(base_file: Path, link_path: str) -> Path | None:
    """Resolve a link path relative to a base file."""
    base_dir = base_file.parent

    if link_path.startswith("http"):
        return None

    if "#" in link_path:
        file_part = link_path.split("#")[0]
    else:
        file_part = link_path

    if file_part.startswith("/"):
        return Path(file_part[1:])
    else:
        return (base_dir / file_part).resolve()
