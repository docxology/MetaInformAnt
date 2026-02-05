"""Symbol indexing and cross-referencing utilities for METAINFORMANT.

Provides functions for indexing functions, classes, and other symbols
across the repository, enabling symbol lookup and reference finding.
"""

from __future__ import annotations

import ast
import difflib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from metainformant.core.io.io import load_json, dump_json


@dataclass
class SymbolDefinition:
    """Information about a symbol definition."""

    name: str
    symbol_type: str  # 'function', 'class', 'variable', etc.
    file_path: Path
    line_number: int
    signature: str | None = None
    docstring: str | None = None
    module: str | None = None


@dataclass
class SymbolReference:
    """Information about a symbol reference."""

    symbol_name: str
    file_path: Path
    line_number: int
    column: int
    context: str | None = None
    is_definition: bool = False


def _get_cache_dir(repo_root: Path) -> Path:
    """Get cache directory for symbol index."""
    cache_dir = Path(repo_root) / "output" / ".discovery_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def _get_cache_file(repo_root: Path, index_type: str) -> Path:
    """Get cache file path for symbol index."""
    cache_dir = _get_cache_dir(repo_root)
    return cache_dir / f"symbol_index_{index_type}.json"


def _parse_function_signature(node: ast.FunctionDef) -> str:
    """Parse function signature from AST node."""
    args = []
    defaults_start = len(node.args.args) - len(node.args.defaults)

    for i, arg in enumerate(node.args.args):
        arg_str = arg.arg
        if arg.annotation:
            try:
                arg_str += f": {ast.unparse(arg.annotation)}"
            except Exception:
                # Fallback for complex annotations
                if isinstance(arg.annotation, ast.Name):
                    arg_str += f": {arg.annotation.id}"
                else:
                    arg_str += ": ..."
        # Add default value if present
        if i >= defaults_start:
            default_idx = i - defaults_start
            try:
                default_str = ast.unparse(node.args.defaults[default_idx])
                arg_str += f" = {default_str}"
            except Exception:
                arg_str += " = ..."
        args.append(arg_str)

    if node.args.vararg:
        vararg_name = node.args.vararg.arg
        if node.args.vararg.annotation:
            try:
                vararg_name += f": {ast.unparse(node.args.vararg.annotation)}"
            except Exception:
                pass
        args.append(f"*{vararg_name}")
    if node.args.kwarg:
        kwarg_name = node.args.kwarg.arg
        if node.args.kwarg.annotation:
            try:
                kwarg_name += f": {ast.unparse(node.args.kwarg.annotation)}"
            except Exception:
                pass
        args.append(f"**{kwarg_name}")

    sig = f"({', '.join(args)})"
    if node.returns:
        try:
            sig += f" -> {ast.unparse(node.returns)}"
        except Exception:
            # Fallback for complex return types
            if isinstance(node.returns, ast.Name):
                sig += f" -> {node.returns.id}"
            else:
                sig += " -> ..."
    return sig


def _parse_class_signature(node: ast.ClassDef) -> str:
    """Parse class signature from AST node."""
    bases = []
    for base in node.bases:
        try:
            bases.append(ast.unparse(base))
        except Exception:
            if isinstance(base, ast.Name):
                bases.append(base.id)
    return f"({', '.join(bases)})" if bases else "()"


def index_functions(repo_root: str | Path, use_cache: bool = True) -> dict[str, list[SymbolDefinition]]:
    """Build function signature index across repository.

    Args:
        repo_root: Root directory of repository
        use_cache: Whether to use cached index if available

    Returns:
        Dictionary mapping function names to lists of SymbolDefinition objects
    """
    repo_root = Path(repo_root)
    cache_file = _get_cache_file(repo_root, "functions")

    # Try to load from cache
    if use_cache and cache_file.exists():
        try:
            cached_data = load_json(cache_file)
            # Check if cache is still valid (simplified - could check mtimes)
            if cached_data:
                index: dict[str, list[SymbolDefinition]] = {}
                for name, defs in cached_data.items():
                    index[name] = [SymbolDefinition(**def_data) for def_data in defs]
                return index
        except Exception:
            # Cache invalid, rebuild
            pass

    index: dict[str, list[SymbolDefinition]] = {}

    # Scan Python files
    for py_file in repo_root.rglob("*.py"):
        if "__pycache__" in str(py_file) or ".pyc" in str(py_file):
            continue

        try:
            with open(py_file, "rt", encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=str(py_file))
        except (SyntaxError, OSError, UnicodeDecodeError):
            continue

        # Extract module name
        rel_path = py_file.relative_to(repo_root)
        module_parts = list(rel_path.parts[:-1]) + [rel_path.stem]
        module = ".".join(module_parts)

        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                # Skip private functions if desired
                if node.name.startswith("_"):
                    continue

                try:
                    signature = _parse_function_signature(node)
                except Exception:
                    signature = None

                docstring = ast.get_docstring(node)

                symbol_def = SymbolDefinition(
                    name=node.name,
                    symbol_type="function",
                    file_path=py_file,
                    line_number=node.lineno,
                    signature=signature,
                    docstring=docstring,
                    module=module,
                )

                if node.name not in index:
                    index[node.name] = []
                index[node.name].append(symbol_def)

    # Cache the index
    if use_cache:
        cache_data = {}
        for name, defs in index.items():
            cache_data[name] = [
                {
                    "name": d.name,
                    "symbol_type": d.symbol_type,
                    "file_path": str(d.file_path),
                    "line_number": d.line_number,
                    "signature": d.signature,
                    "docstring": d.docstring,
                    "module": d.module,
                }
                for d in defs
            ]
        dump_json(cache_data, cache_file)

    return index


def index_classes(repo_root: str | Path, use_cache: bool = True) -> dict[str, list[SymbolDefinition]]:
    """Build class index across repository.

    Args:
        repo_root: Root directory of repository
        use_cache: Whether to use cached index if available

    Returns:
        Dictionary mapping class names to lists of SymbolDefinition objects
    """
    repo_root = Path(repo_root)
    cache_file = _get_cache_file(repo_root, "classes")

    # Try to load from cache
    if use_cache and cache_file.exists():
        try:
            cached_data = load_json(cache_file)
            if cached_data:
                index: dict[str, list[SymbolDefinition]] = {}
                for name, defs in cached_data.items():
                    index[name] = [SymbolDefinition(**def_data) for def_data in defs]
                return index
        except Exception:
            pass

    index: dict[str, list[SymbolDefinition]] = {}

    # Scan Python files
    for py_file in repo_root.rglob("*.py"):
        if "__pycache__" in str(py_file) or ".pyc" in str(py_file):
            continue

        try:
            with open(py_file, "rt", encoding="utf-8") as f:
                tree = ast.parse(f.read(), filename=str(py_file))
        except (SyntaxError, OSError, UnicodeDecodeError):
            continue

        # Extract module name
        rel_path = py_file.relative_to(repo_root)
        module_parts = list(rel_path.parts[:-1]) + [rel_path.stem]
        module = ".".join(module_parts)

        for node in ast.walk(tree):
            if isinstance(node, ast.ClassDef):
                # Skip private classes if desired
                if node.name.startswith("_"):
                    continue

                try:
                    signature = _parse_class_signature(node)
                except Exception:
                    signature = None

                docstring = ast.get_docstring(node)

                symbol_def = SymbolDefinition(
                    name=node.name,
                    symbol_type="class",
                    file_path=py_file,
                    line_number=node.lineno,
                    signature=signature,
                    docstring=docstring,
                    module=module,
                )

                if node.name not in index:
                    index[node.name] = []
                index[node.name].append(symbol_def)

    # Cache the index
    if use_cache:
        cache_data = {}
        for name, defs in index.items():
            cache_data[name] = [
                {
                    "name": d.name,
                    "symbol_type": d.symbol_type,
                    "file_path": str(d.file_path),
                    "line_number": d.line_number,
                    "signature": d.signature,
                    "docstring": d.docstring,
                    "module": d.module,
                }
                for d in defs
            ]
        dump_json(cache_data, cache_file)

    return index


def find_symbol(
    symbol_name: str, symbol_type: str = "function", repo_root: str | Path | None = None
) -> list[SymbolDefinition]:
    """Find symbol definition(s) by name.

    Args:
        symbol_name: Name of symbol to find
        symbol_type: Type of symbol ('function', 'class', or 'all')
        repo_root: Root directory of repository

    Returns:
        List of SymbolDefinition objects matching the symbol name
    """
    if repo_root is None:
        repo_root = Path.cwd()
    else:
        repo_root = Path(repo_root)

    results: list[SymbolDefinition] = []

    if symbol_type in ("function", "all"):
        func_index = index_functions(repo_root)
        if symbol_name in func_index:
            results.extend(func_index[symbol_name])

    if symbol_type in ("class", "all"):
        class_index = index_classes(repo_root)
        if symbol_name in class_index:
            results.extend(class_index[symbol_name])

    return results


def get_symbol_signature(symbol_path: str | Path, symbol_name: str) -> str | None:
    """Get full signature of a symbol from its file.

    Args:
        symbol_path: Path to file containing the symbol
        symbol_name: Name of the symbol

    Returns:
        Signature string or None if not found
    """
    symbol_path = Path(symbol_path)
    if not symbol_path.exists():
        return None

    try:
        with open(symbol_path, "rt", encoding="utf-8") as f:
            tree = ast.parse(f.read(), filename=str(symbol_path))
    except (SyntaxError, OSError):
        return None

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name == symbol_name:
            if isinstance(node, ast.FunctionDef):
                return _parse_function_signature(node)
            elif isinstance(node, ast.ClassDef):
                return _parse_class_signature(node)

    return None


def find_symbol_references(symbol_name: str, repo_root: str | Path) -> list[SymbolReference]:
    """Find all references to a symbol across the repository.

    Args:
        symbol_name: Name of symbol to find references for
        repo_root: Root directory of repository

    Returns:
        List of SymbolReference objects for each occurrence
    """
    repo_root = Path(repo_root)
    references: list[SymbolReference] = []

    # First, find definitions
    definitions = find_symbol(symbol_name, "all", repo_root)
    for defn in definitions:
        ref = SymbolReference(
            symbol_name=symbol_name,
            file_path=defn.file_path,
            line_number=defn.line_number,
            column=0,
            context=None,
            is_definition=True,
        )
        references.append(ref)

    # Then find usages
    for py_file in repo_root.rglob("*.py"):
        if "__pycache__" in str(py_file) or ".pyc" in str(py_file):
            continue

        try:
            with open(py_file, "rt", encoding="utf-8") as f:
                content = f.read()
                lines = content.splitlines()
                tree = ast.parse(content, filename=str(py_file))
        except (SyntaxError, OSError, UnicodeDecodeError):
            continue

        class ReferenceVisitor(ast.NodeVisitor):
            def visit_Name(self, node: ast.Name):
                if node.id == symbol_name:
                    # Skip if this is a definition (already handled above)
                    # We can't easily check parent context without more complex traversal,
                    # so we'll include all Name nodes and let the caller filter if needed
                    line = lines[node.lineno - 1] if node.lineno <= len(lines) else ""
                    ref = SymbolReference(
                        symbol_name=symbol_name,
                        file_path=py_file,
                        line_number=node.lineno,
                        column=node.col_offset,
                        context=line.strip()[:100] if line else None,
                        is_definition=False,
                    )
                    references.append(ref)

        visitor = ReferenceVisitor()
        visitor.visit(tree)

    return references


def get_symbol_metadata(symbol_path: str | Path, symbol_name: str) -> dict[str, Any]:
    """Get metadata (docstring, type hints, etc.) for a symbol.

    Args:
        symbol_path: Path to file containing the symbol
        symbol_name: Name of the symbol

    Returns:
        Dictionary with symbol metadata
    """
    symbol_path = Path(symbol_path)
    if not symbol_path.exists():
        return {}

    try:
        with open(symbol_path, "rt", encoding="utf-8") as f:
            tree = ast.parse(f.read(), filename=str(symbol_path))
    except (SyntaxError, OSError):
        return {}

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name == symbol_name:
            metadata: dict[str, Any] = {
                "name": node.name,
                "type": "function" if isinstance(node, ast.FunctionDef) else "class",
                "line_number": node.lineno,
                "docstring": ast.get_docstring(node),
            }

            if isinstance(node, ast.FunctionDef):
                metadata["signature"] = _parse_function_signature(node)
                metadata["parameters"] = [arg.arg for arg in node.args.args]
                if node.returns:
                    try:
                        metadata["return_type"] = ast.unparse(node.returns)
                    except Exception:
                        pass
            elif isinstance(node, ast.ClassDef):
                metadata["signature"] = _parse_class_signature(node)
                metadata["bases"] = []
                for base in node.bases:
                    try:
                        metadata["bases"].append(ast.unparse(base))
                    except Exception:
                        if isinstance(base, ast.Name):
                            metadata["bases"].append(base.id)

            return metadata

    return {}


def fuzzy_find_symbol(
    symbol_name: str, symbol_type: str = "function", repo_root: str | Path | None = None, threshold: float = 0.6
) -> list[tuple[str, float]]:
    """Find symbols with fuzzy matching.

    Args:
        symbol_name: Name of symbol to find (fuzzy match)
        symbol_type: Type of symbol ('function', 'class', or 'all')
        repo_root: Root directory of repository
        threshold: Similarity threshold (0.0 to 1.0)

    Returns:
        List of tuples (symbol_name, similarity_score) sorted by similarity
    """
    if repo_root is None:
        repo_root = Path.cwd()
    else:
        repo_root = Path(repo_root)

    all_symbols: set[str] = set()

    if symbol_type in ("function", "all"):
        func_index = index_functions(repo_root)
        all_symbols.update(func_index.keys())

    if symbol_type in ("class", "all"):
        class_index = index_classes(repo_root)
        all_symbols.update(class_index.keys())

    # Use difflib for fuzzy matching
    matches = difflib.get_close_matches(symbol_name, all_symbols, n=10, cutoff=threshold)

    # Calculate similarity scores
    results: list[tuple[str, float]] = []
    for match in matches:
        ratio = difflib.SequenceMatcher(None, symbol_name.lower(), match.lower()).ratio()
        results.append((match, ratio))

    # Sort by similarity (descending)
    results.sort(key=lambda x: x[1], reverse=True)

    return results
