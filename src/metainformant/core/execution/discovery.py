"""Symbolic mapping and context discovery utilities for METAINFORMANT.

Provides functions for discovering functions, configs, output patterns,
call graphs, and cross-module relationships across the repository.
"""

from __future__ import annotations

import ast
import re
import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class FunctionInfo:
    """Information about a discovered function."""

    name: str
    signature: str
    file_path: Path
    line_number: int
    docstring: str | None = None
    decorators: list[str] = field(default_factory=list)
    parameters: list[str] = field(default_factory=list)
    return_type: str | None = None


@dataclass
class ConfigInfo:
    """Information about a discovered config file."""

    path: Path
    domain: str | None = None
    format: str = "yaml"  # yaml, toml, json
    size: int = 0
    modified_time: float = 0.0


@dataclass
class OutputPattern:
    """Information about output directory patterns."""

    module: str
    base_pattern: str
    subdirs: list[str] = field(default_factory=list)
    examples: list[str] = field(default_factory=list)


@dataclass
class SymbolUsage:
    """Information about symbol usage."""

    file: Path
    line: int
    column: int
    context: str | None = None


@dataclass
class ModuleDependency:
    """Information about module dependencies."""

    module: str
    imports: list[str] = field(default_factory=list)
    from_imports: dict[str, list[str]] = field(default_factory=dict)


class _DiscoveryCache:
    """Thread-safe mtime-based cache for parsed AST results.

    Stores parsed AST trees keyed by file path.  On cache hit, compares the
    file's current mtime to the cached mtime; if the file has been modified
    the entry is evicted and the caller re-parses.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        # path_str -> (mtime, ast.Module)
        self._entries: dict[str, tuple[float, ast.Module]] = {}

    def get(self, path: Path) -> ast.Module | None:
        """Return cached AST if the file has not been modified since caching."""
        key = str(path)
        with self._lock:
            entry = self._entries.get(key)
            if entry is None:
                return None
            cached_mtime, tree = entry
            try:
                current_mtime = path.stat().st_mtime
            except OSError:
                # File disappeared; evict.
                self._entries.pop(key, None)
                return None
            if current_mtime > cached_mtime:
                # Stale; evict and signal re-parse.
                self._entries.pop(key, None)
                return None
            return tree

    def put(self, path: Path, tree: ast.Module) -> None:
        """Store a parsed AST with the file's current mtime."""
        key = str(path)
        try:
            mtime = path.stat().st_mtime
        except OSError:
            return  # Cannot stat; do not cache.
        with self._lock:
            self._entries[key] = (mtime, tree)

    def invalidate(self, path: Path) -> None:
        """Remove a single path from the cache."""
        with self._lock:
            self._entries.pop(str(path), None)

    def clear(self) -> None:
        """Remove all entries from the cache."""
        with self._lock:
            self._entries.clear()


# Module-level singleton used by discover_functions / find_symbol_usage.
_discovery_cache = _DiscoveryCache()


def invalidate_discovery_cache(path: Path | str | None = None) -> None:
    """Public helper to invalidate the discovery AST cache.

    Args:
        path: If given, only that file is invalidated.  If ``None`` the
              entire cache is cleared.
    """
    if path is None:
        _discovery_cache.clear()
    else:
        _discovery_cache.invalidate(Path(path))


def _cached_parse(path: Path) -> ast.Module:
    """Parse a Python file, using the mtime-based cache when possible.

    Raises:
        SyntaxError: If the file cannot be parsed.
        OSError: If the file cannot be read.
    """
    tree = _discovery_cache.get(path)
    if tree is not None:
        return tree
    with open(path, "rt", encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=str(path))
    _discovery_cache.put(path, tree)
    return tree


def _get_cache_dir(repo_root: Path) -> Path:
    """Get cache directory for discovery results."""
    cache_dir = Path(repo_root) / "output" / ".discovery_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


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

    # Handle *args and **kwargs
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


def discover_functions(module_path: str | Path, pattern: str | None = None) -> list[FunctionInfo]:
    """Discover all functions in a module with their signatures.

    Args:
        module_path: Path to Python module file
        pattern: Optional regex pattern to filter function names

    Returns:
        List of FunctionInfo objects for each discovered function

    Raises:
        FileNotFoundError: If module_path does not exist
        SyntaxError: If module cannot be parsed
    """
    module_path = Path(module_path)
    if not module_path.exists():
        raise FileNotFoundError(f"Module not found: {module_path}")

    try:
        tree = _cached_parse(module_path)
    except SyntaxError as e:
        raise SyntaxError(f"Failed to parse {module_path}: {e}") from e

    functions: list[FunctionInfo] = []
    pattern_re = re.compile(pattern) if pattern else None

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            if pattern_re and not pattern_re.search(node.name):
                continue

            # Extract decorators
            decorators = []
            for decorator in node.decorator_list:
                try:
                    if isinstance(decorator, ast.Name):
                        decorators.append(decorator.id)
                    elif isinstance(decorator, ast.Attribute):
                        decorators.append(ast.unparse(decorator))
                    else:
                        # Try to unparse complex decorators
                        try:
                            decorators.append(ast.unparse(decorator))
                        except Exception:
                            decorators.append("<complex>")
                except Exception:
                    decorators.append("<unknown>")

            # Extract parameters
            parameters = [arg.arg for arg in node.args.args]

            # Extract return type
            return_type = None
            if node.returns:
                try:
                    return_type = ast.unparse(node.returns)
                except Exception:
                    # Fallback for complex return types
                    if isinstance(node.returns, ast.Name):
                        return_type = node.returns.id
                    else:
                        return_type = "..."

            # Extract docstring
            docstring = ast.get_docstring(node)

            # Get signature
            try:
                signature = _parse_function_signature(node)
            except Exception:
                # Fallback to simple signature
                signature = f"({', '.join(parameters)})"

            func_info = FunctionInfo(
                name=node.name,
                signature=signature,
                file_path=module_path,
                line_number=node.lineno,
                docstring=docstring,
                decorators=decorators,
                parameters=parameters,
                return_type=return_type,
            )
            functions.append(func_info)

    return functions


def discover_configs(repo_root: str | Path, domain: str | None = None) -> list[ConfigInfo]:
    """Discover all config files in the repository.

    Args:
        repo_root: Root directory of the repository
        domain: Optional domain filter (e.g., 'rna', 'gwas')

    Returns:
        List of ConfigInfo objects for each discovered config file
    """
    repo_root = Path(repo_root)
    config_dir = repo_root / "config"
    if not config_dir.exists():
        return []

    configs: list[ConfigInfo] = []

    # Search in config directory
    for config_file in config_dir.rglob("*"):
        if not config_file.is_file():
            continue

        suffix = config_file.suffix.lower()
        if suffix not in {".yaml", ".yml", ".toml", ".json"}:
            continue

        # Extract domain from path
        rel_path = config_file.relative_to(config_dir)
        path_parts = rel_path.parts
        file_domain = path_parts[0] if len(path_parts) > 1 else None

        # Filter by domain if specified
        if domain and file_domain != domain:
            continue

        try:
            stat = config_file.stat()
            config_info = ConfigInfo(
                path=config_file,
                domain=file_domain,
                format=suffix.lstrip("."),
                size=stat.st_size,
                modified_time=stat.st_mtime,
            )
            configs.append(config_info)
        except OSError:
            # Skip files that can't be accessed
            continue

    return sorted(configs, key=lambda x: x.path)


def discover_output_patterns(module_name: str) -> OutputPattern:
    """Discover output directory patterns for a module.

    Args:
        module_name: Name of the module (e.g., 'rna', 'gwas', 'dna')

    Returns:
        OutputPattern with base pattern and subdirectories
    """
    # Map from .cursorrules patterns
    patterns_map = {
        "rna": OutputPattern(
            module="rna",
            base_pattern="output/amalgkit/<species>/<step>/",
            subdirs=["quant/", "work/", "logs/"],
            examples=["output/amalgkit/Apis_mellifera/quant/", "output/amalgkit/Apis_mellifera/work/"],
        ),
        "gwas": OutputPattern(
            module="gwas",
            base_pattern="output/gwas/<analysis_type>/",
            subdirs=["association/", "plots/", "qc/"],
            examples=["output/gwas/association/", "output/gwas/plots/"],
        ),
        "life_events": OutputPattern(
            module="life_events",
            base_pattern="output/life_events/<workflow>/",
            subdirs=["embeddings/", "models/", "plots/"],
            examples=["output/life_events/embeddings/", "output/life_events/models/"],
        ),
        "multiomics": OutputPattern(
            module="multiomics",
            base_pattern="output/multiomics/<integration>/",
            subdirs=["integrated/", "plots/"],
            examples=["output/multiomics/integrated/", "output/multiomics/plots/"],
        ),
        "singlecell": OutputPattern(
            module="singlecell",
            base_pattern="output/singlecell/<analysis>/",
            subdirs=["preprocessing/", "clustering/"],
            examples=["output/singlecell/preprocessing/", "output/singlecell/clustering/"],
        ),
        "networks": OutputPattern(
            module="networks",
            base_pattern="output/networks/<network_type>/",
            subdirs=["ppi/", "regulatory/", "pathways/"],
            examples=["output/networks/ppi/", "output/networks/regulatory/"],
        ),
        "information": OutputPattern(
            module="information",
            base_pattern="output/information/<analysis_type>/",
            subdirs=["entropy/", "complexity/"],
            examples=["output/information/entropy/", "output/information/complexity/"],
        ),
        "dna": OutputPattern(
            module="dna",
            base_pattern="output/dna/<analysis_type>/",
            subdirs=["phylogeny/", "population/", "variants/"],
            examples=["output/dna/phylogeny/", "output/dna/population/"],
        ),
        "protein": OutputPattern(
            module="protein",
            base_pattern="output/protein/<analysis_type>/",
            subdirs=["structures/", "alignments/"],
            examples=["output/protein/structures/", "output/protein/alignments/"],
        ),
        "quality": OutputPattern(
            module="quality",
            base_pattern="output/quality/<dataset>/",
            subdirs=["fastq/", "reports/"],
            examples=["output/quality/fastq/", "output/quality/reports/"],
        ),
        "math": OutputPattern(
            module="math",
            base_pattern="output/math/<type>/",
            subdirs=["simulations/", "models/", "plots/"],
            examples=["output/math/simulations/", "output/math/models/"],
        ),
        "simulation": OutputPattern(
            module="simulation",
            base_pattern="output/simulation/<type>/",
            subdirs=["sequences/", "ecosystems/"],
            examples=["output/simulation/sequences/", "output/simulation/ecosystems/"],
        ),
    }

    # Return known pattern or default
    if module_name.lower() in patterns_map:
        return patterns_map[module_name.lower()]

    # Default pattern
    return OutputPattern(
        module=module_name,
        base_pattern=f"output/{module_name}/",
        subdirs=[],
        examples=[f"output/{module_name}/"],
    )


def build_call_graph(entry_point: str | Path, repo_root: str | Path | None = None) -> dict[str, list[str]]:
    """Build function call graph from an entry point.

    Args:
        entry_point: Path to entry point module or function name
        repo_root: Root directory of repository (for resolving imports)

    Returns:
        Dictionary mapping function names to lists of called functions
    """
    # This is a simplified implementation
    # Full implementation would require more sophisticated AST traversal
    entry_path = Path(entry_point)
    if not entry_path.exists():
        return {}

    try:
        tree = _cached_parse(entry_path)
    except (SyntaxError, OSError):
        return {}

    call_graph: dict[str, list[str]] = {}

    class CallVisitor(ast.NodeVisitor):
        current_function: str | None = None

        def visit_FunctionDef(self, node: ast.FunctionDef):
            self.current_function = node.name
            call_graph[self.current_function] = []
            self.generic_visit(node)
            self.current_function = None

        def visit_Call(self, node: ast.Call):
            if self.current_function:
                if isinstance(node.func, ast.Name):
                    call_graph[self.current_function].append(node.func.id)
                elif isinstance(node.func, ast.Attribute):
                    # Handle method calls like obj.method()
                    call_graph[self.current_function].append(node.func.attr)
            self.generic_visit(node)

    visitor = CallVisitor()
    visitor.visit(tree)

    return call_graph


def find_symbol_usage(symbol_name: str, repo_root: str | Path) -> list[SymbolUsage]:
    """Find all usages of a symbol across the repository.

    Args:
        symbol_name: Name of symbol to find
        repo_root: Root directory of repository

    Returns:
        List of SymbolUsage objects for each occurrence
    """
    repo_root = Path(repo_root)
    usages: list[SymbolUsage] = []

    # Search Python files
    for py_file in repo_root.rglob("*.py"):
        if "__pycache__" in str(py_file) or ".pyc" in str(py_file):
            continue

        try:
            # Read source lines (needed for context snippets regardless of cache).
            with open(py_file, "rt", encoding="utf-8") as f:
                content = f.read()
                lines = content.splitlines()

            try:
                tree = _cached_parse(py_file)
            except SyntaxError:
                # Skip files with syntax errors
                continue

            class UsageVisitor(ast.NodeVisitor):
                def visit_Name(self, node: ast.Name) -> None:
                    if node.id == symbol_name:
                        # Get context line
                        line = lines[node.lineno - 1] if node.lineno <= len(lines) else ""
                        usage = SymbolUsage(
                            file=py_file,
                            line=node.lineno,
                            column=node.col_offset,
                            context=line.strip()[:100] if line else None,
                        )
                        usages.append(usage)
                    self.generic_visit(node)

            visitor = UsageVisitor()
            visitor.visit(tree)
        except (OSError, UnicodeDecodeError):
            continue

    return usages


def get_module_dependencies(module_path: str | Path) -> ModuleDependency:
    """Extract import dependencies from a module.

    Args:
        module_path: Path to Python module file

    Returns:
        ModuleDependency with imports and from_imports
    """
    module_path = Path(module_path)
    if not module_path.exists():
        return ModuleDependency(module=str(module_path), imports=[], from_imports={})

    try:
        tree = _cached_parse(module_path)
    except (SyntaxError, OSError):
        return ModuleDependency(module=str(module_path), imports=[], from_imports={})

    imports: list[str] = []
    from_imports: dict[str, list[str]] = {}

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imports.append(alias.name)
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            names = [alias.name for alias in node.names]
            if module not in from_imports:
                from_imports[module] = []
            from_imports[module].extend(names)

    return ModuleDependency(module=str(module_path), imports=imports, from_imports=from_imports)


def discover_workflows(repo_root: str | Path | None = None) -> list[dict[str, Any]]:
    """Discover all workflow entry points and their configs.

    Args:
        repo_root: Root directory of repository (defaults to current directory)

    Returns:
        List of dictionaries with workflow information
    """
    if repo_root is None:
        repo_root = Path.cwd()
    else:
        repo_root = Path(repo_root)

    workflows: list[dict[str, Any]] = []

    # Known workflow patterns
    workflow_files = [
        ("rna", "src/metainformant/rna/workflow.py", "execute_workflow", "config/amalgkit/"),
        ("gwas", "src/metainformant/gwas/workflow.py", "run_gwas_workflow", "config/gwas/"),
        ("life_events", "src/metainformant/life_events/workflow.py", "analyze_life_course", "config/"),
        ("simulation", "src/metainformant/simulation/workflow.py", "run_sequence_simulation_workflow", "config/"),
    ]

    for domain, workflow_path, entry_point, config_pattern in workflow_files:
        full_path = repo_root / workflow_path
        if full_path.exists():
            # Find configs for this workflow
            config_dir = repo_root / config_pattern
            configs = []
            if config_dir.exists():
                for config_file in config_dir.rglob("*.yaml"):
                    if config_file.is_file() and "template" not in config_file.name:
                        configs.append(str(config_file.relative_to(repo_root)))

            workflows.append(
                {
                    "domain": domain,
                    "workflow_file": str(workflow_path),
                    "entry_point": entry_point,
                    "configs": configs,
                    "exists": True,
                }
            )

    return workflows
