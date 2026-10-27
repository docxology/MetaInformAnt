"""
Systematic Cross-Code Verification for METAINFORMANT Documentation

Verifies that all code examples, imports, function calls, and CLI commands
in documentation files match the actual source code.

The thin orchestrator ``scripts/verify_documentation_code.py`` builds the CLI
and delegates here.

Target: Check all 506+ documentation files for broken code references.
"""

from __future__ import annotations

import argparse
import ast
import importlib
import importlib.util
import inspect
import logging
import os
import re
import sys
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Dict, Iterator, List, Optional, Set
from urllib.parse import unquote

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Repo layout: src/metainformant/quality/doc_verification.py → repo root
REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_REPORT_PATH = Path("output") / "cross_code_verification_report.md"
KNOWN_OPTIONAL_IMPORT_ROOTS = {
    "anndata",
    "boto3",
    "botocore",
    "cloudscraper",
    "dask",
    "cryptography",
    "dotenv",
    "fastapi",
    "h5py",
    "hypothesis",
    "imblearn",
    "polars",
    "plotly",
    "prometheus_client",
    "psycopg2",
    "pysam",
    "pydantic",
    "scanpy",
    "tenacity",
    "torch",
    "umap",
    "uvicorn",
    "vcfpy",
}
HISTORICAL_REPORT_NAMES = {
    "AUDIT_INDEX.md",
    "DOCUMENTATION_AUDIT_REPORT.md",
    "DOCUMENTATION_EXECUTIVE_SUMMARY.md",
    "DOCUMENTATION_REVIEW_REPORT.md",
    "GWAS_VALIDATION_REPORT.md",
    "LINK_VALIDATION_REPORT.md",
    "PHARMACOGENOMICS_VALIDATION_REPORT.md",
    "SINGLECELL_VALIDATION_REPORT.md",
    "VALIDATION_REPORT.md",
    "cross_code_verification_report.md",
}


def resolve_repo_path(path: Path, repo_root: Path = REPO_ROOT) -> Path:
    """Resolve an absolute or repo-relative path."""
    return path.resolve() if path.is_absolute() else (repo_root / path).resolve()


def collect_code_examples(
    parser: "DocumentationParser", doc_files: List[Path], verbose: bool = False
) -> List[CodeExample]:
    """Extract code examples from documentation files."""
    all_examples = []
    for doc_file in doc_files:
        examples = parser.extract_code_examples(doc_file)
        all_examples.extend(examples)
        if verbose:
            logger.debug(f"  {doc_file.name}: {len(examples)} code examples")
    return all_examples


@dataclass
class Symbol:
    """Represents a Python symbol (module, class, function, method)."""

    name: str
    module_path: str
    full_name: str
    line: int
    kind: str  # 'module', 'class', 'function', 'method'
    docstring: str = ""
    signature: str = ""


@dataclass
class CodeExample:
    """Represents a code example found in documentation."""

    file_path: str
    line_number: int
    code_type: str  # 'python', 'bash', 'inline'
    code: str
    language: str = ""


@dataclass
class Violation:
    """Represents a broken code reference."""

    doc_file: str
    line_number: int
    code_example: str
    issue_type: str  # 'ImportError', 'AttributeError', 'ModuleNotFound', 'WrongParameter', 'NonExistentCommand'
    details: str
    severity: str = "error"  # 'error', 'warning', 'info'

    def __str__(self) -> str:
        return f"{self.doc_file}:{self.line_number} - {self.issue_type}: {self.details}"


class PythonSymbolIndex:
    """Builds and maintains an index of all Python symbols in the source code."""

    def __init__(self, src_dir: Path):
        self.src_dir = src_dir.resolve()
        self.symbols: Dict[str, Symbol] = {}  # full_name -> Symbol
        self.modules: Dict[str, Path] = {}  # module_name -> file_path
        self.class_methods: Dict[str, Set[str]] = defaultdict(set)  # class_name -> {method_names}
        self.imports_cache: Dict[str, Set[str]] = defaultdict(set)  # module -> {imported_names}
        self.module_exports: Dict[str, Set[str]] = defaultdict(set)  # module -> {public exported names}

    def _module_name_for_file(self, file_path: Path) -> str:
        """Return the importable module name for a Python file under src_dir."""
        rel_parts = list(file_path.resolve().relative_to(self.src_dir).parts)
        if rel_parts[-1] == "__init__.py":
            rel_parts = rel_parts[:-1]
        elif rel_parts[-1].endswith(".py"):
            rel_parts[-1] = rel_parts[-1][:-3]

        src_parts = list(self.src_dir.parts)
        if "metainformant" in src_parts:
            package_start = src_parts.index("metainformant")
            module_parts = src_parts[package_start:] + rel_parts
        else:
            module_parts = rel_parts

        return ".".join(module_parts)

    def build_index(self) -> PythonSymbolIndex:
        """Scan all Python files and build symbol index."""
        logger.info(f"Building Python symbol index from {self.src_dir}...")

        python_files = list(self.src_dir.rglob("*.py"))
        logger.info(f"Found {len(python_files)} Python files")

        for py_file in python_files:
            self._parse_file(py_file)

        logger.info(f"Indexed {len(self.symbols)} symbols across {len(self.modules)} modules")
        return self

    def _parse_file(self, file_path: Path) -> None:
        """Parse a single Python file and extract symbols."""
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            tree = ast.parse(content, filename=str(file_path))

            module_name = self._module_name_for_file(file_path)

            # Register module
            self.modules[module_name] = file_path

            # Extract imports (for reference)
            self._extract_imports(tree, module_name)

            # Walk AST and collect symbols
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    self._add_symbol(node, "class", module_name, file_path)
                    if not node.name.startswith("_"):
                        self.module_exports[module_name].add(node.name)
                    # Track methods for this class
                    for item in node.body:
                        if isinstance(item, ast.FunctionDef):
                            self.class_methods[node.name].add(item.name)

                elif isinstance(node, ast.FunctionDef):
                    # Top-level function
                    if isinstance(node.parent, ast.Module) if hasattr(node, "parent") else True:
                        self._add_symbol(node, "function", module_name, file_path)
                        if not node.name.startswith("_"):
                            self.module_exports[module_name].add(node.name)

        except Exception as e:
            logger.warning(f"Failed to parse {file_path}: {e}")

    def _add_symbol(
        self, node: ast.ClassDef | ast.AsyncFunctionDef | ast.FunctionDef, kind: str, module_name: str, file_path: Path
    ) -> None:
        """Add a symbol to the index."""
        name = node.name
        full_name = f"{module_name}.{name}"
        signature = self._get_signature(node)

        docstring = ast.get_docstring(node) or ""

        symbol = Symbol(
            name=name,
            module_path=str(file_path),
            full_name=full_name,
            line=node.lineno,
            kind=kind,
            docstring=docstring,
            signature=signature,
        )
        self.symbols[full_name] = symbol

    def _get_signature(self, node: ast.AST) -> str:
        """Extract function/method signature."""
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = []
            for arg in node.args.args:
                args.append(arg.arg)
            if node.args.vararg:
                args.append(f"*{node.args.vararg.arg}")
            if node.args.kwarg:
                args.append(f"**{node.args.kwarg.arg}")
            return f"({', '.join(args)})"
        elif isinstance(node, ast.ClassDef):
            bases = []
            for base in node.bases:
                if isinstance(base, ast.Name):
                    bases.append(base.id)
            if bases:
                return f"({', '.join(bases)})"
        return ""

    def _extract_imports(self, tree: ast.Module, module_name: str) -> None:
        """Extract imports from a module."""
        imported = set()
        explicit_all = self._extract_all_exports(tree)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imported.add(alias.name)
                    exported_name = alias.asname or alias.name.split(".", 1)[0]
                    if not exported_name.startswith("_"):
                        self.module_exports[module_name].add(exported_name)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    for alias in node.names:
                        imported.add(f"{node.module}.{alias.name}" if alias.name != "*" else f"{node.module}.*")
                        if alias.name != "*":
                            exported_name = alias.asname or alias.name
                            if not exported_name.startswith("_"):
                                self.module_exports[module_name].add(exported_name)
        if explicit_all:
            self.module_exports[module_name].update(explicit_all)
        self.imports_cache[module_name] = imported

    @staticmethod
    def _extract_all_exports(tree: ast.AST) -> Set[str]:
        """Extract simple string names from a module-level __all__ assignment."""
        exports: Set[str] = set()
        for node in getattr(tree, "body", []):
            if not isinstance(node, (ast.Assign, ast.AnnAssign)):
                continue

            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if not any(isinstance(target, ast.Name) and target.id == "__all__" for target in targets):
                continue

            value = node.value
            if isinstance(value, (ast.List, ast.Tuple, ast.Set)):
                for element in value.elts:
                    if isinstance(element, ast.Constant) and isinstance(element.value, str):
                        exports.add(element.value)
            elif isinstance(value, ast.ListComp):
                continue
        return exports

    def get_symbol(self, full_name: str) -> Optional[Symbol]:
        """Look up a symbol by its full name."""
        return self.symbols.get(full_name)

    def get_module_members(self, module_name: str) -> Set[str]:
        """Get all public members exported by a module."""
        members = set(self.module_exports.get(module_name, set()))
        for sym_full_name, sym in self.symbols.items():
            if sym_full_name.startswith(f"{module_name}.") and sym_full_name.count(".") == module_name.count(".") + 1:
                if not sym.name.startswith("_"):
                    members.add(sym.name)
        return members

    def module_exists(self, module_name: str) -> bool:
        """Check if a module exists."""
        return module_name in self.modules


class DocumentationParser:
    """Parses documentation files and extracts code examples."""

    # Regex patterns for inline code
    INLINE_CODE_PATTERN = r"`([^`]+)`"

    # Code block patterns: ```python, ```python-snippet, ```bash, ```, etc.
    CODE_BLOCK_PATTERN = r"```([A-Za-z0-9_-]*)\n(.*?)```"

    def __init__(self, docs_dir: Path, include_historical: bool = False):
        self.docs_dir = docs_dir
        self.include_historical = include_historical

    def find_all_docs(self) -> List[Path]:
        """Find all markdown documentation files."""
        docs = [
            doc_path
            for doc_path in self.docs_dir.rglob("*.md")
            if self.include_historical or not self._is_historical_report(doc_path)
        ]
        logger.info(f"Found {len(docs)} documentation files")
        return docs

    def _is_historical_report(self, file_path: Path) -> bool:
        """Return True for retained snapshot reports that are not current docs."""
        try:
            head = "\n".join(file_path.read_text(encoding="utf-8").splitlines()[:12])
        except OSError:
            return False
        if "Historical snapshot" in head:
            return True

        try:
            rel = file_path.resolve().relative_to(REPO_ROOT)
        except ValueError:
            rel = file_path
        return len(rel.parts) == 1 and file_path.name in HISTORICAL_REPORT_NAMES

    def extract_code_examples(self, file_path: Path) -> List[CodeExample]:
        """Extract all code examples from a markdown file."""
        examples = []

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            content = "".join(lines)

            # Find all code blocks
            for match in re.finditer(self.CODE_BLOCK_PATTERN, content, re.DOTALL):
                lang = match.group(1).strip()
                code = textwrap.dedent(match.group(2)).strip()
                lang = self._classify_block_language(lang, code)
                # Find line number by counting newlines before the match
                line_num = content[: match.start()].count("\n") + 1

                examples.append(
                    CodeExample(
                        file_path=str(file_path),
                        line_number=line_num,
                        code_type="block",
                        language=lang,
                        code=code.strip(),
                    )
                )

            # Find inline code (but not inside code blocks)
            # Simple approach: split by code blocks first, then search remainder
            parts = re.split(r"```.*?```", content, flags=re.DOTALL)
            current_line = 0
            for part in parts:
                for match in re.finditer(r"`([^`]+)`", part):
                    inline_code = match.group(1)
                    line_num = part[: match.start()].count("\n") + current_line + 1

                    # Filter: only inline code that looks like Python or shell
                    if self._is_code_candidate(inline_code):
                        examples.append(
                            CodeExample(
                                file_path=str(file_path),
                                line_number=line_num,
                                code_type="inline",
                                language=self._guess_language(inline_code),
                                code=inline_code,
                            )
                        )
                current_line += part.count("\n")

        except Exception as e:
            logger.warning(f"Failed to parse {file_path}: {e}")

        return examples

    def _is_code_candidate(self, code: str) -> bool:
        """Check if inline code looks like actual code (not just a filename)."""
        code = code.strip()
        if "\n" in code or len(code) > 160:
            return False
        if "..." in code:
            return False
        # Skip single words, filenames, paths
        if len(code.split()) == 1 and ("." in code or "/" in code or code.endswith(".md")):
            return False
        return bool(
            re.match(r"^(from\s+[A-Za-z_][\w.]*\s+import\s+|import\s+[A-Za-z_][\w.]*)", code)
            or code.startswith(("metainformant ", "uv ", "pytest ", "python ", "python3 "))
        )

    def _guess_language(self, code: str) -> str:
        """Guess the language of a code snippet."""
        code = code.strip()
        # Python indicators
        if any(kw in code for kw in ["import ", "def ", "class ", "from ", "print("]):
            return "python"
        # Shell/bash indicators
        if any(kw in code for kw in ["$", "bash", "sh ", "&&", "||", ";"]):
            return "bash"
        return "unknown"

    def _classify_block_language(self, lang: str, code: str) -> str:
        """Return a verification language for a fenced code block."""
        normalized = (lang or "").lower()
        if normalized != "python":
            return lang

        stripped = code.strip()
        if not stripped:
            return "text"
        if stripped.startswith((">>>", "... ")) or "Traceback (most recent call last):" in stripped:
            return "pycon"
        if re.match(r"^[A-Za-z_][\w.]*Error:", stripped):
            return "text"
        if stripped.startswith("rule "):
            return "snakemake"
        if any(line.lstrip().startswith(("!", "%")) for line in stripped.splitlines()):
            return "python-notebook"
        if re.search(r"\b(import|from)\s+\.\.\.", stripped) or re.search(r"\bimport\s+\.\.\.", stripped):
            return "python-snippet"
        try:
            ast.parse(stripped)
        except SyntaxError:
            return "python-snippet"
        return "python"


class CodeValidator:
    """Validates code examples against the source code index."""

    def __init__(
        self,
        symbol_index: PythonSymbolIndex,
        pyproject_toml: Path,
        *,
        strict_optional_imports: bool = False,
    ):
        self.symbol_index = symbol_index
        self.violations: List[Violation] = []
        self.strict_optional_imports = strict_optional_imports
        self.entry_points = self._load_entry_points(pyproject_toml)

    def _load_entry_points(self, pyproject: Path) -> Set[str]:
        """Load CLI entry points from pyproject.toml."""
        entry_points = set()
        if sys.version_info >= (3, 11):
            import tomllib
        else:
            import tomli as tomllib

        try:
            with open(pyproject, "rb") as f:
                data = tomllib.load(f)
            scripts = data.get("project", {}).get("scripts", {})
            entry_points.update(scripts.keys())
            logger.info(f"Loaded {len(entry_points)} CLI entry points: {entry_points}")
        except Exception as e:
            logger.warning(f"Failed to parse pyproject.toml: {e}")
        return entry_points

    def validate_examples(self, examples: List[CodeExample]) -> List[Violation]:
        """Validate all code examples and return violations."""
        for ex in examples:
            try:
                if ex.language == "python":
                    self._validate_python_example(ex)
                elif ex.language == "bash":
                    self._validate_bash_example(ex)
                # else: skip unknown languages
            except Exception as e:
                logger.debug(f"Validation error for {ex.file_path}:{ex.line_number}: {e}")

        return self.violations

    def _validate_python_example(self, example: CodeExample) -> None:
        """Validate a Python code example."""
        code = example.code

        try:
            tree = ast.parse(code)
        except SyntaxError as e:
            self.violations.append(
                Violation(
                    doc_file=example.file_path,
                    line_number=example.line_number,
                    code_example=code[:100],
                    issue_type="SyntaxError",
                    details=f"Invalid Python syntax: {e.msg} at line {e.lineno}",
                    severity="error",
                )
            )
            return

        # Walk the AST and check imports, attribute accesses, function calls
        for node in ast.walk(tree):
            # Check imports
            if isinstance(node, ast.Import):
                for alias in node.names:
                    module_name = alias.name
                    # Check if module exists (try to resolve)
                    if not self._module_exists(module_name):
                        self.violations.append(
                            Violation(
                                doc_file=example.file_path,
                                line_number=example.line_number,
                                code_example=f"import {module_name}",
                                issue_type="ModuleNotFoundError",
                                details=f"Module '{module_name}' not found in project or standard library",
                                severity="error",
                            )
                        )

            elif isinstance(node, ast.ImportFrom):
                module_name = node.module or ""
                for alias in node.names:
                    imported_name = alias.name

                    # Check if this specific import exists
                    if not self._import_exists(module_name, imported_name):
                        self.violations.append(
                            Violation(
                                doc_file=example.file_path,
                                line_number=example.line_number,
                                code_example=f"from {module_name} import {imported_name}",
                                issue_type="ImportError",
                                details=f"Cannot import '{imported_name}' from module '{module_name}'",
                                severity="error",
                            )
                        )

            # Check attribute access (e.g., module.function() or obj.method())
            elif isinstance(node, ast.Attribute):
                # For attribute chains like a.b.c, we get the Attribute nodes in reverse order
                # We need to reconstruct the full chain
                pass  # Complex - handle in a separate pass

        # Additional check: Try to extract all dotted names and verify
        self._check_dotted_access(code, example)

    def _module_exists(self, module_name: str) -> bool:
        """Check if a module exists in the project or standard library."""
        if self._should_skip_optional_import(module_name):
            return True

        # Check project modules
        if self.symbol_index.module_exists(f"metainformant.{module_name}"):
            return True
        if module_name in self.symbol_index.modules:
            return True
        # Check standard library/common third-party
        try:
            spec = importlib.util.find_spec(module_name)
            return spec is not None
        except (ImportError, ValueError):
            return False

    def _import_exists(self, module_name: str, item_name: str) -> bool:
        """Check if an imported item exists."""
        if not module_name:
            return True

        full_module = f"metainformant.{module_name}" if not module_name.startswith("metainformant") else module_name

        if self._should_skip_optional_import(module_name):
            return True

        # If full_module exists in modules, check if item is a submodule or symbol
        if full_module in self.symbol_index.modules:
            # First: check if the item is a submodule (package or module)
            full_item_module = f"{full_module}.{item_name}"
            if full_item_module in self.symbol_index.modules:
                return True
            if item_name in self.symbol_index.get_module_members(full_module):
                return True
            # Second: check if item is a symbol (class/function) defined in that module
            for sym_full_name in self.symbol_index.symbols:
                if sym_full_name.startswith(full_module + "."):
                    if sym_full_name.endswith("." + item_name):
                        return True
            try:
                module = importlib.import_module(full_module)
                return hasattr(module, item_name)
            except ImportError:
                return False

        # Standard library / third-party check via importlib
        try:
            module = importlib.import_module(module_name)
            try:
                if importlib.util.find_spec(f"{module_name}.{item_name}") is not None:
                    return True
            except (ImportError, ModuleNotFoundError, ValueError):
                pass
            return hasattr(module, item_name)
        except ImportError:
            return False

    def _should_skip_optional_import(self, module_name: str) -> bool:
        """Treat known optional dependencies as non-fatal unless strict mode is enabled."""
        root_name = module_name.split(".", 1)[0]
        return root_name in KNOWN_OPTIONAL_IMPORT_ROOTS and not self.strict_optional_imports

    def _check_dotted_access(self, code: str, example: CodeExample) -> None:
        """Check dotted attribute access patterns (e.g., metainformant.core.io.read_file)."""
        # Extract all dotted names from code
        # Pattern: word(.word)* where not inside string
        dotted_pattern = r"\b([a-zA-Z_][a-zA-Z0-9_]*(\.[a-zA-Z_][a-zA-Z0-9_]*)+)\b"

        for match in re.finditer(dotted_pattern, code):
            dotted_name = match.group(1)
            context_start = max(0, match.start() - 20)
            context = code[context_start : match.end() + 20]

            # Skip if it looks like an attribute access after parentheses
            if "(" in context and context.index("(") < context.index(dotted_name.split(".")[0]):
                continue

            # Check if the chain is valid
            self._validate_attribute_chain(dotted_name, code, example)

    def _validate_attribute_chain(self, chain: str, full_code: str, example: CodeExample) -> None:
        """Validate an attribute access chain."""
        parts = chain.split(".")
        if len(parts) < 2:
            return

        # Try to resolve from known symbols
        # First part might be an import alias or a known module/symbol
        resolved = self._resolve_name(parts[0])

        if resolved is None:
            # Check if it's a standard library module
            try:
                importlib.import_module(parts[0])
                resolved = parts[0]
            except ImportError:
                pass

        if resolved is None:
            return  # Can't resolve, skip further checks

        # Walk the chain
        current: str | Symbol | ModuleType = resolved
        for i, part in enumerate(parts[1:]):
            try:
                if isinstance(current, str):
                    # Try to import; bind the name first so the failure
                    # handler can still use it for symbol lookup.
                    name = current
                    try:
                        current = importlib.import_module(name)
                    except ImportError:
                        # Maybe it's a class/function we know
                        symbol = self.symbol_index.get_symbol(name)
                        if symbol:
                            current = symbol
                        else:
                            # Unknown, give up
                            if i == 0:  # First attribute access failed
                                self.violations.append(
                                    Violation(
                                        doc_file=example.file_path,
                                        line_number=example.line_number,
                                        code_example=chain,
                                        issue_type="AttributeError",
                                        details=f"'{parts[0]}' has no attribute '{parts[1]}' (module not found)",
                                        severity="error",
                                    )
                                )
                            return
                else:
                    # It's a Symbol or module object
                    if hasattr(current, part):
                        current = getattr(current, part)
                    else:
                        # Check if it's a method/class we should know about
                        if isinstance(current, Symbol):
                            # Check sub-symbols
                            full_name = f"{current.full_name}.{part}"
                            sub = self.symbol_index.get_symbol(full_name)
                            if sub:
                                current = sub
                            else:
                                # Not found - report error
                                self.violations.append(
                                    Violation(
                                        doc_file=example.file_path,
                                        line_number=example.line_number,
                                        code_example=chain,
                                        issue_type="AttributeError",
                                        details=f"'{parts[i]}' object has no attribute '{part}'",
                                        severity="error",
                                    )
                                )
                                return
                        else:
                            return  # Can't verify further
            except Exception:
                # Could not resolve - possibly dynamic or external dependency
                break

    def _resolve_name(self, name: str) -> str | Symbol | None:
        """Resolve a simple name to a module or symbol."""
        # Check if it's a project module
        full_name = f"metainformant.{name}"
        if self.symbol_index.module_exists(name):
            return name
        if self.symbol_index.module_exists(full_name):
            return full_name
        # Check symbols
        symbol = self.symbol_index.get_symbol(name)
        if symbol:
            return symbol
        return None

    def _validate_bash_example(self, example: CodeExample) -> None:
        """Validate a bash/shell command example."""
        code = example.code.strip()
        lines = code.split("\n")

        for line in lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Get the command (first word, handling $ and other prefixes)
            cmd_match = re.match(r"^\$?\s*([a-zA-Z0-9_-]+)", line)
            if cmd_match:
                cmd = cmd_match.group(1)

                # Check if it's a known entry point
                if cmd == "metainformant":
                    if "metainformant" not in self.entry_points:
                        self.violations.append(
                            Violation(
                                doc_file=example.file_path,
                                line_number=example.line_number,
                                code_example=line[:80],
                                issue_type="NonExistentCommand",
                                details="'metainformant' command not found in entry points",
                                severity="error",
                            )
                        )
                # Check other common CLI tools (uv, pytest, etc.)
                elif cmd in ["uv", "pytest", "git", "python", "python3", "pip"]:
                    # These are external, assume they exist. Optionally verify with which
                    pass


class ReportGenerator:
    """Generates a comprehensive report of violations."""

    def __init__(self, violations: List[Violation], repo_root: Path):
        self.violations = violations
        self.repo_root = repo_root

    def _format_doc_file(self, doc_file: str) -> str:
        """Return a stable repo-relative path when possible."""
        path = Path(doc_file)
        try:
            return str(path.resolve().relative_to(self.repo_root))
        except ValueError:
            return str(path)

    @staticmethod
    def _table_text(value: str, limit: Optional[int] = None) -> str:
        """Make text safe for a Markdown table cell."""
        text = value.replace("\n", " ").replace("|", "\\|")
        return text[:limit] if limit is not None else text

    def generate_report(self, output_path: Path) -> None:
        """Write report to file."""
        # Group by file
        by_file = defaultdict(list)
        for v in self.violations:
            by_file[v.doc_file].append(v)

        total = len(self.violations)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write("# Cross-Code Verification Report\n\n")
            f.write(f"**Total violations found:** {total}\n\n")

            if total == 0:
                f.write("✓ All code examples validated successfully!\n\n")
                return

            f.write("## Summary by Issue Type\n\n")
            by_type: defaultdict[str, int] = defaultdict(int)
            for v in self.violations:
                by_type[v.issue_type] += 1

            for issue_type, count in sorted(by_type.items()):
                f.write(f"- **{issue_type}**: {count} occurrence(s)\n")

            f.write("\n")

            f.write("## Detailed Violations\n\n")
            f.write("| File | Line | Type | Issue | Code Example |\n")
            f.write("|------|------|------|-------|--------------|\n")

            for v in sorted(self.violations, key=lambda x: (x.doc_file, x.line_number)):
                file_rel = self._format_doc_file(v.doc_file)
                details = self._table_text(v.details)
                code_preview = self._table_text(v.code_example, limit=60)
                f.write(f"| {file_rel} | {v.line_number} | {v.issue_type} | {details} | `{code_preview}` |\n")

        logger.info(f"Report written to {output_path}")

    def print_summary(self) -> None:
        """Print a concise summary."""
        print(f"\n{'='*70}")
        print("CROSS-CODE VERIFICATION SUMMARY")
        print(f"{'='*70}")
        print(f"Total issues found: {len(self.violations)}")

        if not self.violations:
            print("✓ All code examples are valid!")
            return

        by_type: defaultdict[str, int] = defaultdict(int)
        for v in self.violations:
            by_type[v.issue_type] += 1

        print("\nIssues by category:")
        for issue_type, count in sorted(by_type.items()):
            print(f"  {issue_type}: {count}")

        # Show top 10 files with most issues
        by_file: defaultdict[str, int] = defaultdict(int)
        for v in self.violations:
            by_file[v.doc_file] += 1

        print("\nTop files with issues:")
        for file, count in sorted(by_file.items(), key=lambda x: -x[1])[:10]:
            print(f"  {file}: {count}")


# --- Runtime import truth -------------------------------------------------
#
# The example validator above reasons from a static symbol index and treats many
# fences as illustrative snippets. Every quoted import is still a claim about the
# API, so this pass resolves each complete ``metainformant`` import line exactly as
# Python would (attribute first, then submodule) without executing document text.

IMPORT_STATEMENT = re.compile(
    r"^\s*((?:from\s+metainformant(?:\.\w+)*\s+import\s+[^#\n]+?)"
    r"|(?:import\s+metainformant(?:\.\w+)*(?:\s+as\s+\w+)?))\s*(?:#.*)?$"
)
IMPORT_SCAN_SKIP_PARTS = frozenset(
    {
        ".cursor",
        ".git",
        ".tmp",
        ".venv",
        "__pycache__",
        "archive",
        "historical_reports",
        "node_modules",
        "output",
        "projects",
    }
)


def _is_historical_snapshot(lines: List[str]) -> bool:
    return any("Historical snapshot" in line for line in lines[:12])


def _active_markdown(paths: List[Path], include_historical: bool = False) -> Iterator[tuple[Path, List[str]]]:
    """Yield ``(path, lines)`` for active Markdown under ``paths`` (directories or files).

    Archived and generated trees (``IMPORT_SCAN_SKIP_PARTS``) are skipped, and so
    are historical-snapshot documents unless ``include_historical`` is set.
    """
    for root in paths:
        files = [root] if root.is_file() else sorted(root.rglob("*.md"))
        for md in files:
            parts = md.relative_to(root).parts if md != root else (md.name,)
            if any(part in IMPORT_SCAN_SKIP_PARTS for part in parts):
                continue
            try:
                lines = md.read_text(encoding="utf-8").splitlines()
            except (OSError, UnicodeDecodeError):
                continue
            if not include_historical and _is_historical_snapshot(lines):
                continue
            yield md, lines


def collect_doc_import_statements(paths: List[Path], include_historical: bool = False) -> Dict[str, List[str]]:
    """Map each complete ``metainformant`` import quoted in Markdown to its locations.

    ``paths`` may mix directories (scanned recursively) and individual files.
    Single-line imports are read from any line; parenthesized multi-line imports
    are read from Python-fenced blocks that parse. Placeholders such as
    ``from metainformant.dna import ...`` are skipped, and so are archived and
    historical-snapshot documents.
    """
    statements: Dict[str, List[str]] = defaultdict(list)
    for md, lines in _active_markdown(paths, include_historical):
        for first_line, code in _python_fenced_blocks(lines):
            try:
                tree = ast.parse(code)
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, (ast.Import, ast.ImportFrom)) and node.end_lineno != node.lineno:
                    names = [alias.name for alias in node.names]
                    module = node.module if isinstance(node, ast.ImportFrom) else names[0]
                    if module and module.split(".")[0] == "metainformant":
                        statements[ast.unparse(node)].append(f"{md}:{first_line + node.lineno - 1}")
        for number, line in enumerate(lines, 1):
            match = IMPORT_STATEMENT.match(line)
            if not match:
                continue
            statement = match.group(1).strip()
            if statement.endswith(("(", "\\", ",", "import ...")):
                continue
            statements[statement].append(f"{md}:{number}")
    return dict(statements)


def resolve_import_statement(statement: str) -> Optional[str]:
    """Return why ``statement`` cannot be satisfied by the installed package, or None.

    Only a single ``import``/``from ... import`` of a ``metainformant`` module is
    accepted. A missing third-party optional dependency is not a documentation
    error and resolves to None.
    """
    try:
        tree = ast.parse(statement)
    except SyntaxError as exc:
        return f"SyntaxError: {exc.msg}"
    node = tree.body[0] if len(tree.body) == 1 else None
    try:
        if isinstance(node, ast.Import):
            for alias in node.names:
                importlib.import_module(alias.name)
            return None
        if not isinstance(node, ast.ImportFrom) or not node.module or node.level:
            return "not a single absolute import statement"
        module_name = node.module
        module = importlib.import_module(module_name)
        for alias in node.names:
            if alias.name == "*" or hasattr(module, alias.name):
                continue
            try:
                importlib.import_module(f"{module_name}.{alias.name}")
            except ModuleNotFoundError as exc:
                if exc.name == f"{module_name}.{alias.name}":
                    return f"ImportError: cannot import name {alias.name!r} from {module_name!r}"
                raise
    except ModuleNotFoundError as exc:
        if exc.name and exc.name.split(".")[0] != "metainformant":
            return None
        return f"ModuleNotFoundError: {exc}"
    except Exception as exc:  # an import-time failure inside a real module is still a broken example
        return f"{type(exc).__name__}: {exc}"
    return None


def find_stale_doc_imports(paths: List[Path], include_historical: bool = False) -> List[Violation]:
    """Return a violation for every quoted import the current code cannot satisfy."""
    violations = []
    for statement, locations in sorted(collect_doc_import_statements(paths, include_historical).items()):
        reason = resolve_import_statement(statement)
        if reason is None:
            continue
        for location in locations:
            doc_file, _, line = location.rpartition(":")
            violations.append(Violation(doc_file, int(line), statement, "StaleImport", reason))
    return violations


PYTHON_FENCES = frozenset({"python", "python-snippet", "py", "python3"})
FENCE_LINE = re.compile(r"^\s*```\s*([\w+-]*)\s*$")


def _python_fenced_blocks(lines: List[str]) -> List[tuple[int, str]]:
    """Return ``(first_line_number, dedented_code)`` for Python-fenced blocks."""
    blocks: List[tuple[int, str]] = []
    start: Optional[int] = None
    body: List[str] = []
    for number, line in enumerate(lines, 1):
        match = FENCE_LINE.match(line)
        if start is None:
            if match and match.group(1).lower() in PYTHON_FENCES:
                start, body = number + 1, []
        elif match and not match.group(1):
            blocks.append((start, textwrap.dedent("\n".join(body))))
            start = None
        else:
            body.append(line)
    return blocks


def _bound_import_targets(tree: ast.Module) -> Dict[str, object]:
    """Map names bound by a block's resolvable ``metainformant`` imports to the imported objects."""
    bound: Dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module and not node.level:
            if node.module.split(".")[0] != "metainformant":
                continue
            try:
                module = importlib.import_module(node.module)
            except Exception:
                continue
            for alias in node.names:
                if alias.name == "*":
                    continue
                target = getattr(module, alias.name, None)
                if target is None:
                    try:
                        target = importlib.import_module(f"{node.module}.{alias.name}")
                    except Exception:
                        continue
                bound[alias.asname or alias.name] = target
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] != "metainformant" or not alias.asname:
                    continue
                try:
                    bound[alias.asname] = importlib.import_module(alias.name)
                except Exception:
                    continue
    return bound


def _rebound_names(tree: ast.Module) -> Set[str]:
    """Names a block assigns, deletes, or binds as parameters/loop targets after importing them."""
    names: Set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name) and isinstance(node.ctx, (ast.Store, ast.Del)):
            names.add(node.id)
        elif isinstance(node, ast.arg):
            names.add(node.arg)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(node.name)
    return names


def _attribute_chain(node: ast.Attribute) -> Optional[tuple[str, List[str]]]:
    attrs: List[str] = []
    current: ast.AST = node
    while isinstance(current, ast.Attribute):
        attrs.append(current.attr)
        current = current.value
    return (current.id, attrs[::-1]) if isinstance(current, ast.Name) else None


def _missing_attribute(target: object, attrs: List[str]) -> Optional[tuple[int, str]]:
    """Walk ``attrs`` through modules and classes; return ``(index, reason)`` for the first missing step."""
    current = target
    for index, attr in enumerate(attrs):
        if isinstance(current, ModuleType):
            if hasattr(current, attr):
                current = getattr(current, attr)
                continue
            submodule = f"{current.__name__}.{attr}"
            try:
                if importlib.util.find_spec(submodule) is not None:
                    current = importlib.import_module(submodule)
                    continue
            except Exception:
                pass
            return index, f"module {current.__name__!r} has no attribute {attr!r}"
        if isinstance(current, type):
            if not hasattr(current, attr):
                return index, f"class {current.__qualname__!r} has no attribute {attr!r}"
            current = getattr(current, attr)
            continue
        return None  # instances, functions and values may carry dynamic attributes
    return None


def find_stale_doc_attributes(paths: List[Path], include_historical: bool = False) -> List[Violation]:
    """Report attribute chains on imported modules/classes that the current code lacks.

    Only Python-fenced blocks that parse are checked, and only names bound by the
    block's own ``metainformant`` imports and never reassigned in that block.
    """
    violations: List[Violation] = []
    for md, lines in _active_markdown(paths, include_historical):
        for first_line, code in _python_fenced_blocks(lines):
            try:
                tree = ast.parse(code)
            except SyntaxError:
                continue
            bound = _bound_import_targets(tree)
            if not bound:
                continue
            rebound = _rebound_names(tree)
            reported: Set[str] = set()
            for node in ast.walk(tree):
                chain = _attribute_chain(node) if isinstance(node, ast.Attribute) else None
                if not chain or chain[0] not in bound or chain[0] in rebound:
                    continue
                name, attrs = chain
                missing = _missing_attribute(bound[name], attrs)
                if missing is None:
                    continue
                index, reason = missing
                failing = ".".join([name, *attrs[: index + 1]])
                if failing in reported:  # report each missing step once per block
                    continue
                reported.add(failing)
                line = first_line + getattr(node, "lineno", 1) - 1
                violations.append(Violation(str(md), line, failing, "StaleAttribute", reason))
    return violations


# --- Call signatures -------------------------------------------------------
#
# A quoted call is a claim about the callee's parameters. Calls on names a block
# imports from ``metainformant``, directly or through module and class attribute
# chains, are bound against ``inspect.signature`` with placeholder arguments, so a
# renamed keyword or a removed parameter is reported without running the example.


def _resolve_attribute_chain(target: object, attrs: List[str]) -> Optional[object]:
    """Follow ``attrs`` through modules and classes; None when a step is missing or dynamic."""
    current = target
    for attr in attrs:
        if isinstance(current, ModuleType):
            if hasattr(current, attr):
                current = getattr(current, attr)
                continue
            try:
                current = importlib.import_module(f"{current.__name__}.{attr}")
            except Exception:
                return None
        elif isinstance(current, type) and hasattr(current, attr):
            current = getattr(current, attr)
        else:
            return None
    return current


def call_signature_mismatch(callee: object, call: ast.Call) -> Optional[str]:
    """Return why ``call``'s arguments cannot bind to ``callee``'s signature, or None.

    Calls that unpack ``*args`` or ``**kwargs``, and callables without an
    introspectable signature (some builtins and extension types), are not judged.
    """
    if any(isinstance(arg, ast.Starred) for arg in call.args) or any(kw.arg is None for kw in call.keywords):
        return None
    if not callable(callee):
        return f"{type(callee).__name__!r} object is not callable"
    try:
        signature = inspect.signature(callee)
    except (TypeError, ValueError):
        return None
    placeholder = object()
    keywords = {kw.arg: placeholder for kw in call.keywords if kw.arg is not None}
    try:
        signature.bind(*[placeholder] * len(call.args), **keywords)
    except TypeError as exc:
        return str(exc)
    return None


def find_stale_doc_calls(paths: List[Path], include_historical: bool = False) -> List[Violation]:
    """Report quoted calls whose arguments no longer fit the imported callee's signature."""
    violations: List[Violation] = []
    for md, lines in _active_markdown(paths, include_historical):
        for first_line, code in _python_fenced_blocks(lines):
            try:
                tree = ast.parse(code)
            except SyntaxError:
                continue
            bound = _bound_import_targets(tree)
            if not bound:
                continue
            rebound = _rebound_names(tree)
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                chain: Optional[tuple[str, List[str]]] = None
                if isinstance(node.func, ast.Name):
                    chain = (node.func.id, [])
                elif isinstance(node.func, ast.Attribute):
                    chain = _attribute_chain(node.func)
                if not chain or chain[0] not in bound or chain[0] in rebound:
                    continue
                name, attrs = chain
                callee = _resolve_attribute_chain(bound[name], attrs)
                if callee is None:
                    continue
                reason = call_signature_mismatch(callee, node)
                if reason is None:
                    continue
                example = ".".join([name, *attrs]) + "(...)"
                violations.append(Violation(str(md), first_line + node.lineno - 1, example, "StaleCall", reason))
    return violations


# --- Relative links --------------------------------------------------------
#
# A relative link or image is a claim that its target exists in the checkout.
# Paths are matched one component at a time with exact case, so a link that only
# resolves on a case-insensitive filesystem is still reported, and an anchor into
# a Markdown file must name one of its headings or an explicit HTML id.

INLINE_LINK = re.compile(r"\[(?:[^\[\]]|\[[^\[\]]*\])*\]\(\s*(<[^<>\n]*>|[^\s()<>]+(?:\([^\s()]*\)[^\s()<>]*)*)")
REFERENCE_DEFINITION = re.compile(r"^\s{0,3}\[(?!\^)[^\]]+\]:\s*(<[^<>\n]*>|\S+)")
HTML_LINK = re.compile(r"""<(?:a|img)\b[^>]*?\b(?:href|src)\s*=\s*["']([^"']+)["']""", re.IGNORECASE)
INLINE_CODE = re.compile(r"(`+)(?:(?!\1).)+?\1")
URL_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
ATX_HEADING = re.compile(r"^\s{0,3}#{1,6}\s+(.*?)\s*#*\s*$")
EXPLICIT_ANCHOR = re.compile(r"""\b(?:name|id)\s*=\s*["']([^"']+)["']""", re.IGNORECASE)
FENCE_OPEN = re.compile(r"^\s{0,3}(`{3,}|~{3,})")


def _lines_outside_fences(lines: List[str]) -> Iterator[tuple[int, str]]:
    """Yield ``(line_number, line)`` for lines outside fenced code blocks."""
    fence: Optional[str] = None
    for number, line in enumerate(lines, 1):
        if fence is None:
            match = FENCE_OPEN.match(line)
            if match:
                fence = match.group(1)
            else:
                yield number, line
        elif re.fullmatch(rf"\s{{0,3}}{re.escape(fence[0])}{{{len(fence)},}}\s*", line):
            fence = None


def _heading_slug(text: str) -> str:
    """GitHub's anchor for a heading: rendered text, lowercased, punctuation dropped, spaces to hyphens."""
    text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"<[^>]+>", "", text).replace("`", "").strip().lower()
    return re.sub(r"[^\w\- ]", "", text).replace(" ", "-")


def markdown_anchors(lines: List[str]) -> Set[str]:
    """Return the anchors a Markdown document defines (heading slugs and explicit ids)."""
    anchors: Set[str] = set()
    seen: Dict[str, int] = defaultdict(int)
    for _, line in _lines_outside_fences(lines):
        heading = ATX_HEADING.match(line)
        if heading:
            slug = _heading_slug(heading.group(1))
            anchors.add(f"{slug}-{seen[slug]}" if seen[slug] else slug)
            seen[slug] += 1
            # Underscores used as emphasis markers vanish from the rendered heading.
            anchors.add(_heading_slug(re.sub(r"(?<!\w)_+|_+(?!\w)", "", heading.group(1))))
        anchors.update(match.group(1).lower() for match in EXPLICIT_ANCHOR.finditer(line))
    return anchors


def _submodule_paths(repo_root: Path) -> List[Path]:
    try:
        text = (repo_root / ".gitmodules").read_text(encoding="utf-8")
    except OSError:
        return []
    return [repo_root / m.group(1).strip() for m in re.finditer(r"^\s*path\s*=\s*(.+)$", text, re.MULTILINE)]


def find_broken_doc_links(
    paths: List[Path], include_historical: bool = False, repo_root: Path = REPO_ROOT
) -> List[Violation]:
    """Report relative links and images whose file, directory or Markdown anchor is missing.

    External URLs and templated targets are not checked. Targets inside a
    registered submodule that is not checked out cannot be verified and are skipped.
    """
    repo_root = repo_root.resolve()
    submodules = _submodule_paths(repo_root)
    entries: Dict[Path, Set[str]] = {}
    anchors: Dict[Path, Set[str]] = {}

    def listing(directory: Path) -> Set[str]:
        if directory not in entries:
            try:
                entries[directory] = set(os.listdir(directory))
            except OSError:
                entries[directory] = set()
        return entries[directory]

    def anchors_of(document: Path) -> Set[str]:
        if document not in anchors:
            try:
                anchors[document] = markdown_anchors(document.read_text(encoding="utf-8").splitlines())
            except (OSError, UnicodeDecodeError):
                anchors[document] = set()
        return anchors[document]

    def problem(md: Path, raw: str) -> Optional[str]:
        target = raw[1:-1].strip() if raw.startswith("<") and raw.endswith(">") else raw.strip()
        if not target or target.startswith("//") or URL_SCHEME.match(target) or any(c in target for c in "{}$"):
            return None
        path_part, _, anchor = target.partition("#")
        path_part = unquote(path_part.split("?", 1)[0])
        destination = md
        if path_part:
            base = repo_root if path_part.startswith("/") else md.parent
            destination = Path(os.path.normpath(base / path_part.lstrip("/")))
            try:
                parts = destination.relative_to(repo_root).parts
            except ValueError:
                return "link target is outside the repository"
            for submodule in submodules:
                if (destination == submodule or submodule in destination.parents) and not listing(submodule):
                    return None
            current = repo_root
            for part in parts:
                if part not in listing(current):
                    near = [name for name in listing(current) if name.lower() == part.lower()]
                    hint = f" (on disk: {near[0]!r})" if near else ""
                    return f"{current / part} does not exist{hint}"
                current = current / part
        if anchor and destination.suffix.lower() == ".md" and destination.is_file():
            if unquote(anchor).lower() not in anchors_of(destination):
                return f"no heading or id {anchor!r} in {destination.name}"
        return None

    violations: List[Violation] = []
    for md, lines in _active_markdown(paths, include_historical):
        md = Path(os.path.normpath(md.absolute()))
        for number, line in _lines_outside_fences(lines):
            prose = INLINE_CODE.sub(lambda m: " " * len(m.group(0)), line)
            targets = [
                m.group(1)
                for pattern in (INLINE_LINK, REFERENCE_DEFINITION, HTML_LINK)
                for m in pattern.finditer(prose)
            ]
            for raw in targets:
                reason = problem(md, raw)
                if reason:
                    violations.append(Violation(str(md), number, raw, "BrokenLink", reason))
    return violations


def run(args: argparse.Namespace) -> int:
    """Execute the verification pipeline and return the process exit code."""

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    repo_root = REPO_ROOT
    docs_dir = resolve_repo_path(args.docs_dir, repo_root)
    src_dir = resolve_repo_path(args.src_dir, repo_root)
    output_path = resolve_repo_path(args.output, repo_root)
    pyproject = repo_root / "pyproject.toml"

    if not docs_dir.exists():
        logger.error(f"Docs directory not found: {docs_dir}")
        sys.exit(1)
    if not src_dir.exists():
        logger.error(f"Source directory not found: {src_dir}")
        sys.exit(1)

    logger.info(f"Repository root: {repo_root}")
    logger.info(f"Docs directory: {docs_dir}")
    logger.info(f"Source directory: {src_dir}")
    logger.info(f"Report output: {output_path}")

    # Step 1: Build Python symbol index
    index = PythonSymbolIndex(src_dir)
    index.build_index()

    # Step 2: Parse all documentation files
    parser = DocumentationParser(docs_dir, include_historical=args.include_historical)
    doc_files = parser.find_all_docs()

    all_examples = collect_code_examples(parser, doc_files, args.verbose)

    logger.info(f"Extracted {len(all_examples)} total code examples from {len(doc_files)} files")

    # Step 3: Validate all examples
    validator = CodeValidator(index, pyproject, strict_optional_imports=args.strict_optional_imports)
    violations = validator.validate_examples(all_examples)
    violations.extend(find_stale_doc_imports([docs_dir], include_historical=args.include_historical))
    violations.extend(find_stale_doc_attributes([docs_dir], include_historical=args.include_historical))
    violations.extend(find_stale_doc_calls([docs_dir], include_historical=args.include_historical))
    violations.extend(
        find_broken_doc_links([docs_dir], include_historical=args.include_historical, repo_root=repo_root)
    )

    # Step 4: Generate report
    report_gen = ReportGenerator(violations, repo_root)
    report_gen.generate_report(output_path)
    report_gen.print_summary()

    # Exit code: 0 if no violations, 1 if violations found
    return 0 if not violations else 1
