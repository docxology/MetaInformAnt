#!/usr/bin/env python3
"""
Systematic Cross-Code Verification for METAINFORMANT Documentation

Verifies that all code examples, imports, function calls, and CLI commands
in documentation files match the actual source code.

Usage:
    python scripts/verify_documentation_code.py [--docs-dir DIR] [--src-dir DIR] [--output REPORT.md]

Target: Check all 506+ documentation files for broken code references.
"""

import argparse
import ast
import importlib.util
import logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Add src directory to sys.path to enable imports of metainformant package
REPO_ROOT = Path(__file__).parent.parent
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


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

    def __str__(self):
        return f"{self.doc_file}:{self.line_number} - {self.issue_type}: {self.details}"


class PythonSymbolIndex:
    """Builds and maintains an index of all Python symbols in the source code."""

    def __init__(self, src_dir: Path):
        self.src_dir = src_dir
        self.symbols: Dict[str, Symbol] = {}  # full_name -> Symbol
        self.modules: Dict[str, Path] = {}  # module_name -> file_path
        self.class_methods: Dict[str, Set[str]] = defaultdict(set)  # class_name -> {method_names}
        self.imports_cache: Dict[str, Set[str]] = defaultdict(set)  # module -> {imported_names}

    def build_index(self):
        """Scan all Python files and build symbol index."""
        logger.info(f"Building Python symbol index from {self.src_dir}...")

        python_files = list(self.src_dir.rglob("*.py"))
        logger.info(f"Found {len(python_files)} Python files")

        for py_file in python_files:
            self._parse_file(py_file)

        logger.info(f"Indexed {len(self.symbols)} symbols across {len(self.modules)} modules")
        return self

    def _parse_file(self, file_path: Path):
        """Parse a single Python file and extract symbols."""
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                content = f.read()
            tree = ast.parse(content, filename=str(file_path))

            # Determine module path relative to src_dir
            rel_path = file_path.relative_to(self.src_dir)
            module_parts = list(rel_path.parts)
            if module_parts[-1] == "__init__.py":
                module_parts = module_parts[:-1]
            elif module_parts[-1].endswith(".py"):
                module_parts[-1] = module_parts[-1][:-3]

            module_name = ".".join(["metainformant"] + module_parts)

            # Register module
            self.modules[module_name] = file_path

            # Extract imports (for reference)
            self._extract_imports(tree, module_name)

            # Walk AST and collect symbols
            for node in ast.walk(tree):
                if isinstance(node, ast.ClassDef):
                    self._add_symbol(node, "class", module_name, file_path)
                    # Track methods for this class
                    for item in node.body:
                        if isinstance(item, ast.FunctionDef):
                            self.class_methods[node.name].add(item.name)

                elif isinstance(node, ast.FunctionDef):
                    # Top-level function
                    if isinstance(node.parent, ast.Module) if hasattr(node, "parent") else True:
                        self._add_symbol(node, "function", module_name, file_path)

        except Exception as e:
            logger.warning(f"Failed to parse {file_path}: {e}")

    def _add_symbol(self, node, kind: str, module_name: str, file_path: Path):
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

    def _get_signature(self, node) -> str:
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

    def _extract_imports(self, tree, module_name: str):
        """Extract imports from a module."""
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    imported.add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    for alias in node.names:
                        imported.add(f"{node.module}.{alias.name}" if alias.name != "*" else f"{node.module}.*")
        self.imports_cache[module_name] = imported

    def get_symbol(self, full_name: str) -> Optional[Symbol]:
        """Look up a symbol by its full name."""
        return self.symbols.get(full_name)

    def get_module_members(self, module_name: str) -> Set[str]:
        """Get all public members exported by a module."""
        members = set()
        for sym in self.symbols.values():
            if sym.module_path.endswith(module_name.replace(".", "/")):
                # Check if it's 'module.submodule.symbol' - complex
                pass
        # Also check __all__ if defined
        return members

    def module_exists(self, module_name: str) -> bool:
        """Check if a module exists."""
        return module_name in self.modules


class DocumentationParser:
    """Parses documentation files and extracts code examples."""

    # Regex patterns for inline code
    INLINE_CODE_PATTERN = r"`([^`]+)`"

    # Code block patterns: ```python, ```bash, ```, etc.
    CODE_BLOCK_PATTERN = r"```(\w*)\n(.*?)```"

    def __init__(self, docs_dir: Path):
        self.docs_dir = docs_dir

    def find_all_docs(self) -> List[Path]:
        """Find all markdown documentation files."""
        docs = list(self.docs_dir.rglob("*.md"))
        # Exclude certain directories?
        logger.info(f"Found {len(docs)} documentation files")
        return docs

    def extract_code_examples(self, file_path: Path) -> List[CodeExample]:
        """Extract all code examples from a markdown file."""
        examples = []

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()

            content = "".join(lines)

            # Find all code blocks
            for match in re.finditer(r"```(\w*)\n(.*?)```", content, re.DOTALL):
                lang = match.group(1).strip()
                code = match.group(2)
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
        # Skip single words, filenames, paths
        if len(code.split()) == 1 and ("." in code or "/" in code or code.endswith(".md")):
            return False
        # Skip just commands that are clearly not importable
        return True

    def _guess_language(self, code: str) -> str:
        """Guess the language of a code snippet."""
        # Python indicators
        if any(kw in code for kw in ["import ", "def ", "class ", "from ", "print("]):
            return "python"
        # Shell/bash indicators
        if any(kw in code for kw in ["$", "bash", "sh ", "&&", "||", ";"]):
            return "bash"
        return "unknown"


class CodeValidator:
    """Validates code examples against the source code index."""

    def __init__(self, symbol_index: PythonSymbolIndex, pyproject_toml: Path):
        self.symbol_index = symbol_index
        self.violations: List[Violation] = []
        self.entry_points = self._load_entry_points(pyproject_toml)

    def _load_entry_points(self, pyproject: Path) -> Set[str]:
        """Load CLI entry points from pyproject.toml."""
        entry_points = set()
        try:
            import tomllib
        except ImportError:
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

    def _validate_python_example(self, example: CodeExample):
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

        # If full_module exists in modules, check if item is a submodule or symbol
        if full_module in self.symbol_index.modules:
            # First: check if the item is a submodule (package or module)
            full_item_module = f"{full_module}.{item_name}"
            if full_item_module in self.symbol_index.modules:
                return True
            # Second: check if item is a symbol (class/function) defined in that module
            for sym_full_name in self.symbol_index.keys():
                if sym_full_name.startswith(full_module + "."):
                    if sym_full_name.endswith("." + item_name):
                        return True
            return False

        # Standard library / third-party check via importlib
        try:
            module = importlib.import_module(module_name)
            return hasattr(module, item_name)
        except ImportError:
            return False

    def _check_dotted_access(self, code: str, example: CodeExample):
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

    def _validate_attribute_chain(self, chain: str, full_code: str, example: CodeExample):
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
        current = resolved
        for i, part in enumerate(parts[1:]):
            try:
                if isinstance(current, str):
                    # Try to import
                    try:
                        module = importlib.import_module(current)
                        current = module
                    except ImportError:
                        # Maybe it's a class/function we know
                        symbol = self.symbol_index.get_symbol(current)
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

    def _resolve_name(self, name: str):
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

    def _validate_bash_example(self, example: CodeExample):
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

    def __init__(self, violations: List[Violation]):
        self.violations = violations

    def generate_report(self, output_path: Path):
        """Write report to file."""
        # Group by file
        by_file = defaultdict(list)
        for v in self.violations:
            by_file[v.doc_file].append(v)

        total = len(self.violations)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write("# Cross-Code Verification Report\n\n")
            f.write(f"**Total violations found:** {total}\n\n")

            if total == 0:
                f.write("✓ All code examples validated successfully!\n\n")
                return

            f.write("## Summary by Issue Type\n\n")
            by_type = defaultdict(int)
            for v in self.violations:
                by_type[v.issue_type] += 1

            for issue_type, count in sorted(by_type.items()):
                f.write(f"- **{issue_type}**: {count} occurrence(s)\n")

            f.write("\n")

            f.write("## Detailed Violations\n\n")
            f.write("| File | Line | Type | Issue | Code Example |\n")
            f.write("|------|------|------|-------|--------------|\n")

            for v in sorted(self.violations, key=lambda x: (x.doc_file, x.line_number)):
                file_rel = v.doc_file.replace("/home/trim/Documents/Git/MetaInformAnt/", "")
                code_preview = v.code_example.replace("\n", " ")[:60]
                f.write(f"| {file_rel} | {v.line_number} | {v.issue_type} | {v.details} | `{code_preview}` |\n")

        logger.info(f"Report written to {output_path}")

    def print_summary(self):
        """Print a concise summary."""
        print(f"\n{'='*70}")
        print("CROSS-CODE VERIFICATION SUMMARY")
        print(f"{'='*70}")
        print(f"Total issues found: {len(self.violations)}")

        if not self.violations:
            print("✓ All code examples are valid!")
            return

        by_type = defaultdict(int)
        for v in self.violations:
            by_type[v.issue_type] += 1

        print("\nIssues by category:")
        for issue_type, count in sorted(by_type.items()):
            print(f"  {issue_type}: {count}")

        # Show top 10 files with most issues
        by_file = defaultdict(int)
        for v in self.violations:
            by_file[v.doc_file] += 1

        print("\nTop files with issues:")
        for file, count in sorted(by_file.items(), key=lambda x: -x[1])[:10]:
            print(f"  {file}: {count}")


def main():
    parser = argparse.ArgumentParser(description="Verify code examples in documentation")
    parser.add_argument("--docs-dir", type=Path, default=Path("docs"), help="Documentation directory (default: docs/)")
    parser.add_argument("--src-dir", type=Path, default=Path("src"), help="Source code directory (default: src/)")
    parser.add_argument(
        "--output", type=Path, default=Path("cross_code_verification_report.md"), help="Output report file"
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Resolve paths
    repo_root = Path(__file__).parent.parent  # Assuming script is in scripts/
    docs_dir = (repo_root / args.docs_dir).resolve()
    src_dir = (repo_root / args.src_dir).resolve()
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

    # Step 1: Build Python symbol index
    index = PythonSymbolIndex(src_dir)
    index.build_index()

    # Step 2: Parse all documentation files
    parser = DocumentationParser(docs_dir)
    doc_files = parser.find_all_docs()

    all_examples = []
    for doc_file in doc_files:
        examples = parser.extract_code_examples(doc_file)
        all_examples.extend(examples)
        if args.verbose:
            logger.debug(f"  {doc_file.name}: {len(examples)} code examples")

    logger.info(f"Extracted {len(all_examples)} total code examples from {len(doc_files)} files")

    # Step 3: Validate all examples
    validator = CodeValidator(index, pyproject)
    violations = validator.validate_examples(all_examples)

    # Step 4: Generate report
    report_gen = ReportGenerator(violations)
    report_gen.generate_report(args.output)
    report_gen.print_summary()

    # Exit code: 0 if no violations, 1 if violations found
    sys.exit(0 if not violations else 1)


if __name__ == "__main__":
    main()
