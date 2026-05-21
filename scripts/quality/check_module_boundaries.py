#!/usr/bin/env python3
"""Check METAINFORMANT package dependency boundaries.

Domain modules may depend on ``metainformant.core``.  Cross-domain imports are
allowed only from explicit adapter surfaces such as ``integration`` and
``workflow`` modules, or from CLI entry points.
"""

from __future__ import annotations

import ast
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BoundaryViolation:
    """A cross-domain import that is not on an allowed adapter surface."""

    path: Path
    line: int
    source_domain: str
    target_domain: str
    imported: str


def repo_root_from_script() -> Path:
    """Resolve the repository root from this script location."""
    return Path(__file__).resolve().parents[2]


def package_domains(src_root: Path) -> set[str]:
    """Return top-level package domains under ``src/metainformant``."""
    return {p.name for p in src_root.iterdir() if p.is_dir() and not p.name.startswith("_")}


def is_allowed_adapter(path: Path, src_root: Path) -> bool:
    """Return True if ``path`` is an approved cross-domain adapter surface."""
    rel_parts = path.relative_to(src_root).parts
    if rel_parts[0] == "__main__.py":
        return True
    return any(part in {"integration", "workflow"} for part in rel_parts)


def imported_domain(module_name: str, domains: set[str]) -> tuple[str | None, str | None]:
    """Return the target domain and full import when importing metainformant.*."""
    if not module_name.startswith("metainformant."):
        return None, None
    parts = module_name.split(".")
    if len(parts) < 2 or parts[1] not in domains:
        return None, None
    return parts[1], module_name


def check_file(path: Path, src_root: Path, domains: set[str]) -> list[BoundaryViolation]:
    """Check one Python file for disallowed cross-domain imports."""
    rel = path.relative_to(src_root)
    source_domain = rel.parts[0]
    if source_domain.endswith(".py") or source_domain not in domains:
        return []

    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as exc:
        return [BoundaryViolation(path, exc.lineno or 0, source_domain, "syntax", str(exc))]

    violations: list[BoundaryViolation] = []
    for node in ast.walk(tree):
        imports: list[tuple[int, str]] = []
        if isinstance(node, ast.Import):
            imports.extend((node.lineno, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.append((node.lineno, node.module))

        for line, module_name in imports:
            target_domain, imported = imported_domain(module_name, domains)
            if target_domain is None:
                continue
            if target_domain in {source_domain, "core"}:
                continue
            if is_allowed_adapter(path, src_root):
                continue
            violations.append(
                BoundaryViolation(
                    path=path,
                    line=line,
                    source_domain=source_domain,
                    target_domain=target_domain,
                    imported=imported or module_name,
                )
            )
    return violations


def collect_violations(repo_root: Path) -> list[BoundaryViolation]:
    """Collect all package boundary violations."""
    src_root = repo_root / "src" / "metainformant"
    domains = package_domains(src_root)
    violations: list[BoundaryViolation] = []
    for path in sorted(src_root.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        violations.extend(check_file(path, src_root, domains))
    return violations


def main() -> int:
    """Run the boundary check."""
    repo_root = repo_root_from_script()
    violations = collect_violations(repo_root)
    if violations:
        print("Disallowed cross-domain imports:")
        for violation in violations:
            rel = violation.path.relative_to(repo_root)
            print(
                f"- {rel}:{violation.line}: {violation.source_domain} -> "
                f"{violation.target_domain} via {violation.imported}"
            )
        return 1
    print("Module boundary check passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
