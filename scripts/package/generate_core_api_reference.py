#!/usr/bin/env python3
"""Generate the canonical core API inventory from the source AST.

The generated reference intentionally contains signatures and one-line
docstring summaries rather than hand-maintained prose.  This keeps the public
core API inventory synchronized with source changes while the topic guides
remain the place for tutorials, caveats, and design discussion.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import dataclass
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src" / "metainformant" / "core"
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "core" / "API_REFERENCE.md"

# These modules are compatibility facades for the canonical package layout.
# Their re-exported names are documented through the implementation modules,
# avoiding duplicate entries in the generated reference.
EXCLUDED_SOURCE_FILES = frozenset({"__init__.py", "config.py", "db.py", "paths.py"})


@dataclass(frozen=True)
class Symbol:
    """A public source symbol rendered in the API inventory."""

    module: str
    qualified_name: str
    kind: str
    signature: str
    summary: str


def _annotation(node: ast.AST | None) -> str:
    if node is None:
        return ""
    try:
        return ast.unparse(node)
    except Exception:
        return "..."


def _default(node: ast.AST | None) -> str:
    if node is None:
        return ""
    try:
        return ast.unparse(node)
    except Exception:
        return "..."


def _argument_text(arg: ast.arg, default: ast.AST | None = None) -> str:
    text = arg.arg
    annotation = _annotation(arg.annotation)
    if annotation:
        text += f": {annotation}"
    if default is not None:
        text += f" = {_default(default)}"
    return text


def _signature(node: ast.FunctionDef | ast.AsyncFunctionDef, *, method: bool = False) -> str:
    """Render a stable, source-derived function signature."""

    args = node.args
    positional = [*args.posonlyargs, *args.args]
    if method and positional and positional[0].arg in {"self", "cls"}:
        positional = positional[1:]
    all_positional = [*args.posonlyargs, *args.args]
    defaults = [None] * (len(all_positional) - len(args.defaults)) + list(args.defaults)
    if method and all_positional and all_positional[0].arg in {"self", "cls"}:
        defaults = defaults[1:]

    parts = [
        _argument_text(argument, default)
        for argument, default in zip(positional, defaults, strict=True)
    ]
    if args.vararg is not None:
        parts.append(f"*{_argument_text(args.vararg)}")
    elif args.kwonlyargs:
        parts.append("*")
    for argument, default in zip(args.kwonlyargs, args.kw_defaults, strict=True):
        parts.append(_argument_text(argument, default))
    if args.kwarg is not None:
        parts.append(f"**{_argument_text(args.kwarg)}")

    separator = ", ".join(parts)
    if args.posonlyargs:
        # Keep the Python 3.8+ positional-only boundary visible.
        posonly_count = len(args.posonlyargs)
        rendered = parts[:posonly_count]
        remaining = parts[posonly_count:]
        separator = ", ".join([*rendered, "/", *remaining])
    result = f"({separator})"
    return_annotation = _annotation(node.returns)
    if return_annotation:
        result += f" -> {return_annotation}"
    return result


def _is_staticmethod(node: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    return any(
        isinstance(decorator, ast.Name) and decorator.id == "staticmethod"
        for decorator in node.decorator_list
    )


def _summary(node: ast.AST) -> str:
    docstring = ast.get_docstring(node, clean=True)
    if not docstring:
        return "No public docstring recorded; add one before treating this as stable API."
    paragraph = docstring.split("\n\n", 1)[0]
    return " ".join(line.strip() for line in paragraph.splitlines())


def _public_symbols(module_path: Path, module_name: str) -> list[Symbol]:
    source = module_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(module_path))
    symbols: list[Symbol] = []

    def visit_class(node: ast.ClassDef) -> None:
        if node.name.startswith("_"):
            return
        initializer = next(
            (
                child
                for child in node.body
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))
                and child.name == "__init__"
            ),
            None,
        )
        class_signature = (
            _signature(initializer, method=True) if initializer is not None else "()"
        )
        symbols.append(
            Symbol(
                module=module_name,
                qualified_name=f"{node.name}",
                kind="class",
                signature=class_signature,
                summary=_summary(node),
            )
        )
        for child in node.body:
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if child.name.startswith("_"):
                    continue
                symbols.append(
                    Symbol(
                        module=module_name,
                        qualified_name=f"{node.name}.{child.name}",
                        kind="method",
                        signature=_signature(
                            child, method=not _is_staticmethod(child)
                        ),
                        summary=_summary(child),
                    )
                )
            elif isinstance(child, ast.ClassDef):
                visit_class(child)

    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if not node.name.startswith("_"):
                symbols.append(
                    Symbol(
                        module=module_name,
                        qualified_name=node.name,
                        kind="function",
                        signature=_signature(node),
                        summary=_summary(node),
                    )
                )
        elif isinstance(node, ast.ClassDef):
            visit_class(node)
    return symbols


def collect_symbols(src_root: Path) -> list[Symbol]:
    """Collect public symbols from canonical core implementation modules."""

    symbols: list[Symbol] = []
    for path in sorted(src_root.rglob("*.py")):
        if path.name in EXCLUDED_SOURCE_FILES or "__pycache__" in path.parts:
            continue
        relative = path.relative_to(src_root).with_suffix("")
        module_name = "metainformant.core." + ".".join(relative.parts)
        symbols.extend(_public_symbols(path, module_name))
    return sorted(symbols, key=lambda item: (item.module, item.qualified_name))


def _cell(value: str) -> str:
    """Escape a value for a Markdown table cell."""

    return value.replace("|", r"\|").replace("`", "'").replace("\n", " ")


def render(symbols: list[Symbol], *, src_root: Path) -> str:
    """Render the generated API reference."""

    modules: dict[str, list[Symbol]] = {}
    for symbol in symbols:
        modules.setdefault(symbol.module, []).append(symbol)

    lines = [
        "# Core API Reference",
        "",
        "<!-- Generated by scripts/package/generate_core_api_reference.py; do not edit manually. -->",
        "",
        "This inventory is generated from the canonical public symbols under",
        f"`{src_root.relative_to(REPO_ROOT).as_posix()}/`. It is intentionally concise:",
        "topic guides contain behavioral guidance, examples, and caveats.",
        "Compatibility facades are excluded from the inventory and remain",
        "supported through the import paths described in [index.md](index.md).",
        "",
        f"**Public symbols:** {len(symbols)}",
        f"**Canonical modules:** {len(modules)}",
        "",
        "## Module index",
        "",
        "| Module | Public symbols |",
        "| --- | ---: |",
    ]
    for module, module_symbols in modules.items():
        anchor = module.replace(".", "-")
        lines.append(f"| [{module}](#{anchor}) | {len(module_symbols)} |")

    for module, module_symbols in modules.items():
        anchor = module.replace(".", "-")
        lines.extend(
            [
                "",
                f"## {module} {{#{anchor}}}",
                "",
                "| Symbol | Kind | Signature | Summary |",
                "| --- | --- | --- | --- |",
            ]
        )
        for symbol in module_symbols:
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"`{_cell(symbol.qualified_name)}`",
                        _cell(symbol.kind),
                        f"`{_cell(symbol.signature)}`",
                        _cell(symbol.summary),
                    ]
                )
                + " |"
            )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src-root", type=Path, default=SRC_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    src_root = args.src_root.expanduser().resolve()
    output = args.output.expanduser().resolve()
    symbols = collect_symbols(src_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(render(symbols, src_root=src_root), encoding="utf-8")
    print(f"Generated {len(symbols)} public core symbols across {len({s.module for s in symbols})} modules")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
