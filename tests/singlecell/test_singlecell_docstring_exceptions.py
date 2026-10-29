"""The ``Raises:`` sections of the singlecell docstrings name the exceptions the code raises.

Argument checks in this package go through ``validation.validate_type`` / ``validate_range`` and
``errors.ValidationError``, which is not a subclass of ``ValueError`` or ``TypeError``. About fifty docstrings
nevertheless promised ``TypeError`` and ``ValueError`` for exactly those checks, so a caller who trusted them and
wrote ``except ValueError`` was not protected. The scan below reads the source with ``ast``: a function whose body
raises ``ValidationError`` (directly or through a validator) and no ``ValueError`` or ``TypeError`` of its own may not
document one of those two, unless a library it calls raises it; those few are listed with the library.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

import metainformant.singlecell as singlecell
from metainformant.core.utils.errors import IOError as PackageIOError, ValidationError
from metainformant.singlecell.data.preprocessing import load_count_matrix

SINGLECELL_ROOT = Path(singlecell.__file__).resolve().parent
VALIDATORS = {
    "validate_type",
    "validate_range",
    "validate_not_none",
    "validate_not_empty",
    "validate_path_exists",
    "validate_path_is_file",
    "validate_path_is_dir",
    "validate_path_within",
}
BUILTINS_FOR_ARGUMENT_CHECKS = {"TypeError", "ValueError"}

# Functions whose documented ValueError is raised by a library they call; all of their own checks raise ValidationError.
RAISED_BY_A_LIBRARY = {
    ("analysis/clustering.py", "hierarchical_clustering"): "scipy rejects an unknown distance metric",
    ("analysis/nonlinear_methods.py", "tsne_reduction"): "scikit-learn rejects an invalid parameter",
}


def _direct_raises(function: ast.FunctionDef | ast.AsyncFunctionDef) -> set[str]:
    """Exception names raised by the function's own ``raise`` statements, plus ValidationError via validators."""
    names: set[str] = set()

    class Visitor(ast.NodeVisitor):
        def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
            if node is function:  # a nested function raises for itself
                self.generic_visit(node)

        visit_AsyncFunctionDef = visit_FunctionDef  # type: ignore[assignment]

        def visit_Raise(self, node: ast.Raise) -> None:
            if node.exc is not None:
                target = node.exc.func if isinstance(node.exc, ast.Call) else node.exc
                if isinstance(target, ast.Attribute):
                    names.add(target.attr)
                elif isinstance(target, ast.Name):
                    names.add(target.id)
            self.generic_visit(node)

        def visit_Call(self, node: ast.Call) -> None:
            callee = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", None)
            if callee in VALIDATORS:
                names.add("ValidationError")
            self.generic_visit(node)

    Visitor().visit(function)
    return names


def _documented_raises(docstring: str | None) -> set[str]:
    """Exception names listed in the docstring's Google-style ``Raises:`` section."""
    documented: set[str] = set()
    if not docstring:
        return documented
    in_section, base = False, 0
    for line in docstring.splitlines():
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        if re.fullmatch(r"Raises:", stripped):
            in_section, base = True, indent
        elif in_section and stripped and indent <= base:
            in_section = False
        elif in_section and stripped and indent == base + 4:
            entry = re.match(r"([A-Za-z_][\w.]*)\s*:", stripped)
            if entry:
                documented.add(entry.group(1).split(".")[-1])
    return documented


def _offenders(source: str, relative: str) -> tuple[list[str], int, int]:
    """(offences, functions with a Raises section, of those how many raise ValidationError)."""
    offences: list[str] = []
    with_section = raising_validation = 0
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        documented = _documented_raises(ast.get_docstring(node, clean=False))
        if not documented:
            continue
        with_section += 1
        raised = _direct_raises(node)
        if "ValidationError" not in raised:
            continue
        raising_validation += 1
        promised = (documented & BUILTINS_FOR_ARGUMENT_CHECKS) - raised
        if promised and (relative, node.name) not in RAISED_BY_A_LIBRARY:
            offences.append(
                f"{relative}:{node.lineno} {node.name} documents {sorted(promised)} "
                "but its checks raise ValidationError"
            )
    return offences, with_section, raising_validation


class TestRaisesSectionsMatchTheCode:
    def test_no_function_documents_a_builtin_where_it_raises_validation_error(self) -> None:
        offences: list[str] = []
        with_section = raising_validation = 0
        for path in sorted(SINGLECELL_ROOT.rglob("*.py")):
            found, sections, raising = _offenders(path.read_text(), path.relative_to(SINGLECELL_ROOT).as_posix())
            offences.extend(found)
            with_section += sections
            raising_validation += raising

        # 70 and 45 when written: a scan that sees far fewer is not looking at the package any more
        assert with_section >= 60 and raising_validation >= 40, "the scan no longer sees the package it guards"
        assert offences == []

    def test_the_scan_flags_the_old_pattern(self) -> None:
        old_style = '''
def check(data):
    """Check.

    Raises:
        TypeError: If data is not SingleCellData
        ValueError: If a column is not found
    """
    validation.validate_type(data, SingleCellData, "data")
    if "x" not in data.obs:
        raise errors.ValidationError("column not found")
'''
        corrected = old_style.replace(
            "TypeError: If data is not SingleCellData\n        ValueError:", "ValidationError:"
        )

        flagged, _, _ = _offenders(old_style, "demo.py")
        clean, _, _ = _offenders(corrected, "demo.py")

        assert flagged and "['TypeError', 'ValueError']" in flagged[0]
        assert clean == []

    def test_a_builtin_that_the_function_itself_raises_may_still_be_documented(self) -> None:
        mixed = '''
def check(method):
    """Check.

    Raises:
        ValueError: If method is unknown
    """
    validation.validate_range(1, min_val=0)
    if method != "x":
        raise ValueError("unknown method")
'''

        assert _offenders(mixed, "demo.py")[0] == []

    def test_the_documented_library_exceptions_are_really_raised_by_the_library(self) -> None:
        import numpy as np

        from metainformant.singlecell.analysis.clustering import hierarchical_clustering
        from metainformant.singlecell.data.preprocessing import SingleCellData

        data = SingleCellData(np.random.default_rng(0).random((30, 6)))

        with pytest.raises(ValueError, match="Unknown Distance Metric") as library_error:
            hierarchical_clustering(data, n_clusters=3, linkage_method="average", metric="not_a_metric")
        assert not isinstance(library_error.value, ValidationError)
        with pytest.raises(ValidationError, match="Invalid linkage method"):
            hierarchical_clustering(data, n_clusters=3, linkage_method="not_a_linkage")


class TestLoadCountMatrixRaisesWhatItDocuments:
    """``load_count_matrix`` documented FileNotFoundError and ValueError; it raises neither for these cases."""

    def test_a_missing_file_raises_validation_error_and_the_docstring_says_so(self, tmp_path: Path) -> None:
        with pytest.raises(ValidationError, match="does not exist") as raised:
            load_count_matrix(tmp_path / "missing.h5ad")

        assert not isinstance(raised.value, FileNotFoundError)
        assert isinstance(raised.value, ValueError)  # ValidationError is also a ValueError
        assert "FileNotFoundError" not in (load_count_matrix.__doc__ or "")

    def test_the_documented_exceptions_are_exactly_the_ones_the_code_raises(self) -> None:
        assert _documented_raises(load_count_matrix.__doc__) == {"ValidationError", "ConfigError", "IOError"}

    def test_an_unsupported_format_raises_validation_error(self, tmp_path: Path) -> None:
        source = tmp_path / "counts.txt"
        source.write_text("gene\tcell1\ng1\t1\n")

        with pytest.raises(ValidationError, match="Unsupported format: xyz") as raised:
            load_count_matrix(source, format="xyz")

        assert isinstance(raised.value, ValueError)  # ValidationError is also a ValueError

    def test_a_matrix_market_load_without_gene_and_barcode_files_raises_the_packages_ioerror(
        self, tmp_path: Path
    ) -> None:
        matrix = tmp_path / "matrix.mtx"
        matrix.write_text("%%MatrixMarket matrix coordinate integer general\n1 1 1\n1 1 1\n")

        with pytest.raises(PackageIOError, match="genes_file and barcodes_file") as raised:
            load_count_matrix(matrix, format="mtx")

        assert isinstance(raised.value.__cause__, ValueError)  # the ValueError is chained, not raised as such
        assert not isinstance(raised.value, (OSError, ValueError))
        documented = load_count_matrix.__doc__ or ""
        assert "ValidationError" in documented and "IOError" in documented and "ConfigError" in documented
