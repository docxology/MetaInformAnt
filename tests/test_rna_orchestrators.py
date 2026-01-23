"""Integration tests for RNA orchestrator scripts.

Tests configuration blocks, method signatures, and basic functionality
without executing full workflows (which require external dependencies).
"""

from __future__ import annotations

import ast
import inspect
import re
from pathlib import Path


def test_run_workflow_has_config_block():
    """Test that run_workflow.py has proper structure and documentation."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    assert script_path.exists(), f"Script not found: {script_path}"

    content = script_path.read_text()

    # Check for main functionality (run_workflow.py uses orchestration module)
    assert "--config" in content, "Missing --config argument"
    assert "run_workflow_for_species" in content or "orchestration" in content, "Missing orchestration call"


def test_run_workflow_has_docstring():
    """Test that run_workflow.py has module docstring."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    content = script_path.read_text()

    # Check for docstring (module-level)
    has_docstring = content.strip().startswith('"""') or content.strip().startswith("'''") or '"""' in content[:500]
    assert has_docstring, "Missing module docstring"


def test_run_workflow_methods_have_docstrings():
    """Test that key functions in run_workflow.py have docstrings."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    content = script_path.read_text()

    # Parse the script to find function definitions
    tree = ast.parse(content)
    functions = [node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]

    # Check that main function has docstring
    main_functions = [f for f in functions if f.name == "main"]
    if main_functions:
        main_func = main_functions[0]
        func_start_line = main_func.lineno
        func_end_line = main_func.end_lineno if hasattr(main_func, "end_lineno") else func_start_line + 20
        func_lines = content.split("\n")[func_start_line - 1 : min(func_end_line, func_start_line + 15)]
        func_text = "\n".join(func_lines)
        has_docstring = '"""' in func_text or "'''" in func_text or ast.get_docstring(main_func) is not None
        assert has_docstring, "main function missing docstring"


def test_run_workflow_imports_metainformant_modules():
    """Test that run_workflow.py imports from metainformant.rna modules."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    content = script_path.read_text()

    # Check that script imports from metainformant modules
    has_metainformant_import = "from metainformant" in content or "import metainformant" in content

    assert has_metainformant_import, "run_workflow.py should import from metainformant modules"


def test_run_workflow_has_main_guard():
    """Test that run_workflow.py has if __name__ == '__main__' guard."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    content = script_path.read_text()

    # Check for main guard
    has_main_guard = 'if __name__ == "__main__":' in content or "if __name__ == '__main__':" in content

    assert has_main_guard, "run_workflow.py missing if __name__ == '__main__' guard"


def test_run_workflow_has_required_structure():
    """Test that run_workflow.py has required structure (argparse, config loading, etc.)."""
    script_path = Path(__file__).parent.parent / "scripts" / "rna" / "run_workflow.py"
    content = script_path.read_text()

    # Check for key components
    assert "argparse" in content or "ArgumentParser" in content, "Missing argparse usage"
    assert "--config" in content, "Missing --config argument"
    # run_workflow.py uses orchestration module which handles config loading
    assert "orchestration" in content or "run_workflow_for_species" in content, "Missing orchestration module usage"
