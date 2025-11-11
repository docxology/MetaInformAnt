from __future__ import annotations

import sys
from pathlib import Path


def run_all_tests(pytest_args: list[str] | None = None) -> int:
    """Run the repository's pytest suite programmatically.

    Ensures the metainformant package is importable by adding src to Python path
    if not already available. This allows tests to run whether the package is
    installed in editable mode or not.

    Returns the pytest exit code (0 = success).
    """
    try:
        import pytest  # type: ignore
    except Exception as exc:  # pragma: no cover
        print(f"pytest not available: {exc}")
        return 2

    repo_root = Path(__file__).resolve().parents[3]
    test_dir = repo_root / "tests"
    if not test_dir.is_dir():
        print("No tests directory found.")
        return 5

    # Ensure metainformant package is importable
    # Check if package is already importable (installed in editable mode)
    try:
        import metainformant  # noqa: F401
    except ImportError:
        # Package not installed, add src to path
        src_dir = repo_root / "src"
        if src_dir.exists() and str(src_dir) not in sys.path:
            sys.path.insert(0, str(src_dir))

    args = [str(test_dir)]
    if pytest_args:
        args.extend(pytest_args)
    return pytest.main(args)
