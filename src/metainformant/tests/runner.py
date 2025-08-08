from __future__ import annotations

import sys
from pathlib import Path


def run_all_tests(pytest_args: list[str] | None = None) -> int:
    """Run the repository's pytest suite programmatically.

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

    args = [str(test_dir)]
    if pytest_args:
        args.extend(pytest_args)
    return pytest.main(args)


