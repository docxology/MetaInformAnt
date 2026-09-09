#!/usr/bin/env python3
"""
Systematic Cross-Code Verification for METAINFORMANT Documentation

Thin orchestrator: builds the CLI, bootstraps ``src/``, and delegates all
AST-walking verification logic to ``metainformant.quality.doc_verification``.

Usage:
    python scripts/verify_documentation_code.py [--docs-dir DIR] [--src-dir DIR] [--output REPORT.md]
"""

import argparse
import sys
from pathlib import Path

# Path bootstrap: make the project's src/ library importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from metainformant.quality.doc_verification import run  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify code examples in documentation")
    parser.add_argument("--docs-dir", type=Path, default=Path("docs"), help="Documentation directory (default: docs/)")
    parser.add_argument("--src-dir", type=Path, default=Path("src"), help="Source code directory (default: src/)")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output") / "cross_code_verification_report.md",
        help="Output report file (default: output/cross_code_verification_report.md)",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--include-historical",
        action="store_true",
        help="Include historical audit and validation report snapshots",
    )
    parser.add_argument(
        "--strict-optional-imports",
        action="store_true",
        help="Treat optional third-party dependency imports as validation failures",
    )

    args = parser.parse_args()

    # Exit code: 0 if no violations, 1 if violations found
    sys.exit(run(args))


if __name__ == "__main__":
    main()
