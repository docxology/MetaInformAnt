#!/usr/bin/env python3
"""Example discovery functionality for METAINFORMANT test runner."""

from __future__ import annotations

from pathlib import Path
from typing import List


def discover_examples(examples_dir: Path, domain_filter: str | None = None) -> List[Path]:
    """Discover all example Python files.

    Args:
        examples_dir: Path to examples directory
        domain_filter: Optional domain filter

    Returns:
        Sorted list of example file paths

    Raises:
        FileNotFoundError: If examples directory doesn't exist
    """
    if not examples_dir.exists():
        raise FileNotFoundError(f"Examples directory not found: {examples_dir}")

    all_examples = list(examples_dir.rglob("example_*.py"))

    if domain_filter:
        # Filter examples by domain
        filtered = []
        for example in all_examples:
            domain = example.parent.name
            if domain == domain_filter:
                filtered.append(example)
        all_examples = filtered

    return sorted(all_examples)





