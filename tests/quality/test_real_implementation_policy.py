"""Tests for the real-implementation policy scanner."""

from __future__ import annotations

from scripts.quality.verify_real_implementation_policy import scan_repo


def test_repo_real_implementation_policy_scan_passes() -> None:
    """Repository text should avoid old policy names and test-double APIs."""
    assert scan_repo() == []
