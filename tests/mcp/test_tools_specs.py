"""Spec-contract tests for the MCP tools package (real imports, no test doubles)."""

from __future__ import annotations

import json

import pytest

from metainformant.mcp.tools.catalog import iter_tool_specs


def _all_specs() -> list[dict]:
    return list(iter_tool_specs())


def test_tool_count_in_expected_range() -> None:
    n = len(_all_specs())
    assert 10 <= n <= 25, f"expected 10-25 tools, found {n}"


def test_names_unique_and_snake_case() -> None:
    names = [s["name"] for s in _all_specs()]
    assert len(names) == len(set(names))
    for name in names:
        assert name.replace("_", "").isalnum(), f"non-identifier tool name: {name}"


@pytest.mark.parametrize("key", ["name", "description", "input_schema", "handler", "writes"])
def test_every_spec_has_required_keys(key: str) -> None:
    for spec in _all_specs():
        assert key in spec, f"{spec.get('name')} missing {key}"


@pytest.mark.parametrize("writes", ["read-only", "output-dir-only"])
def test_every_spec_declares_allowed_writes_class(writes: str) -> None:
    classes = {s["writes"] for s in _all_specs()}
    assert classes <= {"read-only", "output-dir-only"}


def test_input_schemas_are_json_serializable_objects() -> None:
    for spec in _all_specs():
        schema = spec["input_schema"]
        assert schema.get("type") == "object"
        json.dumps(schema)  # raises if not serializable


def test_handlers_are_callable() -> None:
    for spec in _all_specs():
        assert callable(spec["handler"])
