"""Tests for the MCP tool registry (real code, real data; zero mocks)."""

from __future__ import annotations

from typing import Any

import pytest

from metainformant.mcp import registry as mcp_registry
from metainformant.mcp.registry import SchemaError, Tool, ToolRegistry, validate_arguments


def _identity_tool(name: str = "echo", **schema_overrides: Any) -> Tool:
    schema: dict[str, Any] = {
        "type": "object",
        "properties": {"message": {"type": "string"}},
        "required": ["message"],
    }
    schema.update(schema_overrides)
    return Tool(
        name=name,
        description="echoes its message argument",
        input_schema=schema,
        handler=lambda args: {"echo": args["message"]},
    )


def test_register_and_list_tools_round_trip() -> None:
    reg = ToolRegistry()
    reg.register(_identity_tool())

    assert reg.names() == ["echo"]
    descriptors = reg.list_tools()
    assert descriptors == [
        {
            "name": "echo",
            "description": "echoes its message argument",
            "inputSchema": {
                "type": "object",
                "properties": {"message": {"type": "string"}},
                "required": ["message"],
            },
        }
    ]


def test_call_validates_then_executes_handler() -> None:
    reg = ToolRegistry()
    reg.register(_identity_tool())

    assert reg.call("echo", {"message": "hello"}) == {"echo": "hello"}
    with pytest.raises(SchemaError, match="missing required argument"):
        reg.call("echo", {})


def test_duplicate_registration_rejected() -> None:
    reg = ToolRegistry()
    reg.register(_identity_tool())
    with pytest.raises(SchemaError, match="already registered"):
        reg.register(_identity_tool())


@pytest.mark.parametrize(
    "schema",
    [
        [],  # not an object
        {"type": "string"},  # unsupported top-level type
        {"type": "object", "properties": {"a": "not-a-schema"}},
        {"type": "object", "properties": {"a": {"type": "nonsense"}}},
        {"type": "object", "properties": {}, "required": "a"},  # required not a list
        {"type": "object", "properties": {}, "required": ["ghost"]},  # undeclared required
        {"type": "object", "properties": {"e": {"type": "string", "enum": []}}},
    ],
)
def test_invalid_schemas_rejected_at_registration(schema: Any) -> None:
    reg = ToolRegistry()
    with pytest.raises(SchemaError):
        reg.register(
            Tool(
                name="bad",
                description="invalid schema tool",
                input_schema=schema,
                handler=lambda args: None,
            )
        )


def test_argument_type_and_enum_validation() -> None:
    schema = {
        "type": "object",
        "properties": {
            "count": {"type": "integer"},
            "mode": {"type": "string", "enum": ["fast", "slow"]},
        },
        "required": ["count"],
    }
    validate_arguments({"count": 3, "mode": "fast"}, schema)
    with pytest.raises(SchemaError, match="must be of type integer"):
        validate_arguments({"count": "3"}, schema)
    # bool is not an integer for schema purposes
    with pytest.raises(SchemaError, match="must be of type integer"):
        validate_arguments({"count": True}, schema)
    with pytest.raises(SchemaError, match="must be one of"):
        validate_arguments({"count": 1, "mode": "turbo"}, schema)


def test_additional_properties_false_rejects_unknown_arguments() -> None:
    schema = {
        "type": "object",
        "properties": {"a": {"type": "string"}},
        "required": ["a"],
        "additionalProperties": False,
    }
    validate_arguments({"a": "ok"}, schema)
    with pytest.raises(SchemaError, match="unexpected argument"):
        validate_arguments({"a": "ok", "b": 1}, schema)


def test_register_module_tools_wires_declared_tool_constants() -> None:
    # The real bundled adapter module declares TOOL/TOOLS and must import cleanly.
    reg = ToolRegistry()
    registered = reg.register_module_tools("metainformant.mcp.tool_adapters")
    assert registered == ["amalgkit_monitor"]
    tool = reg.get("amalgkit_monitor")
    assert "Amalgkit" in tool.description
    assert tool.input_schema["type"] == "object"


def test_register_module_tools_fails_loudly_on_missing_module() -> None:
    reg = ToolRegistry()
    with pytest.raises(ImportError):
        reg.register_module_tools("metainformant.mcp.definitely_not_a_module")


def test_default_registry_builds_bundled_tools() -> None:
    from metainformant.mcp.server import build_default_registry

    reg = build_default_registry()
    assert "amalgkit_monitor" in reg.names()


def test_public_registry_aliases_exist() -> None:
    assert mcp_registry.SchemaError is SchemaError
    assert mcp_registry.ToolRegistry is ToolRegistry
