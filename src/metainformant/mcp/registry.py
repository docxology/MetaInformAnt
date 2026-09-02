"""Declarative tool registry for the METAINFORMANT MCP server.

Each tool is described by a :class:`Tool` record with a name, description, a
JSON-Schema-style input schema, and a plain callable handler.  The registry
validates schemas at registration time and validates tool call arguments
against the schema at dispatch time, so the server layer can rely on every
registered tool being well-formed.

Only the subset of JSON Schema that the bundled tools need is enforced:
``type`` (object/string/integer/number/boolean/array), ``required``,
``properties``, and ``enum``.  This keeps the server dependency-free while
remaining deterministic and fully testable.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

JsonDict = dict[str, Any]

#: Mapping from JSON Schema type names to Python types used for validation.
_JSON_TYPE_MAP: dict[str, tuple[type, ...]] = {
    "object": (dict,),
    "array": (list,),
    "string": (str,),
    "boolean": (bool,),
    "integer": (int,),
    "number": (int, float),
}


class SchemaError(ValueError):
    """Raised when a tool schema is malformed or arguments fail validation."""


def validate_json_schema(schema: Mapping[str, Any]) -> None:
    """Validate that ``schema`` is a supported JSON Schema object definition.

    Raises :class:`SchemaError` for unsupported or contradictory schemas so
    that registration failures surface at startup rather than at call time.
    """

    if not isinstance(schema, dict):
        raise SchemaError("input schema must be a JSON object")
    schema_type = schema.get("type", "object")
    if schema_type != "object":
        raise SchemaError(f"unsupported top-level schema type: {schema_type!r}")
    properties = schema.get("properties", {})
    if not isinstance(properties, dict):
        raise SchemaError("schema 'properties' must be an object")
    for name, prop in properties.items():
        if not isinstance(prop, dict):
            raise SchemaError(f"property {name!r} must be a schema object")
        prop_type = prop.get("type")
        if prop_type not in _JSON_TYPE_MAP:
            raise SchemaError(f"property {name!r} has unsupported type {prop_type!r}")
        if "enum" in prop:
            enum_values = prop["enum"]
            if not isinstance(enum_values, list) or not enum_values:
                raise SchemaError(f"property {name!r} enum must be a non-empty array")
    required = schema.get("required", [])
    if not isinstance(required, list) or not all(isinstance(item, str) for item in required):
        raise SchemaError("schema 'required' must be an array of property names")
    missing = set(required) - set(properties)
    if missing:
        raise SchemaError(f"required properties not declared: {sorted(missing)}")


def validate_arguments(arguments: Mapping[str, Any], schema: Mapping[str, Any]) -> None:
    """Validate tool call ``arguments`` against ``schema``.

    Raises :class:`SchemaError` with a human-readable message on the first
    violation.  Unknown properties are rejected only when the schema sets
    ``"additionalProperties": false`` (permissive JSON Schema semantics by
    default).
    """

    validate_json_schema(schema)
    for name in schema.get("required", []):
        if name not in arguments:
            raise SchemaError(f"missing required argument: {name!r}")
    properties = schema.get("properties", {})
    for name, value in arguments.items():
        prop = properties.get(name)
        if prop is None:
            if schema.get("additionalProperties") is False:
                raise SchemaError(f"unexpected argument: {name!r}")
            continue
        expected = prop["type"]
        python_types = _JSON_TYPE_MAP[expected]
        # bool is a subclass of int; exclude it explicitly from numeric slots.
        if expected in ("integer", "number") and isinstance(value, bool):
            raise SchemaError(f"argument {name!r} must be of type {expected}")
        if not isinstance(value, python_types):
            raise SchemaError(f"argument {name!r} must be of type {expected}")
        if "enum" in prop and value not in prop["enum"]:
            raise SchemaError(f"argument {name!r} must be one of {prop['enum']!r}, got {value!r}")


@dataclass(frozen=True)
class Tool:
    """A single MCP-exposed tool."""

    name: str
    description: str
    input_schema: JsonDict
    handler: Callable[[Mapping[str, Any]], Any]
    metadata: JsonDict = field(default_factory=dict)

    def validate_call(self, arguments: Mapping[str, Any]) -> None:
        """Validate ``arguments`` against this tool's input schema."""

        validate_arguments(arguments, self.input_schema)


class ToolRegistry:
    """Registry of MCP tools with schema validation and safe dispatch."""

    def __init__(self) -> None:
        self._tools: dict[str, Tool] = {}

    # -- registration ------------------------------------------------------

    def register(self, tool: Tool) -> Tool:
        """Register ``tool``; returns the validated tool record."""

        if not tool.name or not isinstance(tool.name, str):
            raise SchemaError("tool name must be a non-empty string")
        if not isinstance(tool.description, str) or not tool.description.strip():
            raise SchemaError(f"tool {tool.name!r} must have a description")
        validate_json_schema(tool.input_schema)
        if not callable(tool.handler):
            raise SchemaError(f"tool {tool.name!r} handler must be callable")
        if tool.name in self._tools:
            raise SchemaError(f"tool already registered: {tool.name!r}")
        self._tools[tool.name] = tool
        return tool

    def register_module_tools(self, module_name: str) -> list[str]:
        """Register every ``TOOL``/``TOOLS`` constant exported by ``module_name``.

        Tool modules declare a module-level ``TOOL`` (a :class:`Tool`) or a
        ``TOOLS`` sequence of them.  Returns the names that were registered.
        Unknown or empty modules raise so wiring mistakes fail loudly.
        """

        module = importlib.import_module(module_name)
        candidates: Sequence[Any]
        if hasattr(module, "TOOLS"):
            declared = module.TOOLS
            if not isinstance(declared, (list, tuple)) or not declared:
                raise ValueError(f"module {module_name!r} TOOLS must be a non-empty sequence")
            candidates = list(declared)
        elif hasattr(module, "TOOL"):
            candidates = [module.TOOL]
        else:
            raise ValueError(f"module {module_name!r} declares no TOOL/TOOLS constant")
        registered: list[str] = []
        for candidate in candidates:
            if not isinstance(candidate, Tool):
                raise TypeError(f"module {module_name!r} declared a non-Tool entry: {candidate!r}")
            self.register(candidate)
            registered.append(candidate.name)
        return registered

    # -- lookup / dispatch -------------------------------------------------

    def get(self, name: str) -> Tool:
        """Return the tool registered under ``name`` or raise ``KeyError``."""

        return self._tools[name]

    def names(self) -> list[str]:
        """Sorted list of registered tool names (deterministic ordering)."""

        return sorted(self._tools)

    def list_tools(self) -> list[JsonDict]:
        """MCP ``tools/list`` descriptors for all registered tools."""

        return [
            {
                "name": tool.name,
                "description": tool.description,
                "inputSchema": tool.input_schema,
            }
            for tool in (self._tools[name] for name in sorted(self._tools))
        ]

    def call(self, name: str, arguments: Mapping[str, Any] | None) -> Any:
        """Validate arguments and execute the handler for tool ``name``."""

        tool = self._tools[name]
        tool.validate_call(arguments or {})
        return tool.handler(arguments or {})
