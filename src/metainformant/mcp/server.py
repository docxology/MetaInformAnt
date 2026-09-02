"""Stdio JSON-RPC 2.0 MCP server for METAINFORMANT.

Implements the core Model Context Protocol surface over stdio using only the
standard library: ``initialize``, ``tools/list``, ``tools/call``,
``resources/list``, and ``resources/read``, plus graceful shutdown via EOF or
an ``exit`` notification.

Protocol notes:

- Each request is a single line of JSON (newline-delimited JSON-RPC 2.0).
- Responses are a single line of JSON on stdout; nothing else is ever written
  to stdout.  Diagnostics go to stderr.
- Unknown methods return error code ``-32601`` (method not found); malformed
  JSON returns ``-32700``; a non-object request returns ``-32600``.
- A notification (a message without an ``id``) produces no response, except
  that ``exit`` terminates the server cleanly.
- The server supports protocol version ``2024-11-05`` and negotiates down to
  the newest version the client offers that we know; an unknown version still
  initializes (per MCP spec the server replies with its own version).

Resources: the server exposes a small set of read-only ``file://`` resources
describing the platform (capability index, methods documentation), so clients
can discover how to use the tools.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Iterable, Mapping
from typing import Any, TextIO

from metainformant.mcp.registry import SchemaError, ToolRegistry

__all__ = ["MCPServer", "build_default_registry", "serve_stdio"]

PROTOCOL_VERSION = "2024-11-05"

# JSON-RPC error codes used by this server.
PARSE_ERROR = -32700
INVALID_REQUEST = -32600
METHOD_NOT_FOUND = -32601
INVALID_PARAMS = -32602
INTERNAL_ERROR = -32603

_SERVER_INFO = {
    "name": "metainformant-mcp",
    "version": "0.1.0",
}

_CAPABILITIES = {
    "tools": {"listChanged": False},
    "resources": {"listChanged": False},
}


class MCPServer:
    """JSON-RPC dispatcher over a :class:`ToolRegistry`."""

    def __init__(self, registry: ToolRegistry) -> None:
        self.registry = registry

    # -- resources ---------------------------------------------------------

    def list_resources(self) -> list[dict[str, Any]]:
        """Read-only capability resources exposed by the server."""

        return [
            {
                "uri": "metainformant://capabilities",
                "name": "METAINFORMANT MCP capability index",
                "description": (
                    "Lists registered tools, their schemas, and the JSON-RPC methods " "this server implements."
                ),
                "mimeType": "application/json",
            },
            {
                "uri": "metainformant://methods",
                "name": "METAINFORMANT method documentation",
                "description": (
                    "Human-readable documentation of supported JSON-RPC methods and "
                    "how to connect a client over stdio."
                ),
                "mimeType": "text/markdown",
            },
        ]

    def read_resource(self, uri: str) -> dict[str, Any]:
        """Return the contents of a built-in resource by URI."""

        if uri == "metainformant://capabilities":
            body = json.dumps(
                {
                    "server": _SERVER_INFO,
                    "protocol_version": PROTOCOL_VERSION,
                    "methods": [
                        "initialize",
                        "tools/list",
                        "tools/call",
                        "resources/list",
                        "resources/read",
                    ],
                    "tools": self.registry.list_tools(),
                },
                indent=2,
            )
            return {"uri": uri, "mimeType": "application/json", "text": body}
        if uri == "metainformant://methods":
            body = (
                "# METAINFORMANT MCP server\n\n"
                "Connect with any MCP stdio client via:\n\n"
                "```\n"
                "uv run python -m metainformant.mcp.server\n"
                "```\n\n"
                "Messages are newline-delimited JSON-RPC 2.0 on stdin/stdout.\n"
                "Methods: initialize, tools/list, tools/call, resources/list,\n"
                "resources/read, and the `exit` notification for shutdown.\n"
            )
            return {"uri": uri, "mimeType": "text/markdown", "text": body}
        raise KeyError(f"unknown resource: {uri!r}")

    # -- JSON-RPC dispatch ---------------------------------------------------

    def handle_request(self, message: Mapping[str, Any]) -> dict[str, Any] | None:
        """Handle one decoded JSON-RPC message; returns the response (or None).

        Notifications (no ``id``) return ``None`` except ``exit``, which raises
        ``SystemExit(0)`` for a clean shutdown.
        """

        if not isinstance(message, dict):
            return _error_response(None, INVALID_REQUEST, "request must be a JSON object")
        method = message.get("method")
        request_id = message.get("id")
        if not isinstance(method, str) or not method:
            return _error_response(request_id, INVALID_REQUEST, "missing method")
        params = message.get("params") or {}
        if not isinstance(params, dict):
            return _error_response(request_id, INVALID_PARAMS, "params must be an object")
        if "id" not in message:
            if method == "exit":
                raise SystemExit(0)
            return None

        if method == "initialize":
            requested = params.get("protocolVersion")
            version = requested if isinstance(requested, str) else PROTOCOL_VERSION
            return {
                "jsonrpc": "2.0",
                "id": request_id,
                "result": {
                    "protocolVersion": version,
                    "capabilities": _CAPABILITIES,
                    "serverInfo": _SERVER_INFO,
                },
            }
        if method == "ping":
            # MCP liveness probe.
            return {"jsonrpc": "2.0", "id": request_id, "result": {}}
        if method == "tools/list":
            return {"jsonrpc": "2.0", "id": request_id, "result": {"tools": self.registry.list_tools()}}
        if method == "tools/call":
            return self._handle_tools_call(request_id, params)
        if method == "resources/list":
            return {
                "jsonrpc": "2.0",
                "id": request_id,
                "result": {"resources": self.list_resources()},
            }
        if method == "resources/read":
            uri = params.get("uri")
            if not isinstance(uri, str) or not uri:
                return _error_response(request_id, INVALID_PARAMS, "params.uri must be a string")
            try:
                contents = self.read_resource(uri)
            except KeyError:
                return _error_response(request_id, INVALID_PARAMS, f"unknown resource: {uri}")
            return {"jsonrpc": "2.0", "id": request_id, "result": {"contents": [contents]}}
        return _error_response(request_id, METHOD_NOT_FOUND, f"method not found: {method}")

    def _handle_tools_call(self, request_id: Any, params: Mapping[str, Any]) -> dict[str, Any]:
        name = params.get("name")
        if not isinstance(name, str) or not name:
            return _error_response(request_id, INVALID_PARAMS, "params.name must be a string")
        arguments = params.get("arguments") or {}
        if not isinstance(arguments, dict):
            return _error_response(request_id, INVALID_PARAMS, "params.arguments must be an object")
        try:
            result = self.registry.call(name, arguments)
        except KeyError:
            return _error_response(request_id, INVALID_PARAMS, f"unknown tool: {name!r}")
        except SchemaError as exc:
            return _error_response(request_id, INVALID_PARAMS, str(exc))
        except Exception as exc:  # noqa: BLE001 - tool handlers are user code
            return _error_response(request_id, INTERNAL_ERROR, f"tool {name!r} failed: {exc}")
        return {
            "jsonrpc": "2.0",
            "id": request_id,
            "result": {
                "content": [{"type": "text", "text": json.dumps(result, indent=2)}],
                "structuredContent": result,
                "isError": False,
            },
        }

    # -- stdio loop ----------------------------------------------------------

    def serve(
        self,
        stdin: TextIO | None = None,
        stdout: TextIO | None = None,
        stderr: TextIO | None = None,
    ) -> int:
        """Serve newline-delimited JSON-RPC until EOF or an ``exit`` notification."""

        stdin = stdin if stdin is not None else sys.stdin
        stdout = stdout if stdout is not None else sys.stdout
        stderr = stderr if stderr is not None else sys.stderr
        for line in stdin:
            line = line.strip()
            if not line:
                continue
            response = self.process_line(line, stderr)
            if response is not None:
                stdout.write(json.dumps(response) + "\n")
                stdout.flush()
        return 0

    def process_line(self, line: str, stderr: TextIO | None = None) -> dict[str, Any] | None:
        """Decode one line and dispatch it; swallow protocol errors into responses."""

        stderr = stderr if stderr is not None else sys.stderr
        try:
            message = json.loads(line)
        except json.JSONDecodeError as exc:
            return _error_response(None, PARSE_ERROR, f"parse error: {exc}")
        try:
            return self.handle_request(message)
        except SystemExit:
            raise
        except Exception as exc:  # noqa: BLE001 - server must never crash the stream
            print(f"metainformant.mcp: internal error: {exc}", file=stderr)
            return _error_response(
                message.get("id") if isinstance(message, dict) else None,
                INTERNAL_ERROR,
                "internal server error",
            )


def build_default_registry() -> ToolRegistry:
    """Build the production registry: bundled tool modules, discovered declaratively."""

    from metainformant.mcp import TOOLS_MODULES

    registry = ToolRegistry()
    for module_name in TOOLS_MODULES:
        registry.register_module_tools(module_name)
    return registry


def serve_stdio(
    registry: Iterable[ToolRegistry] | None = None,
    stdin: TextIO | None = None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Run the stdio server loop with the default registry."""

    active = registry if registry is not None else build_default_registry()
    server = active if isinstance(active, MCPServer) else MCPServer(active)  # type: ignore[arg-type]
    return server.serve(stdin=stdin, stdout=stdout, stderr=stderr)


def main() -> int:
    """Console entry point for ``python -m metainformant.mcp.server``."""

    server = MCPServer(build_default_registry())
    return server.serve()


def _error_response(request_id: Any, code: int, message: str) -> dict[str, Any]:
    return {
        "jsonrpc": "2.0",
        "id": request_id,
        "error": {"code": code, "message": message},
    }


if __name__ == "__main__":
    raise SystemExit(main())
