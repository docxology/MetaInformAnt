"""Tests for the stdio JSON-RPC MCP server.

Real implementations only: in-process dispatch plus a genuine subprocess
round-trip against ``python -m metainformant.mcp.server`` (real implementations).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.mcp.registry import SchemaError, Tool, ToolRegistry
from metainformant.mcp.server import PROTOCOL_VERSION, MCPServer


def _make_registry() -> ToolRegistry:
    reg = ToolRegistry()
    reg.register(
        Tool(
            name="add",
            description="adds two integers",
            input_schema={
                "type": "object",
                "properties": {
                    "a": {"type": "integer"},
                    "b": {"type": "integer"},
                },
                "required": ["a", "b"],
            },
            handler=lambda args: {"sum": args["a"] + args["b"]},
        )
    )
    return reg


def _request(method: str, **kwargs: object) -> dict[str, object]:
    payload: dict[str, object] = {"jsonrpc": "2.0", "id": 1, "method": method}
    payload.update(kwargs)
    return payload


def test_initialize_handshake_negotiates_protocol_version() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(
        _request("initialize", params={"protocolVersion": PROTOCOL_VERSION, "clientInfo": {"name": "t"}})
    )
    assert response is not None and "result" in response
    assert response["result"]["protocolVersion"] == PROTOCOL_VERSION
    assert response["result"]["serverInfo"]["name"] == "metainformant-mcp"
    assert "tools" in response["result"]["capabilities"]
    assert "resources" in response["result"]["capabilities"]


def test_tools_list_reports_complete_registry() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(_request("tools/list"))
    tools = response["result"]["tools"]
    assert [t["name"] for t in tools] == ["add"]
    assert tools[0]["inputSchema"]["required"] == ["a", "b"]


def test_tools_call_round_trip_executes_handler() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(_request("tools/call", params={"name": "add", "arguments": {"a": 2, "b": 3}}))
    assert response["result"]["isError"] is False
    assert response["result"]["structuredContent"] == {"sum": 5}
    assert json.loads(response["result"]["content"][0]["text"]) == {"sum": 5}


def test_tools_call_schema_violation_returns_invalid_params() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(
        _request("tools/call", params={"name": "add", "arguments": {"a": "not-an-int", "b": 1}})
    )
    assert response["error"]["code"] == -32602
    assert "must be of type integer" in response["error"]["message"]


def test_tools_call_unknown_tool_returns_invalid_params() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(_request("tools/call", params={"name": "nope", "arguments": {}}))
    assert response["error"]["code"] == -32602
    assert "unknown tool" in response["error"]["message"]


def test_unknown_method_returns_method_not_found() -> None:
    server = MCPServer(_make_registry())
    response = server.handle_request(_request("wat/nope"))
    assert response["error"]["code"] == -32601


def test_malformed_json_line_returns_parse_error() -> None:
    import io

    server = MCPServer(_make_registry())
    response = server.process_line("{not json", stderr=io.StringIO())
    assert response is not None and response["error"]["code"] == -32700


def test_notification_without_id_produces_no_response() -> None:
    server = MCPServer(_make_registry())
    assert server.handle_request({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None


def test_resources_list_and_read_round_trip() -> None:
    server = MCPServer(_make_registry())
    listing = server.handle_request(_request("resources/list"))
    uris = [r["uri"] for r in listing["result"]["resources"]]
    assert "metainformant://capabilities" in uris
    read = server.handle_request(_request("resources/read", params={"uri": "metainformant://capabilities"}))
    body = json.loads(read["result"]["contents"][0]["text"])
    assert body["protocol_version"] == PROTOCOL_VERSION
    assert any(t["name"] == "add" for t in body["tools"])
    missing = server.handle_request(_request("resources/read", params={"uri": "metainformant://ghost"}))
    assert missing["error"]["code"] == -32602


def test_server_loop_serves_lines_until_eof(tmp_path: Path) -> None:
    import io

    inbox = io.StringIO(
        json.dumps(_request("initialize", params={"protocolVersion": PROTOCOL_VERSION}))
        + "\n"
        + json.dumps(_request("tools/list"))
        + "\n"
    )
    out = io.StringIO()
    server = MCPServer(_make_registry())
    server.serve(stdin=inbox, stdout=out, stderr=io.StringIO())
    lines = [json.loads(line) for line in out.getvalue().splitlines()]
    assert len(lines) == 2
    assert lines[0]["result"]["protocolVersion"] == PROTOCOL_VERSION
    assert lines[1]["result"]["tools"][0]["name"] == "add"


def test_subprocess_round_trip_initialize_tools_list_and_call(tmp_path: Path) -> None:
    """Full stdio round-trip against the real server process."""

    repo_src = Path(__file__).resolve().parents[2] / "src"
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(repo_src) + (os.pathsep + existing if existing else "")

    messages = [
        _request("initialize", params={"protocolVersion": PROTOCOL_VERSION}),
        {"jsonrpc": "2.0", "method": "notifications/initialized"},
        _request("tools/list"),
        _request("tools/call", params={"name": "amalgkit_monitor", "arguments": {"inspect_processes": False}}),
    ]
    stdin_payload = "".join(json.dumps(m) + "\n" for m in messages)

    proc = subprocess.run(
        [sys.executable, "-m", "metainformant.mcp.server"],
        input=stdin_payload,
        capture_output=True,
        text=True,
        env=env,
        timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    lines = [json.loads(line) for line in proc.stdout.splitlines()]
    assert len(lines) == 3  # notification produced no response

    assert lines[0]["result"]["protocolVersion"] == PROTOCOL_VERSION
    tool_names = [t["name"] for t in lines[1]["result"]["tools"]]
    assert "amalgkit_monitor" in tool_names
    call_result = lines[2]["result"]
    assert call_result["isError"] is False
    snapshot = call_result["structuredContent"]
    assert snapshot["status"] in {"running", "stopped"}
    assert snapshot["readiness"]["biological_inference"] == "withheld"


def test_subprocess_unknown_method_surfaces_jsonrpc_error(tmp_path: Path) -> None:
    repo_src = Path(__file__).resolve().parents[2] / "src"
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(repo_src) + (os.pathsep + existing if existing else "")

    proc = subprocess.run(
        [sys.executable, "-m", "metainformant.mcp.server"],
        input=json.dumps(_request("bogus/method")) + "\n",
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
    )
    assert proc.returncode == 0
    response = json.loads(proc.stdout.splitlines()[0])
    assert response["error"]["code"] == -32601


def test_exit_notification_shuts_down_cleanly() -> None:
    server = MCPServer(_make_registry())
    with pytest.raises(SystemExit) as excinfo:
        server.handle_request({"jsonrpc": "2.0", "method": "exit"})
    assert excinfo.value.code == 0


def test_schema_error_is_exported_for_clients() -> None:
    from metainformant.mcp import registry

    assert registry.SchemaError is SchemaError
