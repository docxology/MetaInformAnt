"""Tests for OllamaClient HTTP behavior against a real local HTTP server.

These tests exercise the real urllib request/streaming/error paths of the
client without requiring a live Ollama installation. A threaded
http.server on a random localhost port serves canned Ollama-shaped
responses (real implementations; no mocks).
"""

from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from metainformant.ml.llm.ollama.client import ChatMessage, OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig


class _OllamaStubHandler(BaseHTTPRequestHandler):
    """Minimal Ollama-compatible API stub backed by a real HTTP server."""

    def log_message(self, *args):  # silence request logging
        pass

    def _respond(self, payload: dict, ndjson: bool = False) -> None:
        if ndjson:
            body = ("\n".join(json.dumps(chunk) for chunk in payload) + "\n").encode("utf-8")
        else:
            body = json.dumps(payload).encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        data = json.loads(self.rfile.read(length) or b"{}")

        if self.path.endswith("/error"):
            body = json.dumps({"error": "model not found"}).encode("utf-8")
            self.send_response(500)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        if self.path.endswith("/generate"):
            if data.get("stream"):
                self._respond(
                    [
                        {"response": "Hel", "done": False},
                        {"response": "lo", "done": False},
                        {
                            "response": "",
                            "done": True,
                            "total_duration": 2_000_000_000,
                            "eval_count": 4,
                            "prompt_eval_count": 2,
                            "context": [1, 2, 3],
                        },
                    ],
                    ndjson=True,
                )
            else:
                self._respond(
                    {
                        "response": "hello world",
                        "model": data.get("model", ""),
                        "done": True,
                        "eval_count": 4,
                    }
                )
            return

        if self.path.endswith("/chat"):
            if data.get("stream"):
                self._respond(
                    [
                        {"message": {"role": "assistant", "content": "Hi"}, "done": False},
                        {"message": {"role": "assistant", "content": " there"}, "done": True},
                    ],
                    ndjson=True,
                )
            else:
                self._respond(
                    {
                        "message": {"role": "assistant", "content": "hi there"},
                        "model": data.get("model", ""),
                        "done": True,
                        "eval_count": 3,
                    }
                )
            return

        self._respond({})

    def do_GET(self):
        if self.path.endswith("/tags"):
            self._respond({"models": [{"name": "stub:latest", "size": 1024**3}]})
            return
        self._respond({})


@pytest.fixture()
def local_client():
    """OllamaClient pointed at a real local HTTP server stubbing the API."""
    server = ThreadingHTTPServer(("127.0.0.1", 0), _OllamaStubHandler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    config = OllamaConfig(
        host=f"http://127.0.0.1:{server.server_address[1]}",
        model="stub:latest",
        max_retries=1,
        timeout=10.0,
    )
    yield OllamaClient(config)

    server.shutdown()
    server.server_close()


class TestClientAgainstLocalServer:
    def test_is_available_and_list_models(self, local_client):
        assert local_client.is_available() is True
        models = local_client.list_models()
        assert [m.name for m in models] == ["stub:latest"]
        assert models[0].size_gb == pytest.approx(1.0)

    def test_generate_non_streaming(self, local_client):
        response = local_client.generate("Say hello.")

        assert response.text == "hello world"
        assert response.model == "stub:latest"
        assert response.done is True
        assert response.eval_count == 4

    def test_generate_streaming_with_callback(self, local_client):
        chunks: list[str] = []

        response = local_client.generate("Say hello.", stream=True, stream_callback=chunks.append)

        assert response.text == "Hello"
        assert "".join(chunks) == "Hello"
        assert response.eval_count == 4
        assert response.duration_seconds == pytest.approx(2.0)
        assert response.context == [1, 2, 3]

    def test_chat_non_streaming(self, local_client):
        response = local_client.chat([ChatMessage("user", "hello")])

        assert response.text == "hi there"
        assert response.message.role == "assistant"

    def test_chat_streaming(self, local_client):
        response = local_client.chat([ChatMessage("user", "hello")], stream=True)

        assert response.text == "Hi there"

    def test_http_error_raises_runtime_error_not_connection_error(self, local_client):
        """HTTPError must surface immediately as RuntimeError (not retried
        as a connection failure) — regression test for exception ordering."""
        url = f"{local_client.config.generate_url}/error"
        with pytest.raises(RuntimeError, match="Ollama API error 500"):
            local_client._request(url, {"model": "stub"})

    def test_connection_error_on_unreachable_server(self):
        config = OllamaConfig(host="http://127.0.0.1:1", max_retries=1, timeout=2.0)
        client = OllamaClient(config)

        assert client.is_available() is False
        with pytest.raises(ConnectionError, match="Failed to connect"):
            client.generate("Hi")
