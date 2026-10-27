"""Regression tests for ``max_retries`` in the Ollama client, exercised without a server.

``OllamaConfig`` accepted ``max_retries=0`` (and negative values). ``OllamaClient._request`` loops
``range(max_retries)``, so with 0 it made no attempt at all and fell off the end returning ``None``; ``generate()``
then crashed with ``AttributeError: 'NoneType' object has no attribute 'get'`` instead of saying what was wrong.

``max_retries`` is the maximum number of attempts per request (the first one included, as the retry loop always
treated it), so the config now rejects values below 1, and ``_request`` refuses a value that was lowered after
construction before it touches the network. The retry loop itself is checked against ``file://`` URLs, which go
through the same ``urllib`` code path as HTTP but need no server: a missing file is a ``URLError`` (a failed
connection), an existing JSON file is a reply. No test doubles and no network.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import pytest

from metainformant.ml.llm.ollama.client import OllamaClient
from metainformant.ml.llm.ollama.config import OllamaConfig


@pytest.mark.parametrize("bad", [0, -1, -10])
def test_config_rejects_max_retries_below_one(bad: int) -> None:
    with pytest.raises(ValueError, match="max_retries"):
        OllamaConfig(max_retries=bad)


@pytest.mark.parametrize("good", [1, 2, 3, 10])
def test_config_accepts_one_or_more_attempts(good: int) -> None:
    config = OllamaConfig(max_retries=good)

    assert config.max_retries == good
    assert config.with_model("other:latest").max_retries == good


def test_default_config_is_unchanged() -> None:
    assert OllamaConfig().max_retries == 3


def test_documentation_states_that_max_retries_counts_attempts_including_the_first() -> None:
    documented = " ".join((OllamaConfig.__doc__ or "").split("max_retries:")[1].split("stream:")[0].split())

    assert "number of attempts" in documented and "first one included" in documented


def test_request_with_max_retries_lowered_after_construction_fails_clearly_not_with_none() -> None:
    client = OllamaClient(OllamaConfig())
    client.config.max_retries = 0  # the dataclass is mutable: this bypasses __post_init__

    with pytest.raises(ValueError, match="max_retries must be >= 1"):
        client._request("file:///no/such/ollama/endpoint", {"model": "x"})
    with pytest.raises(ValueError, match="max_retries must be >= 1"):
        client.generate("hello")  # used to crash with AttributeError: 'NoneType' object has no attribute 'get'


def _retry_messages(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records if "retrying" in record.getMessage()]


def test_a_single_attempt_fails_with_a_connection_error_and_never_retries(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    client = OllamaClient(OllamaConfig(max_retries=1))
    missing = (tmp_path / "missing.json").as_uri()

    with caplog.at_level(logging.WARNING):
        with pytest.raises(ConnectionError, match="Failed to connect"):
            client._request(missing, {"model": "x"})

    assert _retry_messages(caplog) == []


def test_two_attempts_retry_exactly_once_after_a_pause(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    client = OllamaClient(OllamaConfig(max_retries=2))
    missing = (tmp_path / "missing.json").as_uri()

    started = time.monotonic()
    with caplog.at_level(logging.WARNING):
        with pytest.raises(ConnectionError, match="Failed to connect"):
            client._request(missing, {"model": "x"})

    retries = _retry_messages(caplog)
    assert len(retries) == 1 and "retrying in 1s" in retries[0]
    assert time.monotonic() - started >= 0.9, "the second attempt follows a 1 s back-off"


def test_a_reply_on_the_first_attempt_is_returned(tmp_path: Path) -> None:
    reply = tmp_path / "reply.json"
    reply.write_text(json.dumps({"response": "hello", "model": "stub", "done": True}))
    client = OllamaClient(OllamaConfig(max_retries=1))

    assert client._request(reply.as_uri(), {"model": "x"}) == {"response": "hello", "model": "stub", "done": True}
