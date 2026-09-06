"""Tests for core logging utilities."""

from __future__ import annotations

import logging
import os
from contextlib import contextmanager
from typing import Iterator

from metainformant.core.utils import logging as core_logging


@contextmanager
def _temporary_env(name: str, value: str | None) -> Iterator[None]:
    """Temporarily set or unset an environment variable."""
    original = os.environ.get(name)
    if value is None:
        os.environ.pop(name, None)
    else:
        os.environ[name] = value
    try:
        yield
    finally:
        if original is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = original


def test_get_logger_returns_configured_logger() -> None:
    """Test that get_logger returns a configured logger."""
    logger = core_logging.get_logger("metainformant.test")
    assert isinstance(logger, logging.Logger)
    assert logger.handlers, "Expected at least one handler"


def test_get_logger_with_level() -> None:
    """Test get_logger_with_level function."""
    logger = core_logging.get_logger_with_level("metainformant.test.debug", "DEBUG")
    assert isinstance(logger, logging.Logger)
    assert logger.level == logging.DEBUG

    logger_info = core_logging.get_logger_with_level("metainformant.test.info", "INFO")
    assert logger_info.level == logging.INFO


def test_get_logger_with_level_from_env() -> None:
    """Test that get_logger_with_level reads from CORE_LOG_LEVEL env var."""
    with _temporary_env("CORE_LOG_LEVEL", "WARNING"):
        logger = core_logging.get_logger_with_level("metainformant.test.env")
    assert logger.level == logging.WARNING


def test_configure_logging_from_env() -> None:
    """Test configure_logging_from_env function."""
    with _temporary_env("CORE_LOG_LEVEL", "ERROR"):
        core_logging.configure_logging_from_env()
        assert logging.root.level == logging.ERROR

    with _temporary_env("CORE_LOG_LEVEL", None):
        core_logging.configure_logging_from_env(default_level="INFO")
        assert logging.root.level == logging.INFO


def test_log_with_metadata_basic() -> None:
    """Test log_with_metadata with basic usage."""
    import io

    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)

    logger = logging.getLogger("test_metadata")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    core_logging.log_with_metadata(logger, "Test message", {"key1": "value1", "key2": 42})

    log_output = log_stream.getvalue()
    assert "Test message" in log_output
    assert "key1" in log_output or "value1" in log_output


def test_log_with_metadata_structured() -> None:
    """Test log_with_metadata with structured=True."""
    import io

    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)

    logger = logging.getLogger("test_structured")
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

    core_logging.log_with_metadata(logger, "Structured message", {"batch": 1}, structured=True)

    log_output = log_stream.getvalue()
    assert "Structured message" in log_output
    assert "METADATA:" in log_output


def test_log_with_metadata_level() -> None:
    """Test log_with_metadata with different log levels."""
    import io

    log_stream = io.StringIO()
    handler = logging.StreamHandler(log_stream)

    logger = logging.getLogger("test_level")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)

    core_logging.log_with_metadata(logger, "Debug message", {"data": "test"}, level="DEBUG")

    log_output = log_stream.getvalue()
    assert "Debug message" in log_output


def test_setup_logger_closes_previous_handlers(tmp_path) -> None:
    """Re-running setup_logger must close old handlers, not leak open files."""

    log_file = tmp_path / "logs" / "run.log"

    core_logging.setup_logger("leak_check", log_file=str(log_file))
    logger = logging.getLogger("leak_check")
    first_file_handlers = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    assert first_file_handlers, "expected a FileHandler on first setup"

    core_logging.setup_logger("leak_check", log_file=str(log_file))

    # FileHandler.close() closes the stream and drops it (stream becomes None).
    assert first_file_handlers[0].stream is None
    active = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    assert len(active) == 1
