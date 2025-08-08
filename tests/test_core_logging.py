from __future__ import annotations

import logging
from metainformant.core import logging as core_logging


def test_get_logger_returns_configured_logger() -> None:
    logger = core_logging.get_logger("metainformant.test")
    assert isinstance(logger, logging.Logger)
    assert logger.handlers, "Expected at least one handler"


