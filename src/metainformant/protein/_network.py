"""Shared protein network-client configuration helpers."""

from __future__ import annotations

from metainformant.core.utils import config, logging

logger = logging.get_logger(__name__)


def get_protein_api_timeout(default: float = 30.0) -> float:
    """Return the protein API request timeout from ``PROT_TIMEOUT`` or a default."""
    raw_timeout = config.get_env_or_default("PROT_TIMEOUT", str(default))
    try:
        timeout = float(raw_timeout)
    except ValueError:
        logger.warning("Invalid PROT_TIMEOUT=%r; using default %.1fs", raw_timeout, default)
        return float(default)

    if timeout <= 0:
        logger.warning("PROT_TIMEOUT must be positive; using default %.1fs", default)
        return float(default)

    return timeout
