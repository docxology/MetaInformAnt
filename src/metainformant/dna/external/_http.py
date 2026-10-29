"""Shared HTTP plumbing for the NCBI clients in this package.

The E-utilities clients (:mod:`.entrez`, :mod:`.ncbi`) and the Datasets v2
client (:mod:`.ncbi_datasets`) used to catch ``requests.RequestException`` and
return an empty list, dict or ``None``.  An outage, a rate limit and a query
with no hits then looked identical to the caller.  They now share this module:
a failed request raises :class:`NCBIRequestError`, and the message never
contains the ``api_key`` or ``email`` values that E-utilities requires in the
query string.
"""

from __future__ import annotations

import json
import os
import re
from typing import Any, Dict, Mapping, Optional, Type

import requests

__all__ = [
    "NCBIRequestError",
    "checked_get",
    "decode_json_object",
    "describe_http_error",
    "describe_transport_error",
    "get_json_object",
    "get_text",
    "redact_secrets",
    "resolve_api_key",
]

_SECRET_QUERY_PARAM = re.compile(r"(?P<name>\b(?:api[_-]?key|email)=)[^&\s'\"<>)]+", re.IGNORECASE)
_MAX_DETAIL_CHARS = 300


class NCBIRequestError(requests.RequestException):
    """An NCBI request failed or returned an unusable response.

    Subclasses :class:`requests.RequestException`, so callers that already
    handle ``requests`` failures keep working.  ``status_code`` is set when the
    server answered with an HTTP error status.
    """

    def __init__(self, message: str, *, status_code: Optional[int] = None) -> None:
        super().__init__(message)
        self.status_code = status_code


def resolve_api_key(api_key: Optional[str] = None) -> Optional[str]:
    """The NCBI API key to use: an explicit ``api_key``, else ``NCBI_API_KEY``, else ``None``.

    The documentation has long told users to set ``NCBI_API_KEY``; the E-utilities clients ignored it.
    """
    return (api_key or os.environ.get("NCBI_API_KEY", "")).strip() or None


def redact_secrets(text: str) -> str:
    """Mask ``api_key=`` and ``email=`` values embedded in URLs or error text."""
    return _SECRET_QUERY_PARAM.sub(lambda match: match.group("name") + "***", text)


def _squash(text: str) -> str:
    """Collapse whitespace, redact secrets and bound the length of ``text``."""
    flattened = " ".join(redact_secrets(text).split())
    if len(flattened) > _MAX_DETAIL_CHARS:
        return flattened[: _MAX_DETAIL_CHARS - 3] + "..."
    return flattened


def _error_detail(body: str) -> str:
    """Return the human-readable part of an error body (JSON ``message`` or plain text)."""
    text = body.strip()
    if not text:
        return ""
    try:
        payload = json.loads(text)
    except ValueError:
        payload = None
    if isinstance(payload, Mapping):
        for key in ("message", "error", "detail", "reason"):
            value = payload.get(key)
            if isinstance(value, str) and value.strip():
                return _squash(value)
    return _squash(text)


def describe_http_error(status_code: int, reason: Optional[str], body: str, operation: str) -> str:
    """Describe an HTTP error response without echoing the request URL.

    NCBI documents Datasets errors as ``rpcStatus`` objects (``code``,
    ``message``, ``details``) served as ``text/plain``; E-utilities reports
    errors as JSON or plain text.  Both shapes are handled.
    """
    head = f"{operation} failed: HTTP {status_code}"
    if reason:
        head += f" {reason}"
    detail = _error_detail(body)
    return f"{head}: {detail}" if detail else head


def describe_transport_error(exc: BaseException) -> str:
    """Describe a ``requests`` failure (DNS, connection, timeout) without leaking credentials."""
    return f"{type(exc).__name__}: {_squash(str(exc))}"


def checked_get(
    session: requests.Session,
    url: str,
    *,
    operation: str,
    timeout: float,
    params: Optional[Mapping[str, Any]] = None,
    headers: Optional[Mapping[str, str]] = None,
    stream: bool = False,
    error_class: Type[NCBIRequestError] = NCBIRequestError,
) -> requests.Response:
    """GET ``url`` and raise ``error_class`` unless the HTTP status is below 400.

    The original ``requests`` exception is deliberately not chained: its text
    contains the full URL, including any ``api_key`` query parameter.  With
    ``stream=True`` the caller must close the returned response.
    """
    try:
        response = session.get(url, params=params, timeout=timeout, headers=headers, stream=stream)
    except requests.RequestException as exc:
        raise error_class(f"{operation} failed: {describe_transport_error(exc)}") from None

    status_code = int(response.status_code)
    if status_code >= 400:
        try:
            body = response.text
        except requests.RequestException:
            body = ""
        finally:
            response.close()
        raise error_class(
            describe_http_error(status_code, response.reason, body, operation),
            status_code=status_code,
        )
    return response


def decode_json_object(
    response: requests.Response,
    operation: str,
    *,
    error_class: Type[NCBIRequestError] = NCBIRequestError,
) -> Dict[str, Any]:
    """Decode a JSON object from ``response`` or raise ``error_class``."""
    try:
        payload = response.json()
    except ValueError as exc:  # requests' JSONDecodeError is a ValueError
        raise error_class(f"{operation} failed: response is not valid JSON ({type(exc).__name__})") from None
    if not isinstance(payload, dict):
        raise error_class(f"{operation} failed: expected a JSON object, got {type(payload).__name__}")
    return payload


def get_json_object(
    session: requests.Session,
    url: str,
    *,
    operation: str,
    timeout: float,
    params: Optional[Mapping[str, Any]] = None,
    headers: Optional[Mapping[str, str]] = None,
    error_class: Type[NCBIRequestError] = NCBIRequestError,
) -> Dict[str, Any]:
    """GET ``url`` and return its JSON object; every failure raises ``error_class``."""
    with checked_get(
        session, url, operation=operation, timeout=timeout, params=params, headers=headers, error_class=error_class
    ) as response:
        return decode_json_object(response, operation, error_class=error_class)


def get_text(
    session: requests.Session,
    url: str,
    *,
    operation: str,
    timeout: float,
    params: Optional[Mapping[str, Any]] = None,
    headers: Optional[Mapping[str, str]] = None,
    error_class: Type[NCBIRequestError] = NCBIRequestError,
) -> str:
    """GET ``url`` and return the body text; every failure raises ``error_class``."""
    with checked_get(
        session, url, operation=operation, timeout=timeout, params=params, headers=headers, error_class=error_class
    ) as response:
        text: str = response.text
        return text
