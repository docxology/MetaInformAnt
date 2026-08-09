"""Tests for explicit NCBI contact resolution."""

from __future__ import annotations

import pytest

from metainformant.core.ncbi import NCBIContactError, resolve_ncbi_contact


def test_explicit_contact_takes_precedence(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NCBI_EMAIL", "environment@example.org")
    contact = resolve_ncbi_contact("explicit@example.org")
    assert contact.email == "explicit@example.org"
    assert contact.mode == "explicit"


def test_environment_contact_is_supported(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("NCBI_EMAIL", "environment@example.org")
    contact = resolve_ncbi_contact()
    assert contact.email == "environment@example.org"


def test_missing_contact_requires_anonymous_opt_in(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.raises(NCBIContactError):
        resolve_ncbi_contact()


def test_anonymous_contact_is_unset_and_warns(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    with pytest.warns(RuntimeWarning, match="anonymous NCBI"):
        contact = resolve_ncbi_contact(allow_anonymous=True)
    assert contact.email is None
    assert contact.mode == "anonymous"
