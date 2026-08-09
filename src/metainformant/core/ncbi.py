"""Explicit contact policy for NCBI-facing operations."""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass


class NCBIContactError(ValueError):
    """Raised when an NCBI operation lacks an explicit contact policy."""


@dataclass(frozen=True)
class NCBIContact:
    """Resolved NCBI contact without storing a fabricated identity."""

    email: str | None
    mode: str


def resolve_ncbi_contact(
    email: str | None = None,
    *,
    allow_anonymous: bool = False,
) -> NCBIContact:
    """Resolve explicit or opt-in anonymous NCBI contact settings.

    Explicit ``email`` takes precedence over ``NCBI_EMAIL``. If neither is
    supplied, callers must opt into anonymous mode. Anonymous mode never
    invents, sets, or persists an email address.
    """

    resolved = (email or os.environ.get("NCBI_EMAIL", "")).strip()
    if resolved:
        return NCBIContact(email=resolved, mode="explicit")
    if not allow_anonymous:
        raise NCBIContactError(
            "NCBI contact is required; pass email=... or set NCBI_EMAIL. "
            "Use allow_anonymous=True only for intentionally anonymous requests."
        )
    warnings.warn(
        "Proceeding with an anonymous NCBI request; no email identity will be set or persisted.",
        RuntimeWarning,
        stacklevel=2,
    )
    return NCBIContact(email=None, mode="anonymous")


__all__ = ["NCBIContact", "NCBIContactError", "resolve_ncbi_contact"]
