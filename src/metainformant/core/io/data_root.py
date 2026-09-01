"""Canonical, deterministic data-root and Amalgkit-prefix path resolution.

The campaign data root may live on a separate volume from the repository and is
selected via the ``AMALGKIT_DATA_ROOT`` environment variable.  Several modules
previously re-implemented this resolution with slightly different semantics
(some resolved to absolute, some did not, and the prefix remap lived only in
the RNA workflow layer).  This module is the single shared implementation so
that launchers, monitors, reports, and the clade-replication recipe all agree.

Resolution order (deterministic, no implicit cwd dependence beyond the
documented relative-path defaults):

1. Explicit ``data_root`` argument (highest priority).
2. ``AMALGKIT_DATA_ROOT`` environment variable.
3. ``default`` fallback (typically ``output/amalgkit``).

All results are ``expanduser``-ed and resolved to absolute paths.
"""

from __future__ import annotations

import os
from pathlib import Path

AMALGKIT_DATA_ROOT_ENV = "AMALGKIT_DATA_ROOT"
AMALGKIT_PREFIX = "output/amalgkit"


def resolve_data_root(
    data_root: str | Path | None = None,
    *,
    default: str | Path = AMALGKIT_PREFIX,
    environ: dict[str, str] | None = None,
) -> Path:
    """Resolve the campaign data root to an absolute path.

    Args:
        data_root: Explicit override; wins over the environment when given.
        default: Fallback when neither explicit value nor env var is set.
        environ: Environment mapping to consult (defaults to ``os.environ``);
            injectable so callers can resolve against a non-process context.

    Returns:
        Absolute, ``expanduser``-ed path.
    """
    if data_root is not None:
        return Path(data_root).expanduser().resolve()
    env = os.environ if environ is None else environ
    value = env.get(AMALGKIT_DATA_ROOT_ENV)
    if value:
        return Path(value).expanduser().resolve()
    return Path(default).expanduser().resolve()


def remap_amalgkit_prefix(
    path: str,
    data_root: str | Path | None = None,
    *,
    environ: dict[str, str] | None = None,
) -> str:
    """Remap the repository-relative ``output/amalgkit`` prefix to the data root.

    Species YAML files use paths such as ``output/amalgkit/<species>/work``.
    When the production tree lives on a separate volume, this remaps that
    prefix onto the resolved data root while leaving repository-relative paths
    such as ``config/config_base`` untouched.  ``~/`` and ``$VAR`` expansions
    are applied first, matching the workflow-layer behavior this module
    replaces.

    Args:
        path: Candidate path string from a config value.
        data_root: Explicit data-root override (see :func:`resolve_data_root`).
        environ: Environment mapping to consult (defaults to ``os.environ``).

    Returns:
        The possibly remapped path string.
    """
    expanded = os.path.expanduser(os.path.expandvars(path))
    relative = expanded[2:] if expanded.startswith("./") else expanded
    if relative == AMALGKIT_PREFIX:
        return str(resolve_data_root(data_root, environ=environ))
    if relative.startswith(AMALGKIT_PREFIX + "/"):
        root = resolve_data_root(data_root, environ=environ)
        return str(root / relative[len(AMALGKIT_PREFIX) + 1 :])
    return expanded
