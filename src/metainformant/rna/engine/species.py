"""Shared species/configuration discovery helpers for RNA workflows.

Species inventories are configuration- and data-root-dependent.  Keeping the
discovery rules in one small module prevents launchers, monitors, and reports
from drifting apart or silently omitting newly added configurations.
"""

from __future__ import annotations


from pathlib import Path

_NON_SPECIES_MARKERS = ("template", "test", "cross_species")


def species_name_from_config(config_name: str | Path) -> str:
    """Return the species identifier encoded in an Amalgkit config filename."""

    return Path(config_name).stem.removeprefix("amalgkit_")


def discover_species_config_names(config_dir: str | Path) -> list[str]:
    """Return sorted species config filenames, excluding non-species configs.

    Only ``amalgkit_*.yaml`` files are considered.  Template, test, and
    cross-species orchestration files are intentionally excluded because they
    are not runnable single-species workflows.
    """

    directory = Path(config_dir).expanduser()
    if not directory.is_dir():
        return []

    names: list[str] = []
    for config_path in sorted(directory.glob("amalgkit_*.yaml")):
        species = species_name_from_config(config_path)
        if any(marker in species.lower() for marker in _NON_SPECIES_MARKERS):
            continue
        names.append(config_path.name)
    return names


def discover_species_names(config_dir: str | Path) -> list[str]:
    """Return sorted species identifiers from a configuration directory."""

    return [species_name_from_config(name) for name in discover_species_config_names(config_dir)]


def configured_data_root(default: str | Path = "output/amalgkit") -> Path:
    """Return the configured data root, honoring ``AMALGKIT_DATA_ROOT``.

    Delegates to the canonical :func:`metainformant.core.io.data_root.resolve_data_root`.
    """
    from metainformant.core.io.data_root import resolve_data_root

    return resolve_data_root(None, default=default)
