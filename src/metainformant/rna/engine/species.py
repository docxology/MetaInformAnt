"""Shared species/configuration discovery helpers for RNA workflows.

Species inventories are configuration- and data-root-dependent.  Keeping the
discovery rules in one small module prevents launchers, monitors, and reports
from drifting apart or silently omitting newly added configurations.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

_NON_SPECIES_MARKERS = ("template", "test", "cross_species")

SPECIES_DIRECTORY_ALIASES: dict[str, tuple[str, ...]] = {
    "apis_mellifera": ("apis_mellifera_all", "amellifera"),
    "pogonomyrmex_barbatus": ("pbarbatus",),
}
"""Historical data-root directory names accepted for a canonical species (canonical name first in lookups)."""

_WORK_EVIDENCE_MARKERS = ("metadata", "index", "quant", "merge")


@dataclass(frozen=True)
class ConfiguredSpecies:
    """One runnable species configuration: its slug, file and the identity it declares."""

    slug: str
    config_name: str
    config_path: Path
    display_name: str
    taxon_id: str | None


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


def configured_species(config_dir: str | Path) -> list[ConfiguredSpecies]:
    """Return the runnable species of a configuration directory, validated against their YAML.

    The set is :func:`discover_species_config_names` (template, test and
    cross-species files excluded). Each YAML must declare exactly one
    ``species_list`` entry whose lower-cased form equals the slug encoded in
    its filename, so the filename, the data-root directory and the declared
    species can never drift apart; duplicate slugs or taxon ids are refused.
    Raises :class:`ValueError` on any inconsistency and returns ``[]`` for a
    missing directory.
    """

    import yaml

    directory = Path(config_dir).expanduser()
    records: list[ConfiguredSpecies] = []
    for config_name in discover_species_config_names(directory):
        path = directory / config_name
        slug = species_name_from_config(config_name)
        try:
            config = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError) as error:
            raise ValueError(f"{path} is not a readable YAML configuration: {error}") from None
        if not isinstance(config, dict):
            raise ValueError(f"{path} must contain a mapping")
        declared = config.get("species_list")
        if not isinstance(declared, list) or len(declared) != 1 or not str(declared[0]).strip():
            raise ValueError(f"{path} must define exactly one species_list entry")
        display_name = str(declared[0]).strip()
        if display_name.lower() != slug:
            raise ValueError(f"{path} declares species_list {display_name!r}, which differs from its filename slug")
        taxon = config.get("taxon_id")
        records.append(
            ConfiguredSpecies(
                slug=slug,
                config_name=config_name,
                config_path=path,
                display_name=display_name,
                taxon_id=str(taxon) if taxon is not None else None,
            )
        )
    slugs = [record.slug for record in records]
    if len(set(slugs)) != len(slugs):
        raise ValueError(f"duplicate species slugs in {directory}")
    taxa = [record.taxon_id for record in records if record.taxon_id is not None]
    if len(set(taxa)) != len(taxa):
        raise ValueError(f"duplicate taxon_id values in {directory}")
    return records


def species_directory_candidates(species: str) -> tuple[str, ...]:
    """The canonical directory name of ``species`` followed by its accepted historical aliases."""

    return (species, *(alias for alias in SPECIES_DIRECTORY_ALIASES.get(species, ()) if alias != species))


def resolve_species_work_dir(data_root: str | Path, species: str) -> Path:
    """Return the active ``<data_root>/<dir>/work`` directory for ``species``.

    The canonical directory name is tried first, then the historical aliases
    (:data:`SPECIES_DIRECTORY_ALIASES`). A work tree that already holds
    workflow evidence (``metadata``, ``index``, ``quant`` or ``merge``) wins
    over a merely existing directory, so an empty canonical directory left by
    an interrupted preparation never hides the populated alias. Without any
    existing directory the canonical path is returned (it is not created).
    """

    root = Path(data_root)
    candidates = species_directory_candidates(species)
    for name in candidates:
        work_dir = root / name / "work"
        if work_dir.is_dir() and any((work_dir / marker).exists() for marker in _WORK_EVIDENCE_MARKERS):
            return work_dir
    for name in candidates:
        if (root / name).is_dir():
            return root / name / "work"
    return root / species / "work"
