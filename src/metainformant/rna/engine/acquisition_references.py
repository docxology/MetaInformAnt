"""Prepare explicit reference aliases for frozen local or AWS acquisition."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.streaming_orchestrator import _ensure_reference_alias_indexes, _species_work_dir

def prepare_reference_inputs(
    *, data_root: Path, config_dir: Path, tasks: list[dict[str, Any]],
    snapshot: dict[str, Any] | None = None
) -> list[dict[str, Any]]:
    """Create and record explicit reference aliases before task execution.

    The cloud worker receives a static manifest and therefore does not run the
    full local discovery phase.  Reference aliases must still be materialized
    before Amalgkit quantification; otherwise a metadata target such as a
    subspecies name can fail with a misleading missing-index error even when
    its declared species-level index is present in the input bundle.
    """

    import pandas as pd
    import yaml

    species_configs = {str(task["species"]): str(task["config_name"]) for task in tasks}
    records: list[dict[str, Any]] = []
    for species, config_name in sorted(species_configs.items()):
        config_path = config_dir / config_name
        if not config_path.is_file():
            raise FileNotFoundError(f"config missing for {species}: {config_path}")
        config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        work_dir = _species_work_dir(species)
        metadata_path = work_dir / "metadata" / "metadata_selected.tsv"
        if not metadata_path.is_file():
            raise FileNotFoundError(
                f"selected metadata missing for {species}: {metadata_path}"
            )
        if snapshot is not None:
            hashes = {record["sha256"] for record in snapshot.get("input_files", [])
                      if Path(record["path"]).parts[-4:] == (species, "work", "metadata", "metadata_selected.tsv")}
            if len(hashes) != 1 or sha256_file(metadata_path) not in hashes:
                raise ValueError(f"worker metadata differs from frozen acquisition envelope: {metadata_path}")
        metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)
        target_column = next(
            (
                column
                for column in ("scientific_name", "organism", "species")
                if column in metadata.columns
            ),
            None,
        )
        if target_column is None:
            target_names = [
                str(value).strip()
                for value in config.get("species_list", [])
                if str(value).strip()
            ]
        else:
            target_names = sorted(
                str(value).strip()
                for value in metadata[target_column].dropna().unique()
                if str(value).strip()
            )
        if not target_names:
            raise ValueError(f"no reference target names available for {species}")
        ready, missing, manifest = _ensure_reference_alias_indexes(
            config, species, target_names
        )
        if not ready or manifest is None:
            raise RuntimeError(
                f"reference preflight failed for {species}: {', '.join(missing) or 'unknown error'}"
            )
        expected_indexes = {task["reference_index_sha256"] for task in tasks
                            if task["species"] == species and task.get("reference_index_sha256")}
        if len(expected_indexes) > 1:
            raise ValueError("species tasks disagree about the frozen reference index")
        if expected_indexes:
            aliases = json.loads(manifest.read_text())["aliases"]
            for alias in aliases:
                paths = [alias["index"]] if "index" in alias else alias["index_aliases"]
                for value in paths:
                    index = Path(value)
                    if not index.resolve().is_relative_to(data_root.resolve()) or sha256_file(index) not in expected_indexes:
                        raise ValueError(f"worker reference differs from frozen acquisition envelope: {index}")
        records.append(
            {
                "species": species,
                "targets": target_names,
                "manifest": str(manifest.relative_to(data_root)),
                "manifest_sha256": sha256_file(manifest),
            }
        )
    return records
