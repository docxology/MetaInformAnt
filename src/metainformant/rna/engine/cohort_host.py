"""Prepare a fresh data root to quantify frozen-cohort species on a new host.

Cohort-mode discovery (``run_all_species.py --cohort``) refuses a species until
``<data-root>/<species>/work/reference/reference_aliases.json`` binds the
cohort's kallisto index bytes and the configured index directory holds those
exact bytes (:func:`metainformant.rna.engine.cohort.verify_species_reference`).
It never builds or downloads an index, and a rebuilt index is not a substitute:
``kallisto index`` output is not byte-reproducible, and the cohort pins bytes.

For each requested species this module, idempotently:

1. checks the species YAML hashes to the cohort configuration digest;
2. places the cohort's own index file (``<cohort>/inputs/<species>/work/index``)
   in the index directory the configuration resolves to under the data root,
   copying through a re-hashed temporary file; an existing file must already
   hold the cohort bytes and is never replaced;
3. writes the reference manifest with the producer's own
   ``_ensure_reference_alias_indexes`` (metadata targets from the cohort's
   selected metadata, aliases from the configuration);
4. re-reads everything with ``verify_species_reference``.

The default is a read-only plan; ``--execute`` writes. Imports of already
finished outputs (``python -m metainformant.rna.engine.quant_import``) need
step 3 first, because the importer requires the destination to quantify each
species against the cohort index bytes.

Exit codes: 0 every species ready (or ready to prepare, in a plan), 65 a
species cannot be prepared (existing different index bytes, unresolved
reference targets), 2 usage or cohort errors.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from pathlib import Path
from typing import Any, Iterable

import yaml

from metainformant.core.utils.hash import sha256_bytes, sha256_file
from metainformant.rna.engine.cohort import CohortError, FrozenCohort, load_frozen_cohort, verify_species_reference
from metainformant.rna.engine.reference_digest import IndexDigestCache
from metainformant.rna.engine.species import configured_data_root
from metainformant.rna.engine.streaming_orchestrator import (
    _ensure_reference_alias_indexes,
    _quant_index_candidates,
    _resolve_index_dir,
)

SCHEMA = "metainformant.rna.cohort_host.v1"
EXIT_OK = 0
EXIT_BLOCKED = 65
EXIT_USAGE = 2
TARGET_COLUMNS = ("scientific_name", "organism", "species")


def metadata_targets(cohort: FrozenCohort, species: str, config: dict[str, Any]) -> list[str]:
    """Reference target names of a species: its cohort metadata's organism names, else the config species list."""

    with cohort.metadata_path(species).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        column = next((name for name in TARGET_COLUMNS if name in (reader.fieldnames or [])), None)
        names = {str(row.get(column) or "").strip() for row in reader} if column else set()
    names.discard("")
    if names:
        return sorted(names)
    return [str(value).strip() for value in config.get("species_list", []) if str(value).strip()]


def index_destination(config: dict[str, Any], species: str, index_name: str) -> Path:
    """Where cohort discovery will look for this species' index under the current data root.

    The directory discovery resolves when one already holds an index, else the
    first configured candidate (the quant ``index_dir``).
    """

    resolved = _resolve_index_dir(config, species)
    if resolved:
        return Path(resolved) / index_name
    candidates = _quant_index_candidates(config, species)
    if not candidates:
        raise CohortError(f"{species}: the configuration resolves no kallisto index directory")
    return Path(candidates[0]) / index_name


def _copy_index(source: Path, destination: Path, expected_sha256: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    try:
        shutil.copyfile(source, temporary)
        if sha256_file(temporary) != expected_sha256:
            raise OSError(f"copy of {source} does not hash to the cohort index digest")
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        if destination.exists():
            raise FileExistsError(f"{destination} appeared while the index was copied")
        os.rename(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def prepare_species(
    cohort: FrozenCohort,
    species: str,
    config_dir: Path,
    *,
    execute: bool,
    digests: IndexDigestCache,
) -> dict[str, Any]:
    """Plan or prepare one species; returns a JSON-ready record with ``state`` and ``actions``."""

    entry = cohort.entry(species)
    config_path = config_dir / entry.config_name
    record: dict[str, Any] = {"species": species, "actions": [], "state": "ready"}
    config_bytes = config_path.read_bytes()
    if sha256_bytes(config_bytes) != entry.config_sha256:
        raise CohortError(f"{config_path} does not hash to the cohort configuration digest for {species}")
    config = yaml.safe_load(config_bytes) or {}
    data_root = configured_data_root().expanduser().resolve()

    destination = index_destination(config, species, entry.index_name)
    record["index"] = str(destination)
    if destination.is_file():
        if digests.digest(destination) != entry.index_sha256:
            record["state"] = "blocked"
            record["reason"] = f"{destination} holds different index bytes; it is never replaced"
            return record
    else:
        record["actions"].append("copy_index")
        if execute:
            _copy_index(cohort.index_path(species), destination, entry.index_sha256)

    def verify() -> Any:
        index_dir = _resolve_index_dir(config, species) or str(destination.parent)
        return verify_species_reference(cohort, data_root, species, index_digests=digests, index_dir=index_dir)

    if (execute or destination.is_file()) and verify().ok:
        return record
    record["actions"].append("write_reference_manifest")
    if not execute:
        record["state"] = "to_prepare"
        return record
    targets = metadata_targets(cohort, species, config)
    ready, missing, _ = _ensure_reference_alias_indexes(config, species, targets)
    if not ready:
        record["state"] = "blocked"
        record["reason"] = f"unresolved reference targets: {', '.join(missing)}"
        return record
    check = verify()
    if not check.ok:
        record["state"] = "blocked"
        record["reason"] = check.reason
    return record


def prepare_host(
    cohort: FrozenCohort, config_dir: Path, species: Iterable[str], *, execute: bool
) -> list[dict[str, Any]]:
    digests = IndexDigestCache()
    records = []
    for name in species:
        records.append(prepare_species(cohort, name, config_dir, execute=execute, digests=digests))
    return records


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, help="Data root to prepare (default: AMALGKIT_DATA_ROOT)")
    parser.add_argument("--config-dir", type=Path, required=True, help="Species YAML directory bound by the cohort")
    parser.add_argument("--cohort", type=Path, required=True, help="Installed frozen cohort directory (with indexes)")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless COHORT.json has this sha256")
    parser.add_argument("--species", action="append", help="Species to prepare (default: every cohort species)")
    parser.add_argument("--execute", action="store_true", help="Copy indexes and write reference manifests")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the host preparation CLI; ``--data-root`` is exported only for the duration of the call."""

    args = _parser().parse_args(argv)
    previous = os.environ.get("AMALGKIT_DATA_ROOT")
    try:
        if args.data_root:
            os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
        cohort = load_frozen_cohort(args.cohort.expanduser(), expected_definition_sha256=args.cohort_definition_sha256)
        requested = list(dict.fromkeys(args.species or cohort.species_names))
        unknown = sorted(set(requested) - set(cohort.species_names))
        if unknown:
            print(f"unknown species: {', '.join(unknown)}")
            return EXIT_USAGE
        records = prepare_host(cohort, args.config_dir.expanduser().resolve(), requested, execute=args.execute)
    except CohortError as error:
        print(f"cohort error: {error}")
        return EXIT_USAGE
    finally:
        if previous is None:
            os.environ.pop("AMALGKIT_DATA_ROOT", None)
        else:
            os.environ["AMALGKIT_DATA_ROOT"] = previous
    report = {
        "schema": SCHEMA,
        "mode": "execute" if args.execute else "plan",
        "cohort_id": cohort.cohort_id,
        "species": records,
        "states": {
            state: sum(1 for r in records if r["state"] == state) for state in sorted({r["state"] for r in records})
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return EXIT_BLOCKED if any(record["state"] == "blocked" for record in records) else EXIT_OK
