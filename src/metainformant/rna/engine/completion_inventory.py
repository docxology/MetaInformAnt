"""Freeze the configured species' public RNA-seq inventory and lock existing outputs."""

from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import io
import json
import os
import shutil
import sqlite3
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import requests
import yaml

from metainformant.rna.engine.durable_quant import (
    DirectoryStore,
    S3Store,
    lock_quantification,
)
from metainformant.rna.engine.species import (
    discover_species_config_names,
    species_name_from_config,
)

ENA_FIELDS = "run_accession,scientific_name,tax_id,library_layout,library_strategy,library_source,library_selection,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count,study_accession,sample_accession,experiment_accession,instrument_platform,first_public,last_updated"


def freeze_inventory(
    data_root: Path, config_dir: Path, output_dir: Path
) -> dict[str, Any]:
    """Append newly discovered runs to frozen metadata, leaving canonical inputs untouched."""
    if (data_root / ".full_campaign.lock").exists():
        raise RuntimeError(
            "cannot freeze inputs while the canonical producer lock exists"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "inventory.json"
    if destination.exists():
        return json.loads(destination.read_text())
    with sqlite3.connect(
        f"file:{data_root / 'pipeline_progress.db'}?mode=ro", uri=True
    ) as db:
        exclusions = {
            (s, r): reason
            for s, r, reason in db.execute(
                "SELECT species,srr_id,reason_code FROM sample_exclusions"
            )
        }
        states = {
            (s, r): state
            for s, r, state in db.execute("SELECT species,srr_id,state FROM samples")
        }
    names = discover_species_config_names(config_dir)
    if len(names) != 27:
        raise ValueError(f"expected 27 configured species, found {len(names)}")

    def species_inventory(name: str) -> dict[str, Any]:
        species = species_name_from_config(name)
        config = yaml.safe_load((config_dir / name).read_text())
        taxid = int(config["taxon_id"])
        source_work = data_root / species / "work"
        metadata_path = source_work / "metadata" / "metadata_selected.tsv"
        with metadata_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = list(reader.fieldnames or [])
            rows = list(reader)
        if "run" not in columns:
            raise ValueError(f"selected metadata has no run column: {species}")
        ids = [r["run"] for r in rows]
        if len(ids) != len(set(ids)):
            raise ValueError(f"duplicate selected run accessions: {species}")
        response = requests.get(
            "https://www.ebi.ac.uk/ena/portal/api/search",
            params={
                "result": "read_run",
                "query": f'tax_tree({taxid}) AND library_strategy="RNA-Seq"',
                "fields": ENA_FIELDS,
                "format": "tsv",
                "limit": "0",
            },
            timeout=(15, 180),
        )
        response.raise_for_status()
        ena_path = output_dir / f"ena_{species}.tsv"
        ena_path.write_text(response.text)
        ena = list(csv.DictReader(io.StringIO(response.text), delimiter="\t"))
        if not ena or "run_accession" not in ena[0]:
            raise ValueError(f"ENA returned no usable inventory for {species}")
        by_accession = {r["run_accession"]: r for r in ena}
        if len(by_accession) != len(ena):
            raise ValueError(f"ENA returned duplicated accessions for {species}")
        additions = sorted(set(by_accession) - set(ids))
        for accession in additions:
            r = by_accession[accession]
            base_count, read_count = (
                int(r.get("base_count") or 0),
                int(r.get("read_count") or 0),
            )
            row = dict.fromkeys(columns, "")
            row.update(
                run=accession,
                scientific_name=r["scientific_name"],
                lib_layout=r["library_layout"].lower(),
                lib_strategy="RNA-Seq",
                lib_source=r["library_source"],
                lib_selection=r["library_selection"],
                total_bases=str(base_count),
                total_spots=str(read_count),
                size=str(sum(int(v) for v in r.get("fastq_bytes", "").split(";") if v)),
                bioproject=r["study_accession"],
                biosample=r["sample_accession"],
                experiment=r["experiment_accession"],
                platform=r["instrument_platform"],
                taxid=r["tax_id"],
                published_date=r["first_public"],
                is_sampled="yes",
                is_qualified="yes",
                exclusion="no",
                sample_group="unannotated",
                tissue="unannotated",
            )
            rows.append(row)
        work = output_dir / "inputs" / species / "work"
        (work / "metadata").mkdir(parents=True, exist_ok=True)
        selected = work / "metadata" / "metadata_selected.tsv"
        with selected.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        index_dir = source_work / "index"
        if not index_dir.is_dir():
            index_dir = data_root / species / "genome" / "index"
        indexes = sorted(
            p
            for p in index_dir.glob("*.idx")
            if p.is_file() and not p.name.startswith("._")
        )
        expected_stem = species.casefold()
        exact = [p for p in indexes if p.stem.casefold() == expected_stem]
        if len(exact) != 1:
            raise ValueError(
                f"no unique exact species index for {species}: {[p.name for p in indexes]}"
            )
        index = work / "index" / exact[0].name
        index.parent.mkdir(parents=True, exist_ok=True)
        if not index.exists():
            shutil.copyfile(exact[0], index)
        index_hash = hashlib.sha256(index.read_bytes()).hexdigest()
        tasks, excluded = [], []
        for number, row in enumerate(rows, start=1):
            accession = row["run"]
            if (species, accession) in exclusions:
                excluded.append(
                    {"accession": accession, "reason": exclusions[(species, accession)]}
                )
                continue
            remote = by_accession.get(accession, {})
            size = sum(int(v) for v in remote.get("fastq_bytes", "").split(";") if v)
            if not size:
                size = int(float(row.get("size") or 0))
            tasks.append(
                {
                    "schema": "metainformant.hymenoptera.gcp_task_manifest.v1",
                    "task_id": f"{species}/{accession}",
                    "species": species,
                    "accession": accession,
                    "config_name": name,
                    "batch_index": number,
                    "expected_paired": str(row.get("lib_layout", "")).lower()
                    == "paired",
                    "total_bases": float(row.get("total_bases") or 0),
                    "fastq_bytes": size,
                    "existing_state": states.get((species, accession), "new"),
                    "reference_index_sha256": index_hash,
                }
            )
        return {
            "species": species,
            "config_name": name,
            "taxid": taxid,
            "config_sha256": hashlib.sha256(
                (config_dir / name).read_bytes()
            ).hexdigest(),
            "metadata_sha256": hashlib.sha256(selected.read_bytes()).hexdigest(),
            "index_name": index.name,
            "index_sha256": index_hash,
            "previous_selected": len(ids),
            "new_runs": len(additions),
            "ena_runs": len(ena),
            "tasks": tasks,
            "excluded": excluded,
        }

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        species_rows = list(pool.map(species_inventory, names))
    task_ids = [t["accession"] for s in species_rows for t in s["tasks"]]
    if len(task_ids) != len(set(task_ids)):
        raise ValueError("accession assigned to more than one configured species")
    inventory = {
        "schema": "metainformant.hymenoptera.completion_inventory.v1",
        "frozen_at": datetime.now(UTC).isoformat(),
        "species": species_rows,
        "species_count": len(species_rows),
        "task_count": len(task_ids),
        "new_runs": sum(s["new_runs"] for s in species_rows),
        "excluded_count": sum(len(s["excluded"]) for s in species_rows),
    }
    temporary = destination.with_suffix(".tmp")
    temporary.write_text(json.dumps(inventory, indent=2))
    temporary.rename(destination)
    return inventory


def seal_existing_outputs(
    inventory: dict[str, Any],
    data_root: Path,
    config_dir: Path,
    output_dir: Path,
    *,
    bucket: str,
    cohort: str,
    profile: str,
    region: str,
    workers: int = 12,
) -> dict[str, Any]:
    """Validate and back up local outputs; leave canonical DB states untouched."""
    output_dir.mkdir(parents=True, exist_ok=True)
    local = DirectoryStore(output_dir / "locked")
    remote = S3Store(bucket, "locked-quant-v1", profile=profile, region=region)
    candidates = [
        (s, t)
        for s in inventory["species"]
        for t in s["tasks"]
        if (data_root / s["species"] / "work" / "quant" / t["accession"]).is_dir()
    ]
    results = []

    def seal(item: tuple[dict[str, Any], dict[str, Any]]) -> dict[str, Any]:
        s, task = item
        species, accession = s["species"], task["accession"]
        sample = data_root / species / "work" / "quant" / accession
        try:
            receipt = lock_quantification(
                local,
                cohort,
                species,
                accession,
                sample,
                expected_config_sha256=s["config_sha256"],
            )
            lock_quantification(
                remote,
                cohort,
                species,
                accession,
                sample,
                expected_config_sha256=s["config_sha256"],
            )
        except (OSError, ValueError, KeyError) as exc:
            return {"task_id": task["task_id"], "status": "unlocked", "error": str(exc)}
        return {
            "task_id": task["task_id"],
            "status": "locked",
            "contract_id": receipt["contract_id"],
        }

    with (
        concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor,
        (output_dir / "local_seal_results.jsonl").open("a") as handle,
    ):
        for result in executor.map(seal, candidates):
            results.append(result)
            handle.write(json.dumps(result, sort_keys=True) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
            if len(results) % 100 == 0:
                print(
                    f"Sealed {len(results)}/{len(candidates)} local candidates; locked={sum(r['status'] == 'locked' for r in results)}",
                    flush=True,
                )
    summary = {
        "candidates": len(candidates),
        "locked": sum(r["status"] == "locked" for r in results),
        "unlocked": [r for r in results if r["status"] != "locked"],
    }
    (output_dir / "local_seal_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
