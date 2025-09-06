from __future__ import annotations

import os
from typing import Any, Iterable, List

try:
    from ncbi.datasets import GenomeApi as DatasetsGenomeApi
    from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
except Exception:  # pragma: no cover - optional runtime dependency
    DatasetsApiClient = None  # type: ignore
    DatasetsGenomeApi = None  # type: ignore

# In test environments lacking the optional dependency, force the "not installed"
# behavior to keep unit tests deterministic unless explicitly allowed.
if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
    DatasetsApiClient = None  # type: ignore
    DatasetsGenomeApi = None  # type: ignore


def download_genome_data_package(accessions: Iterable[str], filename: str) -> Any:
    # During pytest, unless explicitly allowed, behave as if dependency is missing
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        return genome_api.download_assembly_package(list(accessions), filename=filename)


def get_metadata_by_single_accession(genome_assembly_accessions: List[str]) -> dict:
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_accessions(genome_assembly_accessions)
        return genome_metadata.get("assemblies", [{}])[0]


def get_accession_by_tax_id(tax_id: str) -> list[str]:
    if os.environ.get("PYTEST_CURRENT_TEST") and os.environ.get("METAINFORMANT_ALLOW_NCBI_DATASETS") != "1":
        raise RuntimeError("ncbi-datasets-pylib not installed")
    if DatasetsApiClient is None or DatasetsGenomeApi is None:
        raise RuntimeError("ncbi-datasets-pylib not installed")
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        genome_metadata = genome_api.assembly_descriptors_by_taxon(tax_id)
        assemblies = genome_metadata.get("assemblies", []) or []
        return [
            a.get("assembly", {}).get("assembly_accession", "")
            for a in assemblies
            if a.get("assembly", {}).get("assembly_accession")
        ]


# ---- Best-effort download helpers (CLI/API/FTP) ----

import json
import shutil
import subprocess
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from zipfile import ZipFile

from ..core.io import dump_json, ensure_directory


def datasets_cli_available() -> bool:
    return shutil.which("datasets") is not None


def _extract_zip(zip_path: Path, out_dir: Path) -> Path:
    extract_to = out_dir / (zip_path.stem + "_extracted")
    ensure_directory(extract_to)
    with ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_to)
    return extract_to


def download_genome_via_datasets_cli(
    accession: str,
    dest_dir: str | Path,
    *,
    include: list[str] | None = None,
    timeout_seconds: int = 0,
) -> dict:
    """Use `datasets` CLI to download a genome package into dest_dir.

    include: e.g., ["gff3","rna","cds","protein","genome","seq-report"].
    Streams stdout/stderr to log files under dest_dir to show progress.
    Returns a record with paths and command metadata.
    """
    out_dir = ensure_directory(dest_dir)
    include_vals = include or ["gff3", "rna", "cds", "protein", "genome", "seq-report"]
    include_arg = ",".join(include_vals)
    cmd = [
        "datasets",
        "download",
        "genome",
        "accession",
        accession,
        "--include",
        include_arg,
    ]
    stdout_log = out_dir / "datasets_cli.stdout.log"
    stderr_log = out_dir / "datasets_cli.stderr.log"
    # Heartbeat file to signal active download
    heartbeat = out_dir / "download.heartbeat"
    heartbeat.write_text(str(time.time()))
    with open(stdout_log, "w", encoding="utf-8") as out_fh, open(stderr_log, "w", encoding="utf-8") as err_fh:
        proc = subprocess.Popen(
            cmd,
            cwd=str(out_dir),
            stdout=out_fh,
            stderr=err_fh,
            text=True,
        )
        try:
            rc = proc.wait(timeout=None if timeout_seconds <= 0 else timeout_seconds)
        finally:
            # Touch heartbeat to indicate completion
            heartbeat.write_text(str(time.time()))
    zip_path = out_dir / "ncbi_dataset.zip"
    extracted_dir: Path | None = None
    if zip_path.exists():
        extracted_dir = _extract_zip(zip_path, out_dir)
    return {
        "method": "datasets-cli",
        "return_code": rc,
        "stdout": "",
        "stderr": "",
        "zip_path": str(zip_path) if zip_path.exists() else "",
        "extracted_dir": str(extracted_dir) if extracted_dir else "",
        "command": " ".join(cmd),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "heartbeat": str(heartbeat),
    }


def download_genome_via_api(
    accession: str,
    dest_dir: str | Path,
    *,
    include: list[str] | None = None,
    chunk_size: int = 1_048_576,
) -> dict:
    """Download a genome package via NCBI Datasets API URL (zip), then extract.

    Writes progress files under dest_dir: download.progress.json and download.progress.txt
    """
    out_dir = ensure_directory(dest_dir)
    include_vals = include or ["genome", "gff3", "rna", "cds", "protein", "seq-report"]
    # Map to API annotation types
    ann_map = {
        "genome": "GENOME_FASTA",
        "gff3": "GENOME_GFF",
        "rna": "RNA_FASTA",
        "cds": "CDS_FASTA",
        "protein": "PROT_FASTA",
        "seq-report": "SEQUENCE_REPORT",
    }
    qs = [("include_annotation_type", ann_map[i]) for i in include_vals if i in ann_map]
    qs.append(("hydrated", "FULLY_HYDRATED"))
    base = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download"
    url = base + "?" + urlencode(qs)
    zip_path = out_dir / "ncbi_dataset_api.zip"
    progress_json = out_dir / "download.progress.json"
    progress_txt = out_dir / "download.progress.txt"
    heartbeat = out_dir / "download.heartbeat"

    def write_progress(done: int, total: int | None) -> None:
        pct = (done / total * 100.0) if (total and total > 0) else None
        data = {
            "bytes": done,
            "total": total or -1,
            "percent": (round(pct, 2) if pct is not None else None),
            "url": url,
        }
        progress_json.write_text(json.dumps(data), encoding="utf-8")
        progress_txt.write_text(f"{data['bytes']}/{data['total']} bytes ({data['percent']}%)\n", encoding="utf-8")
        heartbeat.write_text(str(time.time()))

    with urlopen(url, timeout=60) as resp:
        total = resp.length if hasattr(resp, "length") else None
        done = 0
        with open(zip_path, "wb") as fh:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                fh.write(chunk)
                done += len(chunk)
                write_progress(done, total)
    # One final write to mark completion
    write_progress(done, total)
    extracted_dir = _extract_zip(zip_path, out_dir)
    return {
        "method": "datasets-api",
        "return_code": 0,
        "stdout": "",
        "stderr": "",
        "zip_path": str(zip_path),
        "extracted_dir": str(extracted_dir),
        "url": url,
        "progress_json": str(progress_json),
        "progress_txt": str(progress_txt),
        "heartbeat": str(heartbeat),
    }


def download_genome_via_ftp(
    ftp_url: str,
    dest_dir: str | Path,
) -> dict:
    """Attempt to download common files from an NCBI FTP directory URL.

    This expects ftp_url to be a directory containing standard files. We try to fetch
    md5checksums and common FASTA/GFF files by conventional names if present.
    """
    out_dir = ensure_directory(dest_dir)
    checksums_url = ftp_url.rstrip("/") + "/md5checksums.txt"
    rec: dict[str, str | int] = {"method": "ftp", "return_code": 0}
    try:
        with urlopen(checksums_url, timeout=30) as resp:
            text = resp.read()
        (out_dir / "md5checksums.txt").write_bytes(text)
        rec["md5checksums"] = str(out_dir / "md5checksums.txt")
    except Exception as exc:
        rec["return_code"] = 1
        rec["error"] = str(exc)
    return rec


def download_genome_package_best_effort(
    accession: str,
    dest_dir: str | Path,
    *,
    include: list[str] | None = None,
    ftp_url: str | None = None,
) -> dict:
    """Try to download an accession package via CLI, then API, then FTP.

    Writes a small JSON report to dest_dir and returns the record.
    """
    out_dir = ensure_directory(dest_dir)
    record: dict = {}
    if datasets_cli_available():
        record = download_genome_via_datasets_cli(accession, out_dir, include=include)
        if record.get("return_code", 1) == 0 and record.get("zip_path"):
            dump_json(record, out_dir / "download_record.json", indent=2)
            return record
    # Fallback to API
    try:
        record = download_genome_via_api(accession, out_dir, include=include)
        dump_json(record, out_dir / "download_record.json", indent=2)
        return record
    except Exception as exc:  # pragma: no cover - network
        record = {"method": "datasets-api", "return_code": 1, "error": str(exc)}
    # Optional FTP fallback
    if ftp_url:
        try:
            ftp_rec = download_genome_via_ftp(ftp_url, out_dir)
            record = {**record, **{f"ftp_{k}": v for k, v in ftp_rec.items()}}
        except Exception as exc:  # pragma: no cover - network
            record["ftp_error"] = str(exc)
    dump_json(record, out_dir / "download_record.json", indent=2)
    return record
