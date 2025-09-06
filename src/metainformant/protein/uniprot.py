from __future__ import annotations

import time
from typing import Dict, Iterable

import requests

BASE = "https://rest.uniprot.org"  # documented base for UniProt REST API


def map_ids_uniprot(ids: Iterable[str]) -> Dict[str, str]:
    """Map UniProt accessions/IDs to canonical UniProtKB primary accession.

    Performs an ID mapping job and polls until finished, then returns a mapping
    of input id → primaryAccession for minimal functionality under tests.
    """
    ids_list = list(ids)
    if not ids_list:
        return {}

    # Start mapping job
    run_url = f"{BASE}/idmapping/run"
    data = {
        "from": "UniProtKB_AC-ID",
        "to": "UniProtKB",
        "ids": ",".join(ids_list),
    }
    resp = requests.post(run_url, data=data, headers={"Accept": "application/json"}, timeout=30)
    resp.raise_for_status()
    job = resp.json()
    job_id = job.get("jobId") or job.get("job_id")
    if not job_id:
        return {}

    # Poll for completion
    status_url = f"{BASE}/idmapping/status/{job_id}"
    for _ in range(30):
        s = requests.get(status_url, headers={"Accept": "application/json"}, timeout=30)
        s.raise_for_status()
        js = s.json()
        status = js.get("jobStatus") or js.get("status")
        if status == "FINISHED":
            break
        time.sleep(1)

    # Fetch results
    result_url = f"{BASE}/idmapping/results/{job_id}"
    r = requests.get(result_url, headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    out: Dict[str, str] = {}
    for row in r.json().get("results", []):
        from_id = row.get("from")
        to = row.get("to") or {}
        primary = to.get("primaryAccession") or to.get("primaryAccessionId") or to
        if isinstance(primary, str) and isinstance(from_id, str):
            out[from_id] = primary
    return out


def fetch_uniprot_fasta(accession: str) -> str:
    """Fetch UniProt FASTA sequence for a single accession and return as string.

    Returns the raw FASTA text. Caller can write to file if desired.
    """
    acc = accession.strip().upper()
    url = f"{BASE}/uniprotkb/{acc}.fasta"
    r = requests.get(url, headers={"Accept": "text/plain"}, timeout=60)
    r.raise_for_status()
    return r.text
