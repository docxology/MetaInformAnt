from __future__ import annotations

from typing import List, Dict

import requests


def fetch_interpro_domains(uniprot_acc: str) -> List[Dict]:
    """Fetch InterPro entries for a UniProt accession via REST API.

    Returns a list of result dictionaries (shape depends on API version).
    """
    acc = uniprot_acc.strip().upper()
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{acc}?format=json"
    r = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
    r.raise_for_status()
    js = r.json()
    # Normalize: API may return dict with 'results' or already a list
    if isinstance(js, dict) and "results" in js:
        return list(js.get("results") or [])
    if isinstance(js, list):
        return js
    return []


