from __future__ import annotations

from pathlib import Path


def read_taxon_ids(taxon_id_file: Path) -> list[int]:
    ids: list[int] = []
    for line in taxon_id_file.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            ids.append(int(line))
        except ValueError:
            continue
    return ids


