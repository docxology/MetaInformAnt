from __future__ import annotations

from pathlib import Path


def read_taxon_ids(taxon_id_file: Path) -> list[int]:
    """Read NCBI taxonomy IDs from a text file (one per line).
    
    Args:
        taxon_id_file: Path to file containing taxonomy IDs
        
    Returns:
        List of taxonomy IDs (integers), skipping empty lines and comments
    """
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
