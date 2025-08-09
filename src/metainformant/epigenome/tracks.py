from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd

from metainformant.core.io import open_text_auto


def read_bedgraph(path: str | Path) -> pd.DataFrame:
    """Read a simple bedGraph file into a DataFrame with columns chrom,start,end,value.

    - Supports plain text or gzip-compressed input (".gz").
    - Ignores header/comment lines starting with '#'.
    - start/end are returned as int; value as float.
    """
    rows: List[tuple[str, int, int, float]] = []
    with open_text_auto(path, mode="rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 4:
                continue
            chrom, start, end, value = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            rows.append((chrom, start, end, value))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "value"])



