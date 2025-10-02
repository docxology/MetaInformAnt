from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Mapping


def ensure_directory(path: str | Path) -> Path:
    """Create a directory (and parents) if missing and return it as Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def open_text_auto(path: str | Path, mode: str = "rt", encoding: str = "utf-8") -> io.TextIOBase:
    """Open a text file, handling gzip transparently based on suffix.

    Supports text modes only ("rt", "wt", "at").
    """
    p = Path(path)
    if "b" in mode:
        raise ValueError("open_text_auto supports text modes only; do not include 'b' in mode")
    if p.suffix == ".gz":
        return io.TextIOWrapper(gzip.open(p, mode.replace("t", "")), encoding=encoding)
    return open(p, mode, encoding=encoding)


# JSON utilities
def load_json(path: str | Path) -> Any:
    with open_text_auto(path, mode="rt") as fh:
        return json.load(fh)


def dump_json(obj: Any, path: str | Path, *, indent: int | None = None) -> None:
    ensure_directory(Path(path).parent)
    with open_text_auto(path, mode="wt") as fh:
        json.dump(obj, fh, indent=indent, sort_keys=True)  # Sort keys for consistent output


# JSON Lines utilities
def read_jsonl(path: str | Path) -> Iterator[Dict[str, Any]]:
    with open_text_auto(path, mode="rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def write_jsonl(rows: Iterable[Mapping[str, Any]], path: str | Path) -> None:
    ensure_directory(Path(path).parent)
    with open_text_auto(path, mode="wt") as fh:
        for row in rows:
            fh.write(json.dumps(dict(row)))
            fh.write("\n")


# Delimited text utilities (CSV/TSV)
def read_delimited(path: str | Path, *, delimiter: str = ",") -> Iterator[Dict[str, str]]:
    with open_text_auto(path, mode="rt") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            yield {k: (v if v is not None else "") for k, v in row.items()}


def write_delimited(rows: Iterable[Mapping[str, Any]], path: str | Path, *, delimiter: str = ",") -> None:
    ensure_directory(Path(path).parent)
    rows_iter = iter(rows)
    try:
        first = next(rows_iter)
    except StopIteration:
        # create empty file
        Path(path).write_text("")
        return
    fieldnames = list(first.keys())
    with open_text_auto(path, mode="wt") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        writer.writerow({k: first.get(k, "") for k in fieldnames})
        for row in rows_iter:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


# Pandas-compatible CSV/TSV utilities
def read_csv(path: str | Path):
    """Read CSV file using pandas if available, fallback to native implementation."""
    try:
        import pandas as pd

        return pd.read_csv(path)
    except ImportError:
        # Fallback to native implementation
        rows = list(read_delimited(path, delimiter=","))
        if not rows:
            return None

        # Convert to simple data structure similar to DataFrame
        import collections

        data = collections.defaultdict(list)
        for row in rows:
            for key, value in row.items():
                data[key].append(value)
        return dict(data)


def write_csv(data, path: str | Path) -> None:
    """Write CSV file using pandas if available, fallback to native implementation."""
    try:
        # Assume pandas DataFrame
        data.to_csv(path, index=False)
    except AttributeError:
        # Handle dict-like data
        if isinstance(data, dict):
            rows = []
            keys = list(data.keys())
            if keys:
                num_rows = len(data[keys[0]])
                for i in range(num_rows):
                    row = {key: data[key][i] for key in keys}
                    rows.append(row)
                write_delimited(rows, path, delimiter=",")


def read_tsv(path: str | Path):
    """Read TSV file."""
    with open_text_auto(path, mode="rt") as fh:
        reader = csv.reader(fh, delimiter="\t")
        return list(reader)


def write_tsv(data, path: str | Path) -> None:
    """Write TSV file."""
    ensure_directory(Path(path).parent)
    with open_text_auto(path, mode="wt") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for row in data:
            writer.writerow(row)
