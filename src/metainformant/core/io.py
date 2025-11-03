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
    """Load JSON data from a file.
    
    Args:
        path: Path to JSON file (supports gzip compression automatically)
        
    Returns:
        Parsed JSON data (dict, list, or primitive types)
    """
    with open_text_auto(path, mode="rt") as fh:
        return json.load(fh)


def dump_json(obj: Any, path: str | Path, *, indent: int | None = None) -> None:
    """Write object to JSON file.
    
    Args:
        obj: Object to serialize (must be JSON-serializable)
        path: Output file path
        indent: Number of spaces for indentation (None for compact)
    """
    ensure_directory(Path(path).parent)
    with open_text_auto(path, mode="wt") as fh:
        json.dump(obj, fh, indent=indent, sort_keys=True)  # Sort keys for consistent output


def dump_json_gz(obj: Any, path: str | Path, *, indent: int | None = None) -> None:
    """Dump object as gzipped JSON file."""
    ensure_directory(Path(path).parent)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=indent, sort_keys=True)


def load_json_gz(path: str | Path) -> Any:
    """Load object from gzipped JSON file."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def read_csv(path: str | Path, **kwargs) -> Any:
    """Read CSV file with pandas, handling gzip automatically."""
    try:
        import pandas as pd
        if str(path).endswith('.gz'):
            return pd.read_csv(path, **kwargs)
        else:
            return pd.read_csv(path, **kwargs)
    except ImportError:
        raise ImportError("pandas is required for CSV reading. Install with: uv add pandas")


def write_csv(df: Any, path: str | Path, **kwargs) -> None:
    """Write DataFrame to CSV file."""
    try:
        import pandas as pd
        ensure_directory(Path(path).parent)
        df.to_csv(path, **kwargs)
    except ImportError:
        raise ImportError("pandas is required for CSV writing. Install with: uv add pandas")


def read_parquet(path: str | Path, **kwargs) -> Any:
    """Read Parquet file with pandas."""
    try:
        import pandas as pd
        return pd.read_parquet(path, **kwargs)
    except ImportError:
        raise ImportError("pandas is required for Parquet reading. Install with: uv add pandas")


def write_parquet(df: Any, path: str | Path, **kwargs) -> None:
    """Write DataFrame to Parquet file."""
    try:
        import pandas as pd
        ensure_directory(Path(path).parent)
        df.to_parquet(path, **kwargs)
    except ImportError:
        raise ImportError("pandas is required for Parquet writing. Install with: uv add pandas")


# JSON Lines utilities
def read_jsonl(path: str | Path) -> Iterator[Dict[str, Any]]:
    """Read JSON Lines format (one JSON object per line).
    
    Args:
        path: Path to JSONL file (supports gzip compression automatically)
        
    Yields:
        Dictionary for each line in the file
    """
    with open_text_auto(path, mode="rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def write_jsonl(rows: Iterable[Mapping[str, Any]], path: str | Path) -> None:
    """Write rows as JSON Lines format (one JSON object per line).
    
    Args:
        rows: Iterable of dictionaries to write
        path: Output file path
    """
    ensure_directory(Path(path).parent)
    with open_text_auto(path, mode="wt") as fh:
        for row in rows:
            fh.write(json.dumps(dict(row)))
            fh.write("\n")


# Delimited text utilities (CSV/TSV)
def read_delimited(path: str | Path, *, delimiter: str = ",") -> Iterator[Dict[str, str]]:
    """Read delimited text file (CSV/TSV) as dictionaries.
    
    Args:
        path: Path to delimited file (supports gzip compression automatically)
        delimiter: Field delimiter (default: comma for CSV, use "\\t" for TSV)
        
    Yields:
        Dictionary for each row with column names as keys
    """
    with open_text_auto(path, mode="rt") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            yield {k: (v if v is not None else "") for k, v in row.items()}


def write_delimited(rows: Iterable[Mapping[str, Any]], path: str | Path, *, delimiter: str = ",") -> None:
    """Write rows to delimited text file (CSV/TSV).
    
    Args:
        rows: Iterable of dictionaries to write
        path: Output file path
        delimiter: Field delimiter (default: comma for CSV, use "\\t" for TSV)
    """
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


def download_file(url: str, dest_path: str | Path, *, chunk_size: int = 8192, timeout: int = 30) -> bool:
    """Download a file from a URL to a local path.
    
    Args:
        url: URL to download from
        dest_path: Local path to save the file
        chunk_size: Size of download chunks
        timeout: Request timeout in seconds
        
    Returns:
        True if download successful, False otherwise
    """
    import requests
    from urllib.parse import urlparse
    
    try:
        # Create parent directory if it doesn't exist
        dest_path = Path(dest_path)
        ensure_directory(dest_path.parent)
        
        # Parse URL to get filename if dest_path is a directory
        if dest_path.is_dir():
            parsed_url = urlparse(url)
            filename = Path(parsed_url.path).name
            if not filename:
                filename = "downloaded_file"
            dest_path = dest_path / filename
            
        # Download with progress tracking
        response = requests.get(url, stream=True, timeout=timeout)
        response.raise_for_status()
        
        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    
        return True
        
    except Exception:
        return False


def download_json(url: str, *, timeout: int = 30) -> Any:
    """Download and parse JSON from a URL.
    
    Args:
        url: URL to download JSON from
        timeout: Request timeout in seconds
        
    Returns:
        Parsed JSON data or None if failed
    """
    import requests
    
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except Exception:
        return None


def download_text(url: str, *, timeout: int = 30) -> str | None:
    """Download text content from a URL.
    
    Args:
        url: URL to download text from
        timeout: Request timeout in seconds
        
    Returns:
        Text content or None if failed
    """
    import requests
    
    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except Exception:
        return None


def download_csv(url: str, *, timeout: int = 30, **kwargs) -> Any:
    """Download and parse CSV from a URL.
    
    Args:
        url: URL to download CSV from
        timeout: Request timeout in seconds
        **kwargs: Additional arguments for pandas.read_csv
        
    Returns:
        DataFrame or None if failed
    """
    try:
        import pandas as pd
        import io
        
        text_content = download_text(url, timeout=timeout)
        if text_content:
            return pd.read_csv(io.StringIO(text_content), **kwargs)
    except Exception:
        pass
    
    return None


def batch_download(urls: list[str], dest_dir: str | Path, *, timeout: int = 30) -> dict[str, bool]:
    """Download multiple files in batch.
    
    Args:
        urls: List of URLs to download
        dest_dir: Directory to save files to
        timeout: Request timeout in seconds
        
    Returns:
        Dictionary mapping URLs to success status
    """
    dest_dir = Path(dest_dir)
    ensure_directory(dest_dir)
    
    results = {}
    for i, url in enumerate(urls):
        try:
            from urllib.parse import urlparse
            parsed_url = urlparse(url)
            filename = Path(parsed_url.path).name
            if not filename:
                filename = f"download_{i}.txt"
            
            dest_path = dest_dir / filename
            success = download_file(url, dest_path, timeout=timeout)
            results[url] = success
            
        except Exception:
            results[url] = False
            
    return results
