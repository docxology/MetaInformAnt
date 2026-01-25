from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping


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

    Raises:
        IOError: If file read fails or JSON parsing fails
        FileNotFoundError: If file does not exist
    """
    from .errors import IOError as CoreIOError

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"JSON file not found: {path}")

    try:
        with open_text_auto(p, mode="rt") as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        raise CoreIOError(f"Failed to parse JSON file {path}: {e}") from e
    except Exception as e:
        raise CoreIOError(f"Failed to read JSON file {path}: {e}") from e


def dump_json(obj: Any, path: str | Path, *, indent: int | None = None, atomic: bool = True) -> None:
    """Write object to JSON file with optional atomic write.

    Args:
        obj: Object to serialize (must be JSON-serializable)
        path: Output file path (supports .gz extension for gzip compression)
        indent: Number of spaces for indentation (None for compact)
        atomic: If True, write to temp file then rename (prevents corruption)

    Raises:
        IOError: If file write fails
    """
    from .errors import IOError as CoreIOError

    p = Path(path)
    ensure_directory(p.parent)

    # Check if path ends with .gz to determine compression
    is_gzipped = p.suffix == ".gz"

    try:
        if atomic:
            # Atomic write: write to temp file, then rename
            # For .gz files, keep .gz extension on temp file so compression works
            if is_gzipped:
                # For file.json.gz, create file.json.gz.tmp
                temp_path = p.with_suffix(p.suffix + ".tmp")
            else:
                temp_path = p.with_suffix(p.suffix + ".tmp")

            if is_gzipped:
                with gzip.open(temp_path, "wt", encoding="utf-8") as fh:
                    json.dump(obj, fh, indent=indent, sort_keys=True)
            else:
                with open_text_auto(temp_path, mode="wt") as fh:
                    json.dump(obj, fh, indent=indent, sort_keys=True)
            temp_path.replace(p)  # Atomic rename
        else:
            if is_gzipped:
                with gzip.open(p, "wt", encoding="utf-8") as fh:
                    json.dump(obj, fh, indent=indent, sort_keys=True)
            else:
                with open_text_auto(p, mode="wt") as fh:
                    json.dump(obj, fh, indent=indent, sort_keys=True)
    except Exception as e:
        raise CoreIOError(f"Failed to write JSON file {path}: {e}") from e


def dump_json_gz(obj: Any, path: str | Path, *, indent: int | None = None) -> None:
    """Dump object as gzipped JSON file."""
    ensure_directory(Path(path).parent)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=indent, sort_keys=True)


def load_json_gz(path: str | Path) -> Any:
    """Load object from gzipped JSON file."""
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return json.load(fh)


def dump_yaml(obj: Any, path: str | Path) -> None:
    """Dump object as YAML file.

    Requires PyYAML to be installed. Falls back to JSON if unavailable.

    Args:
        obj: Object to serialize
        path: Output file path

    Raises:
        ImportError: If PyYAML is not available and fallback fails
    """
    ensure_directory(Path(path).parent)

    try:
        import yaml

        with open(path, "w", encoding="utf-8") as fh:
            yaml.dump(obj, fh, default_flow_style=False, sort_keys=False)
    except ImportError:
        # Fallback to JSON with .yaml extension
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2, sort_keys=True)
            fh.write("\n")


def load_yaml(path: str | Path) -> Any:
    """Load YAML data from a file.

    Args:
        path: Path to YAML file

    Returns:
        Parsed YAML data

    Raises:
        ImportError: If PyYAML is not available
        IOError: If file read fails
    """
    from .errors import IOError as CoreIOError

    try:
        import yaml

        with open_text_auto(path, mode="rt") as fh:
            return yaml.safe_load(fh)
    except ImportError as e:
        raise ImportError("PyYAML is required for YAML support. Install with: uv add PyYAML") from e
    except Exception as e:
        raise CoreIOError(f"Failed to read YAML file {path}: {e}") from e


def load_toml(path: str | Path) -> Any:
    """Load TOML data from a file.

    Args:
        path: Path to TOML file

    Returns:
        Parsed TOML data

    Raises:
        IOError: If file read fails
    """
    from .errors import IOError as CoreIOError

    try:
        import tomllib

        with open(path, "rb") as fh:
            return tomllib.load(fh)
    except ImportError as e:
        raise ImportError("tomllib is required for TOML support (Python 3.11+).") from e
    except Exception as e:
        raise CoreIOError(f"Failed to read TOML file {path}: {e}") from e


def read_parquet(path: str | Path, **kwargs) -> Any:
    """Read Parquet file with pandas."""
    try:
        import pandas as pd

        return pd.read_parquet(path, **kwargs)
    except ImportError as e:
        # Preserve original error message if it mentions pyarrow/fastparquet
        error_msg = str(e).lower()
        if "pyarrow" in error_msg or "fastparquet" in error_msg:
            raise ImportError("Parquet support requires pyarrow or fastparquet. Install with: uv add pyarrow") from e
        raise ImportError(f"pandas is required for Parquet reading: {e}. Install with: uv add pandas") from e


def write_parquet(df: Any, path: str | Path, **kwargs) -> None:
    """Write DataFrame to Parquet file."""
    try:
        import pandas as pd

        ensure_directory(Path(path).parent)
        df.to_parquet(path, **kwargs)
    except ImportError as e:
        # Preserve original error message if it mentions pyarrow/fastparquet
        error_msg = str(e).lower()
        if "pyarrow" in error_msg or "fastparquet" in error_msg:
            raise ImportError("Parquet support requires pyarrow or fastparquet. Install with: uv add pyarrow") from e
        raise ImportError(f"pandas is required for Parquet writing: {e}. Install with: uv add pandas") from e


# JSON Lines utilities
def read_jsonl(path: str | Path) -> Iterator[dict[str, Any]]:
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


def write_jsonl(rows: Iterable[Mapping[str, Any]], path: str | Path, *, atomic: bool = True) -> None:
    """Write rows as JSON Lines format (one JSON object per line).

    Args:
        rows: Iterable of dictionaries to write
        path: Output file path
        atomic: If True, write to temp file then rename (prevents corruption)

    Raises:
        IOError: If file write fails
    """
    from .errors import IOError as CoreIOError

    p = Path(path)
    ensure_directory(p.parent)

    try:
        if atomic:
            temp_path = p.with_suffix(p.suffix + ".tmp")
            with open_text_auto(temp_path, mode="wt") as fh:
                for row in rows:
                    fh.write(json.dumps(dict(row)))
                    fh.write("\n")
            temp_path.replace(p)
        else:
            with open_text_auto(p, mode="wt") as fh:
                for row in rows:
                    fh.write(json.dumps(dict(row)))
                    fh.write("\n")
    except Exception as e:
        raise CoreIOError(f"Failed to write JSONL file {path}: {e}") from e


# Delimited text utilities (CSV/TSV)
def read_delimited(path: str | Path, *, delimiter: str = ",") -> Iterator[dict[str, str]]:
    """Read delimited text file (CSV/TSV) as dictionaries.

    Args:
        path: Path to delimited file (supports gzip compression automatically)
        delimiter: Field delimiter (default: comma for CSV, use "\\t" for TSV)

    Yields:
        Dictionary for each row with column names as keys

    Raises:
        FileNotFoundError: If file does not exist
        IOError: If file read fails
    """
    from .errors import IOError as CoreIOError

    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {path}")

    try:
        with open_text_auto(p, mode="rt") as fh:
            reader = csv.DictReader(fh, delimiter=delimiter)
            for row in reader:
                yield {k: (v if v is not None else "") for k, v in row.items()}
    except Exception as e:
        raise CoreIOError(f"Failed to read delimited file {path}: {e}") from e


def write_delimited(
    rows: Iterable[Mapping[str, Any]], path: str | Path, *, delimiter: str = ",", atomic: bool = True
) -> None:
    """Write rows to delimited text file (CSV/TSV).

    Args:
        rows: Iterable of dictionaries to write
        path: Output file path
        delimiter: Field delimiter (default: comma for CSV, use "\\t" for TSV)
        atomic: If True, write to temp file then rename (prevents corruption)

    Raises:
        IOError: If file write fails
    """
    from .errors import IOError as CoreIOError

    p = Path(path)
    ensure_directory(p.parent)

    try:
        rows_iter = iter(rows)
        try:
            first = next(rows_iter)
        except StopIteration:
            # create empty file
            if atomic:
                temp_path = p.with_suffix(p.suffix + ".tmp")
                temp_path.write_text("")
                temp_path.replace(p)
            else:
                p.write_text("")
            return

        fieldnames = list(first.keys())

        if atomic:
            temp_path = p.with_suffix(p.suffix + ".tmp")
            with open_text_auto(temp_path, mode="wt") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
                writer.writeheader()
                writer.writerow({k: first.get(k, "") for k in fieldnames})
                for row in rows_iter:
                    writer.writerow({k: row.get(k, "") for k in fieldnames})
            temp_path.replace(p)
        else:
            with open_text_auto(p, mode="wt") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
                writer.writeheader()
                writer.writerow({k: first.get(k, "") for k in fieldnames})
                for row in rows_iter:
                    writer.writerow({k: row.get(k, "") for k in fieldnames})
    except Exception as e:
        raise CoreIOError(f"Failed to write delimited file {path}: {e}") from e


# Pandas-compatible CSV/TSV utilities
def read_csv(path: str | Path, **kwargs) -> Any:
    """Read CSV file using pandas if available, fallback to native implementation.

    Args:
        path: Path to CSV file (supports gzip compression automatically)
        **kwargs: Additional arguments passed to pandas.read_csv if pandas available

    Returns:
        pandas DataFrame if pandas available, otherwise dict of lists
    """
    try:
        import pandas as pd

        return pd.read_csv(path, **kwargs)
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


def write_csv(data: Any, path: str | Path, **kwargs) -> None:
    """Write CSV file using pandas if available, fallback to native implementation.

    Args:
        data: pandas DataFrame or dict-like data to write
        path: Output file path
        **kwargs: Additional arguments passed to pandas DataFrame.to_csv if pandas available
    """
    try:
        import pandas as pd

        ensure_directory(Path(path).parent)
        # Default to not writing index unless explicitly requested
        if "index" not in kwargs:
            kwargs["index"] = False
        # Assume pandas DataFrame
        if isinstance(data, pd.DataFrame):
            data.to_csv(path, **kwargs)
        else:
            # Try to convert to DataFrame
            df = pd.DataFrame(data)
            df.to_csv(path, **kwargs)
    except (ImportError, AttributeError):
        # Fallback to native implementation
        if isinstance(data, dict):
            rows = []
            keys = list(data.keys())
            if keys:
                num_rows = len(data[keys[0]])
                for i in range(num_rows):
                    row = {key: data[key][i] for key in keys}
                    rows.append(row)
                write_delimited(rows, path, delimiter=",")
        else:
            # Try to convert to list of dicts
            rows = []
            if hasattr(data, "__iter__"):
                for row in data:
                    if isinstance(row, dict):
                        rows.append(row)
                    elif hasattr(row, "_asdict"):  # namedtuple
                        rows.append(row._asdict())
            if rows:
                write_delimited(rows, path, delimiter=",")


def read_tsv(path: str | Path) -> list[list[str]]:
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
        chunk_size: Size of download chunks (in bytes)
        timeout: Request timeout in seconds

    Returns:
        True if download successful, False otherwise
    """
    from metainformant.core import logging as core_logging

    logger = core_logging.get_logger(__name__)

    try:
        from .download import download_with_progress

        # Keep existing signature/behavior: silent + boolean return.
        # Also keep default chunk_size small to match previous behavior.
        result = download_with_progress(
            url,
            dest_path,
            chunk_size=chunk_size,
            timeout=timeout,
            show_progress=False,
            heartbeat_interval=5,
            max_retries=3,
            resume=True,
        )
        return bool(result.success)
    except (OSError, IOError) as e:
        logger.debug(f"Download failed due to I/O error: {e}")
        return False
    except ImportError as e:
        logger.warning(f"Download module unavailable: {e}")
        return False
    except ValueError as e:
        logger.debug(f"Download failed due to invalid value: {e}")
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

    from metainformant.core import logging as core_logging

    logger = core_logging.get_logger(__name__)

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        logger.debug(f"Failed to download JSON from {url}: {e}")
        return None
    except json.JSONDecodeError as e:
        logger.debug(f"Failed to parse JSON from {url}: {e}")
        return None
    except (ValueError, TypeError) as e:
        logger.debug(f"Invalid JSON data from {url}: {e}")
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

    from metainformant.core import logging as core_logging

    logger = core_logging.get_logger(__name__)

    try:
        response = requests.get(url, timeout=timeout)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        logger.debug(f"Failed to download text from {url}: {e}")
        return None
    except (UnicodeDecodeError, ValueError) as e:
        logger.debug(f"Failed to decode text from {url}: {e}")
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
    from metainformant.core import logging as core_logging

    logger = core_logging.get_logger(__name__)

    try:
        import pandas as pd
        import io as std_io

        text_content = download_text(url, timeout=timeout)
        if text_content:
            return pd.read_csv(std_io.StringIO(text_content), **kwargs)
    except ImportError as e:
        logger.debug(f"pandas not available for CSV parsing: {e}")
    except (ValueError, TypeError) as e:
        logger.debug(f"Failed to parse CSV from {url}: {e}")

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
    from metainformant.core import logging as core_logging

    logger = core_logging.get_logger(__name__)

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

        except (ValueError, OSError) as e:
            logger.debug(f"Failed to download {url}: {e}")
            results[url] = False

    return results
