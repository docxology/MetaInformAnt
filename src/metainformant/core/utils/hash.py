from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_bytes(data: bytes) -> str:
    """Return the SHA-256 hex digest of bytes."""
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 hex digest of a file without loading it entirely."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def deterministic_seed(data: str) -> int:
    """Generate deterministic integer seed from string data.

    Args:
        data: String to hash for seed generation

    Returns:
        Integer seed in valid range for random number generators
    """
    hash_digest = sha256_bytes(data.encode("utf-8"))
    # Convert first 8 hex chars to integer and ensure it's in valid range
    seed = int(hash_digest[:8], 16) % (2**31 - 1)
    return seed


def sha256_string(s: str) -> str:
    """Compute SHA256 hash of string content."""
    return sha256_bytes(s.encode("utf-8"))


def file_hash_comparison(file1: str | Path, file2: str | Path) -> bool:
    """Compare two files by their SHA256 hashes."""
    return sha256_file(file1) == sha256_file(file2)


def hash_directory(path: str | Path, pattern: str = "**/*") -> dict[str, str]:
    """Compute hashes for all files in a directory matching a pattern."""
    p = Path(path)
    hashes = {}

    for file_path in p.glob(pattern):
        if file_path.is_file():
            try:
                hashes[str(file_path.relative_to(p))] = sha256_file(file_path)
            except (OSError, IOError):
                # Skip files that can't be read
                continue

    return hashes


def verify_file_integrity(file_path: str | Path, expected_hash: str) -> bool:
    """Verify file integrity against expected hash."""
    try:
        actual_hash = sha256_file(file_path)
        return actual_hash == expected_hash
    except (OSError, IOError):
        return False
