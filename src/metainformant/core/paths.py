from __future__ import annotations

from pathlib import Path


def expand_and_resolve(path: str | Path) -> Path:
    """Expand user (~) and resolve to absolute path without requiring existence."""
    return Path(path).expanduser().resolve(strict=False)


def is_within(path: str | Path, parent: str | Path) -> bool:
    """Return True if path is inside parent directory (after resolving)."""
    p = expand_and_resolve(path)
    r = expand_and_resolve(parent)
    try:
        p.relative_to(r)
        return True
    except ValueError:
        return False


def ensure_directory(path: Path) -> None:
    """Create directory and any missing parent directories.

    Args:
        path: Directory path to create
    """
    path.mkdir(parents=True, exist_ok=True)


def prepare_file_path(file_path: Path) -> None:
    """Ensure parent directories exist for a file path.

    Args:
        file_path: File path to prepare
    """
    if file_path.parent:
        ensure_directory(file_path.parent)


def is_safe_path(path: str) -> bool:
    """Check if path is safe (no path traversal attempts).

    Args:
        path: Path string to check

    Returns:
        True if path appears safe
    """
    # Check for path traversal attempts
    if ".." in path or path.startswith("/etc/") or path.startswith("/root/"):
        return False

    # Check for other potentially dangerous patterns
    dangerous_patterns = ["..//", "..\\", ";", "|", "&", "$"]
    for pattern in dangerous_patterns:
        if pattern in path:
            return False

    return True


def get_file_extension(filename: str) -> str:
    """Get file extension from filename.

    Args:
        filename: Name of file

    Returns:
        File extension including the dot, or empty string
    """
    return Path(filename).suffix


def change_extension(path: str, new_extension: str) -> Path:
    """Change file extension.

    Args:
        path: Original file path
        new_extension: New extension (with or without dot)

    Returns:
        Path with new extension
    """
    p = Path(path)
    if not new_extension.startswith("."):
        new_extension = "." + new_extension

    return p.with_suffix(new_extension)
