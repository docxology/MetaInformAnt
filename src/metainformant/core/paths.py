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


def find_files_by_extension(directory: str | Path, extension: str) -> list[Path]:
    """Find all files with a specific extension in a directory.

    Args:
        directory: Directory to search
        extension: File extension (with or without dot)

    Returns:
        List of matching file paths
    """
    if not extension.startswith("."):
        extension = "." + extension

    dir_path = Path(directory)
    return list(dir_path.rglob(f"*{extension}"))


def get_file_size(path: str | Path) -> int:
    """Get file size in bytes.

    Args:
        path: Path to file

    Returns:
        File size in bytes, or 0 if file doesn't exist
    """
    try:
        return Path(path).stat().st_size
    except (OSError, FileNotFoundError):
        return 0


def get_directory_size(path: str | Path) -> int:
    """Get total size of all files in a directory.

    Args:
        path: Directory path

    Returns:
        Total size in bytes
    """
    total_size = 0
    dir_path = Path(path)

    if not dir_path.exists():
        return 0

    for file_path in dir_path.rglob("*"):
        if file_path.is_file():
            try:
                total_size += file_path.stat().st_size
            except (OSError, FileNotFoundError):
                continue

    return total_size


def sanitize_filename(filename: str) -> str:
    """Sanitize filename for safe filesystem use.

    Args:
        filename: Original filename

    Returns:
        Sanitized filename safe for filesystem
    """
    import re

    # Remove or replace dangerous characters
    sanitized = re.sub(r'[<>:"/\\|?*]', '_', filename)
    # Remove control characters
    sanitized = re.sub(r'[\x00-\x1f\x7f-\x9f]', '', sanitized)
    # Remove leading/trailing dots and spaces
    sanitized = sanitized.strip(' .')

    # Ensure filename is not empty
    if not sanitized:
        sanitized = "untitled"

    # Truncate if too long (255 chars max for most filesystems)
    if len(sanitized) > 255:
        name_part, ext = Path(sanitized).stem, Path(sanitized).suffix
        max_name = 255 - len(ext)
        sanitized = name_part[:max_name] + ext

    return sanitized


def create_temp_file(suffix: str = "", prefix: str = "tmp", directory: str | Path | None = None) -> Path:
    """Create a temporary file path that doesn't exist yet.

    Args:
        suffix: File extension
        prefix: Filename prefix
        directory: Directory for temp file (uses system temp if None)

    Returns:
        Path to non-existent temporary file
    """
    import tempfile

    if directory is None:
        directory = tempfile.gettempdir()

    dir_path = Path(directory)
    ensure_directory(dir_path)

    while True:
        temp_name = f"{prefix}_{hash(str(Path.home()))}_{suffix}"
        temp_path = dir_path / temp_name
        if not temp_path.exists():
            return temp_path
