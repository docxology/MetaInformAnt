"""Core I/O exceptions."""

from metainformant.core.utils.errors import METAINFORMANTError


class IOError(METAINFORMANTError):
    """Base class for I/O errors in metainformant."""

    pass


class FileNotFoundError(IOError):
    """File not found error."""

    pass


class CacheError(IOError):
    """Cache operation error."""

    pass


class DownloadError(IOError):
    """Download operation error."""

    pass
