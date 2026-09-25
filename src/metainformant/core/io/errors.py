"""Core I/O exceptions.

The canonical exception hierarchy lives in
:mod:`metainformant.core.utils.errors`; this module re-exports its I/O classes
so both import paths resolve to the same objects -- a caller catching
``metainformant.core.io.errors.IOError`` also catches what
``metainformant.core.utils.errors.IOError`` (and vice versa).
"""

from metainformant.core.utils.errors import CacheError, DownloadError, IOError

__all__ = ["CacheError", "DownloadError", "IOError", "InputFileMissingError"]


class InputFileMissingError(IOError):
    """An expected input file does not exist.

    Replaces the historical ``FileNotFoundError`` alias, which shadowed the
    builtin exception of the same name.
    """

    pass
