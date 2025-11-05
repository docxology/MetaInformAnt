"""Error handling and resilience utilities for METAINFORMANT.

Provides custom exception hierarchy, retry decorators, and error recovery
utilities for robust operation across all modules.
"""

from __future__ import annotations

import functools
import time
from contextlib import contextmanager
from typing import Any, Callable, TypeVar

T = TypeVar("T")


# Custom exception hierarchy
class METAINFORMANTError(Exception):
    """Base exception for all METAINFORMANT errors."""

    pass


class ConfigError(METAINFORMANTError):
    """Configuration-related errors."""

    pass


class IOError(METAINFORMANTError):
    """Input/output operation errors."""

    pass


class ValidationError(METAINFORMANTError):
    """Data validation errors."""

    pass


class NetworkError(METAINFORMANTError):
    """Network operation errors."""

    pass


class CacheError(METAINFORMANTError):
    """Cache operation errors."""

    pass


def retry_with_backoff(
    max_attempts: int = 3,
    initial_delay: float = 1.0,
    backoff_factor: float = 2.0,
    max_delay: float = 60.0,
    exceptions: tuple[type[Exception], ...] = (Exception,),
):
    """Decorator for retrying functions with exponential backoff.

    Args:
        max_attempts: Maximum number of retry attempts
        initial_delay: Initial delay in seconds before first retry
        backoff_factor: Multiplier for delay between retries
        max_delay: Maximum delay in seconds between retries
        exceptions: Tuple of exception types to catch and retry on

    Returns:
        Decorated function that retries on specified exceptions

    Example:
        @retry_with_backoff(max_attempts=5, initial_delay=2.0)
        def download_file(url):
            # This will retry up to 5 times with exponential backoff
            return requests.get(url)
    """

    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> T:
            delay = initial_delay
            last_exception = None

            for attempt in range(max_attempts):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e
                    if attempt < max_attempts - 1:
                        time.sleep(min(delay, max_delay))
                        delay *= backoff_factor
                    else:
                        raise

            # Should not reach here, but for type safety
            if last_exception:
                raise last_exception
            raise RuntimeError("Retry failed without exception")

        return wrapper

    return decorator


@contextmanager
def error_context(context_msg: str, reraise: bool = True):
    """Context manager for adding error context.

    Args:
        context_msg: Message to add to exception
        reraise: Whether to reraise the exception with context

    Example:
        with error_context("Failed to process file"):
            # Operations that might fail
            process_file(path)
    """
    try:
        yield
    except Exception as e:
        if reraise:
            # Add context to error message
            error_msg = f"{context_msg}: {str(e)}"
            raise type(e)(error_msg).with_traceback(e.__traceback__) from e
        raise


def safe_execute(func: Callable[..., T], *args: Any, default: T | None = None, **kwargs: Any) -> T | None:
    """Execute a function safely, returning default on error.

    Args:
        func: Function to execute
        *args: Positional arguments for function
        default: Default value to return on error
        **kwargs: Keyword arguments for function

    Returns:
        Function result or default value on error

    Example:
        result = safe_execute(risky_function, arg1, arg2, default="fallback")
    """
    try:
        return func(*args, **kwargs)
    except Exception:
        return default


def validate_not_none(value: Any, name: str = "value") -> None:
    """Validate that a value is not None.

    Args:
        value: Value to validate
        name: Name of value for error message

    Raises:
        ValidationError: If value is None
    """
    if value is None:
        raise ValidationError(f"{name} cannot be None")


def validate_type(value: Any, expected_type: type | tuple[type, ...], name: str = "value") -> None:
    """Validate that a value is of expected type.

    Args:
        value: Value to validate
        expected_type: Expected type or tuple of types
        name: Name of value for error message

    Raises:
        ValidationError: If value is not of expected type
    """
    if not isinstance(value, expected_type):
        type_names = expected_type.__name__ if isinstance(expected_type, type) else ", ".join(t.__name__ for t in expected_type)
        raise ValidationError(f"{name} must be of type {type_names}, got {type(value).__name__}")


