"""Comprehensive tests for METAINFORMANT error hierarchy and utilities.

Tests the unified error hierarchy where all errors inherit from METAINFORMANTError,
including io.errors module classes, and validates error handling utilities.
"""

from __future__ import annotations

import time
from typing import Any

import pytest

from metainformant.core.io import errors as io_errors
from metainformant.core.utils.errors import (
    CacheError,
    ConfigError,
    DependencyError,
    DownloadError,
    IOError,
    METAINFORMANTError,
    NetworkError,
    PipelineError,
    ResourceError,
    ValidationError,
    error_context,
    retry_with_backoff,
    safe_execute,
    validate_not_none,
    validate_type,
)

# ============================================================================
# Test Error Hierarchy
# ============================================================================


class TestErrorHierarchy:
    """Test that all error types inherit from METAINFORMANTError correctly."""

    def test_base_error_exists(self) -> None:
        """Test that METAINFORMANTError is the base exception."""
        error = METAINFORMANTError("test error")
        assert isinstance(error, Exception)
        assert str(error) == "test error"

    def test_utils_errors_inherit_from_base(self) -> None:
        """Test that all utils.errors classes inherit from METAINFORMANTError."""
        utils_errors = [
            ConfigError,
            IOError,
            ValidationError,
            NetworkError,
            CacheError,
            DownloadError,
            PipelineError,
            DependencyError,
            ResourceError,
        ]

        for error_class in utils_errors:
            error = error_class("test")
            assert isinstance(
                error, METAINFORMANTError
            ), f"{error_class.__name__} does not inherit from METAINFORMANTError"
            assert isinstance(error, Exception)

    def test_io_errors_inherit_from_base(self) -> None:
        """Test that all io.errors classes inherit from METAINFORMANTError."""
        io_error_classes = [
            io_errors.IOError,
            io_errors.FileNotFoundError,
            io_errors.CacheError,
            io_errors.DownloadError,
        ]

        for error_class in io_error_classes:
            error = error_class("test")
            assert isinstance(
                error, METAINFORMANTError
            ), f"io.errors.{error_class.__name__} does not inherit from METAINFORMANTError"
            assert isinstance(error, Exception)

    def test_download_error_inherits_from_network_error(self) -> None:
        """Test that DownloadError inherits from NetworkError."""
        error = DownloadError("download failed")
        assert isinstance(error, NetworkError)
        assert isinstance(error, METAINFORMANTError)
        assert isinstance(error, Exception)

    def test_io_errors_module_hierarchy(self) -> None:
        """Test io.errors module error hierarchy."""
        # FileNotFoundError inherits from io.errors.IOError
        file_error = io_errors.FileNotFoundError("file not found")
        assert isinstance(file_error, io_errors.IOError)
        assert isinstance(file_error, METAINFORMANTError)

        # CacheError inherits from io.errors.IOError
        cache_error = io_errors.CacheError("cache failed")
        assert isinstance(cache_error, io_errors.IOError)
        assert isinstance(cache_error, METAINFORMANTError)

        # DownloadError inherits from io.errors.IOError
        download_error = io_errors.DownloadError("download failed")
        assert isinstance(download_error, io_errors.IOError)
        assert isinstance(download_error, METAINFORMANTError)

    def test_catching_base_error_catches_all(self) -> None:
        """Test that catching METAINFORMANTError catches all custom errors."""
        errors_to_test = [
            ConfigError("config error"),
            IOError("io error"),
            ValidationError("validation error"),
            NetworkError("network error"),
            CacheError("cache error"),
            DownloadError("download error"),
            PipelineError("pipeline error"),
            DependencyError("dependency error"),
            ResourceError("resource error"),
            io_errors.IOError("io error from io module"),
            io_errors.FileNotFoundError("file not found"),
            io_errors.CacheError("cache error from io module"),
            io_errors.DownloadError("download error from io module"),
        ]

        for error in errors_to_test:
            try:
                raise error
            except METAINFORMANTError as e:
                assert e is error, f"Failed to catch {type(error).__name__} as METAINFORMANTError"
            else:
                pytest.fail(f"Did not catch {type(error).__name__}")


# ============================================================================
# Test retry_with_backoff
# ============================================================================


class TestRetryWithBackoff:
    """Test retry_with_backoff decorator functionality."""

    def test_retry_succeeds_on_first_attempt(self) -> None:
        """Test that successful functions don't retry."""
        call_count = 0

        @retry_with_backoff(max_attempts=3)
        def succeeds_immediately() -> str:
            nonlocal call_count
            call_count += 1
            return "success"

        result = succeeds_immediately()
        assert result == "success"
        assert call_count == 1

    def test_retry_succeeds_on_second_attempt(self) -> None:
        """Test that function retries and succeeds on second attempt."""
        call_count = 0

        @retry_with_backoff(max_attempts=3, initial_delay=0.01, exceptions=(ValueError,))
        def succeeds_on_second() -> str:
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise ValueError("first attempt fails")
            return "success"

        result = succeeds_on_second()
        assert result == "success"
        assert call_count == 2

    def test_retry_raises_on_max_attempts(self) -> None:
        """Test that function raises after max attempts."""
        call_count = 0

        @retry_with_backoff(max_attempts=3, initial_delay=0.01, exceptions=(ValueError,))
        def always_fails() -> str:
            nonlocal call_count
            call_count += 1
            raise ValueError(f"attempt {call_count} failed")

        with pytest.raises(ValueError, match="attempt 3 failed"):
            always_fails()

        assert call_count == 3

    def test_retry_only_on_specified_exceptions(self) -> None:
        """Test that retry only happens for specified exception types."""
        call_count = 0

        @retry_with_backoff(max_attempts=3, initial_delay=0.01, exceptions=(NetworkError,))
        def fails_with_wrong_exception() -> str:
            nonlocal call_count
            call_count += 1
            raise ValueError("not a network error")

        # Should raise immediately without retrying
        with pytest.raises(ValueError, match="not a network error"):
            fails_with_wrong_exception()

        assert call_count == 1

    def test_retry_with_multiple_exception_types(self) -> None:
        """Test retry with multiple exception types."""
        call_count = 0

        @retry_with_backoff(max_attempts=3, initial_delay=0.01, exceptions=(NetworkError, IOError))
        def fails_with_various_errors() -> str:
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise NetworkError("network error")
            elif call_count == 2:
                raise IOError("io error")
            return "success"

        result = fails_with_various_errors()
        assert result == "success"
        assert call_count == 3

    def test_retry_backoff_timing(self) -> None:
        """Test that backoff delays are applied correctly."""
        call_count = 0
        start_times: list[float] = []

        @retry_with_backoff(max_attempts=3, initial_delay=0.05, backoff_factor=2.0, exceptions=(ValueError,))
        def fails_twice() -> str:
            nonlocal call_count
            call_count += 1
            start_times.append(time.time())
            if call_count < 3:
                raise ValueError(f"attempt {call_count}")
            return "success"

        result = fails_twice()
        assert result == "success"
        assert call_count == 3

        # Check that delays increased
        # First delay: ~0.05s, second delay: ~0.1s
        if len(start_times) == 3:
            delay1 = start_times[1] - start_times[0]
            delay2 = start_times[2] - start_times[1]
            assert delay1 >= 0.04, f"First delay too short: {delay1}"
            assert delay2 >= 0.08, f"Second delay too short: {delay2}"
            assert delay2 > delay1, "Backoff not increasing"

    def test_retry_max_delay_cap(self) -> None:
        """Test that delays are capped at max_delay."""
        call_count = 0

        @retry_with_backoff(
            max_attempts=5, initial_delay=1.0, backoff_factor=10.0, max_delay=0.1, exceptions=(ValueError,)
        )
        def fails_multiple_times() -> str:
            nonlocal call_count
            call_count += 1
            if call_count < 5:
                raise ValueError(f"attempt {call_count}")
            return "success"

        start = time.time()
        result = fails_multiple_times()
        elapsed = time.time() - start

        assert result == "success"
        assert call_count == 5
        # With max_delay=0.1, 4 delays should be ~0.4s total (not 10+100+1000+...)
        assert elapsed < 1.0, f"Took too long: {elapsed}s, max_delay not working"


# ============================================================================
# Test error_context
# ============================================================================


class TestErrorContext:
    """Test error_context context manager."""

    def test_error_context_adds_message(self) -> None:
        """Test that error_context adds context to error message."""
        with pytest.raises(ValueError, match="Failed to process: original error"):
            with error_context("Failed to process"):
                raise ValueError("original error")

    def test_error_context_preserves_traceback(self) -> None:
        """Test that error_context preserves the original traceback."""

        def inner_function() -> None:
            raise ValueError("inner error")

        try:
            with error_context("Outer context"):
                inner_function()
        except ValueError as e:
            # Check that traceback includes inner_function
            import traceback

            tb_str = "".join(traceback.format_tb(e.__traceback__))
            assert "inner_function" in tb_str
            assert "inner error" in str(e)
            assert "Outer context" in str(e)

    def test_error_context_preserves_exception_type(self) -> None:
        """Test that error_context preserves the original exception type."""
        with pytest.raises(NetworkError):
            with error_context("Network operation failed"):
                raise NetworkError("connection timeout")

        with pytest.raises(ConfigError):
            with error_context("Config loading failed"):
                raise ConfigError("invalid yaml")

    def test_error_context_no_error(self) -> None:
        """Test that error_context doesn't interfere when no error occurs."""
        result = None
        with error_context("This should not interfere"):
            result = "success"

        assert result == "success"

    def test_error_context_nested(self) -> None:
        """Test nested error_context managers."""
        with pytest.raises(IOError, match="Outer context.*Inner context.*actual error"):
            with error_context("Outer context"):
                with error_context("Inner context"):
                    raise IOError("actual error")

    def test_error_context_no_reraise(self) -> None:
        """Test error_context with reraise=False."""
        with pytest.raises(ValueError, match="^original error$"):
            with error_context("This context should not be added", reraise=False):
                raise ValueError("original error")


# ============================================================================
# Test safe_execute
# ============================================================================


class TestSafeExecute:
    """Test safe_execute utility function."""

    def test_safe_execute_returns_result_on_success(self) -> None:
        """Test that safe_execute returns the function result on success."""

        def succeeds(x: int, y: int) -> int:
            return x + y

        result = safe_execute(succeeds, 2, 3)
        assert result == 5

    def test_safe_execute_returns_default_on_error(self) -> None:
        """Test that safe_execute returns default value on error."""

        def fails() -> str:
            raise ValueError("error")

        result = safe_execute(fails, default="fallback")
        assert result == "fallback"

    def test_safe_execute_with_none_default(self) -> None:
        """Test that safe_execute returns None by default."""

        def fails() -> str:
            raise ValueError("error")

        result = safe_execute(fails)
        assert result is None

    def test_safe_execute_with_kwargs(self) -> None:
        """Test that safe_execute works with keyword arguments."""

        def func_with_kwargs(a: int, b: int = 5) -> int:
            return a * b

        result = safe_execute(func_with_kwargs, 3, b=7)
        assert result == 21

    def test_safe_execute_catches_all_exceptions(self) -> None:
        """Test that safe_execute catches various exception types."""

        def raises_network_error() -> str:
            raise NetworkError("network error")

        def raises_io_error() -> str:
            raise IOError("io error")

        def raises_validation_error() -> str:
            raise ValidationError("validation error")

        assert safe_execute(raises_network_error, default="default") == "default"
        assert safe_execute(raises_io_error, default="default") == "default"
        assert safe_execute(raises_validation_error, default="default") == "default"

    def test_safe_execute_with_complex_return_types(self) -> None:
        """Test safe_execute with complex return types."""

        def returns_dict() -> dict[str, Any]:
            return {"key": "value", "count": 42}

        def fails_dict() -> dict[str, Any]:
            raise ValueError("error")

        result_success = safe_execute(returns_dict)
        assert result_success == {"key": "value", "count": 42}

        result_fail = safe_execute(fails_dict, default={})
        assert result_fail == {}


# ============================================================================
# Test Validation Functions
# ============================================================================


class TestValidationFunctions:
    """Test validate_not_none and validate_type functions."""

    def test_validate_not_none_success(self) -> None:
        """Test that validate_not_none passes for non-None values."""
        validate_not_none(42)
        validate_not_none("string")
        validate_not_none([])
        validate_not_none({})
        validate_not_none(0)
        validate_not_none(False)

    def test_validate_not_none_failure(self) -> None:
        """Test that validate_not_none raises for None."""
        with pytest.raises(ValidationError, match="value cannot be None"):
            validate_not_none(None)

    def test_validate_not_none_custom_name(self) -> None:
        """Test validate_not_none with custom name."""
        with pytest.raises(ValidationError, match="my_parameter cannot be None"):
            validate_not_none(None, name="my_parameter")

    def test_validate_type_success_single_type(self) -> None:
        """Test that validate_type passes for correct types."""
        validate_type(42, int)
        validate_type("string", str)
        validate_type([1, 2, 3], list)
        validate_type({"key": "value"}, dict)

    def test_validate_type_success_multiple_types(self) -> None:
        """Test that validate_type passes for tuple of types."""
        validate_type(42, (int, str))
        validate_type("string", (int, str))
        validate_type([1, 2], (list, tuple))

    def test_validate_type_failure(self) -> None:
        """Test that validate_type raises for incorrect types."""
        with pytest.raises(ValidationError, match="value must be of type int, got str"):
            validate_type("string", int)

        with pytest.raises(ValidationError, match="value must be of type list, got dict"):
            validate_type({}, list)

    def test_validate_type_failure_multiple_types(self) -> None:
        """Test that validate_type raises with multiple expected types."""
        with pytest.raises(ValidationError, match="value must be of type int, str, got list"):
            validate_type([1, 2], (int, str))

    def test_validate_type_custom_name(self) -> None:
        """Test validate_type with custom name."""
        with pytest.raises(ValidationError, match="my_param must be of type str, got int"):
            validate_type(42, str, name="my_param")


# ============================================================================
# Integration Tests
# ============================================================================


class TestErrorIntegration:
    """Integration tests combining multiple error utilities."""

    def test_retry_with_validation(self) -> None:
        """Test retry with validation errors."""
        call_count = 0

        @retry_with_backoff(max_attempts=3, initial_delay=0.01, exceptions=(ValidationError,))
        def validate_and_process(value: Any) -> str:
            nonlocal call_count
            call_count += 1
            validate_not_none(value)
            validate_type(value, str)
            if call_count < 2:
                raise ValidationError("temporary validation error")
            return f"processed: {value}"

        result = validate_and_process("test")
        assert result == "processed: test"
        assert call_count == 2

    def test_error_context_with_safe_execute(self) -> None:
        """Test combining error_context with safe_execute."""

        def risky_operation() -> str:
            with error_context("Risky operation"):
                raise ValueError("something went wrong")

        result = safe_execute(risky_operation, default="safe fallback")
        assert result == "safe fallback"

    def test_hierarchical_error_catching(self) -> None:
        """Test catching errors at different hierarchy levels."""
        # Catch at NetworkError level
        try:
            raise DownloadError("download failed")
        except NetworkError as e:
            assert isinstance(e, DownloadError)
            assert isinstance(e, NetworkError)
            assert isinstance(e, METAINFORMANTError)

        # Catch at METAINFORMANTError level
        try:
            raise DownloadError("download failed")
        except METAINFORMANTError as e:
            assert isinstance(e, DownloadError)

    def test_cross_module_error_catching(self) -> None:
        """Test catching errors from different modules using base class."""
        errors_from_different_modules = [
            NetworkError("from utils"),
            IOError("from utils"),
            io_errors.IOError("from io module"),
            io_errors.FileNotFoundError("from io module"),
        ]

        caught_count = 0
        for error in errors_from_different_modules:
            try:
                raise error
            except METAINFORMANTError:
                caught_count += 1

        assert caught_count == len(errors_from_different_modules)

    def test_error_context_preserves_hierarchy(self) -> None:
        """Test that error_context preserves exception hierarchy."""
        try:
            with error_context("Download operation"):
                raise DownloadError("network timeout")
        except DownloadError as e:
            assert isinstance(e, NetworkError)
            assert isinstance(e, METAINFORMANTError)
            assert "Download operation" in str(e)
            assert "network timeout" in str(e)
        except Exception:
            pytest.fail("Should have caught as DownloadError")
