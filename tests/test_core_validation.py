"""Tests for core validation utilities.

Tests validation functions following NO_MOCKING policy.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.core.data.validation import (
    validate_json_schema,
    validate_not_empty,
    validate_not_none,
    validate_path_exists,
    validate_path_is_dir,
    validate_path_is_file,
    validate_path_within,
    validate_range,
    validate_schema,
    validate_type,
    validator,
)


class TestValidateType:
    """Tests for validate_type function."""

    def test_validate_type_int(self):
        """Test validating integer type."""
        validate_type(42, int, "age")
        # Should not raise

    def test_validate_type_str(self):
        """Test validating string type."""
        validate_type("hello", str, "name")
        # Should not raise

    def test_validate_type_wrong_type(self):
        """Test that wrong type raises ValidationError."""
        with pytest.raises(ValidationError, match="age must be of type int"):
            validate_type("42", int, "age")

    def test_validate_type_tuple(self):
        """Test validating against tuple of types."""
        validate_type(42, (int, float), "number")
        validate_type(3.14, (int, float), "number")
        # Should not raise

    def test_validate_type_tuple_wrong(self):
        """Test that type not in tuple raises ValidationError."""
        with pytest.raises(ValidationError, match="value must be of type"):
            validate_type("hello", (int, float), "value")


class TestValidateRange:
    """Tests for validate_range function."""

    def test_validate_range_within(self):
        """Test validating value within range."""
        validate_range(5, min_val=0, max_val=10, name="value")
        # Should not raise

    def test_validate_range_min_only(self):
        """Test validating with only minimum."""
        validate_range(5, min_val=0, name="value")
        # Should not raise

    def test_validate_range_max_only(self):
        """Test validating with only maximum."""
        validate_range(5, max_val=10, name="value")
        # Should not raise

    def test_validate_range_below_min(self):
        """Test that value below minimum raises ValidationError."""
        with pytest.raises(ValidationError, match="value must be >= 10"):
            validate_range(5, min_val=10, name="value")

    def test_validate_range_above_max(self):
        """Test that value above maximum raises ValidationError."""
        with pytest.raises(ValidationError, match="value must be <= 5"):
            validate_range(10, max_val=5, name="value")


class TestValidatePath:
    """Tests for path validation functions."""

    def test_validate_path_exists(self, tmp_path):
        """Test validating existing path."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")
        result = validate_path_exists(test_file, "test_file")
        assert isinstance(result, Path)
        assert result.exists()

    def test_validate_path_exists_nonexistent(self, tmp_path):
        """Test that nonexistent path raises ValidationError."""
        nonexistent = tmp_path / "nonexistent.txt"
        with pytest.raises(ValidationError, match="path does not exist"):
            validate_path_exists(nonexistent, "path")

    def test_validate_path_is_file(self, tmp_path):
        """Test validating file path."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")
        result = validate_path_is_file(test_file, "test_file")
        assert result.is_file()

    def test_validate_path_is_file_directory(self, tmp_path):
        """Test that directory raises ValidationError when expecting file."""
        with pytest.raises(ValidationError, match="path is not a file"):
            validate_path_is_file(tmp_path, "path")

    def test_validate_path_is_dir(self, tmp_path):
        """Test validating directory path."""
        result = validate_path_is_dir(tmp_path, "test_dir")
        assert result.is_dir()

    def test_validate_path_is_dir_file(self, tmp_path):
        """Test that file raises ValidationError when expecting directory."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")
        with pytest.raises(ValidationError, match="path is not a directory"):
            validate_path_is_dir(test_file, "path")

    def test_validate_path_within(self, tmp_path):
        """Test validating path within parent directory."""
        parent = tmp_path / "parent"
        parent.mkdir()
        child = parent / "child.txt"
        child.write_text("test")
        result = validate_path_within(parent, child, "child")
        assert isinstance(result, Path)

    def test_validate_path_within_outside(self, tmp_path):
        """Test that path outside parent raises ValidationError."""
        parent = tmp_path / "parent"
        parent.mkdir()
        outside = tmp_path / "outside.txt"
        outside.write_text("test")
        with pytest.raises(ValidationError, match="path must be within"):
            validate_path_within(parent, outside, "path")


class TestValidateNotNone:
    """Tests for validate_not_none function."""

    def test_validate_not_none_valid(self):
        """Test validating non-None value."""
        validate_not_none("value", "name")
        validate_not_none(42, "number")
        validate_not_none([1, 2, 3], "list")
        # Should not raise

    def test_validate_not_none_none(self):
        """Test that None raises ValidationError."""
        with pytest.raises(ValidationError, match="value cannot be None"):
            validate_not_none(None, "value")


class TestValidateNotEmpty:
    """Tests for validate_not_empty function."""

    def test_validate_not_empty_valid(self):
        """Test validating non-empty values."""
        validate_not_empty("hello", "string")
        validate_not_empty([1, 2, 3], "list")
        validate_not_empty({"key": "value"}, "dict")
        # Should not raise

    def test_validate_not_empty_empty_string(self):
        """Test that empty string raises ValidationError."""
        with pytest.raises(ValidationError, match="value cannot be empty"):
            validate_not_empty("", "value")

    def test_validate_not_empty_empty_list(self):
        """Test that empty list raises ValidationError."""
        with pytest.raises(ValidationError, match="value cannot be empty"):
            validate_not_empty([], "value")

    def test_validate_not_empty_empty_dict(self):
        """Test that empty dict raises ValidationError."""
        with pytest.raises(ValidationError, match="value cannot be empty"):
            validate_not_empty({}, "value")


class TestValidateSchema:
    """Tests for validate_schema function."""

    def test_validate_schema_valid(self):
        """Test validating data against schema."""
        schema = {"name": str, "age": int, "email": str}
        data = {"name": "John", "age": 30, "email": "john@example.com"}
        validate_schema(data, schema, "person")
        # Should not raise

    def test_validate_schema_missing_field(self):
        """Test that missing field raises ValidationError."""
        schema = {"name": str, "age": int}
        data = {"name": "John"}
        with pytest.raises(ValidationError, match="missing required field"):
            validate_schema(data, schema, "person")

    def test_validate_schema_wrong_type(self):
        """Test that wrong type raises ValidationError."""
        schema = {"name": str, "age": int}
        data = {"name": "John", "age": "30"}
        with pytest.raises(ValidationError, match="must be of type"):
            validate_schema(data, schema, "person")


class TestValidateJsonSchema:
    """Tests for validate_json_schema function."""

    def test_validate_json_schema_missing_package(self, tmp_path):
        """Test that missing jsonschema package raises ValidationError."""
        schema_file = tmp_path / "schema.json"
        schema_file.write_text('{"type": "object"}')
        data = {"key": "value"}

        try:
            import jsonschema  # noqa: F401

            # If jsonschema is available, test should work
            validate_json_schema(data, schema_file)
        except ImportError:
            # If not available, should raise ValidationError
            with pytest.raises(ValidationError, match="jsonschema package required"):
                validate_json_schema(data, schema_file)


class TestValidatorDecorator:
    """Tests for validator decorator."""

    def test_validator_decorator(self):
        """Test validator decorator creates validator function."""

        @validator
        def is_positive(x):
            return x > 0

        # Valid value should not raise
        is_positive(5, "value")

        # Invalid value should raise ValidationError
        with pytest.raises(ValidationError, match="failed validation"):
            is_positive(-1, "value")

    def test_validator_decorator_custom_name(self):
        """Test validator with custom name parameter."""

        @validator
        def is_even(x):
            return x % 2 == 0

        is_even(4, "number")
        with pytest.raises(ValidationError, match="number failed validation"):
            is_even(3, "number")
