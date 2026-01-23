"""Validation utilities for METAINFORMANT.

Provides type validators, range validators, path validators, and schema
validation for configuration and data validation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from metainformant.core.utils.errors import ValidationError
from metainformant.core.io.paths import expand_and_resolve, is_within
from metainformant.core.io.io import load_json


def validate_type(value: Any, expected_type: type | tuple[type, ...], name: str = "value") -> None:
    """Validate that a value is of expected type.

    Args:
        value: Value to validate
        expected_type: Expected type or tuple of types
        name: Name of value for error message

    Raises:
        ValidationError: If value is not of expected type

    Example:
        validate_type(42, int, "age")
    """
    if not isinstance(value, expected_type):
        if isinstance(expected_type, tuple):
            type_names = ", ".join(t.__name__ for t in expected_type)
        else:
            type_names = expected_type.__name__
        raise ValidationError(f"{name} must be of type {type_names}, got {type(value).__name__}")


def validate_range(
    value: float | int,
    min_val: float | int | None = None,
    max_val: float | int | None = None,
    name: str = "value",
) -> None:
    """Validate that a numeric value is within a range.

    Args:
        value: Numeric value to validate
        min_val: Minimum allowed value (None for no minimum)
        max_val: Maximum allowed value (None for no maximum)
        name: Name of value for error message

    Raises:
        ValidationError: If value is outside range

    Example:
        validate_range(0.5, min_val=0.0, max_val=1.0, name="probability")
    """
    if min_val is not None and value < min_val:
        raise ValidationError(f"{name} must be >= {min_val}, got {value}")
    if max_val is not None and value > max_val:
        raise ValidationError(f"{name} must be <= {max_val}, got {value}")


def validate_path_exists(path: str | Path, name: str = "path") -> Path:
    """Validate that a path exists.

    Args:
        path: Path to validate
        name: Name of path for error message

    Returns:
        Resolved Path object

    Raises:
        ValidationError: If path does not exist

    Example:
        file_path = validate_path_exists("data/file.txt", "input_file")
    """
    p = Path(path)
    if not p.exists():
        raise ValidationError(f"{name} does not exist: {path}")
    return p.expanduser().resolve()


def validate_path_is_file(path: str | Path, name: str = "path") -> Path:
    """Validate that a path exists and is a file.

    Args:
        path: Path to validate
        name: Name of path for error message

    Returns:
        Resolved Path object

    Raises:
        ValidationError: If path does not exist or is not a file

    Example:
        file_path = validate_path_is_file("data/file.txt", "input_file")
    """
    p = validate_path_exists(path, name)
    if not p.is_file():
        raise ValidationError(f"{name} is not a file: {path}")
    return p


def validate_path_is_dir(path: str | Path, name: str = "path") -> Path:
    """Validate that a path exists and is a directory.

    Args:
        path: Path to validate
        name: Name of path for error message

    Returns:
        Resolved Path object

    Raises:
        ValidationError: If path does not exist or is not a directory

    Example:
        dir_path = validate_path_is_dir("output/", "output_dir")
    """
    p = validate_path_exists(path, name)
    if not p.is_dir():
        raise ValidationError(f"{name} is not a directory: {path}")
    return p


def validate_path_within(parent: str | Path, path: str | Path, name: str = "path") -> Path:
    """Validate that a path is within a parent directory (security check).

    Args:
        parent: Parent directory path
        path: Path to validate
        name: Name of path for error message

    Returns:
        Resolved Path object

    Raises:
        ValidationError: If path is not within parent directory

    Example:
        safe_path = validate_path_within("/allowed/dir", user_path, "user_path")
    """
    p = expand_and_resolve(path)
    parent_path = expand_and_resolve(parent)

    if not is_within(p, parent_path):
        raise ValidationError(f"{name} must be within {parent}, got {path}")
    return p


def validate_not_none(value: Any, name: str = "value") -> None:
    """Validate that a value is not None.

    Args:
        value: Value to validate
        name: Name of value for error message

    Raises:
        ValidationError: If value is None

    Example:
        validate_not_none(config, "config")
    """
    if value is None:
        raise ValidationError(f"{name} cannot be None")


def validate_not_empty(value: str | list | dict, name: str = "value") -> None:
    """Validate that a value is not empty.

    Args:
        value: Value to validate (string, list, or dict)
        name: Name of value for error message

    Raises:
        ValidationError: If value is empty

    Example:
        validate_not_empty(items, "items")
    """
    if not value:
        raise ValidationError(f"{name} cannot be empty")


def validate_schema(data: dict[str, Any], schema: dict[str, Any], name: str = "data") -> None:
    """Validate data against a simple schema.

    Args:
        data: Data dictionary to validate
        schema: Schema dictionary defining expected structure
        name: Name of data for error message

    Raises:
        ValidationError: If data does not match schema

    Example:
        schema = {"name": str, "age": int, "email": str}
        validate_schema({"name": "John", "age": 30, "email": "john@example.com"}, schema)
    """
    validate_type(data, dict, name)

    # Check required fields
    for key, expected_type in schema.items():
        if key not in data:
            raise ValidationError(f"{name} missing required field: {key}")
        validate_type(data[key], expected_type, f"{name}.{key}")


def validate_json_schema(data: dict[str, Any], schema_path: str | Path) -> None:
    """Validate data against a JSON Schema file.

    Args:
        data: Data dictionary to validate
        schema_path: Path to JSON Schema file

    Raises:
        ValidationError: If data does not match schema
        FileNotFoundError: If schema file does not exist

    Example:
        validate_json_schema(config, "config/schema.json")
    """
    try:
        import jsonschema

        schema = load_json(schema_path)
        jsonschema.validate(instance=data, schema=schema)
    except ImportError:
        raise ValidationError("jsonschema package required for JSON Schema validation. Install with: uv add jsonschema")
    except json.JSONDecodeError as e:
        raise ValidationError(f"Invalid JSON schema file {schema_path}: {e}")
    except jsonschema.ValidationError as e:
        raise ValidationError(f"Data validation failed: {e.message}")


def validator(func: Callable[[Any], bool]) -> Callable[[Any], None]:
    """Decorator to create a validator function from a predicate.

    Args:
        func: Function that returns True if value is valid

    Returns:
        Validator function that raises ValidationError on invalid input

    Example:
        @validator
        def is_positive(x):
            return x > 0

        is_positive(5)  # OK
        is_positive(-1)  # Raises ValidationError
    """

    def wrapper(value: Any, name: str = "value") -> None:
        if not func(value):
            raise ValidationError(f"{name} failed validation: {value}")

    return wrapper
