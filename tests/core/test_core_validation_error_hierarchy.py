"""ValidationError is a METAINFORMANT error and a ValueError.

Callers that follow the standard library convention catch invalid arguments with
``except ValueError``; package code raises ValidationError for the same mistakes.
"""

from __future__ import annotations

import pytest

from metainformant.core.utils.errors import METAINFORMANTError, ValidationError


def test_validation_error_is_both_a_package_error_and_a_value_error() -> None:
    assert issubclass(ValidationError, METAINFORMANTError)
    assert issubclass(ValidationError, ValueError)


def test_except_value_error_catches_validation_error() -> None:
    with pytest.raises(ValueError, match="bad flavor"):
        raise ValidationError("bad flavor")


def test_other_package_errors_are_not_value_errors() -> None:
    from metainformant.core.utils.errors import ConfigError, NetworkError

    assert not issubclass(ConfigError, ValueError)
    assert not issubclass(NetworkError, ValueError)
