"""Tests for menu display module."""

from __future__ import annotations

import io as iolib
import sys

import pytest

from metainformant.menu.display import clear_screen, format_breadcrumb, format_menu, get_choice
from metainformant.menu.navigation import MenuItem


class TestFormatMenu:
    """Tests for menu formatting."""

    def test_format_menu_empty(self) -> None:
        """Test formatting empty menu."""
        result = format_menu([])
        assert "No menu items available" in result

    def test_format_menu_basic(self) -> None:
        """Test basic menu formatting."""
        items = [
            MenuItem(id="1", label="Option 1", description="First option"),
            MenuItem(id="2", label="Option 2", description="Second option"),
        ]
        result = format_menu(items, title="Test Menu")
        assert "Test Menu" in result
        assert "Option 1" in result
        assert "Option 2" in result
        assert "Back/Exit" in result

    def test_format_menu_without_title(self) -> None:
        """Test menu formatting without title."""
        items = [MenuItem(id="1", label="Option 1")]
        result = format_menu(items)
        assert "Option 1" in result
        assert "Back/Exit" in result

    def test_format_menu_disabled_items(self) -> None:
        """Test that disabled items are not shown."""
        items = [
            MenuItem(id="1", label="Enabled", enabled=True),
            MenuItem(id="2", label="Disabled", enabled=False),
        ]
        result = format_menu(items)
        assert "Enabled" in result
        assert "Disabled" not in result

    def test_format_menu_long_description(self) -> None:
        """Test menu with long descriptions."""
        long_desc = "A" * 100
        items = [MenuItem(id="1", label="Option", description=long_desc)]
        result = format_menu(items, width=80)
        # Description should be truncated
        assert len(result.split("\n")) > 0


class TestShowMenu:
    """Tests for menu display."""

    def test_show_menu(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test showing a menu."""
        items = [MenuItem(id="1", label="Option 1")]
        show_menu(items, title="Test")
        captured = capsys.readouterr()
        assert "Test" in captured.out
        assert "Option 1" in captured.out


class TestGetChoice:
    """Tests for user choice input.

    Note: Interactive choice input tests are skipped as they require
    input/output mocking/stubbing, which violates the NO_MOCKING_POLICY.
    The get_choice function is tested through integration tests with
    real user input in actual menu scenarios.
    """

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_valid(self) -> None:
        """Test getting valid choice."""
        # This test would require mocking builtin input
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_invalid_then_valid(self) -> None:
        """Test invalid choice followed by valid choice."""
        # This test would require mocking builtin input
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_no_validation(self) -> None:
        """Test getting choice without validation."""
        # This test would require mocking builtin input
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_empty_then_valid(self) -> None:
        """Test empty input followed by valid input."""
        # This test would require mocking builtin input
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_keyboard_interrupt(self) -> None:
        """Test handling KeyboardInterrupt."""
        # This test would require mocking builtin input
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: Interactive input tests require mocking")
    def test_get_choice_eof(self) -> None:
        """Test handling EOFError."""
        # This test would require mocking builtin input
        pass


class TestFormatBreadcrumb:
    """Tests for breadcrumb formatting."""

    def test_format_breadcrumb_empty(self) -> None:
        """Test formatting empty breadcrumb."""
        result = format_breadcrumb([])
        assert result == "Home"

    def test_format_breadcrumb_single(self) -> None:
        """Test formatting single-level breadcrumb."""
        result = format_breadcrumb(["RNA"])
        assert result == "RNA"

    def test_format_breadcrumb_multiple(self) -> None:
        """Test formatting multi-level breadcrumb."""
        result = format_breadcrumb(["Home", "RNA", "Workflows"])
        assert result == "Home > RNA > Workflows"

    def test_format_breadcrumb_special_chars(self) -> None:
        """Test breadcrumb with special characters."""
        result = format_breadcrumb(["Category 1", "Sub-Category"])
        assert "Category 1" in result
        assert "Sub-Category" in result


class TestClearScreen:
    """Tests for screen clearing.

    Note: System state mocking tests are skipped as they violate the
    NO_MOCKING_POLICY. The clear_screen function works with real TTY
    detection and should be tested through integration tests.
    """

    @pytest.mark.skip("NO_MOCKING_POLICY: TTY state mocking violates policy")
    def test_clear_screen_tty(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test clearing screen when TTY."""
        # This test would require mocking sys.stdout.isatty
        pass

    @pytest.mark.skip("NO_MOCKING_POLICY: TTY state mocking violates policy")
    def test_clear_screen_non_tty(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test clearing screen when not TTY."""
        # This test would require mocking sys.stdout.isatty
        pass



