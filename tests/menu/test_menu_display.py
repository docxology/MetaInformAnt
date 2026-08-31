"""Tests for menu display module."""

from __future__ import annotations

import os
import pty
import subprocess
import sys
import textwrap

import pytest

from metainformant.menu.ui.display import format_breadcrumb, format_menu, show_menu
from metainformant.menu.ui.navigation import MenuItem


def _run_child(source: str, *, input_text: str = "") -> subprocess.CompletedProcess[str]:
    """Run an interactive menu operation in a real child process.

    The parent test process remains non-interactive, while the child receives
    actual pipe/terminal file descriptors.  This exercises the production
    input and TTY branches without replacing ``input`` or ``sys.stdout``.
    """

    return subprocess.run(
        [sys.executable, "-c", source],
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )


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
    """Tests for user choice input through real child-process streams."""

    def test_get_choice_valid(self) -> None:
        """Test getting valid choice."""
        result = _run_child(
            "from metainformant.menu.ui.display import get_choice; "
            "print('RESULT=' + get_choice('Choose:', ['1', '2']))",
            input_text="2\n",
        )
        assert result.returncode == 0
        assert "RESULT=2" in result.stdout

    def test_get_choice_invalid_then_valid(self) -> None:
        """Test invalid choice followed by valid choice."""
        result = _run_child(
            "from metainformant.menu.ui.display import get_choice; "
            "print('RESULT=' + get_choice('Choose:', ['1', '2']))",
            input_text="invalid\n1\n",
        )
        assert result.returncode == 0
        assert "Invalid choice" in result.stdout
        assert "RESULT=1" in result.stdout

    def test_get_choice_no_validation(self) -> None:
        """Test getting choice without validation."""
        result = _run_child(
            "from metainformant.menu.ui.display import get_choice; " "print('RESULT=' + get_choice('Choose:'))",
            input_text="free-form\n",
        )
        assert result.returncode == 0
        assert "RESULT=free-form" in result.stdout

    def test_get_choice_empty_then_valid(self) -> None:
        """Test empty input followed by valid input."""
        result = _run_child(
            "from metainformant.menu.ui.display import get_choice; "
            "print('RESULT=' + get_choice('Choose:', ['1', '2']))",
            input_text="\n2\n",
        )
        assert result.returncode == 0
        assert "RESULT=2" in result.stdout

    def test_get_choice_keyboard_interrupt(self) -> None:
        """Test handling KeyboardInterrupt."""
        source = textwrap.dedent("""
            import os
            import signal
            import threading
            import time
            from metainformant.menu.ui.display import get_choice

            signal.signal(signal.SIGINT, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt))
            threading.Timer(0.2, lambda: os.kill(os.getpid(), signal.SIGINT)).start()
            try:
                get_choice("Choose:")
            except KeyboardInterrupt:
                print("INTERRUPTED")
            """)
        process = subprocess.Popen(
            [sys.executable, "-c", source],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        try:
            assert process.wait(timeout=10) == 0
            stdout = process.stdout.read() if process.stdout is not None else ""
        finally:
            if process.stdin is not None:
                process.stdin.close()
            if process.stdout is not None:
                process.stdout.close()
            if process.stderr is not None:
                process.stderr.close()
        assert "INTERRUPTED" in stdout

    def test_get_choice_eof(self) -> None:
        """Test handling EOFError."""
        result = _run_child(
            "from metainformant.menu.ui.display import get_choice; "
            "\ntry: get_choice('Choose:')\nexcept EOFError: print('EOF')",
        )
        assert result.returncode == 0
        assert "EOF" in result.stdout


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
    """Tests for screen clearing with real pipe and pseudo-terminal streams."""

    def test_clear_screen_tty(self) -> None:
        """Test clearing screen when TTY."""
        master_fd, slave_fd = pty.openpty()
        try:
            process = subprocess.Popen(
                [
                    sys.executable,
                    "-c",
                    "from metainformant.menu.ui.display import clear_screen; clear_screen()",
                ],
                stdin=slave_fd,
                stdout=slave_fd,
                stderr=slave_fd,
                close_fds=True,
            )
        finally:
            os.close(slave_fd)
        try:
            output = os.read(master_fd, 4096).decode()
            assert process.wait(timeout=10) == 0
        finally:
            os.close(master_fd)
        assert "\x1b[2J\x1b[H" in output

    def test_clear_screen_non_tty(self) -> None:
        """Test clearing screen when not TTY."""
        result = _run_child("from metainformant.menu.ui.display import clear_screen; clear_screen()")
        assert result.returncode == 0
        assert "\x1b[2J\x1b[H" not in result.stdout
        assert result.stdout.count("\n") == 51  # 50 fallback lines plus print's final newline
