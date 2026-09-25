"""Integration tests for the interactive menu application loop.

These tests exercise the real action protocol end to end: submenu navigation,
script execution via a real child interpreter, callable actions, and the
Back/Exit loop. Interactive runs are driven through a real child process with
piped stdin (no input patching), matching the policy in
``tests/REAL_IMPLEMENTATION_TESTING_POLICY.md``.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Callable

from metainformant.menu.app import MenuApp
from metainformant.menu.ui.navigation import Menu, MenuItem, MenuSystem

REPO_ROOT = Path(__file__).resolve().parents[2]

# Child processes must import this repository's src tree, not any installed
# copy of the package (see tests/conftest.py for the in-process equivalent).
CHILD_ENV = {**os.environ, "PYTHONPATH": str(REPO_ROOT / "src")}


def _menu_system(
    root_items: list[MenuItem], extra_menus: dict[str, Menu] | None = None
) -> MenuSystem:
    """Build a MenuSystem rooted at 'root'."""
    menus: dict[str, Menu] = {"root": Menu(id="root", title="Root", items=root_items)}
    menus.update(extra_menus or {})
    return MenuSystem(menus=menus, current_menu_id="root")


class TestSubmenuAction:
    """'submenu:<id>' actions drive navigation state."""

    def test_navigates_into_submenu(self) -> None:
        child = Menu(id="menu_cat", title="Category", items=[])
        item = MenuItem(id="open", label="Open category", action="submenu:menu_cat")
        system = _menu_system([item], {"menu_cat": child})
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) is None
        assert system.current_menu_id == "menu_cat"
        assert system.history.menu_ids == ["root", "menu_cat"]
        assert system.history.get_path() == ["Root", "Category"]

    def test_unknown_submenu_keeps_current_menu(self, capsys: Callable) -> None:
        item = MenuItem(id="open", label="Open", action="submenu:missing")
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        app.handle_action(item)

        assert system.current_menu_id == "root"
        assert system.history.menu_ids == ["root"]
        assert "Submenu not found: missing" in capsys.readouterr().out


class TestScriptAction:
    """'script:<path>' actions execute the script in a child interpreter."""

    def test_runs_python_script_with_side_effect(self, tmp_path: Path) -> None:
        marker = tmp_path / "marker.txt"
        script = tmp_path / "ok.py"
        script.write_text(
            f"from pathlib import Path\nPath({str(marker)!r}).write_text('ran')\n"
        )
        item = MenuItem(id="run", label="Run script", action=f"script:{script}")
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) == 0
        assert marker.read_text() == "ran"

    def test_reports_failing_script_exit_code(self, tmp_path: Path) -> None:
        script = tmp_path / "fail.py"
        script.write_text("import sys\nsys.exit(3)\n")
        item = MenuItem(id="run", label="Run script", action=f"script:{script}")
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) == 3

    def test_missing_script_returns_error_code(
        self, tmp_path: Path, capsys: Callable
    ) -> None:
        item = MenuItem(
            id="run", label="Run", action=f"script:{tmp_path / 'missing.py'}"
        )
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) == 1
        assert "Script not found" in capsys.readouterr().out


class TestCallableAction:
    """Callable actions are invoked and keep the loop alive on errors."""

    def test_invokes_callable_and_prints_result(
        self, tmp_path: Path, capsys: Callable
    ) -> None:
        marker = tmp_path / "called.txt"

        def action() -> str:
            marker.write_text("called")
            return "callable-ok"

        item = MenuItem(id="call", label="Call action", action=action)
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) is None
        assert marker.read_text() == "called"
        assert "callable-ok" in capsys.readouterr().out

    def test_callable_exception_does_not_crash(self, capsys: Callable) -> None:
        def action() -> str:
            raise ValueError("boom")

        item = MenuItem(id="call", label="Call action", action=action)
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) is None
        assert "Error executing action: boom" in capsys.readouterr().out

    def test_empty_action_reports_no_action(self, capsys: Callable) -> None:
        item = MenuItem(id="idle", label="Idle entry")
        system = _menu_system([item])
        app = MenuApp(menu_system=system, clear_display=False)

        assert app.handle_action(item) is None
        assert "No action defined for: Idle entry" in capsys.readouterr().out


def _write_driver(tmp_path: Path, body: str) -> Path:
    """Write a child-process driver that builds menus and runs MenuApp."""
    driver = tmp_path / "menu_driver.py"
    driver.write_text(textwrap.dedent(body))
    return driver


LOOP_DRIVER = """
    import sys
    from pathlib import Path

    from metainformant.menu.app import MenuApp
    from metainformant.menu.ui.navigation import Menu, MenuItem, MenuSystem

    def callable_action():
        Path(sys.argv[2]).write_text("callable")
        return "callable-action-ok"

    root = Menu(
        id="root",
        title="Root",
        items=[
            MenuItem(id="sub", label="Open submenu", action="submenu:menu_sub"),
            MenuItem(id="run", label="Run script", action="script:{script}"),
            MenuItem(id="call", label="Call action", action=callable_action),
        ],
    )
    sub = Menu(id="menu_sub", title="Sub", items=[])
    system = MenuSystem(menus={{"root": root, "menu_sub": sub}}, current_menu_id="root")
    app = MenuApp(menu_system=system, clear_display=False)
    app.run()
    print("ENDED_AT", system.current_menu_id)
    print("HISTORY", system.history.menu_ids)
"""


class TestInteractiveLoop:
    """End-to-end loop driven through a real child process with piped stdin."""

    def test_loop_dispatches_all_action_kinds_with_back_and_exit(
        self, tmp_path: Path
    ) -> None:
        script_marker = tmp_path / "script_marker.txt"
        callable_marker = tmp_path / "callable_marker.txt"
        script = tmp_path / "ok.py"
        script.write_text(
            f"from pathlib import Path\nPath({str(script_marker)!r}).write_text('ran')\n"
        )
        driver = _write_driver(tmp_path, LOOP_DRIVER.format(script=str(script)))

        result = subprocess.run(
            [sys.executable, str(driver), str(script), str(callable_marker)],
            input="1\n0\n2\n3\n0\n",
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(REPO_ROOT),
            env=CHILD_ENV,
        )

        assert result.returncode == 0, result.stderr
        # "1" navigated into the submenu, "0" navigated back to the root.
        assert "ENDED_AT root" in result.stdout
        assert "HISTORY ['root']" in result.stdout
        # The script action really ran in a child interpreter...
        assert script_marker.read_text() == "ran"
        assert "Script completed successfully" in result.stdout
        # ...and the callable action really ran.
        assert callable_marker.read_text() == "callable"
        assert "callable-action-ok" in result.stdout

    def test_invalid_choice_is_reported_then_exit_works(self, tmp_path: Path) -> None:
        driver = _write_driver(
            tmp_path,
            """
            import sys

            from metainformant.menu.app import MenuApp
            from metainformant.menu.ui.navigation import Menu, MenuItem, MenuSystem

            root = Menu(id="root", title="Root", items=[MenuItem(id="idle", label="Idle")])
            system = MenuSystem(menus={"root": root}, current_menu_id="root")
            app = MenuApp(menu_system=system, clear_display=False)
            app.run()
            print("ENDED_AT", system.current_menu_id)
            """,
        )

        result = subprocess.run(
            [sys.executable, str(driver)],
            input="not-a-number\n0\n",
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(REPO_ROOT),
            env=CHILD_ENV,
        )

        assert result.returncode == 0, result.stderr
        assert "Invalid choice: 'not-a-number'" in result.stdout
        assert "ENDED_AT root" in result.stdout

    def test_eof_exits_gracefully(self, tmp_path: Path) -> None:
        driver = _write_driver(
            tmp_path,
            """
            from metainformant.menu.app import MenuApp
            from metainformant.menu.ui.navigation import Menu, MenuItem, MenuSystem

            root = Menu(id="root", title="Root", items=[MenuItem(id="idle", label="Idle")])
            system = MenuSystem(menus={"root": root}, current_menu_id="root")
            app = MenuApp(menu_system=system, clear_display=False)
            app.run()
            print("LOOP_ENDED")
            """,
        )

        result = subprocess.run(
            [sys.executable, str(driver)],
            input="",
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(REPO_ROOT),
            env=CHILD_ENV,
        )

        assert result.returncode == 0, result.stderr
        assert "Goodbye!" in result.stdout
        assert "LOOP_ENDED" in result.stdout

    def test_run_menu_app_requires_root_menu(self, tmp_path: Path) -> None:
        """Convenience builder validates the root menu id."""
        import pytest

        from metainformant.menu.app import run_menu_app

        with pytest.raises(KeyError):
            run_menu_app({}, root_id="missing")
