"""Interactive menu application loop.

Consumes the action protocol produced by
:mod:`metainformant.menu.core.discovery` (``"submenu:<id>"``, ``"script:<path>"``
or a callable) and drives end-to-end user interaction: menu rendering, choice
handling, submenu navigation, script execution and back/exit handling.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from metainformant.menu.core.executor import execute_script
from metainformant.menu.ui.display import (
    clear_screen,
    format_breadcrumb,
    get_choice,
    show_menu,
)
from metainformant.menu.ui.navigation import MenuItem, MenuSystem


@dataclass
class MenuApp:
    """Run menus from :mod:`metainformant.menu.core.discovery` and dispatch actions.

    Attributes:
        menu_system: Navigation state (menus, current menu, history).
        clear_display: Clear the terminal before each menu render.
    """

    menu_system: MenuSystem
    clear_display: bool = True

    def run(self) -> None:
        """Run the interactive loop until Back/Exit at the root menu or EOF.

        Choice ``0`` navigates back one level; at the root menu it exits the
        application. EOF (``Ctrl+D``) or ``Ctrl+C`` exits gracefully.
        """
        try:
            while True:
                self.display()
                try:
                    choice = get_choice("Select option: ")
                except (KeyboardInterrupt, EOFError):
                    print("\nGoodbye!")
                    return

                if choice == "0":
                    if len(self.menu_system.history.menu_ids) <= 1:
                        print("Goodbye!")
                        return
                    self.menu_system.go_back()
                    continue

                menu = self.menu_system.get_current_menu()
                try:
                    index = int(choice)
                except ValueError:
                    print(
                        f"Invalid choice: {choice!r}. Enter a menu number or 0 to go back."
                    )
                    continue
                if not 1 <= index <= len(menu.items):
                    print(
                        f"Invalid choice: {choice}. Enter a number between 1 and {len(menu.items)}."
                    )
                    continue

                item = menu.items[index - 1]
                if not item.enabled:
                    print(f"{item.label} is currently disabled.")
                    continue
                self.handle_action(item)
        except KeyboardInterrupt:
            print("\nGoodbye!")

    def display(self) -> None:
        """Render the breadcrumb and the current menu."""
        if self.clear_display:
            clear_screen()
        breadcrumb = format_breadcrumb(self.menu_system.history.get_path())
        if breadcrumb:
            print(breadcrumb)
        menu = self.menu_system.get_current_menu()
        show_menu(menu.items, menu.title)

    def handle_action(self, item: MenuItem) -> int | None:
        """Resolve and dispatch a single menu item action.

        Args:
            item: The selected menu item.

        Returns:
            The script exit code for ``script:`` actions, otherwise ``None``.
        """
        action = item.action
        if isinstance(action, str):
            if not action.strip():
                print(f"No action defined for: {item.label}")
                return None
            if action.startswith("submenu:"):
                submenu_id = action.split(":", 1)[1]
                if not self.menu_system.navigate_to(submenu_id):
                    print(f"Submenu not found: {submenu_id}")
                return None
            if action.startswith("script:"):
                script_path = Path(action.split(":", 1)[1])
                print(f"Executing: {script_path}")
                exit_code = execute_script(script_path)
                if exit_code == 0:
                    print("Script completed successfully")
                else:
                    print(f"Script exited with code: {exit_code}")
                return exit_code
            print(f"Unknown action: {action!r}")
            return None

        if callable(action):
            try:
                result = action()
            except Exception as exc:  # keep the loop alive on action errors
                print(f"Error executing action: {exc}")
                return None
            if isinstance(result, str):
                print(result)
            return None

        print(f"No action defined for: {item.label}")
        return None


def run_menu_app(
    menus: dict, root_id: str = "root", *, clear_display: bool = True
) -> None:
    """Build a :class:`MenuSystem` and run the interactive loop.

    Args:
        menus: Menu ID to :class:`~metainformant.menu.ui.navigation.Menu` mapping,
            typically from ``generate_menu_from_scripts``.
        root_id: ID of the root menu.
        clear_display: Clear the terminal before each menu render.
    """
    if root_id not in menus:
        raise KeyError(f"Root menu '{root_id}' not found")
    app = MenuApp(
        menu_system=MenuSystem(menus=menus, current_menu_id=root_id),
        clear_display=clear_display,
    )
    app.run()
