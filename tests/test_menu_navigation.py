"""Tests for menu navigation module."""

from __future__ import annotations

import pytest

from metainformant.menu.navigation import Menu, MenuHistory, MenuItem, MenuSystem, get_current_menu, go_back, navigate_to_submenu


class TestMenuItem:
    """Tests for MenuItem class."""

    def test_menu_item_creation(self) -> None:
        """Test creating a menu item."""
        item = MenuItem(id="test", label="Test", description="Test item")
        assert item.id == "test"
        assert item.label == "Test"
        assert item.description == "Test item"
        assert item.enabled is True

    def test_menu_item_disabled(self) -> None:
        """Test creating disabled menu item."""
        item = MenuItem(id="test", label="Test", enabled=False)
        assert item.enabled is False


class TestMenu:
    """Tests for Menu class."""

    def test_menu_creation(self) -> None:
        """Test creating a menu."""
        items = [MenuItem(id="1", label="Option 1")]
        menu = Menu(id="test", title="Test Menu", items=items)
        assert menu.id == "test"
        assert menu.title == "Test Menu"
        assert len(menu.items) == 1
        assert menu.parent_id is None

    def test_menu_with_parent(self) -> None:
        """Test menu with parent ID."""
        menu = Menu(id="child", title="Child", items=[], parent_id="parent")
        assert menu.parent_id == "parent"


class TestMenuHistory:
    """Tests for MenuHistory class."""

    def test_history_push(self) -> None:
        """Test pushing to history."""
        history = MenuHistory()
        history.push("menu1", "Menu 1")
        assert len(history.path) == 1
        assert len(history.menu_ids) == 1
        assert history.path[0] == "Menu 1"
        assert history.menu_ids[0] == "menu1"

    def test_history_pop(self) -> None:
        """Test popping from history."""
        history = MenuHistory()
        history.push("menu1", "Menu 1")
        history.push("menu2", "Menu 2")
        
        result = history.pop()
        assert result == ("menu2", "Menu 2")
        assert len(history.path) == 1
        assert len(history.menu_ids) == 1

    def test_history_pop_empty(self) -> None:
        """Test popping from empty history."""
        history = MenuHistory()
        result = history.pop()
        assert result is None

    def test_history_get_path(self) -> None:
        """Test getting breadcrumb path."""
        history = MenuHistory()
        history.push("menu1", "Menu 1")
        history.push("menu2", "Menu 2")
        
        path = history.get_path()
        assert path == ["Menu 1", "Menu 2"]
        # Should be a copy, not reference
        path.append("Menu 3")
        assert len(history.path) == 2

    def test_history_clear(self) -> None:
        """Test clearing history."""
        history = MenuHistory()
        history.push("menu1", "Menu 1")
        history.clear()
        assert len(history.path) == 0
        assert len(history.menu_ids) == 0


class TestMenuSystem:
    """Tests for MenuSystem class."""

    def test_menu_system_creation(self) -> None:
        """Test creating menu system."""
        menus = {
            "root": Menu(
                id="root",
                title="Root",
                items=[MenuItem(id="1", label="Option 1")],
            )
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        assert system.current_menu_id == "root"
        assert len(system.history.path) == 1

    def test_get_current_menu(self) -> None:
        """Test getting current menu."""
        menus = {
            "root": Menu(id="root", title="Root", items=[]),
            "sub": Menu(id="sub", title="Sub", items=[]),
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        menu = system.get_current_menu()
        assert menu.id == "root"

    def test_get_current_menu_invalid(self) -> None:
        """Test getting invalid current menu."""
        menus = {"root": Menu(id="root", title="Root", items=[])}
        system = MenuSystem(menus=menus, current_menu_id="invalid")
        with pytest.raises(KeyError):
            system.get_current_menu()

    def test_navigate_to(self) -> None:
        """Test navigating to a menu."""
        menus = {
            "root": Menu(id="root", title="Root", items=[]),
            "sub": Menu(id="sub", title="Sub", items=[]),
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        success = system.navigate_to("sub")
        assert success is True
        assert system.current_menu_id == "sub"
        assert len(system.history.path) == 2

    def test_navigate_to_invalid(self) -> None:
        """Test navigating to invalid menu."""
        menus = {"root": Menu(id="root", title="Root", items=[])}
        system = MenuSystem(menus=menus, current_menu_id="root")
        success = system.navigate_to("invalid")
        assert success is False
        assert system.current_menu_id == "root"

    def test_go_back(self) -> None:
        """Test going back in navigation."""
        menus = {
            "root": Menu(id="root", title="Root", items=[]),
            "sub": Menu(id="sub", title="Sub", items=[]),
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        system.navigate_to("sub")
        success = system.go_back()
        assert success is True
        assert system.current_menu_id == "root"

    def test_go_back_at_root(self) -> None:
        """Test going back when at root."""
        menus = {"root": Menu(id="root", title="Root", items=[])}
        system = MenuSystem(menus=menus, current_menu_id="root")
        success = system.go_back()
        assert success is False
        assert system.current_menu_id == "root"


class TestNavigationFunctions:
    """Tests for navigation module functions."""

    def test_navigate_to_submenu(self) -> None:
        """Test navigate_to_submenu function."""
        menus = {
            "root": Menu(id="root", title="Root", items=[]),
            "sub": Menu(id="sub", title="Sub", items=[]),
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        navigate_to_submenu(system, "sub")
        assert system.current_menu_id == "sub"

    def test_navigate_to_submenu_invalid(self) -> None:
        """Test navigating to invalid submenu."""
        menus = {"root": Menu(id="root", title="Root", items=[])}
        system = MenuSystem(menus=menus, current_menu_id="root")
        with pytest.raises(KeyError):
            navigate_to_submenu(system, "invalid")

    def test_go_back_function(self) -> None:
        """Test go_back function."""
        menus = {
            "root": Menu(id="root", title="Root", items=[]),
            "sub": Menu(id="sub", title="Sub", items=[]),
        }
        system = MenuSystem(menus=menus, current_menu_id="root")
        system.navigate_to("sub")
        success = go_back(system)
        assert success is True
        assert system.current_menu_id == "root"

    def test_get_current_menu_function(self) -> None:
        """Test get_current_menu function."""
        menus = {"root": Menu(id="root", title="Root", items=[])}
        system = MenuSystem(menus=menus, current_menu_id="root")
        menu = get_current_menu(system)
        assert menu.id == "root"



