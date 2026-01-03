# AI Agents in Menu Module Development

This document outlines AI assistance in developing METAINFORMANT's menu system for script orchestration.

## AI Contributions

### Module Architecture Design
**Code Assistant Agent** (grok-code-fast-1) designed:
- Modular menu system with separate concerns (display, navigation, discovery, execution)
- Hierarchical menu structure with breadcrumb navigation
- Automatic script discovery from repository structure
- Script execution with argument prompting
- Integration with orchestrator script

### Core Components

**Code Assistant Agent** implemented:

#### Display Module (`display.py`)
- `format_menu()`: Menu formatting with consistent styling and width control
- `show_menu()`: Display menus to users with proper formatting
- `get_choice()`: User input validation with retry logic
- `format_breadcrumb()`: Breadcrumb path formatting for navigation
- `clear_screen()`: Cross-platform screen clearing with ANSI codes

#### Navigation Module (`navigation.py`)
- `MenuItem`: Dataclass for menu items with actions and metadata
- `Menu`: Menu container with items and parent relationships
- `MenuSystem`: Navigation state management with history tracking
- `MenuHistory`: Breadcrumb tracking with push/pop operations
- `navigate_to_submenu()`: Submenu navigation function
- `go_back()`: Back navigation with history management
- `get_current_menu()`: Current menu retrieval

#### Discovery Module (`discovery.py`)
- `discover_scripts()`: Recursive script discovery from scripts directory
- `extract_script_metadata()`: Metadata extraction from Python and bash scripts
- `categorize_script()`: Script categorization by directory structure
- `generate_menu_from_scripts()`: Menu structure generation from discovered scripts
- `ScriptInfo`: Script metadata dataclass
- Python AST parsing for docstring and argument extraction
- Bash script comment parsing for description extraction

#### Executor Module (`executor.py`)
- `execute_script()`: Generic script execution (Python or bash)
- `execute_python_script()`: Python script execution via subprocess
- `execute_bash_script()`: Bash script execution with proper environment
- `prompt_for_args()`: Interactive argument prompting for scripts
- `validate_script_executable()`: Script validation and existence checking

### Orchestrator Script

**Code Assistant Agent** implemented:
- `metainformant.sh`: Top-level bash orchestrator script
- Interactive menu launch functionality
- Direct script execution shortcuts
- Script discovery and listing
- Category-based menu display
- Environment setup and validation
- Help system and usage documentation

### Testing Framework

**Code Assistant Agent** implemented comprehensive test suites:

#### `test_menu_display.py`
- Menu formatting tests (empty, basic, with/without title)
- Disabled item handling
- Long description truncation
- User input validation
- Breadcrumb formatting
- Screen clearing (TTY and non-TTY)

#### `test_menu_navigation.py`
- MenuItem and Menu creation
- MenuHistory push/pop operations
- MenuSystem navigation state management
- Submenu navigation
- Back navigation
- Invalid menu handling

#### `test_menu_discovery.py`
- Script categorization by directory
- Python script metadata extraction
- Bash script metadata extraction
- Script discovery (Python and bash)
- Test file skipping
- Cache directory skipping
- Menu generation from scripts

#### `test_menu_executor.py`
- Script validation
- Python script execution
- Bash script execution
- Argument handling
- Error handling
- Exit code propagation

### Documentation

**Documentation Agent** created:
- `README.md`: Comprehensive module documentation with usage examples
- `AGENTS.md`: This file documenting AI contributions
- Integration guides and architecture diagrams
- API reference documentation

## Development Approach

- **Modular Design**: Separated concerns into display, navigation, discovery, and execution
- **Type Safety**: Comprehensive type hints throughout
- **Error Handling**: Graceful error handling with user-friendly messages
- **Testing**: Comprehensive test coverage for all components
- **Documentation**: Clear documentation with examples

## Function Signatures

### Display Functions (`display.py`)
- `format_menu(items: list[MenuItem], title: str = "", width: int = 80) -> str`
- `show_menu(items: list[MenuItem], title: str = "") -> None`
- `get_choice(prompt: str, valid_choices: list[str] | None = None) -> str`
- `format_breadcrumb(path: list[str]) -> str`
- `clear_screen() -> None`

### Navigation Functions (`navigation.py`)
- `navigate_to_submenu(menu_system: MenuSystem, submenu_id: str) -> None`
- `go_back(menu_system: MenuSystem) -> bool`
- `get_current_menu(menu_system: MenuSystem) -> Menu`

### Discovery Functions (`discovery.py`)
- `discover_scripts(repo_root: Path) -> dict[str, list[ScriptInfo]]`
- `extract_script_metadata(script_path: Path) -> ScriptInfo`
- `categorize_script(script_path: Path) -> str`
- `generate_menu_from_scripts(scripts: dict[str, list[ScriptInfo]]) -> dict[str, Menu]`

### Execution Functions (`executor.py`)
- `execute_script(script_path: Path, args: list[str] | None = None) -> int`
- `execute_python_script(script_path: Path, args: list[str] | None = None) -> int`
- `execute_bash_script(script_path: Path, args: list[str] | None = None) -> int`
- `prompt_for_args(script_info: ScriptInfo) -> list[str]`
- `validate_script_executable(script_path: Path) -> bool`

## Data Structures

### MenuItem
```python
@dataclass
class MenuItem:
    id: str
    label: str
    description: str = ""
    action: str | Callable = ""
    enabled: bool = True
```

### Menu
```python
@dataclass
class Menu:
    id: str
    title: str
    items: list[MenuItem]
    parent_id: str | None = None
```

### ScriptInfo
```python
@dataclass
class ScriptInfo:
    path: Path
    name: str
    description: str
    category: str
    script_type: str
    required_args: list[str] = field(default_factory=list)
    optional_args: list[str] = field(default_factory=list)
```

## Quality Assurance

- **Comprehensive Testing**: Full test coverage for all modules
- **Type Checking**: Complete type hints with mypy compatibility
- **Error Handling**: Robust error handling throughout
- **Documentation**: Complete API documentation and usage examples
- **Integration**: Seamless integration with orchestrator script

## Future Enhancements

Potential future improvements:
- Menu customization and theming
- Script argument auto-completion
- Menu search functionality
- Favorites/recent scripts tracking
- Menu configuration files
- Multi-language support

This menu system provides a user-friendly interface for accessing METAINFORMANT's extensive script collection.



