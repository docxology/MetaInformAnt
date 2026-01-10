# Terminal UI (TUI) Module

## Overview
The `metainformant.core.ui` module provides a lightweight, dependency-free framework for creating "oldschool" terminal visualizations. It is designed to handle multi-threaded progress tracking, colorized output, and cursor management using standard ANSI escape codes.

## Design Philosophy
- **Zero Dependencies**: Uses only Python standard library.
- **Thread Safe**: Can safely receive updates from multiple background threads.
- **Elegant Aesthetic**: Uses ANSI colors and block characters for a clean, retro feel.

## Key Components
### `TerminalInterface`
The main controller that manages the rendering loop and screen updates.
- `add_bar()`: Registers a new tracked task.
- `update()`: Updates the state of a specific task.
- `start()`: Begins the background rendering thread (10fps).
- `stop()`: Cleanly shuts down the renderer and restores the cursor.

### `ProgressState`
Data class holding the current state of a task (progress, speed, status, color).
