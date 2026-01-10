# UI Module Specification

## Class: TerminalInterface

### Methods
- `__init__(self)`: Initialize locks and state.
- `add_bar(self, task_id: str, label: str, total: float, unit: str = "B")`: add a new progress bar.
- `update(self, task_id: str, ...)`: Update task state.
- `start(self)`: Start the valid rendering loop.
- `stop(self)`: Stop rendering.

## ANSI Standards
- Uses standard CSI codes.
- Colors: 3-bit/4-bit standard colors (30-37).
- Cursor: `?25l` (hide), `?25h` (show), `A` (up), `2K` (clear line).
