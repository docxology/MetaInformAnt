# UI

Terminal user interface for displaying workflow progress with live-updating bars, phase tracking, and elapsed time formatting.

## Contents

| File | Purpose |
|------|---------|
| `tui.py` | Terminal progress display with multi-phase tracking and ANSI rendering |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `TerminalInterface` | Full-screen terminal UI with phase progress, sample counts, and ETA |
| `ProgressState` | Dataclass tracking current/total items, phase name, and timing |

## Usage

```python
from metainformant.core.ui.tui import TerminalInterface

tui = TerminalInterface()
tui.start()
tui.update_progress(phase="download", current=5, total=100)
tui.stop()
```
