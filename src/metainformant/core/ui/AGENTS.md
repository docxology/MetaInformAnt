# UI Module Agents

## TerminalInterface Agent
**Role**: Display Manager
**Responsibility**: 
- Manages the terminal screen buffer.
- Renders the progress dashboard at 10Hz.
- Ensures thread safety for ui updates.
- Cleanup of terminal state (cursor visibility) on exit.

**Interactions**:
- Receives state updates from `DownloadManager` or other orchestrators.
- Writes directly to `sys.stdout`.
