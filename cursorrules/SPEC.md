# Specification: Cursor Rules

## 🎯 Scope

Defines the behavior and structure of module-specific coding rules for AI assistants.

## 🧱 Architecture

- **Component**: Configuration / Documentation
- **Format**: text/markdown-like
- **Location**: `cursorrules/`

## 💾 Data Structures

- **Files**: `*.cursorrules`
- **Content**:
  - Context & Intent
  - Directory Structure (Module specific)
  - Code Style
  - Testing Rules (Zero Mock)

## 🔌 Integration

- **Cursor**: Automatically detected and applied
- **Scripts**: `scripts/maintenance/generate_cursorrules.py`
