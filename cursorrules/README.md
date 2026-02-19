# Cursor Rules

## 🧠 Context & Intent

- **Purpose**: Module-specific AI coding assistant rules
- **Consistency**: Extends root `.cursorrules`
- **Automation**: Generated and maintained to ensure domain specificity

## 📂 Structure

- `*.cursorrules`: Domain-specific rules (e.g., `dna.cursorrules`)
- `README.md`: Index and guide
- `AGENTS.md`: AI attribution and workflow

## 🔄 Workflow

1. **Load**: Cursor loads these based on file context
2. **Override**: Specific rules override or extend general rules
3. **Update**: Keep synchronized with `src/metainformant/<domain>/` structure
