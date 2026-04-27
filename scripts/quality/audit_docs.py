import os
from pathlib import Path

REQUIRED_FILES = ["README.md", "AGENTS.md", "SPEC.md", "PAI.md"]


def is_relevant_dir(path: Path) -> bool:
    # Skip hidden folders, venv, cache, etc.
    if any(part.startswith(".") for part in path.parts):
        return False
    if "output" in path.parts or "__pycache__" in path.parts:
        return False

    # Needs to be a significant folder.
    # Logic: Contains __init__.py OR is a direct child of 'scripts', 'docs', 'tests'
    if (path / "__init__.py").exists():
        return True

    if path.parent.name in ["scripts", "docs", "tests", "config", "examples"]:
        return True

    return False


def audit_directory(root_path: Path):
    print(f"{'Directory':<60} | {'Missing Files'}")
    print("-" * 80)

    for root, dirs, files in os.walk(root_path):
        current_dir = Path(root)

        # Skip checking the root itself again for now, as I just did it (or I can double check)
        # Check relevance
        if not is_relevant_dir(current_dir) and current_dir != root_path:
            continue

        missing = []
        for req in REQUIRED_FILES:
            if not (current_dir / req).exists():
                missing.append(req)

        if missing:
            rel_path = current_dir.relative_to(root_path)
            print(f"{str(rel_path):<60} | {', '.join(missing)}")


if __name__ == "__main__":
    audit_directory(Path("."))
