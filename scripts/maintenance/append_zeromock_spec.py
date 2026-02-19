#!/usr/bin/env python3
"""
Append Zero Mock policy to SPEC.md files.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DOCS_DIR = REPO_ROOT / "docs"

POLICY = """
## 🧪 Testing Policy
- **Zero Mock**: All tests must use real implementations. Mocks are strictly prohibited.
"""

def main():
    print(f"Checking SPEC.md files in {DOCS_DIR}...")
    
    count = 0
    for spec_file in DOCS_DIR.rglob("SPEC.md"):
        content = spec_file.read_text()
        if "Zero Mock" not in content:
            print(f"Appending policy to {spec_file.relative_to(REPO_ROOT)}")
            with open(spec_file, "a") as f:
                f.write(POLICY)
            count += 1
            
    print(f"Updated {count} files.")

if __name__ == "__main__":
    main()
