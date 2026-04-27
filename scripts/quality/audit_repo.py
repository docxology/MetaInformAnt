import os
from pathlib import Path


def audit_repo():
    base_dir = Path("src/metainformant")
    exclude = {"__pycache__", "tests"}

    print(f"{'Module':<20} {'README':<8} {'AGENTS':<8} {'SPEC':<8} {'__init__':<8}")
    print("-" * 60)

    missing_docs = []

    for item in sorted(base_dir.iterdir()):
        if item.is_dir() and item.name not in exclude:
            readme = (item / "README.md").exists()
            agents = (item / "AGENTS.md").exists()
            spec = (item / "SPEC.md").exists()
            init = (item / "__init__.py").exists()

            print(
                f"{item.name:<20} {'✅' if readme else '❌':<8} {'✅' if agents else '❌':<8} {'✅' if spec else '❌':<8} {'✅' if init else '❌':<8}"
            )

            if not (readme and agents and spec):
                missing_docs.append(item.name)

    print("\nSummary:")
    if missing_docs:
        print(f"Modules missing documentation: {', '.join(missing_docs)}")
    else:
        print("All modules have the Triple Play docs!")


if __name__ == "__main__":
    audit_repo()
