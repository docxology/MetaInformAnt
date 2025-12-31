#!/usr/bin/env python3
"""
Simple test to verify METAINFORMANT modules work.
"""

import sys
import traceback
from pathlib import Path

# Add src to path
repo_root = Path(__file__).resolve().parent
src_dir = repo_root / "src"
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

def test_imports():
    """Test that all major modules can be imported."""
    print("Testing module imports...")

    modules_to_test = [
        'metainformant.core',
        'metainformant.core.cache',
        'metainformant.core.config',
        'metainformant.core.io',
        'metainformant.ontology',
        'metainformant.phenotype',
        'metainformant.quality',
        'metainformant.visualization.basic',
        'metainformant.ml.classification',
    ]

    failed = []
    for module in modules_to_test:
        try:
            __import__(module)
            print(f"✅ {module}")
        except Exception as e:
            print(f"❌ {module}: {e}")
            failed.append(module)

    return len(failed) == 0

def test_basic_functionality():
    """Test basic functionality."""
    print("\nTesting basic functionality...")

    try:
        from metainformant.core import io
        from metainformant.core import paths

        # Test path utilities
        test_path = paths.expand_and_resolve(".")
        print(f"✅ Path resolution works: {test_path}")

        # Test basic validation
        from metainformant.core import validation
        validation.validate_type("test", str, "test_string")
        print("✅ Validation works")

        return True
    except Exception as e:
        print(f"❌ Basic functionality failed: {e}")
        traceback.print_exc()
        return False

def main():
    print("🧬 METAINFORMANT Basic Test Suite")
    print("=" * 40)

    imports_ok = test_imports()
    functionality_ok = test_basic_functionality()

    print("\n" + "=" * 40)
    if imports_ok and functionality_ok:
        print("🎉 All basic tests passed!")
        return 0
    else:
        print("❌ Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())
