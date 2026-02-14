#!/usr/bin/env python3
import os
import sys
from pathlib import Path

BLUE_PATH = Path("blue").resolve()
INDEX_PATH = BLUE_PATH / "amalgkit/anoplolepis_gracilipes/work/index"

print(f"Checking permissions for: {BLUE_PATH}")

# 1. Check if blue is readable
try:
    contents = os.listdir(BLUE_PATH)
    print(f"✅ 'blue' is readable. Contents: {contents[:5]}...")
except Exception as e:
    print(f"❌ FAILED to read 'blue': {e}")
    print("\n⚠️  ACTION REQUIRED: Grant 'Full Disk Access' to Terminal.")
    print("   System Settings -> Privacy & Security -> Full Disk Access -> Toggle ON for Terminal")
    sys.exit(1)

# 2. Check if index is readable
print(f"\nChecking index path: {INDEX_PATH}")
if not INDEX_PATH.exists():
    print(f"❌ Index path does not exist: {INDEX_PATH}")
    # Don't exit, might be just missing file
else:
    try:
        files = os.listdir(INDEX_PATH)
        print(f"✅ Index directory is readable. Files: {files[:3]}...")
    except Exception as e:
        print(f"❌ FAILED to read index directory: {e}")
        print("\n⚠️  ACTION REQUIRED: This specific directory seems blocked.")
        print("   If 'blue' was readable but this is not, it might be an extended attribute.")
        print("   Try: xattr -d com.apple.quarantine blue/amalgkit")
        sys.exit(1)

# 3. Check write access (creating a temp file)
test_file = BLUE_PATH / "write_test.tmp"
print(f"\nChecking write access to {test_file}...")
try:
    with open(test_file, "w") as f:
        f.write("test")
    print("✅ Write access confirmed.")
    os.remove(test_file)
except Exception as e:
    print(f"❌ FAILED to write to 'blue': {e}")
    print("\n⚠️  ACTION REQUIRED: The volume might be read-only or restricted.")
    print("   Check 'mount' output or grant Full Disk Access.")
    sys.exit(1)

print("\n✅ All permission checks passed! You can run the pipeline now.")
