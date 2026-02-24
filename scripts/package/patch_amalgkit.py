import sys
import os
import re
import site
from pathlib import Path

def patch_amalgkit_getfastq():
    """Patches amalgkit to prevent NCBI SRA toolkit SQLite locking during parallel extraction."""
    site_packages = site.getsitepackages()
    getfastq_path = None
    
    for sp in site_packages:
        candidate = Path(sp) / "amalgkit" / "getfastq.py"
        if candidate.exists():
            getfastq_path = candidate
            break
            
    if not getfastq_path:
        print("Could not find amalgkit/getfastq.py. Is amalgkit installed?")
        sys.exit(1)
        
    with open(getfastq_path, "r") as f:
        content = f.read()
        
    # Check if already patched
    if "amalgkit_sra_home_" in content:
        print(f"amalgkit is already patched at {getfastq_path}")
        sys.exit(0)
        
    # The target function signature to patch
    target_func = "def execute_fasterq_dump_command(fasterq_dump_command, args, prefix='Command'):\n"
    original_code = "    fqd_out = subprocess.run(fasterq_dump_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)\n"
    
    # We find the target function and replace the subprocess.run call
    patch_code = """    import os, tempfile, shutil
    env = os.environ.copy()
    temp_home = tempfile.mkdtemp(prefix="amalgkit_sra_home_")
    env["HOME"] = temp_home
    env["NCBI_SETTINGS"] = os.path.join(temp_home, "user-settings.mkfg")
    fqd_out = subprocess.run(fasterq_dump_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    shutil.rmtree(temp_home, ignore_errors=True)
"""
    
    if target_func not in content or original_code not in content:
        print(f"Failed to apply patch. The source code in {getfastq_path} does not match expected structure.")
        sys.exit(1)
        
    # Inject our patch instead of the original code
    patched_content = content.replace(original_code, patch_code)
    
    with open(getfastq_path, "w") as f:
        f.write(patched_content)
        
    print(f"Successfully patched {getfastq_path} to enable safe parallel SRA extraction without NCBI database locks.")

if __name__ == "__main__":
    patch_amalgkit_getfastq()
