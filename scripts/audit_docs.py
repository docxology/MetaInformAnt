
import os
from pathlib import Path

def check_triple_play(root_dir):
    root = Path(root_dir)
    src_dir = root / "src" / "metainformant"
    
    # Modules to check (top-level directories in src/metainformant excluding __pycache__ etc)
    modules = [d for d in src_dir.iterdir() if d.is_dir() and not d.name.startswith("__") and d.name != "tests"]
    
    print(f"Checking {len(modules)} top-level modules in {src_dir}...\n")
    
    missing_docs = {}
    
    for module in sorted(modules):
        module_name = module.name
        required_files = ["README.md", "AGENTS.md", "SPEC.md"]
        missing = []
        
        for f in required_files:
            if not (module / f).exists():
                missing.append(f)
                
        if missing:
            missing_docs[module_name] = missing
            print(f"[FAIL] {module_name}: Missing {', '.join(missing)}")
        else:
            print(f"[PASS] {module_name}")
            
    # Also check subdirectories recursively? 
    # The prompt listed top-level modules specifically, but let's see.
    # For now, top-level is the primary requirement.
    
    if not missing_docs:
        print("\nAll top-level modules have Triple Play documentation!")
    else:
        print(f"\nFound missing documentation in {len(missing_docs)} modules.")

if __name__ == "__main__":
    check_triple_play(".")
