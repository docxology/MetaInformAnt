import sys, glob
for path in glob.glob("config/amalgkit/amalgkit_*.yaml"):
    if "template" in path or "test" in path or "cross_species" in path: continue
    with open(path, "r") as f: lines = f.readlines()
    
    in_getfastq, in_quant = False, False
    out_lines = []
    
    for line in lines:
        if line.startswith("  getfastq:"):
            in_getfastq, in_quant = True, False
            out_lines.append(line)
        elif line.startswith("  quant:"):
            in_getfastq, in_quant = False, True
            out_lines.append(line)
        elif line.startswith("  merge:"):
            in_getfastq, in_quant = False, False
            out_lines.append(line)
        elif in_getfastq and "jobs:" in line:
            out_lines.append("    jobs: 8\n")
        elif in_quant and "jobs:" in line:
            pass # We'll check if we need to insert it
        elif line.startswith("  ") and not line.startswith("    "):
            in_getfastq, in_quant = False, False
            out_lines.append(line)
        else:
            out_lines.append(line)
            
    # Re-inject quant jobs if missing
    final_lines = []
    in_q = False
    for line in out_lines:
        if line.startswith("  quant:"):
            in_q = True
            final_lines.append(line)
            final_lines.append("    jobs: 1\n")
        elif in_q and "jobs:" in line:
            continue
        elif line.startswith("  ") and not line.startswith("    "):
            in_q = False
            final_lines.append(line)
        else:
            final_lines.append(line)
            
    with open(path, "w") as f:
        f.writelines(final_lines)
print("Yamls patched")
