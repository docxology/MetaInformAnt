import os, glob

config_dir = "/home/trim/Documents/Git/MetaInformAnt/config/amalgkit"
files = glob.glob(os.path.join(config_dir, "amalgkit_*.yaml"))

count = 0
for f in files:
    if "template" in f or "test" in f or "cross_species" in f:
        continue
    
    with open(f, 'r') as file:
        lines = file.readlines()
        
    in_getfastq = False
    in_quant = False
    
    new_lines = []
    
    for line in lines:
        if line.startswith('threads:') and not line.startswith(' '):
            new_lines.append('threads: 16' + ('  ' + line[line.find('#'):] if '#' in line else '\n'))
            continue

        if line.startswith('  getfastq:'):
            in_getfastq = True
            in_quant = False
            new_lines.append(line)
            continue
        elif line.startswith('  quant:'):
            in_quant = True
            in_getfastq = False
            new_lines.append(line)
            continue
        elif line.startswith('  ') and not line.startswith('    ') and (in_getfastq or in_quant):
            in_getfastq = False
            in_quant = False

        if in_getfastq:
            if line.startswith('    threads:'):
                new_lines.append('    threads: 16\n')
                continue
            elif line.startswith('    jobs:'):
                new_lines.append('    jobs: 16\n')
                continue
                
        if in_quant:
            if line.startswith('    threads:'):
                new_lines.append('    threads: 16\n')
                continue

        new_lines.append(line)
        
    with open(f, 'w') as file:
        file.writelines(new_lines)
        
    count += 1

print(f"Updated {count} config files to 16 threads/jobs.")
