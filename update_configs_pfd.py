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
    
    new_lines = []
    
    for line in lines:

        if line.startswith('  getfastq:'):
            in_getfastq = True
            new_lines.append(line)
            continue
        elif line.startswith('  ') and not line.startswith('    ') and in_getfastq:
            in_getfastq = False

        if in_getfastq:
            if line.startswith('    pfd:'):
                comment = line[line.find('#'):] if '#' in line else ''
                new_lines.append(f'    pfd: yes  {comment}'.rstrip() + '\n')
                continue
            elif line.startswith('    # pfd:'):
                new_lines.append('    pfd: yes  # Automatically enabled for speed\n')
                continue

        new_lines.append(line)
        
    with open(f, 'w') as file:
        file.writelines(new_lines)
        
    count += 1

print(f"Updated {count} config files to enable pfd successfully.")
