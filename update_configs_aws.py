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
            if line.startswith('    aws:'):
                comment = line[line.find('#'):] if '#' in line else ''
                new_lines.append(f'    aws: yes  {comment}'.rstrip() + '\n')
                continue
            elif line.startswith('    # aws:'):
                new_lines.append('    aws: yes  # Automatically enabled for speed and correct sra format\n')
                continue
            if line.startswith('    ncbi:'):
                comment = line[line.find('#'):] if '#' in line else ''
                new_lines.append(f'    ncbi: no  {comment}'.rstrip() + '\n')
                continue
            elif line.startswith('    # ncbi:'):
                new_lines.append('    ncbi: no  # Disabled because sralite lacks quality scores\n')
                continue

        new_lines.append(line)
        
    with open(f, 'w') as file:
        file.writelines(new_lines)
        
    count += 1

print(f"Updated {count} config files to enable aws and disable ncbi successfully.")
