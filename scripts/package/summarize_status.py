import subprocess, ast, sys

def main():
    # Run the status command
    result = subprocess.run(['bash', 'scripts/rna/run_all_species.sh', '--status'], cwd='/home/trim/Documents/Git/MetaInformAnt', capture_output=True, text=True)
    if result.returncode != 0:
        print('Error running status command', file=sys.stderr)
        sys.exit(1)
    lines = result.stdout.splitlines()
    summary = {
        'total_species': 0,
        'completed': 0,
        'failed': 0,
        'undownloaded_samples': 0,
        'quantified_and_deleted': 0,
        'other': 0,
    }
    for line in lines:
        if 'Status:' in line:
            # Extract the dict after 'Status:'
            try:
                dict_str = line.split('Status:')[1].strip()
                data = ast.literal_eval(dict_str)
                summary['total_species'] += 1
                if data.get('completed'):
                    summary['completed'] += 1
                if data.get('failed'):
                    summary['failed'] += 1
                categories = data.get('categories', {})
                summary['undownloaded_samples'] += categories.get('undownloaded', 0)
                summary['quantified_and_deleted'] += categories.get('quantified_and_deleted', 0)
            except Exception as e:
                # ignore parsing errors
                continue
    # Print summary
    print('Sample Status Summary:')
    for k, v in summary.items():
        print(f'{k}: {v}')

if __name__ == '__main__':
    main()
