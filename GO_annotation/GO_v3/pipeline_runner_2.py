import subprocess
import time
import os

scripts = [
    ('1_uniprot_ID_extract.py', 'Script 1', 'uniprot_bees_proteome.fasta'),
    ('2_genetogo.py', 'Script 2', 'bees_uniprot_ids.txt'),
    ('3_genetogotoanno.py', 'Script 3', 'gene2go.txt'),
    ('4_genetogo_summary.py', 'Script 4', 'gene_ontology_annotations.txt')
]

def check_input_file(input_file):
    if not os.path.exists(input_file):
        print(f"Input file '{input_file}' not found. Please ensure it is in the same directory as the script.")
        exit(1)

def run_script(script, name, input_file):
    check_input_file(input_file)
    start_time = time.time()
    try:
        subprocess.run(['python', script], check=True)
        elapsed_time = time.time() - start_time
        print(f"{name} finished successfully in {elapsed_time:.2f} seconds")
    except subprocess.CalledProcessError:
        print(f"{name} failed")
        exit(1)

if __name__ == "__main__":
    for script, name, input_file in scripts:
        run_script(script, name, input_file)
