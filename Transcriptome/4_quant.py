import os
import subprocess

def run_command(command):
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(process.stdout)
    if process.returncode != 0:
        print(f"Error: {process.stderr}")

# Set the desired k-mer length for kallisto index
kmer_length = 31

# Set the desired number of processor threads
number_threads = 14

# Create the index directory if it does not exist
index_dir = "index"
os.makedirs(index_dir, exist_ok=True)

# Create the kallisto index using the specified k-mer length and input fasta file
kallisto_index_cmd = f"kallisto index -k {kmer_length} -i {index_dir}/Apis_mellifera.idx ./seq/Apis_mellifera.fasta"
run_command(kallisto_index_cmd)

# Run amalgkit quant with the specified parameters
amalgkit_quant_cmd = f"amalgkit quant --fasta_dir seq --threads {number_threads} --index_dir {index_dir}"
run_command(amalgkit_quant_cmd)
