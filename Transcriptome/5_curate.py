import subprocess

def run_command(command):
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(process.stdout)
    if process.returncode != 0:
        print(f"Error: {process.stderr}")

# Run amalgkit merge to combine expression data from multiple samples
amalgkit_merge_cmd = "amalgkit merge"
run_command(amalgkit_merge_cmd)

# Set the normalization method for curation
normalization_method = "log2p1-fpkm"

# Run amalgkit curate with the specified normalization method
amalgkit_curate_cmd = f"amalgkit curate --norm {normalization_method}"
run_command(amalgkit_curate_cmd)
