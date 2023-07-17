#!/bin/bash
#SBATCH --job-name=download_genome
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=64  
#SBATCH --mail-user=sacsures@ttu.edu
#SBATCH --mail-type=ALL

# Change the current working directory to the metainformant directory
cd /lustre/research/linksvayer/metainformant/

# This function runs a command and prints its output.
run_command() {
    echo "Running command: $1"
    output=$($1 2>&1)
    exit_status=$?
    echo "$output"
    if [ $exit_status -ne 0 ]; then
        echo "Error: $output"
        exit $exit_status
    fi
}

# Generate metadata using amalgkit
run_command "amalgkit metadata --config_dir config"

# Organize metadata output
mkdir -p intermediates/metadata
mv metadata intermediates/metadata

# 'os.replace' in Python completely replaces the destination, so if the destination exists, it is removed first.
# There is no direct equivalent in bash, so we use 'rm -rf' before 'mv'.
rm -rf metadata
mv intermediates/metadata metadata

echo "Metadata file created in metadata folder"
