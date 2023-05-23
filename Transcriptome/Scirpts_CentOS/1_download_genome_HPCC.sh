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

# Create the seq directory if it doesn't exist and move into it.
mkdir -p seq
cd seq

filename="GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz"

if [ ! -f "$filename" ]; then
    echo "Downloading ref seq"
    # Run the esearch command to get the FTP path.
    cmd="esearch -db assembly -query 'GCF_003254395.2' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
    ftp_path=$(eval $cmd)
    fname=$(basename $ftp_path)_genomic.fna.gz
    url="${ftp_path}/${fname}"
    run_command "wget $url"
fi

# Unzip the reference sequence and rename the file.
if [ -f "$filename" ]; then
    run_command "gzip -d $filename"
    mv "GCA_003254395.2_Amel_HAv3.1_genomic.fna" "Apis_mellifera.fasta"
fi

echo "Ref seq downloaded to seq folder"

#------- end of the script ------------


