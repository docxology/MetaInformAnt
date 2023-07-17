#!/bin/bash

# This script installs required dependencies and sets up the environment for amalgkit.

# Stop on first error
set -e

# Install dependencies
echo "Installing dependencies"
sudo apt update && sudo apt-get dist-upgrade
sudo apt-get -y install fastp
sudo apt-get install pigz
sudo apt-get install kallisto
sudo apt-get install ncbi-entrez-direct

# Install amalgkit from GitHub
pip install git+https://github.com/kfuku52/amalgkit

# Create 'tools' and 'intermediates' directories if they do not exist
mkdir -p tools intermediates

# Change to tools directory
cd tools

# Install Miniconda3 if not already installed
if [ ! -d "$HOME/miniconda3" ]; then
    echo "Installing Miniconda3"
    if [ ! -f "Miniconda3-py39_23.1.0-1-Linux-x86_64.sh" ]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
    fi
    chmod +x Miniconda*.sh
    ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
fi

# Add Conda to PATH in .bashrc if it's not already present
if ! grep -q "export PATH=\"$HOME/miniconda3/bin:\$PATH\"" "$HOME/.bashrc" ; then
    echo "export PATH=\"$HOME/miniconda3/bin:\$PATH\"" >> "$HOME/.bashrc"
    source "$HOME/.bashrc"
fi

# Install conda packages
conda install -c bioconda parallel-fastq-dump entrez-direct seqkit fastp bioconductor-pcamethods
conda install -c r r-colorspace r-rcolorbrewer pcaMethods sva r-NMF r-amap
conda install -c bioconda bioconductor-sva
conda install -c r r-mass r-nmf
conda install -c conda-forge r-dendextend
conda install -c bioconda amap
conda install -c conda-forge r-pvclust r-rtsne

# Change back to the parent directory
cd ../

# Set up amalgkit configurations
amalgkit config --config base --overwrite yes
mv config_base intermediates
cp -R ./intermediates/config_base config

cd config
rm search_term_species.config search_term_keyword.config

# Change this for whichever tissues you want.
printf '# case-insensitive\n# regular expressions not allowed\n\n"Apis mellifera"' > search_term_species.config
cd ../

# Fix numpy compatibility issue in amalgkit/getfastq.py
sed -i 's/numpy.object, numpy.str/object, str/' ~/.local/lib/python3.9/site-packages/amalgkit/getfastq.py

echo ""
echo "Dependencies installed successfully"

