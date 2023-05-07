#!/bin/bash

# This script installs required dependencies and sets up the environment for amalgkit.

# Install dependencies
echo "Installing dependencies"
sudo apt update && sudo apt-get dist-upgrade
sudo apt-get -y install fastp
sudo apt-get install pigz
sudo apt-get install kallisto
sudo apt-get install ncbi-entrez-direct
sudo apt install pigz

# Install amalgkit from GitHub
pip install git+https://github.com/kfuku52/amalgkit

# Create 'tools' and 'intermediates' directories if they do not exist
if [ ! -d "tools" ]; then
    mkdir tools
fi

if [ -d "intermediates" ]; then
    sudo rm -r intermediates
fi
mkdir intermediates

# Change to tools directory
cd tools
echo $PWD

# Install Miniconda3 if not already installed
if [ ! -d "$HOME/miniconda3" ]; then
    echo "Installing Miniconda3"
    if [ ! -f "Miniconda3-py39_23.1.0-1-Linux-x86_64.sh" ]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.1.0-1-Linux-x86_64.sh
    fi

    chmod -v +x Miniconda*.sh
    echo "5dc619babc1d19d6688617966251a38d245cb93d69066ccde9a013e1ebb5bf18 *Miniconda3-py39_23.1.0-1-Linux-x86_64.sh" | shasum --check
    bash ./Miniconda3-py39_23.1.0-1-Linux-x86_64.sh

    if ! [[ $PATH =~ "$HOME/miniconda3/bin" ]]; then
        PATH="$HOME/miniconda3/bin:$PATH"
        conda list
        export PATH="$HOME/miniconda3/bin:$PATH"
    fi
else
    echo "Conda already installed"
fi

# Install conda packages
export PATH="/usr/local/bin:$PATH"

conda install -c bioconda parallel-fastq-dump
conda install -c bioconda entrez-direct
conda install -c bioconda seqkit
conda install -c bioconda fastp
conda install -c bioconda bioconductor-pcamethods
conda install -c r r-colorspace
conda install -c r r-rcolorbrewer pcaMethods sva r-NMF r-amap
conda install -c bioconda bioconductor-sva
conda install -c r r-mass
conda install -c r r-nmf
conda install -c conda-forge r-dendextend
conda install -c bioconda amap
conda install -c conda-forge r-pvclust
conda install -c conda-forge r-rtsne

# Change back to the parent directory
cd ../
echo $PWD

# Set up amalgkit configurations
amalgkit config --config base --overwrite yes
mv config_base intermediates
cp -R ./intermediates/config_base config

cd config
rm search_term_species.config
rm search_term_keyword.config

# Change this for whichever tissues you want.
printf '# case-insensitive\n# regular expressions not allowed\n\n\"Apis mellifera\"' > search_term_species.config
# printf '# case-insensitive\n# regular expressions not allowed\n\n\"brain"\n"liver"\n"kidney"' > search_term_keyword.config
# printf '# case-insensitive\n# regular expressions not allowed\n\n"brain"\n"liver"\n"kidney"' > search_term_tissue.config
cd ../

# Fix numpy compatibility issue in amalgkit/getfastq.py
if [ -d "$HOME/.local/lib/python3.9/site-packages/amalgkit/getfastq.py" ]; then
sudo sed -i 's/numpy.object, numpy.str/object, str/' ~/.local/lib/python3.9/site-packages/amalgkit/getfastq.py
fi
if [ -d "$HOME/.local/lib/python3.9/site-packages/amalgkit/getfastq.py" ]; then
sudo cp -r /usr/local/lib/python3.9/dist-packages/amalgkit $HOME/.local/lib/python3.9/site-packages/
sudo sed -i 's/numpy.object, numpy.str/object, str/' ~/.local/lib/python3.9/site-packages/amalgkit/getfastq.py
fi

if [ -d "/usr/local/lib/python3.9/dist-packages/amalgkit/getfastq.py" ]; then
sudo sed -i 's/numpy.object, numpy.str/object, str/' ~/.local/lib/python3.9/site-packages/amalgkit/getfastq.py
fi

echo ""
echo "Dependencies installed successfully"
