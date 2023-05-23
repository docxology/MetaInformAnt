#!/bin/bash

# This script installs required dependencies and sets up the environment for amalgkit. 
# Note that this documentation is for setting up the environment in the High Performance Computing Core (HPCC) at Texas Tech University (TTU). HPCC at TTU is based on CentOS; more info here: https://www.depts.ttu.edu/hpcc/ 

# Installing dependencies
# All the required dependenicies are installed in the following path: /home/sacsures/
# Note that the dependencies 'fastp', 'pigz', 'kallisto', and 'ncbi-entrez-direct' are downloaded and installed as prebuild binaries since I don't have root access for system wide installation (i.e., cannot use 'sudo yum install *any software*). Alternatively, these can be installed using conda as well

# Change the current working directory to the user's home directory; therefore, this script assumes all permissions to this directory
cd /home/sacsures

# Downloading and Installing fastp
wget http://opengene.org/fastp/fastp # Download the fastp binary from the OpenGene website.
chmod a+x ./fastp # Make it executable
export PATH=$PATH:/home/sacsures/ # Add the directory containing the fastp binary to the PATH environment variable so it can be run from anywhere.

# Downloading and Installing pigz
wget https://zlib.net/pigz/pigz.tar.gz
tar -xzf /path/to/pigz.tar.gz # Extract the downloaded pigz tar.gz file.
cd pigz
make # Compile the source code.
export PATH=$PATH:/home/sacsures/pigz/

# Downloading and Installing kallisto
wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
tar -xzf kallisto_linux-v0.46.1.tar.gz # Direct executable will be downloaded after completing this step
export PATH=$PATH:/home/sacsures/kallisto/

# Downloading and Installing ncbi-entrez-direct
sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=/home/sacsures/edirect:${PATH}


# Install amalgkit from GitHub
pip install --user git+https://github.com/kfuku52/amalgkit # Install the Python package amalgkit from GitHub using pip.
export PATH=$PATH:~/.local/bin

# Install minconda3
# Download the most recent version of Miniconda v3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh 
# Execute the installation script to install to $HOME/conda.
bash ~/miniconda.sh -b -p ~/conda


# Install conda packages
# Create a seperate conda environment for this project to isolate dependencies and avoid conflicts
conda create --name metainformant #Create environment
conda activate metainformant #Activate the environment

# Now install all the necessay packages
conda install -c bioconda parallel-fastq-dump
conda install -c bioconda entrez-direct
conda install -c bioconda seqkit
conda install -c bioconda fastp
conda install -c bioconda bioconductor-pcamethods
conda install -c r r-colorspace
conda install -c r r-rcolorbrewer
conda install -c bioconda bioconductor-sva
conda install -c r r-mass
conda install -c r r-nmf
conda install -c conda-forge r-dendextend
conda install -c bioconda amap
conda install -c conda-forge r-pvclust
conda install -c conda-forge r-rtsne

# All the above packages and software are installed in my user directory, but is executable in any directory which I have access
# Now, move to the parent directory
cd /lustre/research/linksvayer/metainformant/

# Set up amalgkit configurations
mkdir intermediates
amalgkit config --config base --overwrite yes
mv config_base intermediates
cp -R ./intermediates/config_base config

cd config
rm search_term_species.config
rm search_term_keyword.config

# Create a new search term configuration file for the species "Apis mellifera". hange this for whichever tissues you want.
printf '# case-insensitive\n# regular expressions not allowed\n\n\"Apis mellifera\"' > search_term_species.config

#---------- end of this script -------------------