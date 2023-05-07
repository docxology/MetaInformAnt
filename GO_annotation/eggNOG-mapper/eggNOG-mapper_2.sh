#!/bin/bash

# GO annotation based on eggNOG-mapper

# -----------------------
# Install dependencies
# -----------------------

# Install eggNOG-mapper via conda
conda install -c bioconda eggnog-mapper

# Install diamond if not already installed
conda install -c bioconda diamond

# -----------------------
# Download databases
# -----------------------

# Define URLs for database files
EGGNOG_DB_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz"
EGGNOG_PROTEINS_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz"
MMSEQS_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz" #optional
EGGNOG_TAXA_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz" #optional

# Download the necessary database files
wget $EGGNOG_DB_URL
wget $EGGNOG_PROTEINS_URL
wget $MMSEQS_URL #optional
wget $EGGNOG_TAXA_URL #optional

# -----------------------
# Run eggNOG-mapper
# -----------------------

# Define input and output file names
SEQUENCE_FILE="uniprot_bees_proteome.fasta"
OUTPUT_FILE_NAME="output"

# Run eggNOG-mapper
emapper.py --cpu 10 -i $SEQUENCE_FILE --output $OUTPUT_FILE_NAME -d euk -m diamond
