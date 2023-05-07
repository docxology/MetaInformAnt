#!/bin/bash

# This script generates metadata, downloads reference sequence, and organizes the output.

# Generate metadata using amalgkit
amalgkit metadata --config_dir config --overwrite yes --max_sample 10

# Uncomment the following lines if integrating local fastq files into amalgkit metadata output
# amalgkit integrate --fastq_dir /PATH/TO/FASTQ/DIRECTORY/ --out_dir /WORKING/DIRECTORY/
# amalgkit integrate --fastq_dir ./getfastq --out_dir ./ --metadata ./metadata.tsv

# Organize metadata output
cd intermediates
if [ -d "metadata" ]; then
    sudo rm -r metadata
fi
cd ../
mv metadata intermediates
cp -R ./intermediates/metadata ./metadata
echo "Metadata file created in metadata folder"

# Download reference sequence
echo 'Downloading ref seq'
if [ -d "seq" ]; then
    cd seq
    if [ ! -f "GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz" ]; then
        # Download reference sequence if not available
        esearch -db assembly -query 'GCF_003254395.2' \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r line ;
          do
              fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
              wget "$line/$fname" ;
          done
    fi
else
    mkdir seq
    cd seq
    # Download reference sequence
    esearch -db assembly -query 'GCF_003254395.2' \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r line ;
      do
          fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
          wget "$line/$fname" ;
      done
fi

# Unzip reference sequence and rename
gzip -d GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz
mv GCA_003254395.2_Amel_HAv3.1_genomic.fna Apis_mellifera.fasta

echo 'Ref seq downloaded to seq folder'

cd ../
