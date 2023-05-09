#!/bin/bash

# This script downloads fastq files using amalgkit.
# It optionally removes any existing getfastq directory and downloads fastq files with the specified number of threads.

# Set the desired number of processor threads
number_threads=14

# Remove existing getfastq directory if it exists
# if [ -d "getfastq" ]; then
#    sudo rm -r getfastq
# fi

# Create a new getfastq directory
# mkdir getfastq

# Change the working directory to the newly created getfastq directory
# cd getfastq

# Download fastq files using amalgkit with the specified number of threads
amalgkit getfastq --threads $number_threads

# Change back to the original working directory
cd ..
