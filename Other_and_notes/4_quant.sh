#!/bin/bash

# This script creates an index for kallisto using the specified k-mer length,
# and then runs amalgkit quant with the provided parameters.

# Set the desired k-mer length for kallisto index
kmer_length=31

# Set the desired number of processor threads
number_threads=14

# Create the index directory if it does not exist
index_dir="index"
if [ ! -d "$index_dir" ]; then
    mkdir "$index_dir"
fi

# Create the kallisto index using the specified k-mer length and input fasta file
kallisto index -k $kmer_length -i "$index_dir/Apis_mellifera.idx" "./seq/Apis_mellifera.fasta"

# Run amalgkit quant with the specified parameters
amalgkit quant --fasta_dir seq --threads $number_threads --index_dir "$index_dir"
