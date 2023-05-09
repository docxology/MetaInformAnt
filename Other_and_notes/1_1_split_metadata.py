import pandas as pd
import os

def split_tsv(input_file, output_prefix, num_files):
    # Read input TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Calculate the number of rows per output file
    rows_per_file = len(df) // num_files
    remainder = len(df) % num_files

    # Split the dataframe and write to output TSV files
    start = 0
    for i in range(1, num_files + 1):
        end = start + rows_per_file

        # Distribute the remainder across the output files
        if remainder > 0:
            end += 1
            remainder -= 1

        # Write the current chunk to an output TSV file
        output_file = f"{output_prefix}_{i}.tsv"
        df[start:end].to_csv(output_file, index=False, sep='\t')
        print(f"Created {output_file}")

        # Update the start index for the next chunk
        start = end

# Set input file name, output file prefix, and number of output files
input_file = "./metadata/metadata/metadata.tsv"
output_prefix = "metadata"
num_files = 10

# Ensure the input file exists
if os.path.isfile(input_file):
    # Split the input TSV into N output TSV files
    split_tsv(input_file, output_prefix, num_files)
else:
    print(f"Input file '{input_file}' not found.")
