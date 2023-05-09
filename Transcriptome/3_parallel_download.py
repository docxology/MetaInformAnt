# This Python script allows you to download SRA files using a specified number of parallel threads and a fastp option. The script first removes the header row from the metadata.tsv file, and then iterates through the remaining rows to download each entry in parallel. The number of threads and fastp option can be set via command-line arguments.
# python 2_parallel_download.py 4 --fastp no


import argparse
import subprocess
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import random
import time

def download_fastq_file(entry, fastp_option: str) -> None:
    """
    Download fastq files using amalgkit with the specified entry and fastp option.

    :param entry: The metadata entry to be downloaded.
    :param fastp_option: Fastp option to be used with amalgkit, either "yes" or "no".
    """

    # Create a temporary metadata file with the selected entry
    temp_metadata_file = f"temp_metadata_{entry.name}.tsv"
    entry.to_frame().T.to_csv(temp_metadata_file, sep='\t', index=False)

    # Download fastq files using amalgkit with the specified entry file and fastp option
    subprocess.run(["amalgkit", "getfastq", "--threads", "1", "--fastp", fastp_option, "--metadata", temp_metadata_file])

    # Remove the temporary metadata file
    os.remove(temp_metadata_file)

def download_fastq_files(number_threads: int, fastp_option: str) -> None:
    """
    Download fastq files using amalgkit with a specified number of parallel threads and fastp option.

    :param number_threads: Number of parallel threads to use.
    :param fastp_option: Fastp option to be used with amalgkit, either "yes" or "no".
    """

    # Update the path to the metadata.tsv file
    # metadata_file = "./MetaInformAnt/Transcriptome/metadata/metadata/metadata.tsv"
    metadata_file = "/media/tet/1C0842A829B15012/bio/MetaInformAnt/Transcriptome/metadata/metadata/metadata.tsv"

    if os.path.isfile(metadata_file):
        df = pd.read_csv(metadata_file, sep='\t')
    else:
        print("metadata.tsv not found in the specified directory.")
        return

    # Download fastq files using multiple threads
    with ThreadPoolExecutor(max_workers=number_threads) as executor:
        for _, row in df.iterrows():
            if row.name == 0:
                continue
            executor.submit(download_fastq_file, row, fastp_option)

    # Change back to the original working directory
    os.chdir("..")

def main():
    parser = argparse.ArgumentParser(description="Download SRA files using N parallel threads.")
    parser.add_argument("number_threads", type=int, help="Number of parallel threads to use.")
    parser.add_argument("--fastp", default="no", choices=["yes", "no"], help="Fastp option to be used with amalgkit.")
    args = parser.parse_args()

    # Download fastq files using amalgkit with the specified number of threads and fastp option
    download_fastq_files(args.number_threads, args.fastp)

if __name__ == "__main__":
    main()

