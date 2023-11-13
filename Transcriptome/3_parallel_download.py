import argparse
import subprocess
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

def download_fastq_file(entry, fastp_option: str) -> None:
    """
    Download fastq files using amalgkit with the specified entry and fastp option.

    :param entry: The metadata entry to be downloaded.
    :param fastp_option: Fastp option to be used with amalgkit, either "yes" or "no".
    """
    temp_metadata_file = f"temp_metadata_{entry.name}.tsv"
    entry.to_frame().T.to_csv(temp_metadata_file, sep='\t', index=False)

    subprocess.run(["amalgkit", "getfastq", "--threads", "1", "--fastp", fastp_option, "--metadata", temp_metadata_file])

    os.remove(temp_metadata_file)

def download_fastq_files(metadata_file: str, number_threads: int, fastp_option: str) -> None:
    """
    Download fastq files using amalgkit with a specified number of parallel threads and fastp option.

    :param metadata_file: Path to the metadata.tsv file.
    :param number_threads: Number of parallel threads to use.
    :param fastp_option: Fastp option to be used with amalgkit, either "yes" or "no".
    """
    if os.path.isfile(metadata_file):
        df = pd.read_csv(metadata_file, sep='\t')
    else:
        print("metadata.tsv not found in the specified directory.")
        return

    with ThreadPoolExecutor(max_workers=number_threads) as executor:
        for _, row in df.iterrows():
            if row.name == 0:
                continue
            executor.submit(download_fastq_file, row, fastp_option)

def main() -> None:
    parser = argparse.ArgumentParser(description="Download SRA files using N parallel threads.")
    parser.add_argument("number_threads", type=int, help="Number of parallel threads to use.")
    parser.add_argument("--fastp", default="yes", choices=["yes", "no"], help="Fastp option to be used with amalgkit.")
    parser.add_argument("metadata_file", type=str, help="Path to the metadata.tsv file.")
    args = parser.parse_args()

    download_fastq_files(args.metadata_file, args.number_threads, args.fastp)

if __name__ == "__main__":
    main()
