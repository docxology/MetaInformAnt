import os
import subprocess

def run_command(command: str) -> None:
    """
    Run a shell command and print its stdout and stderr.

    :param command: The command to be executed.
    """
    process = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(process.stdout)
    if process.returncode != 0:
        print(f"Error: {process.stderr}")

def create_kallisto_index(kmer_length: int, index_dir: str, input_fasta: str) -> None:
    """
    Create the kallisto index using the specified k-mer length and input fasta file.

    :param kmer_length: The desired k-mer length for kallisto index.
    :param index_dir: The directory to save the kallisto index.
    :param input_fasta: The input fasta file to create the kallisto index.
    """
    os.makedirs(index_dir, exist_ok=True)
    kallisto_index_cmd = f"kallisto index -k {kmer_length} -i {index_dir}/Apis_mellifera.idx {input_fasta}"
    run_command(kallisto_index_cmd)

def run_amalgkit_quant(fasta_dir: str, number_threads: int, index_dir: str) -> None:
    """
    Run amalgkit quant with the specified parameters.

    :param fasta_dir: The directory containing the fasta files.
    :param number_threads: The number of processor threads to use.
    :param index_dir: The directory containing the kallisto index.
    """
    amalgkit_quant_cmd = f"amalgkit quant --fasta_dir {fasta_dir} --threads {number_threads} --index_dir {index_dir}"
    run_command(amalgkit_quant_cmd)

def main() -> None:
    kmer_length = 31
    number_threads = 10
    index_dir = "index"
    input_fasta = "./seq/Apis_mellifera.fasta"

 #   create_kallisto_index(kmer_length, index_dir, input_fasta)
    run_amalgkit_quant("seq", number_threads, index_dir)

if __name__ == "__main__":
    main()
