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

def run_amalgkit_merge() -> None:
    """
    Run amalgkit merge to combine expression data from multiple samples.
    """
    amalgkit_merge_cmd = "amalgkit merge"
    run_command(amalgkit_merge_cmd)

def run_amalgkit_curate(normalization_method: str) -> None:
    """
    Run amalgkit curate with the specified normalization method.

    :param normalization_method: The normalization method to be used for curation.
    """
    amalgkit_curate_cmd = f"amalgkit curate --norm {normalization_method}"
    run_command(amalgkit_curate_cmd)

def main() -> None:
    run_amalgkit_merge()

    normalization_method = "log2p1-fpkm"
    run_amalgkit_curate(normalization_method)

if __name__ == "__main__":
    main()
