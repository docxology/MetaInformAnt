import os
import subprocess
import sys

def run_command(command: str) -> None:
    """
    Run a shell command and print the output.

    :param command: The shell command to run.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    stdout, stderr = process.communicate()
    print(stdout)
    if process.returncode != 0:
        print(f"Error: {stderr}")
        sys.exit(process.returncode)

def main() -> None:
    # Download reference sequence
    print("Downloading ref seq")
    os.makedirs("seq", exist_ok=True)
    os.chdir("seq")

    filename = "GCA_003254395.2_Amel_HAv3.1_genomic.fna.gz"
    if not os.path.exists(filename):
        # Download reference sequence if not available
        cmd = "esearch -db assembly -query 'GCF_003254395.2' | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
        ftp_path = subprocess.check_output(cmd, shell=True, text=True).strip()
        fname = ftp_path.split('/')[-1] + "_genomic.fna.gz"
        url = f"{ftp_path}/{fname}"
        run_command(f"wget {url}")

    # Unzip reference sequence and rename
    if os.path.exists(filename):
        run_command(f"gzip -d {filename}")
        os.rename("GCA_003254395.2_Amel_HAv3.1_genomic.fna", "Apis_mellifera.fasta")

    print("Ref seq downloaded to seq folder")

    os.chdir("../")

if __name__ == "__main__":
    main()
