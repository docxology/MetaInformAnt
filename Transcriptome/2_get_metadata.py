import os
import subprocess
import sys

def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
    stdout, stderr = process.communicate()
    print(stdout)
    if process.returncode != 0:
        print(f"Error: {stderr}")
        sys.exit(process.returncode)

# Generate metadata using amalgkit
run_command("amalgkit metadata --out_dir \"./\" --entrez_email \"DanielAriFriedman@gmail.com\" --search_string '(\"platform illumina\"[Properties]) AND (\"type rnaseq\"[Filter]) AND (\"Apis mellifera\"[Organism]) NOT (\"DNA\"[Source])'")
# amalgkit metadata --out_dir "./" --entrez_email "example@gmail.com" --search_string '("platform illumina"[Properties]) AND ("type rnaseq"[Filter]) AND ("Apis mellifera"[Organism]) NOT ("DNA"[Source])'

# Organize metadata output
os.makedirs("intermediates/metadata", exist_ok=True)
os.rename("metadata", "intermediates/metadata")
os.replace("intermediates/metadata", "metadata")

print("Metadata file created in metadata folder")
