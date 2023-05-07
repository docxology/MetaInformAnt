import os
import requests
from importlib.util import find_spec

# Install required libraries
os.system("pip install requests")

# Check if ete3 is installed and install it if needed
if find_spec("ete3") is None:
    os.system("pip install ete3")

from ete3 import NCBITaxa

# Constants
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/stream"
FORMAT = "fasta"
OUTPUT_DIR = "proteomes"

# Function to download the proteome FASTA file
def download_proteome_fasta(taxon_id, url, params, output_filename):
    response = requests.get(url, params=params)
    if response.status_code == 200:
        with open(output_filename, "wb") as output_file:
            output_file.write(response.content)

        # Count the number of lines starting with '>' and file size
        num_lines = 0
        file_size_mb = os.path.getsize(output_filename) / (1024 * 1024)
        with open(output_filename, "r") as input_file:
            for line in input_file:
                if line.startswith('>'):
                    num_lines += 1

        # Get the Latin name of the species
        ncbi = NCBITaxa()
        species_name = ncbi.get_taxid_translator([taxon_id])[taxon_id]

        print(f"Downloaded proteome FASTA file for TAXON_ID {taxon_id} ({species_name}) and saved as '{output_filename}'")
        print(f"File size: {file_size_mb:.2f} MB")
        print(f"Number of lines starting with '>': {num_lines}")
    else:
        print(f"Error occurred during download for TAXON_ID {taxon_id}:", response.status_code, response.reason)

# Read TAXON_IDs from file
with open("taxon_id_list.txt", "r") as f:
    taxon_ids = [int(line.strip()) for line in f]

# Create output directory if it does not exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Main script
if __name__ == "__main__":
    for taxon_id in taxon_ids:
        query = f"(taxonomy_id:{taxon_id})"
        params = {"query": query, "format": FORMAT}
        output_filename = os.path.join(OUTPUT_DIR, f"uniprot_proteome_{taxon_id}.fasta")
        download_proteome_fasta(taxon_id, UNIPROT_API_URL, params, output_filename)
