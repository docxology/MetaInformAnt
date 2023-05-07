import os
import requests

# Install required libraries
os.system("pip install requests")

# Constants
TAXON_ID = 7460
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/stream"
QUERY = f"(taxonomy_id:{TAXON_ID})"
FORMAT = "fasta"
PARAMS = {"query": QUERY, "format": FORMAT}

# Function to download the proteome FASTA file
def download_proteome_fasta(url, params, output_filename):
    response = requests.get(url, params=params)
    if response.status_code == 200:
        with open(output_filename, "wb") as output_file:
            output_file.write(response.content)
        print(f"Downloaded proteome FASTA file and saved as '{output_filename}'")
    else:
        print("Error occurred during download:", response.status_code, response.reason)
        exit(1)

# Main script
if __name__ == "__main__":
    download_proteome_fasta(UNIPROT_API_URL, PARAMS, "uniprot_bees_proteome.fasta")
