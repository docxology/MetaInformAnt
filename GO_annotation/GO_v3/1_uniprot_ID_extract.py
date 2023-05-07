from Bio import SeqIO

input_fasta = "uniprot_bees_proteome.fasta"  # Downloaded UniProt FASTA file
output_file = "bees_uniprot_ids.txt"  # UniProt ID output file

def extract_uniprot_ids(input_fasta, output_file):
    with open(input_fasta, "r") as f, open(output_file, "w") as out:
        for record in SeqIO.parse(f, "fasta"):
            uniprot_id = record.id.split("|")[1]
            out.write(uniprot_id + "\n")

extract_uniprot_ids(input_fasta, output_file)
