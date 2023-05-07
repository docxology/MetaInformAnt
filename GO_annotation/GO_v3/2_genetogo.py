import os
import urllib.request
from Bio import SeqIO
from goatools.obo_parser import GODag

input_fasta = "uniprot_bees_proteom.fasta" # Proteome file downloaded from UniProt database.
obo_file = "go-basic.obo"
gaf_file = "goa_uniprot_7460.gaf"
id_file = "bees_uniprot_ids.txt" # UniProt ID file
gene2go_file = "gene2go.txt"

def download_annotation_files(obo_file, gaf_file):
    if not os.path.exists(obo_file):
        url = "http://purl.obolibrary.org/obo/go/go-basic.obo"
        urllib.request.urlretrieve(url, obo_file)

    if not os.path.exists(gaf_file):
        url = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gaf.gz"
        urllib.request.urlretrieve(url, gaf_file+".gz")
        os.system("gunzip "+gaf_file+".gz")

def read_uniprot_ids(id_file):
    with open(id_file, 'r') as f:
        return [line.strip() for line in f]

def generate_gene2go_file(gaf_file, uniprot_ids, gene2go_file):
    with open(gene2go_file, 'w') as f:
        for uniprot_id in uniprot_ids:
            with open(gaf_file, 'r') as g:
                go_terms = [line.split('\t')[4] for line in g if not line.startswith('!') and line.split('\t')[1] == uniprot_id]
                if go_terms:
                    go_terms = list(set(go_terms))
                    f.write(uniprot_id+'\t'+';'.join(go_terms)+'\n')

download_annotation_files(obo_file, gaf_file)
uniprot_ids = read_uniprot_ids(id_file)
generate_gene2go_file(gaf_file, uniprot_ids, gene2go_file)
