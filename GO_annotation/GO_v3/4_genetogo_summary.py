import pandas as pd
from goatools.obo_parser import GODag

obo_file = "go-basic.obo"
gene2go_file = "gene2go.txt"
output_file = "gene_ontology_annotations.txt"

def load_obo_file(obo_file):
    return GODag(obo_file)

def generate_gene_ontology_annotations(gene2go_file, obodag, output_file):
    gene2go_df = pd.read_csv(gene2go_file, sep='\t', header=None, names=['GeneID', 'GO_IDs'])
    gene2go_df['GO_IDs'] = gene2go_df['GO_IDs'].apply(lambda x: x.split(';'))
    gene2go_df = gene2go_df.explode('GO_IDs').rename(columns={'GO_IDs': 'GO_ID'})

    gene2go_df['GO_Description'] = gene2go_df['GO_ID'].apply(lambda x: obodag[x].name if x in obodag else 'NA')

    gene2go_df.to_csv(output_file, sep='\t', index=False)

obodag = load_obo_file(obo_file)
generate_gene_ontology_annotations(gene2go_file, obodag, output_file)
