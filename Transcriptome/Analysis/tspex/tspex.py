"""
Tissue-specificity computation using tspex for all available methods

This script computes gene tissue-specificity for all the available methods
using the tspex library and saves the output for each method in a separate
TSV file.

Install tspex library before running the script:
conda install -c conda-forge -c bioconda tspex

Usage:
1. Replace 'input_file' with the path to your expression matrix file
2. Replace 'output_dir' with the desired output directory
3. Run the script
"""

import tspex
import os


def compute_tissue_specificity_for_all_methods(input_file, output_dir):
    """
    Compute tissue-specificity for all available methods in tspex
    and save the output to TSV files in the specified output directory.

    Args:
    input_file (str): Path to the input expression matrix file
    output_dir (str): Path to the output directory
    """

    # Read expression matrix
    expression_matrix = tspex.read_expression_matrix(input_file)

    # List of available methods
    methods = [
        "counts",
        "tau",
        "gini",
        "simpson",
        "shannon_specificity",
        "roku_specificity",
        "tsi",
        "zscore",
        "spm",
        "spm_dpm",
        "js_specificity",
        "js_specificity_dpm",
    ]

    # Compute tissue-specificity for each method and save the output to a TSV file
    for method in methods:
        specificity = tspex.compute_tissue_specificity(expression_matrix, method=method)
        output_file = os.path.join(output_dir, f"tspex_{method}.tsv")
        specificity.to_csv(output_file, sep="\t", index=False)
        print(f"Created {output_file}")


def main():
    # Set input file and output directory
    input_file = "gene_expression.tsv"
    output_dir = "tspex_output"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Compute tissue-specificity for all available methods and save the output
    compute_tissue_specificity_for_all_methods(input_file, output_dir)


if __name__ == "__main__":
    main()
