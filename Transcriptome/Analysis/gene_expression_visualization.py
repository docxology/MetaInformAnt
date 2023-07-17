import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy


def load_data(file_path):
    """Load gene expression data from a file and return a pandas DataFrame."""
    gene_expression_data = pd.read_csv(file_path, index_col=0)
    return gene_expression_data


def heatmap(data, output_file):
    """Create a heatmap of gene expression data."""
    plt.figure(figsize=(10, 10))
    sns.heatmap(data, cmap='viridis')
    plt.xlabel("Samples")
    plt.ylabel("Genes")
    plt.title("Gene Expression Heatmap")
    plt.savefig(output_file)
    plt.close()


def clustermap(data, output_file):
    """Create a clustermap of gene expression data."""
    sns.clustermap(data, method='ward', cmap='viridis')
    plt.xlabel("Samples")
    plt.ylabel("Genes")
    plt.title("Gene Expression Clustermap")
    plt.savefig(output_file)
    plt.close()


def pca_plot(data, output_file):
    """Create a PCA plot of gene expression data."""
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data.T)
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', s=100, alpha=0.8)
    plt.xlabel("Principal Component 1")
    plt.ylabel("Principal Component 2")
    plt.title("PCA Plot of Gene Expression Data")
    plt.savefig(output_file)
    plt.close()


if __name__ == '__main__':
    # Load the gene expression data
    gene_expression_matrix = load_data('gene_expression_data.csv')

    # Create visualizations
    heatmap(gene_expression_matrix, 'heatmap.png')
    clustermap(gene_expression_matrix, 'clustermap.png')
    pca_plot(gene_expression_matrix, 'pca_plot.png')

