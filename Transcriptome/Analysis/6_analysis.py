# Read in the meta data
import pandas as pd

# Load the metadata
metadata = pd.read_csv('/.../data/Apis_mellifera.metadata.tsv', sep='\t')

# Display the first few rows of the metadata
metadata.head()


# Load the gene expression data
gene_expression = pd.read_csv('/.../data/Apis_mellifera.tc.tsv', sep='\t')

# Display the first few rows of the gene expression data
gene_expression.head()

# Subset the new metadata to include only the samples that are present in the gene expression data
metadata_new_subset = metadata_new[metadata_new['unique_ID'].isin(gene_expression_transposed['SampleID'])]

# Save the subset metadata to a new file
metadata_new_subset.to_csv('/mnt/data/reduced_metadata.tsv', sep='\t', index=False)

# Confirm the operation
"Metadata subset saved to 'reduced_metadata.tsv'."

# Merge the curate_group information from the reduced metadata into the gene expression data
merged_expression_data = pd.merge(gene_expression_transposed, reduced_metadata[['unique_ID', 'curate_group']], left_on='SampleID', right_on='unique_ID')

# Set the index to be the SampleID
merged_expression_data.set_index('SampleID', inplace=True)

# Remove the 'unique_ID' column as it is now redundant
merged_expression_data.drop(columns=['unique_ID'], inplace=True)

# Display the first few rows of the merged expression data
merged_expression_data.head()







import matplotlib.pyplot as plt
import seaborn as sns
from pandas.tseries.offsets import YearBegin

# Convert the 'published_date' to datetime format
reduced_metadata['published_date'] = pd.to_datetime(reduced_metadata['published_date'])

# Round down the dates to the start of the year
reduced_metadata['year_published'] = reduced_metadata['published_date'].dt.to_period('Y')

# Group by the 'year_published' and 'curate_group', and count the number of samples in each group
tissue_counts_per_year = reduced_metadata.groupby(['year_published', 'curate_group']).size().reset_index(name='sample_count')

# Create a pivot table for better visualization
pivot_tissue_counts = tissue_counts_per_year.pivot('year_published', 'curate_group', 'sample_count')

# Plot the data
plt.figure(figsize=(15, 8))
sns.heatmap(pivot_tissue_counts, cmap="YlGnBu", annot=True, fmt=".0f")
plt.title('Number of Samples Studied Per Tissue Over Time')
plt.xlabel('Tissue')
plt.ylabel('Year Published')
plt.show()


# Fill NaN values with 0 for proper stacking
pivot_tissue_counts_filled = pivot_tissue_counts.fillna(0)

# Plotting the data as a stacked area plot
plt.figure(figsize=(15, 8))
plt.stackplot(pivot_tissue_counts_filled.index, pivot_tissue_counts_filled.transpose().values, labels=pivot_tissue_counts_filled.columns, alpha=0.7)
plt.xlabel('Year Published')
plt.ylabel('Number of Samples')
plt.title('Differential Focus on Different Tissues Over Time')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.show()









# Convert the 'published_date' to datetime format in the full metadata
metadata_new['published_date'] = pd.to_datetime(metadata_new['published_date'])

# Round down the dates to the start of the year
metadata_new['year_published'] = metadata_new['published_date'].dt.to_period('Y')

# Group by the 'year_published' and 'curate_group', and count the number of samples in each group
tissue_counts_per_year_full = metadata_new.groupby(['year_published', 'curate_group']).size().reset_index(name='sample_count')

# Create a pivot table for better visualization
pivot_tissue_counts_full = tissue_counts_per_year_full.pivot('year_published', 'curate_group', 'sample_count')

# Convert the Period to datetime and extract the year
pivot_tissue_counts_full.index = pivot_tissue_counts_full.index.astype(str).astype('datetime64[ns]').year

# Fill NaN values with 0 for proper stacking
pivot_tissue_counts_full_filled = pivot_tissue_counts_full.fillna(0)

# Plotting the data as a stacked area plot
plt.figure(figsize=(15, 8))
plt.stackplot(pivot_tissue_counts_full_filled.index, pivot_tissue_counts_full_filled.transpose().values, labels=pivot_tissue_counts_full_filled.columns, alpha=0.7)
plt.xlabel('Year Published')
plt.ylabel('Number of Samples')
plt.title('Differential Focus on Different Tissues Over Time (Full Metadata)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.show()





#### Pie-chart for samples with expression.

# Convert the 'published_date' to datetime format in the full metadata
metadata_new['published_date'] = pd.to_datetime(metadata_new['published_date'])

# Round down the dates to the start of the year
metadata_new['year_published'] = metadata_new['published_date'].dt.to_period('Y')

# Group by the 'year_published' and 'curate_group', and count the number of samples in each group
tissue_counts_per_year_full = metadata_new.groupby(['year_published', 'curate_group']).size().reset_index(name='sample_count')

# Create a pivot table for better visualization
pivot_tissue_counts_full = tissue_counts_per_year_full.pivot('year_published', 'curate_group', 'sample_count')

# Convert the Period to datetime and extract the year
pivot_tissue_counts_full.index = pivot_tissue_counts_full.index.astype(str).astype('datetime64[ns]').year

# Fill NaN values with 0 for proper stacking
pivot_tissue_counts_full_filled = pivot_tissue_counts_full.fillna(0)

# Plotting the data as a stacked area plot
plt.figure(figsize=(15, 8))
plt.stackplot(pivot_tissue_counts_full_filled.index, pivot_tissue_counts_full_filled.transpose().values, labels=pivot_tissue_counts_full_filled.columns, alpha=0.7)
plt.xlabel('Year Published')
plt.ylabel('Number of Samples')
plt.title('Differential Focus on Different Tissues Over Time (Full Metadata)')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.show()

#### Pie-chart for all samples.
# Count the number of samples in each curate_group in the full metadata
curate_group_counts_full = metadata_new['curate_group'].value_counts()

# Plotting the data as a pie chart
plt.figure(figsize=(10, 10))
plt.pie(curate_group_counts_full, labels=curate_group_counts_full.index, autopct='%1.1f%%', startangle=140)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title('Distribution of Samples Across Tissues (Full Metadata)')
plt.show()


