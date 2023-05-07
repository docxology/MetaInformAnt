import os
import numpy as np

def create_histogram(data, bins=10):
    # Create logarithmic bins
    log_bins = np.logspace(np.log10(min(data)), np.log10(max(data)), num=bins+1)

    # Calculate histogram with log bins
    hist, bin_edges = np.histogram(data, bins=log_bins)

    # Print histogram
    for i in range(len(hist)):
        bar = "#" * hist[i]
        bin_range = f"{bin_edges[i]:.2f} - {bin_edges[i+1]:.2f}"
        print(f"{bin_range}: {bar}")

# Define the folder path
folder_path = "proteomes/"

# Initialize a list to store the counts
count_list = []

# Loop through all files in the folder
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)

    # Ensure that the file is a regular file
    if os.path.isfile(file_path):
        with open(file_path, 'r') as f:
            count = 0

            # Count lines starting with ">"
            for line in f:
                if line.startswith(">"):
                    count += 1

        count_list.append(count)

# Calculate summary statistics
if count_list:
    min_count = np.min(count_list)
    max_count = np.max(count_list)
    mean_count = np.mean(count_list)
    median_count = np.median(count_list)
    std_dev = np.std(count_list)

    # Print the summary statistics
    print(f"Summary Statistics:")
    print(f"Minimum count: {min_count}")
    print(f"Maximum count: {max_count}")
    print(f"Mean count: {mean_count:.2f}")
    print(f"Median count: {median_count}")
    print(f"Standard deviation: {std_dev:.2f}")

    # Print the histogram
    print("\nHistogram with Logarithmic Bins:")
    create_histogram(count_list, bins=10)
else:
    print("No valid files were found in the folder.")
