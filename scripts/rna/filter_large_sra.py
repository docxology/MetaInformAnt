import pandas as pd
import sys
import os

def filter_large_sra(input_path, max_bases=20000000000):
    print(f"Filtering {input_path}...")
    try:
        df = pd.read_csv(input_path, sep='\t')
        
        # Ensure total_bases is numeric
        df['total_bases'] = pd.to_numeric(df['total_bases'], errors='coerce').fillna(0)
        
        # Filter
        original_count = len(df)
        df_filtered = df[df['total_bases'] <= max_bases]
        filtered_count = len(df_filtered)
        removed_count = original_count - filtered_count
        
        print(f"Original samples: {original_count}")
        print(f"filtered samples: {filtered_count}")
        print(f"Removed samples (> {max_bases/1e9} GB): {removed_count}")
        
        if removed_count > 0:
            df_filtered.to_csv(input_path, sep='\t', index=False)
            print(f"Overwrote {input_path} with filtered data.")
        else:
            print("No samples removed.")
            
    except Exception as e:
        print(f"Error filtering {input_path}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    else:
        input_file = "blue/amalgkit/amellifera/work/metadata/metadata_selected.tsv"
    
    if os.path.exists(input_file):
        filter_large_sra(input_file)
    else:
        print(f"File not found: {input_file}")
