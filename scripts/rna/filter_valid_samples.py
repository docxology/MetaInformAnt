
import csv
import sys
from pathlib import Path

def filter_samples(metadata_path, output_path, n=5):
    print(f"Reading {metadata_path}")
    with open(metadata_path, 'r', newline='', encoding='utf-8') as fin:
        reader = csv.DictReader(fin, delimiter='\t')
        rows = list(reader)
    
    print(f"Total samples: {len(rows)}")
    
    valid_samples = []
    for row in rows:
        run_id = row.get('run', '')
        aws_link = row.get('AWS_Link', '')
        size_str = row.get('size', '0')
        
        try:
            size_bytes = float(size_str)
        except ValueError:
            size_bytes = 0
            
        # Filter criteria:
        # 1. AWS Link exists
        # 2. Not lite
        # 3. Reasonable size (>100MB)
        
        is_lite = '.lite' in aws_link or aws_link.endswith('.lite')
        
        if aws_link and not is_lite and size_bytes > 100 * 1024 * 1024:
            valid_samples.append(row)
            
        if len(valid_samples) >= n:
            break
            
    print(f"Found {len(valid_samples)} valid samples")
    if not valid_samples:
        print("WARNING: No valid samples found!")
        # Fallback: try relaxed criteria (size > 0)
        for row in rows:
            run_id = row.get('run', '')
            aws_link = row.get('AWS_Link', '')
            if aws_link and '.lite' not in aws_link:
                 valid_samples.append(row)
                 if len(valid_samples) >= n:
                     break
        print(f"Relaxed check found {len(valid_samples)} samples")

    if not valid_samples:
        print("ERROR: Could not find any suitable samples.")
        sys.exit(1)
        
    print(f"Writing {len(valid_samples)} samples to {output_path}")
    with open(output_path, 'w', newline='', encoding='utf-8') as fout:
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(valid_samples)
        
    print("Done")

if __name__ == "__main__":
    base_dir = Path("/Users/mini/Documents/GitHub/metainformant/output/amalgkit/pbarbatus_test5/work/metadata")
    filter_samples(base_dir / "metadata.tsv", base_dir / "metadata_selected.tsv")
