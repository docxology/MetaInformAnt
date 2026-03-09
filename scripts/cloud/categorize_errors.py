import subprocess
import re
from collections import Counter

# Connect to VM and just dump the errors
cmd = [
    "gcloud", "compute", "ssh", "metainformant-pipeline", 
    "--zone", "us-central1-a", 
    "--project", "cryptoptera", 
    "--command", 
    'docker exec metainformant-pipeline python3 -c "import sqlite3; db=sqlite3.connect(\\\"/app/output/amalgkit/pipeline_progress.db\\\"); [print(r[0]) for r in db.execute(\\\"SELECT error FROM samples WHERE state=\\\\\'failed\\\\\'\\\").fetchall()]"'
]

result = subprocess.run(cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("Error getting data:", result.stderr)
    exit(1)

counts = Counter()
for line in result.stdout.strip().split('\n'):
    line = line.strip()
    if not line: continue
    
    # Clean up specific IDs to group them together
    clean_line = re.sub(r'(SRR|ERR|DRR)\d+', '{ID}', line)
    
    # Also clean up species paths
    clean_line = re.sub(r'output/amalgkit/[^/]+/work/getfastq', 'output/amalgkit/{SPECIES}/work/getfastq', clean_line)
    
    counts[clean_line] += 1

print("\n--- Error Category Summary ---")
for error, count in counts.most_common():
    print(f"{count:5d} | {error}")
