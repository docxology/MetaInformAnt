import collections
import subprocess

# Run the query
cmd = [
    "gcloud",
    "compute",
    "ssh",
    "metainformant-pipeline",
    "--zone",
    "us-central1-a",
    "--project",
    "cryptoptera",
    "--command",
    'docker exec metainformant-pipeline python3 -c "import sqlite3; db=sqlite3.connect(\\"/app/output/amalgkit/pipeline_progress.db\\"); [print(f\\"{r[0]}|{r[1]}|{r[2]}\\") for r in db.execute(\\"SELECT species, state, count(*) FROM samples GROUP BY species, state ORDER BY species, state;\\").fetchall()]"',
]

result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print("Error:")
    print(result.stderr)
    exit(1)

# Parse output
data = collections.defaultdict(collections.Counter)
all_states = set()

for line in result.stdout.strip().split("\n"):
    if not line:
        continue

    parts = line.split("|")
    if len(parts) == 3:
        sp, st, ct = parts
        data[sp][st] = int(ct)
        data[sp]["total"] += int(ct)
        all_states.add(st)

# Get sorted states (except 'total')
sorted_states = sorted(list(all_states))

# Print header
header_cols = ["Species"] + [st.capitalize() for st in sorted_states] + ["Total"]
print("| " + " | ".join(header_cols) + " |")
print("|" + "|".join(["---"] * len(header_cols)) + "|")

# Print rows
total_counts = collections.Counter()
for sp in sorted(data.keys()):
    row_cols = [sp]
    for st in sorted_states:
        val = data[sp].get(st, 0)
        row_cols.append(str(val))
        total_counts[st] += val
    total_val = data[sp]["total"]
    row_cols.append(str(total_val))
    total_counts["total"] += total_val
    print("| " + " | ".join(row_cols) + " |")

# Print total row
total_cols = ["**TOTAL**"]
for st in sorted_states:
    total_cols.append(f"**{total_counts[st]}**")
total_cols.append(f"**{total_counts['total']}**")
print("| " + " | ".join(total_cols) + " |")
