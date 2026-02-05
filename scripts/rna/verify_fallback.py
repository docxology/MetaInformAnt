import sys
from pathlib import Path

# Add script dir to path
sys.path.insert(0, str(Path(__file__).parent))
from _setup_utils import check_environment_or_exit, ensure_venv_activated

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

# Activate venv
ensure_venv_activated(auto_setup=True)
check_environment_or_exit(auto_setup=True)

from metainformant.core.utils.logging import get_logger
from metainformant.rna.engine.orchestration import run_workflow_for_species

# Setup logging
logger = get_logger("verify_fallback")

config_path = Path("config/amalgkit/amalgkit_pbarbatus_5sample.yaml")

logger.info(f"Running direct verification for {config_path}")

import yaml

# Force metadata to single sample to avoid disk issues
from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig

# Manually read work_dir from yaml because AmalgkitWorkflowConfig might require valid file structure first
with open(config_path) as f:
    raw_config = yaml.safe_load(f)
    work_dir_rel = raw_config.get("work_dir")

# Config loaded relative to CWD
work_dir = Path(work_dir_rel)
metadata_path = work_dir / "metadata" / "metadata.tsv"
logger.info(f"Resolved metadata path from YAML: {metadata_path}")

if metadata_path.exists():
    logger.info(f"Editing metadata to keep only SRR1817176")
    lines = metadata_path.read_text().splitlines()
    header = lines[0]
    # Keep only header and SRR1817176
    new_lines = [header] + [l for l in lines[1:] if "SRR1817176" in l]
    logger.info(f"Filtered to {len(new_lines)-1} samples")
    metadata_path.write_text("\n".join(new_lines) + "\n")
    logger.info("Metadata file updated.")

    # Copy to metadata_selected.tsv since workflow expects it
    selected_path = metadata_path.parent / "metadata_selected.tsv"
    import shutil

    shutil.copy(metadata_path, selected_path)
    logger.info(f"Copied filtered metadata to {selected_path}")
else:
    logger.error(f"Metadata file not found at {metadata_path}!")

# Clean status info to force run
status_file = Path("output/amalgkit/pbarbatus_test5/work/status.info")
if status_file.exists():
    status_file.unlink()
    logger.info("Removed previous status.info to force clean run")

# Remove getfastq dir to ensure fresh start (but keep the SRA file I put there!)
# Wait, if I clean getfastq dir, I lose the SRA file I moved.
# The user wants me to VERIFY FALLBACK.
# Fallback needs SRA file.
# So I should NOT delete the SRA file.
# I will delete everything ELSE.
fastq_dir = Path("output/amalgkit/pbarbatus_test5/fastq/getfastq")
if fastq_dir.exists():
    for p in fastq_dir.iterdir():
        if p.name == "SRR1817176" or p.name == "SRR1817176.sra":
            # Preserve the directory and the file (if at root)
            continue
        if p.is_dir():
            import shutil

            shutil.rmtree(p)
        else:
            p.unlink()
    logger.info("Cleaned getfastq directory (preserved SRR1817176.sra)")

logger.info("Steps: getfastq, integrate, quant, cleanup")

# Run specific steps, bypassing pre-download
# This relies on SRA files being in output/.../fastq/getfastq/ (or sra/)
# Run getfastq first
logger.info("Running getfastq step...")
run_workflow_for_species(config_path, steps=["getfastq"], check=False)

# Check for extracted files
sample_id = "SRR1817176"
fastq_dir_path = fastq_dir / sample_id
work_dir = Path("output/amalgkit/pbarbatus_test5/work")
work_link_dir = work_dir / "getfastq" / sample_id

if not any(fastq_dir_path.glob("*.fastq.gz")):
    logger.error("❌ Fallback extraction failed to produce files!")
    sys.exit(1)

logger.info("✓ Fallback extraction confirmed.")

# Manual Integration (Symlink to work_dir)
logger.info("Manually linking files for quant...")
work_link_dir.mkdir(parents=True, exist_ok=True)
for fq in fastq_dir_path.glob("*.fastq.gz"):
    dest = work_link_dir / fq.name
    if dest.exists():
        dest.unlink()
    dest.symlink_to(fq.resolve())

# Run quant
logger.info("Running quant step...")
results = run_workflow_for_species(config_path, steps=["quant"], check=False)

if results["success"]:
    logger.info("✅ Verification SUCCESS")
    sys.exit(0)
else:
    logger.error("❌ Verification FAILED")
    sys.exit(1)
