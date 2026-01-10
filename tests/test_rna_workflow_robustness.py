from pathlib import Path
from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, plan_workflow

def test_intelligent_redo_override_when_files_exist(tmp_path: Path):
    """Test that plan_workflow overrides redo='yes' to 'no' when SRA files are found."""
    # Setup: Create dummy SRA files in the expected location
    # Note: workflow.py defaults out_dir to work_dir, then looks for getfastq subdir
    work_dir = tmp_path / "work"
    fastq_dir = work_dir / "getfastq"
    fastq_dir.mkdir(parents=True)
    (fastq_dir / "sample1.sra").touch()

    # Config requesting redo: yes
    cfg = AmalgkitWorkflowConfig(
        work_dir=work_dir,
        threads=4,
        per_step={
            "getfastq": {"redo": "yes"}
        }
    )

    # Execute
    steps = plan_workflow(cfg)
    getfastq_params = next(p for s, p in steps if s == "getfastq")

    # Verify: redo should be overridden to 'no' because files exist
    assert getfastq_params["redo"] == "no"


def test_intelligent_redo_respects_user_choice_when_no_files(tmp_path: Path):
    """Test that plan_workflow respects redo='yes' when no SRA files exist."""
    work_dir = tmp_path / "work2"
    
    cfg = AmalgkitWorkflowConfig(
        work_dir=work_dir,
        threads=4,
        per_step={
            "getfastq": {"redo": "yes"}
        }
    )

    steps = plan_workflow(cfg)
    getfastq_params = next(p for s, p in steps if s == "getfastq")
    
    # Verify: redo remains 'yes' because no files found
    assert getfastq_params["redo"] == "yes"


def test_intelligent_redo_keeps_no(tmp_path: Path):
    """Test that redo='no' is kept as 'no' (start state)."""
    work_dir = tmp_path / "work3"
    
    cfg = AmalgkitWorkflowConfig(
        work_dir=work_dir,
        threads=4,
        per_step={
            "getfastq": {"redo": "no"}
        }
    )

    steps = plan_workflow(cfg)
    getfastq_params = next(p for s, p in steps if s == "getfastq")
    
    # Verify: redo remains 'no'
    assert getfastq_params["redo"] == "no"
