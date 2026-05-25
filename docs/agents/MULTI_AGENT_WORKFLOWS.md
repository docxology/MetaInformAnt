# Multi-Agent Workflows: Practical Examples

**Status**: Evolving pattern collection  
**Audience**: Bioinformaticians, workflow engineers, module developers  
**Scope**: Real-world multi-agent workflow templates and compositions

## Introduction

This guide documents concrete examples of multi-agent workflows in METAINFORMANT. Each example demonstrates a coordination pattern (sequential, parallel, conditional) and shows how specific domain agents collaborate through the orchestration layer.

## Table of Example Workflows

| Workflow | Agents Involved | Pattern | Complexity |
|----------|-----------------|---------|------------|
| [1. RNA-seq Amalgkit Pipeline](#1-rna-seq-amalgkit-pipeline) | `rna.amalgkit`, `core.io`, `core.io.paths` | Parallel fan-out + sequential | High |
| [2. GWAS Association + Fine-Mapping](#2-gwas-association--fine-mapping) | `gwas.assoc`, `gwas.finemap`, `gwas.coloc` | Sequential with branching | Medium |
| [3. Multi-Omic Integration](#3-multi-omic-integration) | `multiomics`, `dna`, `rna`, `gwas` | Fan-in aggregation | High |
| [4. Single-Cell Analysis Pipeline](#4-single-cell-analysis-pipeline) | `singlecell.preprocess`, `singlecell.cluster`, `singlecell.trajectory` | Sequential with quality gates | High |
| [5. Batch Download + QC](#5-batch-download--qc) | `core.workflow`, `core.io.download` | Parallel fan-out | Low |
| [6. Cloud VM Orchestration](#6-cloud-vm-orchestration) | `cloud.gcp`, `cloud.docker`, `core.workflow` | Conditional parallel | Medium |
| [7. Metagenomics Taxonomic Profiling](#7-metagenomics-taxonomic-profiling) | `metagenomics.kraken2`, `metagenomics.bracken`, `metagenomics.report` | Sequential + parallel samples | Medium |
| [8. Protein Structure + Function](#8-protein-structure--function) | `protein.alphafold`, `protein.uniprot`, `protein.interpro` | Parallel then merge | High |

---

## 1. RNA-seq Amalgkit Pipeline

**Scale**: 8,300+ samples across 28 Hymenoptera species  
**Agents**:  
- `rna.retrieval` (metadata download from ENA)  
- `rna.amalgkit` (SRA download, FASTQ extraction, quantification)  
- `rna.analysis` (TMM normalization, quality assessment)  

**Pattern**: Parallel fan-out (downloads) → Sequential phases per sample

### Coordination Flow

```mermaid
flowchart TB
    subgraph Orchestrator[RNA Workflow Manager]
        direction LR
        W1[WorkflowManager] --> W2[Phase: download]
        W2 --> W3[Phase: getfastq]
        W3 --> W4[Phase: quantify]
        W4 --> W5[Phase: analyze]
    end

    subgraph Samples[Parallel Samples<br/>8,300 independent items]
        S1[Sample A]
        S2[Sample B]
        S3[Sample C]
    end

    W1 -.->|add_item()| Samples
    W2 -.->|executor.submit()| S1
    W2 -.->|executor.submit()| S2
    W2 -.->|executor.submit()| S3

    S1 -->|robust_download_url| DL1[Download SRA]
    S2 -->|robust_download_url| DL2[Download SRA]
    S3 -->|robust_download_url| DL3[Download SRA]

    DL1 -->|sra_path| EX1[Amalgkit getfastq]
    DL2 -->|sra_path| EX2[Amalgkit getfastq]
    DL3 -->|sra_path| EX3[Amalgkit getfastq]

    EX1 -->|fastq_dir| QT1[Kallisto quant]
    EX2 -->|fastq_dir| QT2[Kallisto quant]
    EX3 -->|fastq_dir| QT3[Kallisto quant]

    QT1 -->|abundance.tsv| FINAL1
    QT2 -->|abundance.tsv| FINAL2
    QT3 -->|abundance.tsv| FINAL3
```

### Code Structure

```python-snippet
# src/metainformant/rna/engine/workflow.py

from metainformant.core.engine.workflow_manager import (
    BasePipelineManager, PipelinePhase, WorkflowManager
)

class RNAWorkflowManager(WorkflowManager):
    """Specialized manager for RNA-seq amalgkit pipeline."""

    def __init__(self, work_dir: Path, config: dict, max_threads: int = 5):
        super().__init__(work_dir, config, max_threads)
        self.phases = [
            PipelinePhase("Download", self._download_phase, filter_fn=self._not_downloaded),
            PipelinePhase("GetFastq", self._getfastq_phase, filter_fn=self._downloaded_not_extracted),
            PipelinePhase("Quantify", self._quantify_phase, filter_fn=self._extracted_not_quantified),
        ]
        self._register_samples()

    def _download_phase(self, manager: BasePipelineManager, items: list[PipelineItem]) -> None:
        # Parallel download via _download_single_sample (uses ThreadPoolExecutor)
        # Updates manager.samples[sid].stage to DOWNLOADED or FAILED
        # ...
```

**Key coordination points**:
- 32 parallel download workers (I/O bound)
- Each sample independent → perfect fan-out
- Metadata (SRA URL, destination path) stored in `item.metadata` and `manager.samples`
- Progress tracked via file size polling (`_update_download_progress()`)

**Result**: 8,300 samples processed across 28 species with per-sample error isolation (one failed download doesn't block others).

### Learning Points

- **Fan-out for independent items**: Each sample can be processed independently → maximal parallelism
- **Progress polling**: TUI updated by checking partial file sizes (`.part` files)
- **Caching**: Already-downloaded files detected and skipped
- **Metadata inheritance**: Manager maintains `self.samples` dict as canonical state alongside `PipelineItem` objects

---

## 2. GWAS Association + Fine-Mapping

**Agents**:
- `gwas.assoc` (PLINK/SAIGE association)
- `gwas.finemap` (FINEMAP, SuSiE)
- `gwas.coloc` (colocalization with eQTL)

**Pattern**: Sequential → conditional branching → parallelize within each stage

### Workflow Description

1. **Association**: Batch-test millions of SNPs across phenotype (CPU-intensive)
2. **Significant subset**: Filter to p < 5e-8; if none found, **early exit**
3. **Fine-mapping**: For each lead SNP, run fine-mapping (parallel over loci)
4. **Coloc** (optional): If GTEx eQTL data available, run coloc per locus

### Coordination Code Outline

```python
class GWASWorkflow(BasePipelineManager):
    def __init__(self, config):
        phases = [
            PipelinePhase("Association", self._association),
            PipelinePhase("Filter", self._filter_significant),
            PipelinePhase("FineMap", self._finemap),
            PipelinePhase("Coloc", self._coloc),
        ]
        super().__init__(phases, max_threads=4)
        # Add items: one GWAS run (single phenotype) is one item
        self.add_item("phenotype_001", metadata={"gwas_config": config})

    def _association(self, manager, items):
        for item in items:
            manager.mark_running(item)
            # Single-item batch: CPU-bound, use executor for parallelization internally
            future = manager.executor.submit(run_association, item.metadata["gwas_config"])
            try:
                sumstats = future.result(timeout=3600)  # 1hr timeout
                item.metadata["sumstats"] = str(sumstats_path)
                manager.mark_done(item)
            except Exception as exc:
                manager.mark_failed(item, str(exc))

    def _filter_significant(self, manager, items):
        """Conditional branching: skip fine-mapping if no hits."""
        for item in items:
            hits = filter_hits(item.metadata["sumstats"])
            if len(hits) == 0:
                manager.mark_done(item, status="No hits — early exit")
                item.metadata["skip_finemap"] = True
            else:
                item.metadata["hits"] = hits
                manager.mark_done(item)

    def _finemap(self, manager, items):
        for item in items:
            if item.metadata.get("skip_finemap"):
                continue  # skip items filtered out

            # Parallel fine-mapping per locus
            loci = group_by_locus(item.metadata["hits"])
            futures = {
                manager.executor.submit(run_finemap_locus, locus): locus
                for locus in loci
            }
            results = {}
            for future in as_completed(futures):
                locus = futures[future]
                try:
                    results[locus] = future.result()
                except Exception as exc:
                    logger.error(f"Fine-map failed for {locus}: {exc}")
            item.metadata["finemap_results"] = results
            manager.mark_done(item)
```

**Key coordination points**:
- `_filter_significant` sets flag on item metadata to skip later phases
- Fine-mapping uses nested parallelism: one item per phenotype, multiple loci per item
- TUI shows overall status per phenotype (not per locus)

---

## 3. Multi-Omic Integration

**Agents**: `multiomics.integrate` (data harmonization), `dna.variants`, `rna.expression`, `gwas.assoc`

**Pattern**: Fan-in aggregation — combine multiple omic layers into joint matrix

### Workflow

```python
class MultiOmicsWorkflow(BasePipelineManager):
    """Integrate DNA variants, RNA expression, and GWAS summary stats."""

    def __init__(self, samples, config):
        phases = [
            PipelinePhase("LoadDNA", self._load_dna),
            PipelinePhase("LoadRNA", self._load_rna),
            PipelinePhase("LoadGWAS", self._load_gwas),
            PipelinePhase("AlignSamples", self._align_samples),  # fan-in
            PipelinePhase("Integrate", self._integrate),
            PipelinePhase("Analyze", self._analyze),
        ]
        super().__init__(phases, max_threads=3)
        for sample_id in samples:
            self.add_item(sample_id)

    def _load_dna(self, manager, items):
        """Load genotype matrices for all samples (parallel)."""
        for item in items:
            future = manager.executor.submit(load_genotype, item.item_id)
            # collect later...

    def _load_rna(self, manager, items):
        """Load expression matrices (parallel)."""
        # similar pattern

    def _load_gwas(self, manager, items):
        """Load GWAS summary stats (single item most likely)."""
        # ...

    def _align_samples(self, manager, items):
        """Fan-in: wait for all preceding phases to complete, then match samples."""
        # Wait for all previous phases to finish
        self.wait_for_phase_completion()  # custom helper

        # Gather all loaded data into shared context
        shared = gather_shared_data(manager.items)
        # Each item now gets aligned data
        for item in items:
            aligned = extract_for_sample(item.item_id, shared)
            item.metadata["aligned"] = aligned
            manager.mark_done(item)

    def _integrate(self, manager, items):
        """Multi-omic integration (parallel per sample)."""
        for item in items:
            future = manager.executor.submit(integrate, item.metadata["aligned"])
            # ...
```

**Key coordination points**:
- First three phases populate interim state (external to items) via futures
- `_align_samples` uses a barrier-like pattern: wait for all prior work, then broadcast aligned data
- Items share a common `manager`-level context: store intermediate data structures at `manager` scope, not per-item

**Implementation tip**: For barrier synchronization, track active phase count:

```python
def track_phase_completion(manager, phase_idx):
    """Custom filter: phase runs only when all previous phases DONE."""
    def filter_fn(item):
        return all(
            get_phase_item(item, i).stage == Stage.DONE
            for i in range(phase_idx)
        )
    return filter_fn
```

---

## 4. Single-Cell Analysis Pipeline

**Agents**: `singlecell.preprocess`, `singlecell.qc`, `singlecell.cluster`, `singlecell.trajectory`

**Pattern**: Sequential with quality gates (early abort on QC failure)

### Workflow

```python
class SingleCellWorkflow(BasePipelineManager):
    """scRNA-seq pipeline with conditional termination."""

    def __init__(self, matrix_path, config):
        phases = [
            PipelinePhase("Preprocess", self._preprocess),      # filter cells/genes
            PipelinePhase("QC", self._qc, filter_fn=self._qc_passable),  # per-item gate
            PipelinePhase("Cluster", self._cluster),
            PipelinePhase("Trajectory", self._trajectory),
        ]
        super().__init__(phases, max_threads=2)
        # Single item: the entire dataset
        self.add_item("dataset", metadata={"matrix": matrix_path})

    def _qc(self, manager, items):
        """Quality gate — can fail item, triggering downstream skip."""
        item = items[0]  # single item
        manager.mark_running(item, status="Running QC")
        qc_report = run_qc(item.metadata["matrix"])
        if qc_report.passed:
            item.metadata["qc"] = qc_report
            manager.mark_done(item)
        else:
            manager.mark_failed(item, f"QC failed: {qc_report.failure_reason}")
            # Subsequent phases' filter_fn will exclude FAILED items automatically

    def _cluster(self, manager, items):
        # Only runs if QC passed (filter excludes FAILED)
        item = items[0]
        clusters = cluster_cells(item.metadata["matrix"], item.metadata["qc"])
        item.metadata["clusters"] = clusters
        manager.mark_done(item)

    def _trajectory(self, manager, items):
        item = items[0]
        # runs only if clustering succeeded (item.stage == DONE)
        trajectory = infer_trajectory(item.metadata["clusters"])
        item.metadata["trajectory"] = trajectory
        manager.mark_done(item)
```

**Key coordination points**:
- Filter functions (`filter_fn` parameter) enable per-item inclusion/exclusion per phase
- Failed QC causes item to be `FAILED` → subsequent phase handlers skip it automatically (default filter = `PENDING` only)
- Quality gates enable fail-fast (don't waste compute on bad data)

---

## 5. Batch Download + QC

**Pattern**: Minimal parallel fan-out, no inter-agent dependencies

### Minimal Config-Driven Workflow

```yaml
# batch_download.yaml
downloads:
  ncbi_sra:
    url: "https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR123/ERR123456/ERR123456.sra"
    filename: "ERR123456.sra"

processing:
  verify:
    type: "checksum"
    algorithm: "md5"
```

```python
from metainformant.core.execution.workflow import download_and_process_data

results = download_and_process_data("batch_download.yaml", output_dir="output/batch1")
# Downloads sequentially (config-driven uses ThreadPoolExecutor internally for parallel downloads if multiple items)

for name, result in results["downloads"].items():
    if result["status"] == "success":
        print(f"{name} → {result['path']}")
```

**Key coordination points**:
- Config-driven approach: no Python code needed, just YAML
- Parallel downloads if multiple entries in `downloads` dict
- Processing step runs after all downloads complete
- Atomic result dictionary captures per-task status

---

## 6. Cloud VM Orchestration

**Agents**: `cloud.gcp` (VM lifecycle), `cloud.docker` (container ops), `rna.amalgkit` (pipeline)

**Pattern**: Conditional parallel (spin up VMs → distribute work → collect results)

### Workflow Description

1. **Provision**: Create N GCP VMs based on `max_samples_per_vm`
2. **Distribute**: Assign sample batches to VMs
3. **Execute**: Each VM runs RNA merger kit locally
4. **Collect**: Copy results from VMs to central storage
5. **Terminate**: Shut down VMs

```python
class CloudRNAWorkflow(BasePipelineManager):
    def __init__(self, samples, vm_config):
        phases = [
            PipelinePhase("ProvisionVMs", self._provision_vms),
            PipelinePhase("Distribute", self._distribute_samples),
            PipelinePhase("ExecuteRemote", self._execute_remote),  # parallel across VMs
            PipelinePhase("Collect", self._collect_results),
            PipelinePhase("Terminate", self._terminate_vms),
        ]
        super().__init__(phases, max_threads=10)
        self.samples = samples
        for s in samples:
            self.add_item(s, metadata={"vm_id": None})

    def _provision_vms(self, manager, items):
        # Create one VM per batch of samples
        batch_size = self.config["samples_per_vm"]
        vm_count = (len(self.samples) + batch_size - 1) // batch_size
        vms = []
        for i in range(vm_count):
            vm = create_gcp_vm(f"rna-worker-{i}")
            vms.append(vm)
            manager.ui.set_footer(f"Provisioned VM {i+1}/{vm_count}")
        manager.config["vms"] = vms
        manager.mark_done(items[0])  # single item representing provisioning phase

    def _distribute_samples(self, manager, items):
        """Assign samples to VMs (round-robin)."""
        vms = manager.config["vms"]
        vm_batches = [[] for _ in vms]
        for i, sample in enumerate(self.samples):
            vm_batch_id = i % len(vms)
            vm_batches[vm_batch_id].append(sample)

        for vm, batch in zip(vms, vm_batches):
            upload_sample_list(vm, batch)
        manager.mark_done(items[0])

    def _execute_remote(self, manager, items):
        # Execute RNA pipeline on each VM in parallel
        vms = manager.config["vms"]
        futures = {
            manager.executor.submit(run_remote_workflow, vm): vm
            for vm in vms
        }
        for future in as_completed(futures):
            vm = futures[future]
            try:
                future.result()
                manager.ui.update(vm.id, status="Completed", color=GREEN)
            except Exception as exc:
                manager.ui.update(vm.id, status="Failed", color=RED)
        manager.mark_done(items[0])

    def _collect_results(self, manager, items):
        # Download output/ directory from each VM
        for vm in manager.config["vms"]:
            download_from_vm(vm, "output/", f"output/vm_{vm.id}/")
        manager.mark_done(items[0])

    def _terminate_vms(self, manager, items):
        for vm in manager.config["vms"]:
            vm.delete()
        manager.mark_done(items[0])
```

**Key coordination points**:
- Single-item phases (`items[0]`) used for infrastructure operations
- Items representing samples ignored during infrastructure phases (filter_fn not needed since only one item exists)
- Manager-level config stores VM references across phase boundaries
- Parallel execution across VMs via executor; TUI shows per-VM progress

---

## 7. Metagenomics Taxonomic Profiling

**Agents**: `metagenomics.kraken2`, `metagenomics.bracken`, `metagenomics.report`

**Pattern**: Sequential sample-level pipeline → parallelize over samples

### Workflow

```python
class MetagenomicsWorkflow(BasePipelineManager):
    def __init__(self, fastq_files, database):
        phases = [
            PipelinePhase("Kraken2", self._kraken2),
            PipelinePhase("Bracken", self._bracken),
            PipelinePhase("Report", self._report),
        ]
        super().__init__(phases, max_threads=4)
        for fastq in fastq_files:
            self.add_item(fastq.stem, metadata={"fastq": str(fastq)})

    def _kraken2(self, manager, items):
        for item in items:
            future = manager.executor.submit(
                run_kraken2,
                item.metadata["fastq"],
                database=manager.config["kraken_db"]
            )
            futures[future] = item
        # Collect kraken reports
        for future in as_completed(futures):
            item = futures[future]
            report = future.result()
            item.metadata["kraken_report"] = report
            manager.mark_done(item)

    def _bracken(self, manager, items):
        # Bracken needs kraken report + database
        for item in items:
            future = manager.executor.submit(
                run_bracken,
                item.metadata["kraken_report"],
                database=manager.config["bracken_db"],
                read_len=150,
            )
            futures[future] = item
        # ...
```

**Key coordination points**:
- Each sample tracked independently across all three tools
- Parallel execution within each phase (4 workers)
- Metadata chain: `fastq → kraken_report → bracken_abundance → report`
- TUI shows progress per sample across all phases (same sample ID appears in each phase listing)

---

## 8. Protein Structure + Function Integration

**Agents**: `protein.alphafold` (structure prediction), `protein.uniprot` (annotation DB), `protein.interpro` (domain detection)

**Pattern**: Parallel prediction → merge annotation → sequential analysis

### Workflow

```python
class ProteinAnnotationWorkflow(BasePipelineManager):
    """Predict structure, annotate function for protein sequences."""

    def __init__(self, sequences: dict[str, str]):
        phases = [
            PipelinePhase("Predict", self._alphafold_predict),  # parallel over sequences
            PipelinePhase("Annotate", self._uniprot_lookup),
            PipelinePhase("Domains", self._interpro_scan),
            PipelinePhase("Integrate", self._integrate_all),
        ]
        super().__init__(phases, max_threads=4)
        for protein_id, seq in sequences.items():
            self.add_item(protein_id, metadata={"sequence": seq})

    def _alphafold_predict(self, manager, items):
        # GPU-bound: limit to 1 worker typically
        gpu_manager = GPUTaskManager()  # separate from ThreadPoolExecutor
        futures = {
            gpu_manager.submit(run_alphafold, item.item_id, item.metadata["sequence"]): item
            for item in items
        }
        for future in as_completed(futures):
            item = futures[future]
            pdb_path = future.result()
            item.metadata["pdb"] = str(pdb_path)
            manager.mark_done(item)

    def _uniprot_lookup(self, manager, items):
        # I/O-bound (API calls): use thread pool
        for item in items:
            future = manager.executor.submit(query_uniprot, item.item_id)
            futures[future] = item
        for future in as_completed(futures):
            item = futures[future]
            annotation = future.result()
            item.metadata["uniprot"] = annotation
            manager.mark_done(item)

    def _interpro_scan(self, manager, items):
        for item in items:
            domains = run_interpro(item.metadata["pdb"])
            item.metadata["domains"] = domains
            manager.mark_done(item)

    def _integrate_all(self, manager, items):
        # Fan-in: all items now have all metadata
        integrated = {}
        for item in items:
            integrated[item.item_id] = {
                "structure": item.metadata["pdb"],
                "annotation": item.metadata["uniprot"],
                "domains": item.metadata["domains"],
            }
        # Save aggregated results
        io.dump_json(integrated, "output/proteins_integrated.json")
        manager.mark_done(items[0])  # single aggregate item
```

**Key coordination points**:
- Mixed resource types: `alphafold` uses GPU (separate executor), others use thread pool
- Metadata aggregation across parallel branches
- Fan-in: individual protein results merged into single JSON file

---

## Workflow Composition Patterns

### Pattern: Sub-Workflow Invocation

Nest workflows by creating a child `BasePipelineManager`:

```python
class ParentWorkflow(BasePipelineManager):
    def __init__(self):
        phases = [PipelinePhase("Child", self._run_child)]
        super().__init__(phases)
        self.add_item("parent_task")

    def _run_child(self, manager, items):
        child = ChildWorkflow(config=self.config)
        child_results = child.run()  # nested execution, independent executor
        items[0].metadata["child_results"] = child_results
        manager.mark_done(items[0])

class ChildWorkflow(BasePipelineManager):
    # ... child phases ...
    pass
```

**Use case**: GWAS pipeline (association → fine-mapping → coloc). Each major stage is its own sub-workflow.

### Pattern: Phase Handlers That Are Themselves Pipelines

Handler runs another `BasePipelineManager` on a subset:

```python
def batch_handler(manager: BasePipelineManager, items: list[PipelineItem]) -> None:
    # Process items as a batch using a sub-pipeline
    batch = BatchWorkflow(items)
    batch.run()
    for item in items:
        item.metadata["batch_result"] = batch.get_item_result(item.item_id)
        manager.mark_done(item)
```

---

## Scaling: From 10 to 10,000 Samples

| Metric | 10 samples | 100 samples | 8,300 samples (RNA) |
|--------|------------|-------------|---------------------|
| Threads | 2 | 8 | 32 (download) / 4 (quant) |
| Est. time | 2 min | 15 min | 3 weeks (distributed) |
| Memory per thread | 2 GB | 2 GB | 2 GB |
| Total RAM | 4 GB | 16 GB | 64-128 GB (phases separate) |
| TUI items | 10 bars | 100 bars | 8,300 bars (one per sample) |

**Tuning for scale**:
- Use phase-specific `max_threads`: downloads use many workers; CPU-intensive phases use few
- TUI still performant with thousands of bars (uses terminal control sequences)
- Monitor disk space: `output/` grows linearly with samples

---

## Monitoring and Debugging Multi-Agent Workflows

### Instrument: Log Milestones

```python
def handler(manager, items):
    logger.info(f"Phase '{phase_name}' starting", extra={"n_items": len(items)})
    for item in items:
        logger.debug(f"Processing {item.item_id}", extra={"metadata": item.metadata})
    logger.info(f"Phase '{phase_name}' completed", extra={"n_done": n_done})
```

**View logs**:

```bash
tail -f output/core/logs/metainformant.log | grep "Phase 'Download'"
```

### Inspect Intermediate State

```python
def debug_handler(manager, items):
    for item in items:
        print(f"Item {item.item_id}:")
        print(f"  stage = {item.stage}")
        print(f"  metadata keys = {list(item.metadata.keys())}")
```

### Checkpoint Serialization

```python
io.dump_json(
    {iid: {"stage": item.stage.value, "metadata": item.metadata}
     for iid, item in manager.items.items()},
    "checkpoint.json"
)
```

Restore:

```python
state = io.load_json("checkpoint.json")
for item_id, item_data in state.items():
    item = manager.items[item_id]
    item.stage = Stage(item_data["stage"])
    item.metadata.update(item_data["metadata"])
```

---

## Next Steps

- **Read [Architecture](ARCHITECTURE.md)** to understand where workflows fit in system
- **Study [Orchestration](ORCHESTRATION.md)** for API reference
- **Review [Safety](SAFETY.md)** for error handling and validation best practices
- **Consult [Communication Protocols](COMMUNICATION_PROTOCOLS.md)** for data-sharing strategies
