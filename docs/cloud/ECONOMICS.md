# Cloud Economics and Cost Management

Deploying large-scale bioinformatics pipelines on Google Cloud Platform requires strategic cost management. This document outlines the economics of the METAINFORMANT cloud module, the tradeoffs of different provisioning types, and lessons learned from past billing incidents.

## Compute Economics: Standard vs. Spot Provisioning

Google Cloud offers two primary provisioning models for Virtual Machines. The choice between them is the single largest factor in pipeline compute costs.

### Spot Provisioning (Preemptible)
- **Cost**: Extremely cheap (~70-90% discount).
- **Behavior**: Google can arbitrarily terminate (preempt) these instances at any time if standard capacity is needed elsewhere.
- **Use Case**: Best for short, fault-tolerant batch workloads or highly redundant distributed systems.
- **Incidents**: During the processing of the massive 8,000+ sample RNA-seq dataset, SPOT instances were repeatedly terminated by Google mid-quantification, causing days of pipeline stalling and forcing the orchestrator to constantly reconstruct its SQLite database.

### Standard Provisioning (On-Demand)
- **Cost**: Standard retail pricing (e.g., ~$0.65/hour for `n2-standard-16`).
- **Behavior**: Guaranteed uptime until manually stopped or deleted.
- **Use Case**: **Required for pipeline stability.** The Amalgkit quantifier and the ENA streaming orchestrator perform heavy, hours-long I/O multiplexing across massive outlier samples (e.g., >20GB). These tasks are critically sensitive to unexpected interruption.
- **Verdict**: METAINFORMANT **mandates** Standard VMs for all production deployments. The higher hourly rate is completely offset by the fact that the pipeline actually finishes without multi-day preemption loops.

## Storage Economics: Persistent Disks

Disk I/O and storage capacity are the second largest cost vector.

- **Disk Type**: `pd-standard` (Standard HDD) is perfectly adequate for the pipeline. While `pd-ssd` offers higher IOPS, the pipeline's ENA streaming bottleneck is almost always network-bound, not disk-bound. `pd-standard` is significantly cheaper per GB.
- **Sizing Strategy**: The pipeline deletes FASTQ files immediately after quantification to minimize disk usage. However, for massive >8,000 sample runs encompassing multiple species, the resulting abundance matrices and intermediate SQLite tracking databases accumulate. 
- **Recommendation**: A 4TB `pd-standard` disk provides ample safety margin for edge cases without breaking the bank (approx. $160/month, but prorated down to ~$5 for the 24 hours the pipeline actually runs).
- **Crucial Rule**: Always set `--auto-delete=no` on your massive data disks if you ever need to reconstruct the VM (e.g., migrating from Spot to Standard). If this is not set, deleting the VM obliviates the data disk instantly.

## The API Billing Incident (March 2026)

During a routine cloud execution, the project incurred a sudden **$149.00** charge originating from the Gemini API (`generativelanguage.googleapis.com`).

### Root Cause Analysis
The AI agents orchestrating the project were making millions of high-context automated calls to the commercial Gemini API utilizing the default GCP project billing account. This bypassed the user's intended `$0.00` Google AI Ultra consumer subscription (which uses a different quota mechanism).

### Remediation & Policy
1. **API Disabled**: The `generativelanguage.googleapis.com` API was administratively disabled on the `cryptoptera` GCP project.
2. **Billing Alarm**: Any future AI orchestration must explicitly identify the billing sink being utilized for its LLM queries before initiating autonomous loops.
3. **Consumer Subscriptions**: AI Agents should fallback to utilizing the owner's AI Premium/Ultra credentials via browser automation or strictly bounded service accounts, never unchecked standard commercial API endpoints.

## Pipeline Optimization for Cost

To minimize the absolute cost of a run (Cost = Hourly Rate * Total Run Time):

1. **Max-GB Ceiling**: The orchestrator natively utilizes `PIPELINE_MAX_GB` (default `50.0`). Biological sample datasets larger than 50 GB are almost always extreme anomalies, whole-body pooled multiplexes, or sequencing artifacts. They take exponentially longer to quantify (hours instead of minutes) and consume 90% of the compute budget for 0.1% of the biological payload. **Filtering them saves immense money.**
2. **0.00 GB FASTQ Protection**: The orchestrator is heavily armored against ENA API "ghost" files. If ENA returns a corrupted `0.00 GB` file, the orchestrator immediately crashes that specific sample gracefully instead of passing it to Kallisto. Passing empty buffers to Kallisto causes silent segmentation faults, stalling the pipeline threads indefinitely while you continue to pay hourly compute rates.

## Summary

- **Never** use SPOT instances for this pipeline.
- **Always** use `pd-standard` 4TB disks.
- **Enforce** the 50GB size limits.
- **Monitor** automated API usage carefully.
