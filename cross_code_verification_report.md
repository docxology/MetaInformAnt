# Cross-Code Verification Report

> Historical snapshot: this generated report is retained for provenance and may
> not describe the current checkout. Regenerate it with
> `uv run python scripts/verify_documentation_code.py`; the default output is
> now `output/cross_code_verification_report.md`.

**Total violations found:** 3141

## Summary by Issue Type

- **AttributeError**: 1442 occurrence(s)
- **ImportError**: 1263 occurrence(s)
- **ModuleNotFoundError**: 7 occurrence(s)
- **SyntaxError**: 429 occurrence(s)

## Detailed Violations

| File | Line | Type | Issue | Code Example |
|------|------|------|-------|--------------|
| docs/COMPARISON_GUIDES.md | 65 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/COMPARISON_GUIDES.md | 110 | ImportError | Cannot import 'coalescent' from module 'metainformant.math' | `from metainformant.math import coalescent` |
| docs/COMPARISON_GUIDES.md | 110 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/COMPARISON_GUIDES.md | 110 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/COMPARISON_GUIDES.md | 202 | ImportError | Cannot import 'ml' from module 'metainformant' | `from metainformant import ml` |
| docs/COMPARISON_GUIDES.md | 319 | ImportError | Cannot import 'syntactic' from module 'metainformant.information' | `from metainformant.information import syntactic` |
| docs/COMPARISON_GUIDES.md | 319 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information` |
| docs/DEVELOPMENT.md | 173 | ImportError | Cannot import 'gc_content' from module 'metainformant.dna.composition' | `from metainformant.dna.composition import gc_content` |
| docs/DEVELOPMENT.md | 173 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.composition` |
| docs/DEVELOPMENT.md | 216 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `"""Brief module docstring."""    from . import submodule1, s` |
| docs/DEVELOPMENT.md | 282 | ImportError | Cannot import 'setup_logging' from module 'metainformant.core.utils.logging' | `from metainformant.core.utils.logging import setup_logging` |
| docs/DEVELOPMENT.md | 282 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/FAQ.md | 74 | ImportError | Cannot import 'run_parallel' from module 'metainformant.core.parallel' | `from metainformant.core.parallel import run_parallel` |
| docs/FAQ.md | 74 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/FAQ.md | 74 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.parallel` |
| docs/FAQ.md | 104 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/FAQ.md | 104 | ImportError | Cannot import 'coalescent' from module 'metainformant.math' | `from metainformant.math import coalescent` |
| docs/FAQ.md | 104 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/FAQ.md | 104 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/FAQ.md | 129 | ImportError | Cannot import 'analyze_sequence_information' from module 'metainformant.information.analysis' | `from metainformant.information.analysis import analyze_seque` |
| docs/FAQ.md | 129 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/FAQ.md | 129 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.analysis` |
| docs/FAQ.md | 143 | ImportError | Cannot import 'codon' from module 'metainformant.dna' | `from metainformant.dna import codon` |
| docs/FAQ.md | 143 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/FAQ.md | 176 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/FAQ.md | 195 | ImportError | Cannot import 'structure' from module 'metainformant.gwas' | `from metainformant.gwas import structure` |
| docs/FAQ.md | 195 | ImportError | Cannot import 'association' from module 'metainformant.gwas' | `from metainformant.gwas import association` |
| docs/FAQ.md | 195 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/FAQ.md | 213 | ImportError | Cannot import 'download' from module 'metainformant.gwas' | `from metainformant.gwas import download` |
| docs/FAQ.md | 213 | ImportError | Cannot import 'parsing' from module 'metainformant.gwas' | `from metainformant.gwas import parsing` |
| docs/FAQ.md | 213 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/FAQ.md | 240 | ImportError | Cannot import 'syntactic' from module 'metainformant.information' | `from metainformant.information import syntactic` |
| docs/FAQ.md | 240 | ImportError | Cannot import 'semantic' from module 'metainformant.information' | `from metainformant.information import semantic` |
| docs/FAQ.md | 240 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information` |
| docs/FAQ.md | 273 | ImportError | Cannot import 'ml' from module 'metainformant' | `from metainformant import ml` |
| docs/FAQ.md | 293 | ImportError | Cannot import 'SMOTE' from module 'imblearn.over_sampling' | `from imblearn.over_sampling import SMOTE` |
| docs/FAQ.md | 293 | ImportError | Cannot import 'ml' from module 'metainformant' | `from metainformant import ml` |
| docs/FAQ.md | 312 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/FAQ.md | 329 | ImportError | Cannot import 'JsonCache' from module 'metainformant.core.cache' | `from metainformant.core.cache import JsonCache` |
| docs/FAQ.md | 329 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.cache` |
| docs/FAQ.md | 464 | ImportError | Cannot import 'run_config_based_workflow' from module 'metainformant.core.workflow' | `from metainformant.core.workflow import run_config_based_wor` |
| docs/FAQ.md | 464 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.workflow` |
| docs/FAQ.md | 464 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.load_json` |
| docs/FAQ.md | 464 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.analysis.batch_entropy_analysis` |
| docs/FAQ.md | 464 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.dump_json` |
| docs/FAQ.md | 497 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `# Install in notebook !curl -LsSf https://astral.sh/uv/insta` |
| docs/GETTING_STARTED.md | 37 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/GETTING_STARTED.md | 128 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/GETTING_STARTED.md | 166 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/GETTING_STARTED.md | 206 | ImportError | Cannot import 'setup_logging' from module 'metainformant.core.utils.logging' | `from metainformant.core.utils.logging import setup_logging` |
| docs/GETTING_STARTED.md | 206 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/GETTING_STARTED.md | 231 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/GETTING_STARTED.md | 231 | ImportError | Cannot import 'setup_logging' from module 'metainformant.core.utils.logging' | `from metainformant.core.utils.logging import setup_logging` |
| docs/GETTING_STARTED.md | 231 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/GETTING_STARTED.md | 231 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/GETTING_STARTED.md | 231 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/INTEGRATION.md | 53 | SyntaxError | Invalid Python syntax: unexpected indent at line 13 | `""" Minimal DNA → GWAS → Visualization pipeline. Assumes VCF` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'preprocessing' from module 'metainformant.dna' | `from metainformant.dna import preprocessing` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'composite' from module 'metainformant.visualization' | `from metainformant.visualization import composite` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'manhattan' from module 'metainformant.visualization.interactive' | `from metainformant.visualization.interactive import manhatta` |
| docs/INTEGRATION.md | 122 | ImportError | Cannot import 'expression' from module 'metainformant.rna' | `from metainformant.rna import expression` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/INTEGRATION.md | 122 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.interactive` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'parse_vcf' from module 'metainformant.dna.variants' | `from metainformant.dna.variants import parse_vcf` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'filter_maf' from module 'metainformant.dna.variants' | `from metainformant.dna.variants import filter_maf` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'to_genotypes' from module 'metainformant.dna.variants' | `from metainformant.dna.variants import to_genotypes` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'linear_regression' from module 'metainformant.gwas.analysis' | `from metainformant.gwas.analysis import linear_regression` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'manhattan' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import manhattan` |
| docs/INTEGRATION.md | 276 | ImportError | Cannot import 'qq' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import qq` |
| docs/INTEGRATION.md | 276 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variants` |
| docs/INTEGRATION.md | 276 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis` |
| docs/INTEGRATION.md | 276 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/INTEGRATION.md | 305 | ImportError | Cannot import 'checkpoint' from module 'metainformant.core' | `from metainformant.core import checkpoint` |
| docs/INTEGRATION.md | 305 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/INTEGRATION.md | 334 | ImportError | Cannot import 'map_over_chromosomes' from module 'metainformant.core.parallel' | `from metainformant.core.parallel import map_over_chromosomes` |
| docs/INTEGRATION.md | 334 | ImportError | Cannot import 'submit_batch' from module 'metainformant.cloud' | `from metainformant.cloud import submit_batch` |
| docs/INTEGRATION.md | 334 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/INTEGRATION.md | 398 | ImportError | Cannot import 'matrix' from module 'metainformant.multiomics' | `from metainformant.multiomics import matrix` |
| docs/INTEGRATION.md | 398 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/INTEGRATION.md | 490 | ImportError | Cannot import 'parse_vcf' from module 'metainformant.dna.variants' | `from metainformant.dna.variants import parse_vcf` |
| docs/INTEGRATION.md | 490 | ImportError | Cannot import 'dosage_matrix' from module 'metainformant.dna.variants' | `from metainformant.dna.variants import dosage_matrix` |
| docs/INTEGRATION.md | 490 | ImportError | Cannot import 'association_map' from module 'metainformant.gwas.analysis' | `from metainformant.gwas.analysis import association_map` |
| docs/INTEGRATION.md | 490 | ImportError | Cannot import 'heatmap_pvalues' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import heatmap_pva` |
| docs/INTEGRATION.md | 490 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variants` |
| docs/INTEGRATION.md | 490 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis` |
| docs/INTEGRATION.md | 490 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/INTEGRATION.md | 520 | ImportError | Cannot import 'meta_forest' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import meta_forest` |
| docs/INTEGRATION.md | 520 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/INTEGRATION.md | 520 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/INTEGRATION.md | 549 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/INTEGRATION.md | 573 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dashboards` |
| docs/INTEGRATION.md | 599 | ImportError | Cannot import 'genomics' from module 'metainformant.visualization.interactive' | `from metainformant.visualization.interactive import genomics` |
| docs/INTEGRATION.md | 599 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.interactive` |
| docs/ORCHESTRATION.md | 112 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/ORCHESTRATION.md | 137 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/ORCHESTRATION.md | 137 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.workflow.workflow_config` |
| docs/ORCHESTRATION.md | 137 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/ORCHESTRATION.md | 159 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/ORCHESTRATION.md | 159 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/TROUBLESHOOTING.md | 14 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/TROUBLESHOOTING.md | 24 | ImportError | Cannot import 'parse_vcf_full' from module 'metainformant.gwas.quality' | `from metainformant.gwas.quality import parse_vcf_full` |
| docs/TROUBLESHOOTING.md | 24 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.quality` |
| docs/TROUBLESHOOTING.md | 43 | ImportError | Cannot import 'autocast' from module 'torch.cuda.amp' | `from torch.cuda.amp import autocast` |
| docs/TROUBLESHOOTING.md | 65 | ImportError | Cannot import 'ParallelProcessor' from module 'metainformant.core.parallel' | `from metainformant.core.parallel import ParallelProcessor` |
| docs/TROUBLESHOOTING.md | 65 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.parallel` |
| docs/TROUBLESHOOTING.md | 176 | ImportError | Cannot import 'validate_vcf' from module 'metainformant.gwas.quality' | `from metainformant.gwas.quality import validate_vcf` |
| docs/TROUBLESHOOTING.md | 176 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.quality` |
| docs/TROUBLESHOOTING.md | 234 | ImportError | Cannot import 'download_sra_run' from module 'metainformant.rna.sra_extraction' | `from metainformant.rna.sra_extraction import download_sra_ru` |
| docs/TROUBLESHOOTING.md | 234 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.sra_extraction` |
| docs/TROUBLESHOOTING.md | 264 | ImportError | Cannot import 'JsonCache' from module 'metainformant.core.cache' | `from metainformant.core.cache import JsonCache` |
| docs/TROUBLESHOOTING.md | 264 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.cache` |
| docs/TROUBLESHOOTING.md | 279 | ImportError | Cannot import 'validate_vcf' from module 'metainformant.gwas.quality' | `from metainformant.gwas.quality import validate_vcf` |
| docs/TROUBLESHOOTING.md | 279 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.quality` |
| docs/TROUBLESHOOTING.md | 304 | ImportError | Cannot import 'detect_quality_encoding' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import detect_quality_encod` |
| docs/TROUBLESHOOTING.md | 304 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/TROUBLESHOOTING.md | 355 | ImportError | Cannot import 'validate_config_schema' from module 'metainformant.core.validation' | `from metainformant.core.validation import validate_config_sc` |
| docs/TROUBLESHOOTING.md | 355 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.validation` |
| docs/TROUBLESHOOTING.md | 427 | ImportError | Cannot import 'get_directory_size' from module 'metainformant.core.paths' | `from metainformant.core.paths import get_directory_size` |
| docs/TROUBLESHOOTING.md | 427 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.paths` |
| docs/TUTORIALS.md | 43 | ImportError | Cannot import 'distances' from module 'metainformant.dna' | `from metainformant.dna import distances` |
| docs/TUTORIALS.md | 43 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/TUTORIALS.md | 83 | ImportError | Cannot import 'coalescent' from module 'metainformant.math' | `from metainformant.math import coalescent` |
| docs/TUTORIALS.md | 83 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/TUTORIALS.md | 103 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/TUTORIALS.md | 146 | ImportError | Cannot import 'association' from module 'metainformant.gwas' | `from metainformant.gwas import association` |
| docs/TUTORIALS.md | 146 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/TUTORIALS.md | 179 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/TUTORIALS.md | 205 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/TUTORIALS.md | 243 | ImportError | Cannot import 'quant' from module 'metainformant.rna.engine.workflow_steps' | `from metainformant.rna.engine.workflow_steps import quant` |
| docs/TUTORIALS.md | 243 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow_steps` |
| docs/TUTORIALS.md | 267 | ImportError | Cannot import 'analysis' from module 'metainformant.information' | `from metainformant.information import analysis` |
| docs/TUTORIALS.md | 267 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information` |
| docs/TUTORIALS.md | 291 | ImportError | Cannot import 'workflows' from module 'metainformant.information' | `from metainformant.information import workflows` |
| docs/TUTORIALS.md | 291 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information` |
| docs/TUTORIALS.md | 323 | ImportError | Cannot import 'multiomics' from module 'metainformant' | `from metainformant import multiomics` |
| docs/TUTORIALS.md | 353 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/TUTORIALS.md | 353 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/TUTORIALS.md | 381 | ImportError | Cannot import 'ml' from module 'metainformant' | `from metainformant import ml` |
| docs/TUTORIALS.md | 449 | ImportError | Cannot import 'networks' from module 'metainformant' | `from metainformant import networks` |
| docs/TUTORIALS.md | 506 | ImportError | Cannot import 'run_config_based_workflow' from module 'metainformant.core.workflow' | `from metainformant.core.workflow import run_config_based_wor` |
| docs/TUTORIALS.md | 506 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.workflow` |
| docs/TUTORIALS.md | 506 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.load_json` |
| docs/TUTORIALS.md | 506 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequences.read_fasta` |
| docs/TUTORIALS.md | 506 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.analysis.batch_entropy_analysis` |
| docs/TUTORIALS.md | 506 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.dump_json` |
| docs/TUTORIALS.md | 558 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/TUTORIALS.md | 558 | ImportError | Cannot import 'analyze_sequence_information' from module 'metainformant.information.analysis' | `from metainformant.information.analysis import analyze_seque` |
| docs/TUTORIALS.md | 558 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/TUTORIALS.md | 558 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.analysis` |
| docs/TUTORIALS.md | 590 | ImportError | Cannot import 'coalescent' from module 'metainformant.math' | `from metainformant.math import coalescent` |
| docs/TUTORIALS.md | 590 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/UV_SETUP.md | 437 | ImportError | Cannot import 'detect_filesystem_type' from module 'metainformant.core.filesystem' | `from metainformant.core.filesystem import detect_filesystem_` |
| docs/UV_SETUP.md | 437 | ImportError | Cannot import 'supports_symlinks' from module 'metainformant.core.filesystem' | `from metainformant.core.filesystem import supports_symlinks` |
| docs/UV_SETUP.md | 437 | ImportError | Cannot import 'get_uv_cache_dir' from module 'metainformant.core.filesystem' | `from metainformant.core.filesystem import get_uv_cache_dir` |
| docs/UV_SETUP.md | 437 | ImportError | Cannot import 'get_venv_location' from module 'metainformant.core.filesystem' | `from metainformant.core.filesystem import get_venv_location` |
| docs/UV_SETUP.md | 437 | ImportError | Cannot import 'is_fat_filesystem' from module 'metainformant.core.filesystem' | `from metainformant.core.filesystem import is_fat_filesystem` |
| docs/UV_SETUP.md | 437 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.filesystem` |
| docs/agents/ARCHITECTURE.md | 463 | SyntaxError | Invalid Python syntax: unexpected indent at line 6 | `phases = [        PipelinePhase("Fetch", fetch_handler),    ` |
| docs/agents/BEST_PRACTICES.md | 83 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.config` |
| docs/agents/BEST_PRACTICES.md | 170 | ImportError | Cannot import 'RNAWorkflowManager' from module 'metainformant.rna.engine.workflow' | `from metainformant.rna.engine.workflow import RNAWorkflowMan` |
| docs/agents/BEST_PRACTICES.md | 170 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.engine.workflow_manager` |
| docs/agents/BEST_PRACTICES.md | 170 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/agents/BEST_PRACTICES.md | 279 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.parallel` |
| docs/agents/BEST_PRACTICES.md | 300 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `del large_dataframe   import gc; gc.collect()` |
| docs/agents/BEST_PRACTICES.md | 331 | ImportError | Cannot import 'setup_logging' from module 'metainformant.core.utils.logging' | `from metainformant.core.utils.logging import setup_logging` |
| docs/agents/BEST_PRACTICES.md | 331 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/agents/BEST_PRACTICES.md | 368 | ImportError | Cannot import 'start_metrics_server' from module 'metainformant.monitoring' | `from metainformant.monitoring import start_metrics_server` |
| docs/agents/BEST_PRACTICES.md | 368 | AttributeError | 'metainformant' has no attribute 'monitoring' (module not found) | `metainformant.monitoring` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 126 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.config` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 355 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Good: ephemeral, in-memory, no cleanup needed    item.meta` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 362 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Good: survives process restart, visible in output/    io.d` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 369 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `item.metadata["produced_by"] = "rna.amalgkit.quantify"    it` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 379 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Good: metadata only has path    item.metadata["matrix_path` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 402 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# BAD: huge in-memory DataFrame increases memory pressure   ` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 409 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# BAD: unclear which phase provides which key    value = ite` |
| docs/agents/COMMUNICATION_PROTOCOLS.md | 416 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# BAD: race condition if multiple items write to same dict  ` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 316 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 317 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 317 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.core.io import ...` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 318 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 321 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 321 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 325 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 392 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 392 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/agents/DOCUMENTATION_REVIEW_REPORT.md | 418 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/MULTI_AGENT_WORKFLOWS.md | 78 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 19 at line 22 | `# src/metainformant/rna/engine/workflow.py  from metainforma` |
| docs/agents/MULTI_AGENT_WORKFLOWS.md | 358 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/agents/ORCHESTRATION.md | 46 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.engine.workflow_manager` |
| docs/agents/ORCHESTRATION.md | 124 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.engine.workflow_manager` |
| docs/agents/ORCHESTRATION.md | 169 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.engine.workflow_manager` |
| docs/agents/ORCHESTRATION.md | 372 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.parallel` |
| docs/agents/ORCHESTRATION.md | 462 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/agents/ORCHESTRATION.md | 475 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/agents/ORCHESTRATION.md | 484 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/agents/ORCHESTRATION.md | 534 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.config` |
| docs/agents/ORCHESTRATION.md | 548 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.engine.workflow_manager` |
| docs/agents/SAFETY.md | 99 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution.workflow` |
| docs/agents/SAFETY.md | 240 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.data` |
| docs/agents/SAFETY.md | 255 | ImportError | Cannot import 'retry' from module 'tenacity' | `from tenacity import retry` |
| docs/agents/SAFETY.md | 255 | ImportError | Cannot import 'stop_after_attempt' from module 'tenacity' | `from tenacity import stop_after_attempt` |
| docs/agents/SAFETY.md | 255 | ImportError | Cannot import 'wait_exponential' from module 'tenacity' | `from tenacity import wait_exponential` |
| docs/agents/SAFETY.md | 493 | ImportError | Cannot import 'Counter' from module 'prometheus_client' | `from prometheus_client import Counter` |
| docs/agents/SAFETY.md | 493 | ImportError | Cannot import 'Histogram' from module 'prometheus_client' | `from prometheus_client import Histogram` |
| docs/agents/SPEC.md | 38 | SyntaxError | Invalid Python syntax: invalid character '→' (U+2192) at line 33 | `from metainformant.core.engine.workflow_manager import (    ` |
| docs/agents/TROUBLESHOOTING.md | 167 | SyntaxError | Invalid Python syntax: unindent does not match any outer indentation level at line 3 | `if not output_path.exists():       raise FileNotFoundError(f` |
| docs/agents/TROUBLESHOOTING.md | 246 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `ImportError: cannot import name 'X'` |
| docs/agents/rules/core.md | 55 | ImportError | Cannot import 'load_mapping_from_file' from module 'metainformant.core.config' | `from metainformant.core.config import load_mapping_from_file` |
| docs/agents/rules/core.md | 55 | ImportError | Cannot import 'apply_env_overrides' from module 'metainformant.core.config' | `from metainformant.core.config import apply_env_overrides` |
| docs/agents/rules/core.md | 55 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.config` |
| docs/agents/rules/core.md | 73 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/core.md | 93 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/agents/rules/core.md | 93 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/agents/rules/core.md | 109 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/core.md | 109 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/core.md | 129 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/agents/rules/core.md | 192 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/core.md | 192 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/agents/rules/core.md | 192 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/core.md | 192 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils` |
| docs/agents/rules/core.md | 213 | ModuleNotFoundError | Module 'optional_library' not found in project or standard library | `import optional_library` |
| docs/agents/rules/dna.md | 47 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence` |
| docs/agents/rules/dna.md | 74 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.phylogeny` |
| docs/agents/rules/dna.md | 87 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/agents/rules/dna.md | 120 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external` |
| docs/agents/rules/dna.md | 204 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/dna.md | 217 | ImportError | Cannot import 'SeqIO' from module 'Bio' | `from Bio import SeqIO` |
| docs/agents/rules/dna.md | 233 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/dna.md | 233 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ecology.md | 63 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ecology.md | 76 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/ecology.md | 76 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/epigenome.md | 35 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome.assays` |
| docs/agents/rules/epigenome.md | 73 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/epigenome.md | 86 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/epigenome.md | 86 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/gwas.md | 46 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis` |
| docs/agents/rules/gwas.md | 64 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis` |
| docs/agents/rules/gwas.md | 103 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization` |
| docs/agents/rules/gwas.md | 146 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/gwas.md | 159 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/gwas.md | 159 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/index.md | 74 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/agents/rules/index.md | 77 | ImportError | Cannot import 'align' from module 'metainformant.dna' | `from metainformant.dna import align` |
| docs/agents/rules/index.md | 77 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/agents/rules/index.md | 87 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/agents/rules/information.md | 35 | ImportError | Cannot import 'syntactic' from module 'metainformant.information.metrics' | `from metainformant.information.metrics import syntactic` |
| docs/agents/rules/information.md | 35 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics` |
| docs/agents/rules/information.md | 114 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/information.md | 127 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/information.md | 127 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/life_events.md | 44 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models` |
| docs/agents/rules/life_events.md | 93 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/life_events.md | 107 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/life_events.md | 107 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/longread.md | 40 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.io` |
| docs/agents/rules/longread.md | 40 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.quality` |
| docs/agents/rules/longread.md | 40 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.analysis` |
| docs/agents/rules/longread.md | 40 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.assembly` |
| docs/agents/rules/longread.md | 40 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.workflow` |
| docs/agents/rules/math.md | 81 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics` |
| docs/agents/rules/math.md | 106 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics` |
| docs/agents/rules/math.md | 120 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics` |
| docs/agents/rules/math.md | 134 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.quantitative_genetics` |
| docs/agents/rules/math.md | 147 | ImportError | Cannot import 'logistic_map' from module 'metainformant.math' | `from metainformant.math import logistic_map` |
| docs/agents/rules/math.md | 147 | ImportError | Cannot import 'lotka_volterra_step' from module 'metainformant.math' | `from metainformant.math import lotka_volterra_step` |
| docs/agents/rules/math.md | 147 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/agents/rules/math.md | 170 | ImportError | Cannot import 'simulate_generations' from module 'metainformant.math.selection_experiments' | `from metainformant.math.selection_experiments import simulat` |
| docs/agents/rules/math.md | 170 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.selection_experiments` |
| docs/agents/rules/math.md | 187 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/math.md | 200 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/math.md | 200 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/menu.md | 21 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core` |
| docs/agents/rules/menu.md | 21 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui` |
| docs/agents/rules/metabolomics.md | 32 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/metabolomics.md | 32 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/agents/rules/metabolomics.md | 42 | ImportError | Cannot import 'load_mapping_from_file' from module 'metainformant.core.config' | `from metainformant.core.config import load_mapping_from_file` |
| docs/agents/rules/metabolomics.md | 42 | ImportError | Cannot import 'apply_env_overrides' from module 'metainformant.core.config' | `from metainformant.core.config import apply_env_overrides` |
| docs/agents/rules/metabolomics.md | 42 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.config` |
| docs/agents/rules/metagenomics.md | 29 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.amplicon` |
| docs/agents/rules/metagenomics.md | 29 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.shotgun` |
| docs/agents/rules/metagenomics.md | 29 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.functional` |
| docs/agents/rules/ml.md | 45 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models` |
| docs/agents/rules/ml.md | 62 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.features` |
| docs/agents/rules/ml.md | 83 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.evaluation` |
| docs/agents/rules/ml.md | 101 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ml.md | 115 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/ml.md | 115 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/multiomics.md | 31 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.analysis` |
| docs/agents/rules/multiomics.md | 50 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/multiomics.md | 63 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/multiomics.md | 63 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/networks.md | 35 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks.analysis` |
| docs/agents/rules/networks.md | 65 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks.analysis` |
| docs/agents/rules/networks.md | 82 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/networks.md | 95 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/networks.md | 95 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ontology.md | 55 | ImportError | Cannot import 'Term' from module 'metainformant.ontology.types' | `from metainformant.ontology.types import Term` |
| docs/agents/rules/ontology.md | 55 | ImportError | Cannot import 'Ontology' from module 'metainformant.ontology.types' | `from metainformant.ontology.types import Ontology` |
| docs/agents/rules/ontology.md | 55 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.types` |
| docs/agents/rules/ontology.md | 88 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ontology.md | 101 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/ontology.md | 101 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/ontology.md | 130 | ImportError | Cannot import 'set_cache_enabled' from module 'metainformant.ontology' | `from metainformant.ontology import set_cache_enabled` |
| docs/agents/rules/ontology.md | 130 | ImportError | Cannot import 'set_cache_ttl' from module 'metainformant.ontology' | `from metainformant.ontology import set_cache_ttl` |
| docs/agents/rules/ontology.md | 130 | ImportError | Cannot import 'clear_cache' from module 'metainformant.ontology' | `from metainformant.ontology import clear_cache` |
| docs/agents/rules/ontology.md | 130 | ImportError | Cannot import 'ancestors' from module 'metainformant.ontology.query' | `from metainformant.ontology.query import ancestors` |
| docs/agents/rules/ontology.md | 130 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology` |
| docs/agents/rules/ontology.md | 130 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.query` |
| docs/agents/rules/pharmacogenomics.md | 30 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles` |
| docs/agents/rules/pharmacogenomics.md | 30 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations` |
| docs/agents/rules/pharmacogenomics.md | 30 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical` |
| docs/agents/rules/phenotype.md | 94 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/phenotype.md | 107 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/phenotype.md | 107 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/protein.md | 52 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein.structure` |
| docs/agents/rules/protein.md | 106 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/protein.md | 119 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/protein.md | 119 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/protein.md | 136 | ImportError | Cannot import 'PDB' from module 'Bio' | `from Bio import PDB` |
| docs/agents/rules/quality.md | 31 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.io` |
| docs/agents/rules/quality.md | 58 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/quality.md | 71 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/quality.md | 71 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/rna.md | 45 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/agents/rules/rna.md | 65 | ImportError | Cannot import 'run_amalgkit_step' from module 'metainformant.rna.amalgkit' | `from metainformant.rna.amalgkit import run_amalgkit_step` |
| docs/agents/rules/rna.md | 65 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.amalgkit` |
| docs/agents/rules/rna.md | 86 | ImportError | Cannot import 'ProgressTracker' from module 'metainformant.rna.progress_tracker' | `from metainformant.rna.progress_tracker import ProgressTrack` |
| docs/agents/rules/rna.md | 86 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.progress_tracker` |
| docs/agents/rules/rna.md | 103 | ImportError | Cannot import 'run_metadata' from module 'metainformant.rna.steps' | `from metainformant.rna.steps import run_metadata` |
| docs/agents/rules/rna.md | 103 | ImportError | Cannot import 'run_quant' from module 'metainformant.rna.steps' | `from metainformant.rna.steps import run_quant` |
| docs/agents/rules/rna.md | 103 | ImportError | Cannot import 'run_download_quant_workflow' from module 'metainformant.rna.steps' | `from metainformant.rna.steps import run_download_quant_workf` |
| docs/agents/rules/rna.md | 103 | ImportError | Cannot import 'quantify_sample' from module 'metainformant.rna.steps' | `from metainformant.rna.steps import quantify_sample` |
| docs/agents/rules/rna.md | 103 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.steps` |
| docs/agents/rules/rna.md | 137 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.core.environment` |
| docs/agents/rules/rna.md | 162 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.monitoring` |
| docs/agents/rules/rna.md | 185 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.orchestration` |
| docs/agents/rules/rna.md | 207 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.core.cleanup` |
| docs/agents/rules/rna.md | 228 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.discovery` |
| docs/agents/rules/rna.md | 259 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/rna.md | 272 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/rna.md | 272 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/rna.md | 325 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/agents/rules/rna.md | 348 | ImportError | Cannot import 'run_workflow_for_species' from module 'metainformant.rna.orchestration' | `from metainformant.rna.orchestration import run_workflow_for` |
| docs/agents/rules/rna.md | 348 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/agents/rules/rna.md | 348 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.orchestration` |
| docs/agents/rules/simulation.md | 35 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.models.sequences` |
| docs/agents/rules/simulation.md | 51 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.models.agents` |
| docs/agents/rules/simulation.md | 68 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.models.rna` |
| docs/agents/rules/simulation.md | 89 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.models.popgen` |
| docs/agents/rules/simulation.md | 111 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.workflow.workflow` |
| docs/agents/rules/simulation.md | 181 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/simulation.md | 194 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/simulation.md | 194 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/singlecell.md | 40 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.analysis` |
| docs/agents/rules/singlecell.md | 40 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.features` |
| docs/agents/rules/singlecell.md | 59 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.analysis` |
| docs/agents/rules/singlecell.md | 74 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.data` |
| docs/agents/rules/singlecell.md | 95 | ImportError | Cannot import 'plot_qc_metrics' from module 'metainformant.singlecell' | `from metainformant.singlecell import plot_qc_metrics` |
| docs/agents/rules/singlecell.md | 95 | ImportError | Cannot import 'plot_dimensionality_reduction' from module 'metainformant.singlecell' | `from metainformant.singlecell import plot_dimensionality_red` |
| docs/agents/rules/singlecell.md | 95 | ImportError | Cannot import 'plot_gene_expression' from module 'metainformant.singlecell' | `from metainformant.singlecell import plot_gene_expression` |
| docs/agents/rules/singlecell.md | 95 | ImportError | Cannot import 'plot_clusters' from module 'metainformant.singlecell' | `from metainformant.singlecell import plot_clusters` |
| docs/agents/rules/singlecell.md | 95 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell` |
| docs/agents/rules/singlecell.md | 117 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/singlecell.md | 130 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/singlecell.md | 130 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/spatial.md | 29 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io` |
| docs/agents/rules/spatial.md | 29 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis` |
| docs/agents/rules/spatial.md | 29 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.integration` |
| docs/agents/rules/structural_variants.md | 28 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection` |
| docs/agents/rules/structural_variants.md | 28 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation` |
| docs/agents/rules/structural_variants.md | 28 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering` |
| docs/agents/rules/visualization.md | 47 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/agents/rules/visualization.md | 136 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/visualization.md | 148 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/agents/rules/visualization.md | 148 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/agents/rules/visualization.md | 166 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/cli.md | 3 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `import metainformant...` |
| docs/cli.md | 44 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/cloud/README.md | 60 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/cloud/README.md | 97 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/cloud/README.md | 141 | ImportError | Cannot import 'prepare_genome_index' from module 'metainformant.cloud.genome_prep' | `from metainformant.cloud.genome_prep import prepare_genome_i` |
| docs/cloud/README.md | 141 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud.genome_prep` |
| docs/cloud/README.md | 175 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/cloud/SPEC.md | 41 | SyntaxError | Invalid Python syntax: parameter without a default follows parameter with a default at line 26 | `class GCPDeployer:     def __init__(self, config: CloudConfi` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | ImportError | Cannot import 'calling' from module 'metainformant.dna' | `from metainformant.dna import calling` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/comparisons/dna_vs_rna_vs_transcriptome.md | 108 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 324 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 324 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 453 | ImportError | Cannot import 'MorphometricProfile' from module 'metainformant.phenotype' | `from metainformant.phenotype import MorphometricProfile` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 453 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 468 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 468 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | ImportError | Cannot import 'run_gwas_workflow' from module 'metainformant.gwas.workflow' | `from metainformant.gwas.workflow import run_gwas_workflow` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | ImportError | Cannot import 'run_eqtl' from module 'metainformant.gwas.finemapping.eqtl' | `from metainformant.gwas.finemapping.eqtl import run_eqtl` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.workflow` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping` |
| docs/comparisons/gwas_vs_phenotype_vs_multiomics.md | 478 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.eqtl` |
| docs/comparisons/methods_matrix.md | 425 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/comparisons/methods_matrix.md | 425 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/comparisons/methods_matrix.md | 538 | ImportError | Cannot import 'shannon_index' from module 'metainformant.ecology.analysis' | `from metainformant.ecology.analysis import shannon_index` |
| docs/comparisons/methods_matrix.md | 538 | ImportError | Cannot import 'bray_curtis' from module 'metainformant.ecology.analysis' | `from metainformant.ecology.analysis import bray_curtis` |
| docs/comparisons/methods_matrix.md | 538 | ImportError | Cannot import 'nmds' from module 'metainformant.ecology.analysis' | `from metainformant.ecology.analysis import nmds` |
| docs/comparisons/methods_matrix.md | 538 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology.analysis` |
| docs/comparisons/methods_matrix.md | 563 | ImportError | Cannot import 'simulate_sequences' from module 'metainformant.simulation.models.coalescent' | `from metainformant.simulation.models.coalescent import simul` |
| docs/comparisons/methods_matrix.md | 563 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.models.coalescent` |
| docs/comparisons/methods_matrix.md | 563 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/comparisons/methods_matrix.md | 688 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics` |
| docs/comparisons/visualization_approaches.md | 185 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import manhattan_p` |
| docs/comparisons/visualization_approaches.md | 185 | ImportError | Cannot import 'volcano_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import volcano_plo` |
| docs/comparisons/visualization_approaches.md | 185 | ImportError | Cannot import 'regional_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import regional_pl` |
| docs/comparisons/visualization_approaches.md | 185 | ImportError | Cannot import 'ideogram' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import ideogram` |
| docs/comparisons/visualization_approaches.md | 185 | ImportError | Cannot import 'circos_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import circos_plot` |
| docs/comparisons/visualization_approaches.md | 185 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/comparisons/visualization_approaches.md | 221 | ImportError | Cannot import 'expression_heatmap' from module 'metainformant.visualization.expression' | `from metainformant.visualization.expression import expressio` |
| docs/comparisons/visualization_approaches.md | 221 | ImportError | Cannot import 'ma_plot' from module 'metainformant.visualization.expression' | `from metainformant.visualization.expression import ma_plot` |
| docs/comparisons/visualization_approaches.md | 221 | ImportError | Cannot import 'enrichment_plot' from module 'metainformant.visualization.expression' | `from metainformant.visualization.expression import enrichmen` |
| docs/comparisons/visualization_approaches.md | 221 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.expression` |
| docs/comparisons/visualization_approaches.md | 253 | ImportError | Cannot import 'plot_phylo_tree' from module 'metainformant.visualization.trees' | `from metainformant.visualization.trees import plot_phylo_tre` |
| docs/comparisons/visualization_approaches.md | 253 | ImportError | Cannot import 'plot_tree_rectangular' from module 'metainformant.visualization.trees' | `from metainformant.visualization.trees import plot_tree_rect` |
| docs/comparisons/visualization_approaches.md | 253 | ImportError | Cannot import 'plot_tree_circular' from module 'metainformant.visualization.trees' | `from metainformant.visualization.trees import plot_tree_circ` |
| docs/comparisons/visualization_approaches.md | 253 | ImportError | Cannot import 'plot_tree_radial' from module 'metainformant.visualization.trees' | `from metainformant.visualization.trees import plot_tree_radi` |
| docs/comparisons/visualization_approaches.md | 253 | ImportError | Cannot import 'Phylo' from module 'Bio' | `from Bio import Phylo` |
| docs/comparisons/visualization_approaches.md | 253 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.trees` |
| docs/comparisons/visualization_approaches.md | 283 | ImportError | Cannot import 'plot_network' from module 'metainformant.visualization.networks' | `from metainformant.visualization.networks import plot_networ` |
| docs/comparisons/visualization_approaches.md | 283 | ImportError | Cannot import 'plot_network_circular' from module 'metainformant.visualization.networks' | `from metainformant.visualization.networks import plot_networ` |
| docs/comparisons/visualization_approaches.md | 283 | ImportError | Cannot import 'plot_network_hierarchical' from module 'metainformant.visualization.networks' | `from metainformant.visualization.networks import plot_networ` |
| docs/comparisons/visualization_approaches.md | 283 | ImportError | Cannot import 'plot_pathway' from module 'metainformant.visualization.networks' | `from metainformant.visualization.networks import plot_pathwa` |
| docs/comparisons/visualization_approaches.md | 283 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.networks` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import histogra` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'boxplot' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import boxplot` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'violin_plot' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import violin_p` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import qq_plot` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'roc_curve' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import roc_curv` |
| docs/comparisons/visualization_approaches.md | 316 | ImportError | Cannot import 'density_plot' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import density_` |
| docs/comparisons/visualization_approaches.md | 316 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.statistical` |
| docs/comparisons/visualization_approaches.md | 358 | ImportError | Cannot import 'pca_plot' from module 'metainformant.visualization.dimred' | `from metainformant.visualization.dimred import pca_plot` |
| docs/comparisons/visualization_approaches.md | 358 | ImportError | Cannot import 'umap_plot' from module 'metainformant.visualization.dimred' | `from metainformant.visualization.dimred import umap_plot` |
| docs/comparisons/visualization_approaches.md | 358 | ImportError | Cannot import 'tsne_plot' from module 'metainformant.visualization.dimred' | `from metainformant.visualization.dimred import tsne_plot` |
| docs/comparisons/visualization_approaches.md | 358 | ImportError | Cannot import 'scree_plot' from module 'metainformant.visualization.dimred' | `from metainformant.visualization.dimred import scree_plot` |
| docs/comparisons/visualization_approaches.md | 358 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dimred` |
| docs/comparisons/visualization_approaches.md | 397 | ImportError | Cannot import 'quality_scores' from module 'metainformant.visualization.quality' | `from metainformant.visualization.quality import quality_scor` |
| docs/comparisons/visualization_approaches.md | 397 | ImportError | Cannot import 'adapter_content' from module 'metainformant.visualization.quality' | `from metainformant.visualization.quality import adapter_cont` |
| docs/comparisons/visualization_approaches.md | 397 | ImportError | Cannot import 'length_distribution' from module 'metainformant.visualization.quality' | `from metainformant.visualization.quality import length_distr` |
| docs/comparisons/visualization_approaches.md | 397 | ImportError | Cannot import 'gc_content' from module 'metainformant.visualization.quality' | `from metainformant.visualization.quality import gc_content` |
| docs/comparisons/visualization_approaches.md | 397 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.quality` |
| docs/comparisons/visualization_approaches.md | 433 | ImportError | Cannot import 'entropy_plot' from module 'metainformant.visualization.information' | `from metainformant.visualization.information import entropy_` |
| docs/comparisons/visualization_approaches.md | 433 | ImportError | Cannot import 'mutual_information_heatmap' from module 'metainformant.visualization.information' | `from metainformant.visualization.information import mutual_i` |
| docs/comparisons/visualization_approaches.md | 433 | ImportError | Cannot import 'information_content_profile' from module 'metainformant.visualization.information' | `from metainformant.visualization.information import informat` |
| docs/comparisons/visualization_approaches.md | 433 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.information` |
| docs/comparisons/visualization_approaches.md | 468 | ImportError | Cannot import 'animate_time_series' from module 'metainformant.visualization.animations' | `from metainformant.visualization.animations import animate_t` |
| docs/comparisons/visualization_approaches.md | 468 | ImportError | Cannot import 'animate_evolution' from module 'metainformant.visualization.animations' | `from metainformant.visualization.animations import animate_e` |
| docs/comparisons/visualization_approaches.md | 468 | ImportError | Cannot import 'animate_clustering' from module 'metainformant.visualization.animations' | `from metainformant.visualization.animations import animate_c` |
| docs/comparisons/visualization_approaches.md | 468 | ImportError | Cannot import 'animate_network' from module 'metainformant.visualization.animations' | `from metainformant.visualization.animations import animate_n` |
| docs/comparisons/visualization_approaches.md | 468 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.animations` |
| docs/comparisons/visualization_approaches.md | 513 | ImportError | Cannot import 'set_publication_style' from module 'metainformant.visualization' | `from metainformant.visualization import set_publication_styl` |
| docs/comparisons/visualization_approaches.md | 513 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/comparisons/visualization_approaches.md | 513 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/comparisons/visualization_approaches.md | 530 | ImportError | Cannot import 'scatter_plot' from module 'metainformant.visualization' | `from metainformant.visualization import scatter_plot` |
| docs/comparisons/visualization_approaches.md | 530 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization' | `from metainformant.visualization import histogram` |
| docs/comparisons/visualization_approaches.md | 530 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/comparisons/visualization_approaches.md | 544 | ImportError | Cannot import 'enable_interactive' from module 'metainformant.visualization' | `from metainformant.visualization import enable_interactive` |
| docs/comparisons/visualization_approaches.md | 544 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/comparisons/visualization_approaches.md | 560 | ImportError | Cannot import 'create_dashboard' from module 'metainformant.visualization.dashboards' | `from metainformant.visualization.dashboards import create_da` |
| docs/comparisons/visualization_approaches.md | 560 | ImportError | Cannot import 'multi_panel_figure' from module 'metainformant.visualization.dashboards' | `from metainformant.visualization.dashboards import multi_pan` |
| docs/comparisons/visualization_approaches.md | 560 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dashboards` |
| docs/comparisons/visualization_approaches.md | 582 | ImportError | Cannot import 'set_presentation_style' from module 'metainformant.visualization' | `from metainformant.visualization import set_presentation_sty` |
| docs/comparisons/visualization_approaches.md | 582 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/comparisons/visualization_approaches.md | 621 | ImportError | Cannot import 'set_theme' from module 'metainformant.visualization.styling' | `from metainformant.visualization.styling import set_theme` |
| docs/comparisons/visualization_approaches.md | 621 | ImportError | Cannot import 'set_color_palette' from module 'metainformant.visualization.styling' | `from metainformant.visualization.styling import set_color_pa` |
| docs/comparisons/visualization_approaches.md | 621 | ImportError | Cannot import 'get_palette' from module 'metainformant.visualization.styling' | `from metainformant.visualization.styling import get_palette` |
| docs/comparisons/visualization_approaches.md | 621 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.styling` |
| docs/comparisons/visualization_approaches.md | 644 | ImportError | Cannot import 'create_figure_grid' from module 'metainformant.visualization.dashboards' | `from metainformant.visualization.dashboards import create_fi` |
| docs/comparisons/visualization_approaches.md | 644 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import manhattan_p` |
| docs/comparisons/visualization_approaches.md | 644 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import qq_plot` |
| docs/comparisons/visualization_approaches.md | 644 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization.statistical' | `from metainformant.visualization.statistical import histogra` |
| docs/comparisons/visualization_approaches.md | 644 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dashboards` |
| docs/comparisons/visualization_approaches.md | 644 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/comparisons/visualization_approaches.md | 644 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.statistical` |
| docs/comparisons/visualization_approaches.md | 670 | ImportError | Cannot import 'plot_manhattan' from module 'metainformant.gwas.visualization' | `from metainformant.gwas.visualization import plot_manhattan` |
| docs/comparisons/visualization_approaches.md | 670 | ImportError | Cannot import 'plot_qq' from module 'metainformant.gwas.visualization' | `from metainformant.gwas.visualization import plot_qq` |
| docs/comparisons/visualization_approaches.md | 670 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import manhattan_p` |
| docs/comparisons/visualization_approaches.md | 670 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization.genomics' | `from metainformant.visualization.genomics import qq_plot` |
| docs/comparisons/visualization_approaches.md | 670 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization` |
| docs/comparisons/visualization_approaches.md | 670 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics` |
| docs/comparisons/visualization_approaches.md | 688 | ImportError | Cannot import 'expression_heatmap' from module 'metainformant.visualization.expression' | `from metainformant.visualization.expression import expressio` |
| docs/comparisons/visualization_approaches.md | 688 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/comparisons/visualization_approaches.md | 707 | ImportError | Cannot import 'BasePlotter' from module 'metainformant.visualization.plots.base' | `from metainformant.visualization.plots.base import BasePlott` |
| docs/comparisons/visualization_approaches.md | 707 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots.base` |
| docs/core/GETTING_STARTED.md | 27 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/GETTING_STARTED.md | 27 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/GETTING_STARTED.md | 27 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/GETTING_STARTED.md | 27 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/core/GETTING_STARTED.md | 27 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/GETTING_STARTED.md | 27 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution` |
| docs/core/GETTING_STARTED.md | 27 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.data` |
| docs/core/GETTING_STARTED.md | 361 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/GETTING_STARTED.md | 361 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/README.md | 50 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.core import ...` |
| docs/core/cache.md | 69 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `CacheEntry(value: Any, ttl_seconds: int | None = None)` |
| docs/core/cache.md | 84 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `JsonCache(cache_dir: str | Path, ttl_seconds: int = 3600)` |
| docs/core/cache.md | 229 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/cache.md | 229 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/cache.md | 259 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/cache.md | 259 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/cache.md | 286 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/cache.md | 286 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/cache.md | 310 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/config.md | 140 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/config.md | 140 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/config.md | 277 | SyntaxError | Invalid Python syntax: invalid character '→' (U+2192) at line 1 | `_bool_map = {"1", "true", "yes", "y", "on"} → True _bool_map` |
| docs/core/config.md | 325 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 7 | `{     "top_level_keys": ["pipeline", "input", "reference"], ` |
| docs/core/config.md | 393 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/config.md | 393 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/config.md | 393 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/config.md | 610 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.data.validation` |
| docs/core/config.md | 623 | ImportError | Cannot import 'BaseModel' from module 'pydantic' | `from pydantic import BaseModel` |
| docs/core/config.md | 623 | ImportError | Cannot import 'Field' from module 'pydantic' | `from pydantic import Field` |
| docs/core/config.md | 635 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/config.md | 635 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/config.md | 635 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/config.md | 635 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/db.md | 96 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `PostgresConnection(     host: str = "localhost",     port: i` |
| docs/core/db.md | 444 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 444 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 481 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 481 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 506 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 506 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 536 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 536 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 578 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 578 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 598 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 598 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/db.md | 598 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 622 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 622 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 645 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 645 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 669 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 669 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/db.md | 669 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 701 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 701 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 701 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/db.md | 734 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 734 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/db.md | 734 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 833 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 833 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 1059 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 1059 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/db.md | 1101 | SyntaxError | Invalid Python syntax: invalid syntax at line 5 | `client = PostgresConnection(     host="db.example.com",     ` |
| docs/core/db.md | 1120 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/db.md | 1174 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/db.md | 1174 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 170 | SyntaxError | Invalid Python syntax: expected ':' at line 13 | `def download_with_progress(     url: str,     dest_path: str` |
| docs/core/download.md | 237 | SyntaxError | Invalid Python syntax: expected ':' at line 11 | `def monitor_subprocess_directory_growth(     *,     process:` |
| docs/core/download.md | 264 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 293 | SyntaxError | Invalid Python syntax: expected ':' at line 11 | `def monitor_subprocess_file_count(     *,     process: Any, ` |
| docs/core/download.md | 330 | SyntaxError | Invalid Python syntax: expected ':' at line 12 | `def monitor_subprocess_sample_progress(     *,     process: ` |
| docs/core/download.md | 374 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `get_file_size(url: str, *, timeout: int = 60) -> int | None ` |
| docs/core/download.md | 435 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 452 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 482 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/download.md | 482 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 482 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 521 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 570 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 601 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 625 | ModuleNotFoundError | Module 'boto3' not found in project or standard library | `import boto3` |
| docs/core/download.md | 625 | ModuleNotFoundError | Module 'boto3' not found in project or standard library | `import boto3` |
| docs/core/download.md | 625 | ImportError | Cannot import 'ClientError' from module 'botocore.exceptions' | `from botocore.exceptions import ClientError` |
| docs/core/download.md | 625 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 782 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 988 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/download.md | 988 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 1031 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/download.md | 1031 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 1031 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 1064 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/core/download.md | 1088 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/download.md | 1088 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/core/download.md | 1088 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 1109 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/download.md | 1109 | ImportError | Cannot import 'download' from module 'metainformant.core' | `from metainformant.core import download` |
| docs/core/download.md | 1109 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/download.md | 1127 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `# In a Snakemake workflow, use heartbeat for resume across r` |
| docs/core/download.md | 1154 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/download.md | 1256 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io.download` |
| docs/core/hash.md | 44 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils` |
| docs/core/hash.md | 209 | ImportError | Cannot import 'hash' from module 'metainformant.core' | `from metainformant.core import hash` |
| docs/core/hash.md | 209 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/hash.md | 235 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils` |
| docs/core/hash.md | 368 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution` |
| docs/core/hash.md | 368 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils` |
| docs/core/index.md | 83 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/index.md | 83 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/index.md | 83 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/index.md | 83 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils` |
| docs/core/io.md | 453 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/core/io.md | 453 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/io.md | 477 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/core/logging.md | 92 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/logging.md | 151 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/logging.md | 250 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/logging.md | 286 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/logging.md | 330 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/logging.md | 439 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/parallel.md | 417 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 417 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 435 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 435 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 459 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 459 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/parallel.md | 459 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 490 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 490 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 564 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 564 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 592 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 592 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 645 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution` |
| docs/core/parallel.md | 822 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 822 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 902 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.execution` |
| docs/core/parallel.md | 933 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/parallel.md | 933 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 1034 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 1034 | ImportError | Cannot import 'download' from module 'metainformant.core' | `from metainformant.core import download` |
| docs/core/parallel.md | 1034 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 1058 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 1058 | ImportError | Cannot import 'db' from module 'metainformant.core' | `from metainformant.core import db` |
| docs/core/parallel.md | 1058 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 1084 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 1084 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/parallel.md | 1108 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/parallel.md | 1108 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/parallel.md | 1108 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/paths.md | 144 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/paths.md | 144 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/paths.md | 172 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/paths.md | 172 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/paths.md | 501 | SyntaxError | Invalid Python syntax: invalid syntax. Perhaps you forgot a comma? at line 9 | `{     "total_dirs": 150,     "total_files": 12000,     "tota` |
| docs/core/paths.md | 558 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/paths.md | 558 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/paths.md | 748 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 23) at line 23 | `import tempfile  def process_with_temp_dir() -> Path:     ""` |
| docs/core/text.md | 556 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 15) at line 15 | `def make_fasta_filenames(accessions: list[str], extension: s` |
| docs/core/text.md | 675 | ImportError | Cannot import 'parallel' from module 'metainformant.core' | `from metainformant.core import parallel` |
| docs/core/text.md | 675 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/text.md | 702 | ImportError | Cannot import 'given' from module 'hypothesis' | `from hypothesis import given` |
| docs/core/text.md | 702 | ImportError | Cannot import 'strategies' from module 'hypothesis' | `from hypothesis import strategies` |
| docs/core/text.md | 751 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 34) at line 34 | `import random import string  def fuzz_text_functions():     ` |
| docs/core/workflow.md | 55 | ImportError | Cannot import 'BaseWorkflowOrchestrator' from module 'metainformant.core.workflow' | `from metainformant.core.workflow import BaseWorkflowOrchestr` |
| docs/core/workflow.md | 55 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.workflow` |
| docs/core/workflow.md | 174 | ImportError | Cannot import 'load_mapping_from_file' from module 'metainformant.core.config' | `from metainformant.core.config import load_mapping_from_file` |
| docs/core/workflow.md | 174 | ImportError | Cannot import 'apply_env_overrides' from module 'metainformant.core.config' | `from metainformant.core.config import apply_env_overrides` |
| docs/core/workflow.md | 174 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.config` |
| docs/core/workflow.md | 250 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/core/workflow.md | 402 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/core/workflow.md | 402 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/workflow.md | 419 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/core/workflow.md | 419 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/workflow.md | 433 | ImportError | Cannot import 'progress' from module 'metainformant.core' | `from metainformant.core import progress` |
| docs/core/workflow.md | 433 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/core/workflow.md | 452 | ImportError | Cannot import 'BaseWorkflowOrchestrator' from module 'metainformant.core.workflow' | `from metainformant.core.workflow import BaseWorkflowOrchestr` |
| docs/core/workflow.md | 452 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.workflow` |
| docs/core/workflow.md | 533 | ImportError | Cannot import 'BaseWorkflowOrchestrator' from module 'metainformant.core.workflow' | `from metainformant.core.workflow import BaseWorkflowOrchestr` |
| docs/core/workflow.md | 533 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.workflow` |
| docs/dna/README.md | 53 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.dna import ...` |
| docs/dna/SPEC.md | 50 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/dna/SPEC.md | 50 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/dna/SPEC.md | 50 | ImportError | Cannot import 'composition' from module 'metainformant.dna' | `from metainformant.dna import composition` |
| docs/dna/SPEC.md | 50 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/SPEC.md | 63 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/dna/SPEC.md | 63 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/accessions.md | 36 | ImportError | Cannot import 'is_valid_assembly_accession' from module 'metainformant.dna.genomes' | `from metainformant.dna.genomes import is_valid_assembly_acce` |
| docs/dna/accessions.md | 36 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.genomes` |
| docs/dna/accessions.md | 57 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/accessions.md | 100 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/alignment.md | 15 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/annotation.md | 30 | ImportError | Cannot import 'predict_orfs' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import predict_orfs` |
| docs/dna/annotation.md | 30 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 41 | ImportError | Cannot import 'annotate_coding_regions' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import annotate_coding_reg` |
| docs/dna/annotation.md | 41 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 56 | ImportError | Cannot import 'find_regulatory_elements' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import find_regulatory_ele` |
| docs/dna/annotation.md | 56 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 67 | ImportError | Cannot import 'annotate_cpg_islands' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import annotate_cpg_island` |
| docs/dna/annotation.md | 67 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 77 | ImportError | Cannot import 'compute_codon_usage' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import compute_codon_usage` |
| docs/dna/annotation.md | 77 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 88 | ImportError | Cannot import 'mask_repeats' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import mask_repeats` |
| docs/dna/annotation.md | 88 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 99 | ImportError | Cannot import 'find_splice_sites' from module 'metainformant.dna.annotation' | `from metainformant.dna.annotation import find_splice_sites` |
| docs/dna/annotation.md | 99 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation` |
| docs/dna/annotation.md | 113 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation.functional` |
| docs/dna/annotation.md | 123 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation.functional` |
| docs/dna/annotation.md | 134 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.annotation.functional` |
| docs/dna/calling.md | 23 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/calling.md | 38 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/calling.md | 49 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/calling.md | 59 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/calling.md | 76 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/calling.md | 87 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.variation.calling` |
| docs/dna/codon.md | 13 | ImportError | Cannot import 'codon' from module 'metainformant.dna' | `from metainformant.dna import codon` |
| docs/dna/codon.md | 13 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/composition.md | 14 | ImportError | Cannot import 'composition' from module 'metainformant.dna' | `from metainformant.dna import composition` |
| docs/dna/composition.md | 14 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/consensus.md | 33 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 46 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 66 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 94 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 122 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 143 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 167 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/consensus.md | 180 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence.consensus` |
| docs/dna/distances.md | 15 | ImportError | Cannot import 'distances' from module 'metainformant.dna' | `from metainformant.dna import distances` |
| docs/dna/distances.md | 15 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/fastq.md | 15 | ImportError | Cannot import 'fastq' from module 'metainformant.dna' | `from metainformant.dna import fastq` |
| docs/dna/fastq.md | 15 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/integration.md | 20 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.integration` |
| docs/dna/motifs.md | 13 | ImportError | Cannot import 'motifs' from module 'metainformant.dna' | `from metainformant.dna import motifs` |
| docs/dna/motifs.md | 13 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/msa.md | 15 | ImportError | Cannot import 'msa' from module 'metainformant.dna' | `from metainformant.dna import msa` |
| docs/dna/msa.md | 15 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/mutations.md | 20 | ImportError | Cannot import 'mutations' from module 'metainformant.dna' | `from metainformant.dna import mutations` |
| docs/dna/mutations.md | 20 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/ncbi.md | 28 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.entrez` |
| docs/dna/ncbi.md | 52 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.entrez` |
| docs/dna/ncbi.md | 123 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.ncbi` |
| docs/dna/ncbi.md | 140 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.ncbi` |
| docs/dna/ncbi.md | 175 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/ncbi.md | 191 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/ncbi.md | 210 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/ncbi.md | 221 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/ncbi.md | 231 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/ncbi.md | 242 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.external.genomes` |
| docs/dna/phylogeny.md | 17 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/dna/phylogeny.md | 17 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/population.md | 37 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/population.md | 165 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/population.md | 226 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/population.md | 226 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics.coalescent` |
| docs/dna/population.md | 243 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/population.md | 243 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/dna/population_visualization.md | 26 | ImportError | Cannot import 'plot_diversity_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_diversity_` |
| docs/dna/population_visualization.md | 26 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 60 | ImportError | Cannot import 'plot_tajimas_d_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_tajimas_d_` |
| docs/dna/population_visualization.md | 60 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 91 | ImportError | Cannot import 'plot_fst_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_fst_compar` |
| docs/dna/population_visualization.md | 91 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 122 | ImportError | Cannot import 'plot_pca_results' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_pca_result` |
| docs/dna/population_visualization.md | 122 | ImportError | Cannot import 'compute_pca' from module 'metainformant.gwas.structure' | `from metainformant.gwas.structure import compute_pca` |
| docs/dna/population_visualization.md | 122 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 122 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.structure` |
| docs/dna/population_visualization.md | 148 | ImportError | Cannot import 'plot_kinship_matrix' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_kinship_ma` |
| docs/dna/population_visualization.md | 148 | ImportError | Cannot import 'compute_kinship_matrix' from module 'metainformant.gwas.structure' | `from metainformant.gwas.structure import compute_kinship_mat` |
| docs/dna/population_visualization.md | 148 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 148 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.structure` |
| docs/dna/population_visualization.md | 173 | ImportError | Cannot import 'plot_site_frequency_spectrum' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_site_frequ` |
| docs/dna/population_visualization.md | 173 | ImportError | Cannot import 'generate_site_frequency_spectrum' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_site_fr` |
| docs/dna/population_visualization.md | 173 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 173 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/dna/population_visualization.md | 207 | ImportError | Cannot import 'plot_neutrality_test_summary' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_neutrality` |
| docs/dna/population_visualization.md | 207 | ImportError | Cannot import 'neutrality_test_suite' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import neutrality` |
| docs/dna/population_visualization.md | 207 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 207 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_visualization.md | 240 | ImportError | Cannot import 'plot_demographic_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_demographi` |
| docs/dna/population_visualization.md | 240 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 281 | ImportError | Cannot import 'plot_summary_statistics_grid' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_summary_st` |
| docs/dna/population_visualization.md | 281 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/dna/population_visualization.md | 281 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 281 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_visualization.md | 311 | ImportError | Cannot import 'plot_linkage_disequilibrium_decay' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_linkage_di` |
| docs/dna/population_visualization.md | 311 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'compare_populations' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import compare_po` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'neutrality_test_suite' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import neutrality` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'plot_diversity_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_diversity_` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'plot_neutrality_test_summary' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_neutrality` |
| docs/dna/population_visualization.md | 328 | ImportError | Cannot import 'plot_summary_statistics_grid' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_summary_st` |
| docs/dna/population_visualization.md | 328 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_visualization.md | 328 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_visualization.md | 358 | ImportError | Cannot import 'plot_diversity_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_diversity_` |
| docs/dna/population_visualization.md | 358 | ImportError | Cannot import 'plot_tajimas_d_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_tajimas_d_` |
| docs/dna/population_visualization.md | 358 | ImportError | Cannot import 'plot_fst_comparison' from module 'metainformant.dna.population_viz' | `from metainformant.dna.population_viz import plot_fst_compar` |
| docs/dna/population_visualization.md | 358 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_viz` |
| docs/dna/population_workflows.md | 11 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/dna/population_workflows.md | 11 | ImportError | Cannot import 'compare_populations' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import compare_po` |
| docs/dna/population_workflows.md | 11 | ImportError | Cannot import 'neutrality_test_suite' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import neutrality` |
| docs/dna/population_workflows.md | 11 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 25 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/dna/population_workflows.md | 25 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 64 | ImportError | Cannot import 'compare_populations' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import compare_po` |
| docs/dna/population_workflows.md | 64 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 95 | ImportError | Cannot import 'neutrality_test_suite' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import neutrality` |
| docs/dna/population_workflows.md | 95 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 117 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/dna/population_workflows.md | 117 | ImportError | Cannot import 'wattersons_theta' from module 'metainformant.math.population_genetics.coalescent' | `from metainformant.math.population_genetics.coalescent impor` |
| docs/dna/population_workflows.md | 117 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 117 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics.coalescent` |
| docs/dna/population_workflows.md | 140 | ImportError | Cannot import 'compare_populations' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import compare_po` |
| docs/dna/population_workflows.md | 140 | ImportError | Cannot import 'compute_pca' from module 'metainformant.gwas.structure' | `from metainformant.gwas.structure import compute_pca` |
| docs/dna/population_workflows.md | 140 | ImportError | Cannot import 'compute_kinship_matrix' from module 'metainformant.gwas.structure' | `from metainformant.gwas.structure import compute_kinship_mat` |
| docs/dna/population_workflows.md | 140 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/dna/population_workflows.md | 140 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.structure` |
| docs/dna/population_workflows.md | 176 | ImportError | Cannot import 'align_sequences' from module 'metainformant.dna.alignment' | `from metainformant.dna.alignment import align_sequences` |
| docs/dna/population_workflows.md | 176 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.alignment` |
| docs/dna/population_workflows.md | 241 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/dna/population_workflows.md | 241 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics.coalescent` |
| docs/dna/restriction.md | 13 | ImportError | Cannot import 'restriction' from module 'metainformant.dna' | `from metainformant.dna import restriction` |
| docs/dna/restriction.md | 13 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/sequences.md | 15 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/dna/sequences.md | 15 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/transcription.md | 13 | ImportError | Cannot import 'transcription' from module 'metainformant.dna' | `from metainformant.dna import transcription` |
| docs/dna/transcription.md | 13 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/translation.md | 16 | ImportError | Cannot import 'translation' from module 'metainformant.dna' | `from metainformant.dna import translation` |
| docs/dna/translation.md | 16 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/dna/variants.md | 12 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/dna/variants.md | 12 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/ecology/README.md | 41 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.ecology import ...` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'shannon_diversity' from module 'metainformant.ecology' | `from metainformant.ecology import shannon_diversity` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'simpson_diversity' from module 'metainformant.ecology' | `from metainformant.ecology import simpson_diversity` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'beta_diversity' from module 'metainformant.ecology' | `from metainformant.ecology import beta_diversity` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'community_metrics' from module 'metainformant.ecology' | `from metainformant.ecology import community_metrics` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'rarefaction_curve' from module 'metainformant.ecology' | `from metainformant.ecology import rarefaction_curve` |
| docs/ecology/community.md | 83 | ImportError | Cannot import 'alpha_beta_gamma_diversity' from module 'metainformant.ecology' | `from metainformant.ecology import alpha_beta_gamma_diversity` |
| docs/ecology/community.md | 83 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'functional_richness' from module 'metainformant.ecology' | `from metainformant.ecology import functional_richness` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'functional_evenness' from module 'metainformant.ecology' | `from metainformant.ecology import functional_evenness` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'functional_divergence' from module 'metainformant.ecology' | `from metainformant.ecology import functional_divergence` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'raos_quadratic_entropy' from module 'metainformant.ecology' | `from metainformant.ecology import raos_quadratic_entropy` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'community_weighted_mean' from module 'metainformant.ecology' | `from metainformant.ecology import community_weighted_mean` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'functional_diversity_suite' from module 'metainformant.ecology' | `from metainformant.ecology import functional_diversity_suite` |
| docs/ecology/functional.md | 63 | ImportError | Cannot import 'trait_distance_matrix' from module 'metainformant.ecology' | `from metainformant.ecology import trait_distance_matrix` |
| docs/ecology/functional.md | 63 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/index.md | 28 | ImportError | Cannot import 'community' from module 'metainformant.ecology' | `from metainformant.ecology import community` |
| docs/ecology/index.md | 28 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/index.md | 45 | ImportError | Cannot import 'community' from module 'metainformant.ecology' | `from metainformant.ecology import community` |
| docs/ecology/index.md | 45 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/index.md | 68 | ImportError | Cannot import 'community' from module 'metainformant.ecology' | `from metainformant.ecology import community` |
| docs/ecology/index.md | 68 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/index.md | 93 | ImportError | Cannot import 'community' from module 'metainformant.ecology' | `from metainformant.ecology import community` |
| docs/ecology/index.md | 93 | ImportError | Cannot import 'antwiki' from module 'metainformant.phenotype' | `from metainformant.phenotype import antwiki` |
| docs/ecology/index.md | 93 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/index.md | 93 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/ecology/index.md | 106 | ImportError | Cannot import 'community' from module 'metainformant.ecology' | `from metainformant.ecology import community` |
| docs/ecology/index.md | 106 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'indval' from module 'metainformant.ecology' | `from metainformant.ecology import indval` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'anosim' from module 'metainformant.ecology' | `from metainformant.ecology import anosim` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'permanova' from module 'metainformant.ecology' | `from metainformant.ecology import permanova` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'simper' from module 'metainformant.ecology' | `from metainformant.ecology import simper` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'cluster_communities' from module 'metainformant.ecology' | `from metainformant.ecology import cluster_communities` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'multivariate_dispersion' from module 'metainformant.ecology' | `from metainformant.ecology import multivariate_dispersion` |
| docs/ecology/indicators.md | 45 | ImportError | Cannot import 'distance_matrix' from module 'metainformant.ecology' | `from metainformant.ecology import distance_matrix` |
| docs/ecology/indicators.md | 45 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/ordination.md | 47 | ImportError | Cannot import 'distance_matrix' from module 'metainformant.ecology' | `from metainformant.ecology import distance_matrix` |
| docs/ecology/ordination.md | 47 | ImportError | Cannot import 'pcoa' from module 'metainformant.ecology' | `from metainformant.ecology import pcoa` |
| docs/ecology/ordination.md | 47 | ImportError | Cannot import 'nmds' from module 'metainformant.ecology' | `from metainformant.ecology import nmds` |
| docs/ecology/ordination.md | 47 | ImportError | Cannot import 'cca' from module 'metainformant.ecology' | `from metainformant.ecology import cca` |
| docs/ecology/ordination.md | 47 | ImportError | Cannot import 'procrustes' from module 'metainformant.ecology' | `from metainformant.ecology import procrustes` |
| docs/ecology/ordination.md | 47 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/ecology/visualization.md | 53 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology.visualization.visualization` |
| docs/epigenome/README.md | 36 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.epigenome import ...` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'ATACPeak' from module 'metainformant.epigenome' | `from metainformant.epigenome import ATACPeak` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'load_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import load_peaks` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'peak_statistics' from module 'metainformant.epigenome' | `from metainformant.epigenome import peak_statistics` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'calculate_nucleosome_fractions' from module 'metainformant.epigenome' | `from metainformant.epigenome import calculate_nucleosome_fra` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'tss_enrichment' from module 'metainformant.epigenome' | `from metainformant.epigenome import tss_enrichment` |
| docs/epigenome/atac_seq.md | 66 | ImportError | Cannot import 'compare_conditions' from module 'metainformant.epigenome' | `from metainformant.epigenome import compare_conditions` |
| docs/epigenome/atac_seq.md | 66 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'ChIPPeak' from module 'metainformant.epigenome' | `from metainformant.epigenome import ChIPPeak` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'load_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import load_peaks` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'filter_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import filter_peaks` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'peak_statistics' from module 'metainformant.epigenome' | `from metainformant.epigenome import peak_statistics` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'find_overlapping_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import find_overlapping_peaks` |
| docs/epigenome/chip_seq.md | 66 | ImportError | Cannot import 'calculate_enrichment' from module 'metainformant.epigenome' | `from metainformant.epigenome import calculate_enrichment` |
| docs/epigenome/chip_seq.md | 66 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/chromatin.md | 61 | ImportError | Cannot import 'learn_chromatin_states' from module 'metainformant.epigenome' | `from metainformant.epigenome import learn_chromatin_states` |
| docs/epigenome/chromatin.md | 61 | ImportError | Cannot import 'assign_states' from module 'metainformant.epigenome' | `from metainformant.epigenome import assign_states` |
| docs/epigenome/chromatin.md | 61 | ImportError | Cannot import 'interpret_states' from module 'metainformant.epigenome' | `from metainformant.epigenome import interpret_states` |
| docs/epigenome/chromatin.md | 61 | ImportError | Cannot import 'compute_state_enrichment' from module 'metainformant.epigenome' | `from metainformant.epigenome import compute_state_enrichment` |
| docs/epigenome/chromatin.md | 61 | ImportError | Cannot import 'segment_genome' from module 'metainformant.epigenome' | `from metainformant.epigenome import segment_genome` |
| docs/epigenome/chromatin.md | 61 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/index.md | 30 | ImportError | Cannot import 'methylation' from module 'metainformant.epigenome' | `from metainformant.epigenome import methylation` |
| docs/epigenome/index.md | 30 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/index.md | 54 | ImportError | Cannot import 'tracks' from module 'metainformant.epigenome' | `from metainformant.epigenome import tracks` |
| docs/epigenome/index.md | 54 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/methylation.md | 60 | ImportError | Cannot import 'MethylationSite' from module 'metainformant.epigenome' | `from metainformant.epigenome import MethylationSite` |
| docs/epigenome/methylation.md | 60 | ImportError | Cannot import 'compute_beta_values' from module 'metainformant.epigenome' | `from metainformant.epigenome import compute_beta_values` |
| docs/epigenome/methylation.md | 60 | ImportError | Cannot import 'find_dmrs' from module 'metainformant.epigenome' | `from metainformant.epigenome import find_dmrs` |
| docs/epigenome/methylation.md | 60 | ImportError | Cannot import 'find_cpg_islands' from module 'metainformant.epigenome' | `from metainformant.epigenome import find_cpg_islands` |
| docs/epigenome/methylation.md | 60 | ImportError | Cannot import 'calculate_methylation_entropy' from module 'metainformant.epigenome' | `from metainformant.epigenome import calculate_methylation_en` |
| docs/epigenome/methylation.md | 60 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/peak_calling.md | 60 | ImportError | Cannot import 'call_peaks_simple' from module 'metainformant.epigenome' | `from metainformant.epigenome import call_peaks_simple` |
| docs/epigenome/peak_calling.md | 60 | ImportError | Cannot import 'call_peaks_broad' from module 'metainformant.epigenome' | `from metainformant.epigenome import call_peaks_broad` |
| docs/epigenome/peak_calling.md | 60 | ImportError | Cannot import 'filter_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import filter_peaks` |
| docs/epigenome/peak_calling.md | 60 | ImportError | Cannot import 'compute_frip' from module 'metainformant.epigenome' | `from metainformant.epigenome import compute_frip` |
| docs/epigenome/peak_calling.md | 60 | ImportError | Cannot import 'differential_peaks' from module 'metainformant.epigenome' | `from metainformant.epigenome import differential_peaks` |
| docs/epigenome/peak_calling.md | 60 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/epigenome/workflow.md | 89 | ImportError | Cannot import 'EpigenomeConfig' from module 'metainformant.epigenome' | `from metainformant.epigenome import EpigenomeConfig` |
| docs/epigenome/workflow.md | 89 | ImportError | Cannot import 'run_methylation_workflow' from module 'metainformant.epigenome' | `from metainformant.epigenome import run_methylation_workflow` |
| docs/epigenome/workflow.md | 89 | ImportError | Cannot import 'run_chipseq_workflow' from module 'metainformant.epigenome' | `from metainformant.epigenome import run_chipseq_workflow` |
| docs/epigenome/workflow.md | 89 | ImportError | Cannot import 'run_atacseq_workflow' from module 'metainformant.epigenome' | `from metainformant.epigenome import run_atacseq_workflow` |
| docs/epigenome/workflow.md | 89 | ImportError | Cannot import 'integrate_epigenome_results' from module 'metainformant.epigenome' | `from metainformant.epigenome import integrate_epigenome_resu` |
| docs/epigenome/workflow.md | 89 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/eqtl/pipeline_guide.md | 76 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.eqtl` |
| docs/gwas/BENCHMARKING.md | 14 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis.benchmarking` |
| docs/gwas/BENCHMARKING.md | 135 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis.benchmarking` |
| docs/gwas/README.md | 35 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.colocalization` |
| docs/gwas/README.md | 84 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.analysis` |
| docs/gwas/README.md | 84 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.colocalization` |
| docs/gwas/amellifera_config.md | 133 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/amellifera_config.md | 133 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/amellifera_data_acquisition.md | 172 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/amellifera_data_acquisition.md | 172 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/index.md | 18 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/index.md | 18 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/index.md | 222 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/gwas/index.md | 222 | ImportError | Cannot import 'apply_qc_filters' from module 'metainformant.gwas' | `from metainformant.gwas import apply_qc_filters` |
| docs/gwas/index.md | 222 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/gwas/index.md | 222 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/index.md | 237 | ImportError | Cannot import 'hardy_weinberg_genotype_freqs' from module 'metainformant.math' | `from metainformant.math import hardy_weinberg_genotype_freqs` |
| docs/gwas/index.md | 237 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/gwas/index.md | 237 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/index.md | 250 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.gwas' | `from metainformant.gwas import manhattan_plot` |
| docs/gwas/index.md | 250 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/pbarbatus_config.md | 82 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/pbarbatus_config.md | 82 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/real_data_acquisition.md | 154 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/real_data_acquisition.md | 154 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/real_data_acquisition.md | 207 | ImportError | Cannot import 'download_variant_data' from module 'metainformant.gwas.download' | `from metainformant.gwas.download import download_variant_dat` |
| docs/gwas/real_data_acquisition.md | 207 | ImportError | Cannot import 'download_sra_run' from module 'metainformant.gwas.sra_download' | `from metainformant.gwas.sra_download import download_sra_run` |
| docs/gwas/real_data_acquisition.md | 207 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.download` |
| docs/gwas/real_data_acquisition.md | 207 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.sra_download` |
| docs/gwas/real_data_acquisition.md | 312 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/real_data_acquisition.md | 312 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 35 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 69 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 100 | ImportError | Cannot import 'estimate_population_structure' from module 'metainformant.gwas' | `from metainformant.gwas import estimate_population_structure` |
| docs/gwas/structure.md | 100 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 121 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 138 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 157 | ImportError | Cannot import 'estimate_population_structure' from module 'metainformant.gwas' | `from metainformant.gwas import estimate_population_structure` |
| docs/gwas/structure.md | 157 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/structure.md | 181 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/structure.md | 181 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/visualization_gallery.md | 7 | ImportError | Cannot import 'generate_all_plots' from module 'metainformant.gwas' | `from metainformant.gwas import generate_all_plots` |
| docs/gwas/visualization_gallery.md | 7 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/visualization_gallery.md | 44 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.gwas.visualization_genome' | `from metainformant.gwas.visualization_genome import manhatta` |
| docs/gwas/visualization_gallery.md | 44 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_genome` |
| docs/gwas/visualization_gallery.md | 96 | ImportError | Cannot import 'qq_plot' from module 'metainformant.gwas.visualization_statistical' | `from metainformant.gwas.visualization_statistical import qq_` |
| docs/gwas/visualization_gallery.md | 96 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_statistical` |
| docs/gwas/visualization_gallery.md | 163 | ImportError | Cannot import 'regional_plot' from module 'metainformant.gwas.visualization_regional' | `from metainformant.gwas.visualization_regional import region` |
| docs/gwas/visualization_gallery.md | 163 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_regional` |
| docs/gwas/visualization_gallery.md | 221 | ImportError | Cannot import 'pca_plot' from module 'metainformant.gwas.visualization_population' | `from metainformant.gwas.visualization_population import pca_` |
| docs/gwas/visualization_gallery.md | 221 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_population` |
| docs/gwas/visualization_gallery.md | 288 | ImportError | Cannot import 'maf_distribution' from module 'metainformant.gwas.visualization_variants' | `from metainformant.gwas.visualization_variants import maf_di` |
| docs/gwas/visualization_gallery.md | 288 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_variants` |
| docs/gwas/visualization_gallery.md | 356 | ImportError | Cannot import 'effect_size_forest_plot' from module 'metainformant.gwas.visualization_effects' | `from metainformant.gwas.visualization_effects import effect_` |
| docs/gwas/visualization_gallery.md | 356 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_effects` |
| docs/gwas/visualization_gallery.md | 407 | ImportError | Cannot import 'miami_plot' from module 'metainformant.gwas.visualization_comparison' | `from metainformant.gwas.visualization_comparison import miam` |
| docs/gwas/visualization_gallery.md | 407 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.visualization_comparison` |
| docs/gwas/visualization_gallery.md | 457 | ImportError | Cannot import 'generate_all_plots' from module 'metainformant.gwas' | `from metainformant.gwas import generate_all_plots` |
| docs/gwas/visualization_gallery.md | 457 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/visualization_gallery.md | 494 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/visualization_gallery.md | 494 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/workflow.md | 38 | ImportError | Cannot import 'execute_gwas_workflow' from module 'metainformant.gwas.workflow' | `from metainformant.gwas.workflow import execute_gwas_workflo` |
| docs/gwas/workflow.md | 38 | ImportError | Cannot import 'GWASWorkflowConfig' from module 'metainformant.gwas.workflow' | `from metainformant.gwas.workflow import GWASWorkflowConfig` |
| docs/gwas/workflow.md | 38 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.workflow` |
| docs/gwas/workflow.md | 63 | ImportError | Cannot import 'run_gwas_analysis' from module 'metainformant.gwas.workflow' | `from metainformant.gwas.workflow import run_gwas_analysis` |
| docs/gwas/workflow.md | 63 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.workflow` |
| docs/gwas/workflow.md | 110 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas.config' | `from metainformant.gwas.config import load_gwas_config` |
| docs/gwas/workflow.md | 110 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.config` |
| docs/gwas/workflow.md | 239 | ImportError | Cannot import 'load_gwas_config' from module 'metainformant.gwas' | `from metainformant.gwas import load_gwas_config` |
| docs/gwas/workflow.md | 239 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/workflow.md | 328 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/workflow.md | 339 | ImportError | Cannot import 'estimate_population_structure' from module 'metainformant.gwas' | `from metainformant.gwas import estimate_population_structure` |
| docs/gwas/workflow.md | 339 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/gwas/workflow.md | 351 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.gwas' | `from metainformant.gwas import manhattan_plot` |
| docs/gwas/workflow.md | 351 | ImportError | Cannot import 'qq_plot' from module 'metainformant.gwas' | `from metainformant.gwas import qq_plot` |
| docs/gwas/workflow.md | 351 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/information/README.md | 38 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.information import ...` |
| docs/information/channel.md | 77 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.advanced.channel` |
| docs/information/decomposition.md | 68 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.advanced.decomposition` |
| docs/information/entropy.md | 66 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.core.syntactic` |
| docs/information/estimation.md | 73 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.core.estimation` |
| docs/information/geometry.md | 88 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.advanced.geometry` |
| docs/information/index.md | 20 | ImportError | Cannot import 'shannon_entropy' from module 'metainformant.information' | `from metainformant.information import shannon_entropy` |
| docs/information/index.md | 20 | ImportError | Cannot import 'mutual_information' from module 'metainformant.information' | `from metainformant.information import mutual_information` |
| docs/information/index.md | 20 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information` |
| docs/information/mutual_information.md | 74 | AttributeError | 'metainformant' has no attribute 'information' (module not found) | `metainformant.information.metrics.core.syntactic` |
| docs/life_events/README.md | 41 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.life_events import ...` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'learn_event_embeddings' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import learn_event_emb` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'biological_embedding' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import biological_embe` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'domain_specific_embeddings' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import domain_specific` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'train_event_predictor' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import train_event_pre` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'predict_outcomes' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import predict_outcome` |
| docs/life_events/VALIDATION_REPORT.md | 71 | ImportError | Cannot import 'save_model' from module 'metainformant.life_events.models' | `from metainformant.life_events.models import save_model` |
| docs/life_events/VALIDATION_REPORT.md | 71 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models` |
| docs/life_events/VALIDATION_REPORT.md | 84 | SyntaxError | Invalid Python syntax: positional argument follows keyword argument at line 2 | `embedding_results = embeddings.learn_event_embeddings(sequen` |
| docs/life_events/VALIDATION_REPORT.md | 95 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 4 | `def train_event_predictor(sequences, outcomes, embedding_res` |
| docs/life_events/VALIDATION_REPORT.md | 103 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `predictor = EventSequencePredictor(**kwargs)    predictor.fi` |
| docs/life_events/events.md | 21 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.events` |
| docs/life_events/events.md | 43 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.events` |
| docs/life_events/events.md | 66 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.events` |
| docs/life_events/events.md | 91 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.events` |
| docs/life_events/index.md | 27 | ImportError | Cannot import 'Event' from module 'metainformant.life_events' | `from metainformant.life_events import Event` |
| docs/life_events/index.md | 27 | ImportError | Cannot import 'EventSequence' from module 'metainformant.life_events' | `from metainformant.life_events import EventSequence` |
| docs/life_events/index.md | 27 | ImportError | Cannot import 'analyze_life_course' from module 'metainformant.life_events' | `from metainformant.life_events import analyze_life_course` |
| docs/life_events/index.md | 27 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events` |
| docs/life_events/models.md | 24 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.predictor` |
| docs/life_events/models.md | 44 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.sequence_models` |
| docs/life_events/models.md | 59 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.statistical_models` |
| docs/life_events/models.md | 74 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.predictor` |
| docs/life_events/models.md | 87 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.statistical_models` |
| docs/life_events/models.md | 106 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.embeddings` |
| docs/life_events/models.md | 117 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.utils` |
| docs/life_events/models.md | 117 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.predictor` |
| docs/life_events/survival.md | 25 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/survival.md | 52 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/survival.md | 77 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/survival.md | 99 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/survival.md | 121 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/survival.md | 145 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.utils` |
| docs/life_events/survival.md | 145 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.survival.time_to_event` |
| docs/life_events/visualization.md | 29 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.visualization.timeline` |
| docs/life_events/visualization.md | 94 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.visualization.timeline` |
| docs/life_events/visualization.md | 94 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.visualization.statistical` |
| docs/life_events/visualization.md | 94 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.visualization.network` |
| docs/life_events/visualization.md | 115 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.models.embeddings` |
| docs/life_events/visualization.md | 115 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.visualization.statistical` |
| docs/life_events/workflow.md | 24 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.workflow.workflow` |
| docs/life_events/workflow.md | 44 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.workflow.workflow` |
| docs/life_events/workflow.md | 61 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.workflow.workflow` |
| docs/life_events/workflow.md | 78 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.workflow.workflow` |
| docs/life_events/workflow.md | 96 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.events` |
| docs/life_events/workflow.md | 96 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.config` |
| docs/life_events/workflow.md | 96 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.workflow.workflow` |
| docs/life_events/workflow.md | 139 | AttributeError | 'metainformant' has no attribute 'life_events' (module not found) | `metainformant.life_events.core.config` |
| docs/longread/README.md | 61 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.workflow.orchestrator` |
| docs/longread/README.md | 61 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.io` |
| docs/longread/README.md | 61 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.quality` |
| docs/longread/README.md | 61 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.analysis` |
| docs/longread/README.md | 61 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.assembly` |
| docs/longread/analysis.md | 23 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def detect_sv_from_long_reads(     alignments: list[Any],   ` |
| docs/longread/analysis.md | 33 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_insertions(alignments: list[Any], min_size: int =` |
| docs/longread/analysis.md | 43 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_methylation(alignments: list[Any]) -> dict[str, A` |
| docs/longread/analysis.md | 60 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def phase_reads(     alignments: list[Any],     heterozygous` |
| docs/longread/analysis.md | 80 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/longread/assembly.md | 76 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def minimizer_sketch(     sequence: str, k: int = 15, w: int` |
| docs/longread/assembly.md | 84 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def find_overlaps(     reads: dict[str, str],     k: int = 1` |
| docs/longread/assembly.md | 93 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def compute_overlap_graph(overlaps: list[Overlap]) -> dict[s` |
| docs/longread/assembly.md | 102 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def generate_consensus(reads: list[str], backbone_idx: int |` |
| docs/longread/assembly.md | 113 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def hybrid_assemble(     long_reads: dict[str, str], short_r` |
| docs/longread/assembly.md | 129 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 7 | `from metainformant.longread import assembly  # Compute minim` |
| docs/longread/index.md | 36 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.io` |
| docs/longread/index.md | 36 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.analysis` |
| docs/longread/index.md | 36 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.assembly` |
| docs/longread/index.md | 36 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread.workflow` |
| docs/longread/io.md | 74 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def read_fast5(filepath: str | Path) -> list[Fast5Read] def ` |
| docs/longread/io.md | 85 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def read_long_read_bam(filepath: str | Path, region: str | N` |
| docs/longread/io.md | 96 | SyntaxError | Invalid Python syntax: expected ':' at line 2 | `def fast5_to_fastq(fast5_path: str | Path, output_path: str ` |
| docs/longread/io.md | 107 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/longread/methylation.md | 34 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def call_methylation_from_signal(     signal_data: Any,     ` |
| docs/longread/methylation.md | 47 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def compute_methylation_stats(     methylation_calls: list[d` |
| docs/longread/methylation.md | 58 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def aggregate_methylation(     per_read_calls: list[dict],  ` |
| docs/longread/methylation.md | 69 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def detect_dmrs(     group1_methylation: dict[str, Any],    ` |
| docs/longread/methylation.md | 83 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def methylation_pattern_analysis(     per_read_calls: list[d` |
| docs/longread/methylation.md | 95 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/longread/phasing.md | 31 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def phase_reads(     alignments: list[Any],     heterozygous` |
| docs/longread/phasing.md | 43 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def build_phase_blocks(     phased_data: dict[str, Any],    ` |
| docs/longread/phasing.md | 54 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def haplotag_reads(     alignments: list[Any],     phase_blo` |
| docs/longread/phasing.md | 65 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def compute_switch_errors(     inferred_phases: dict[str, An` |
| docs/longread/phasing.md | 76 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def allele_specific_analysis(     haplotagged_reads: list[di` |
| docs/longread/phasing.md | 87 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/longread/quality.md | 64 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def calculate_n50(read_lengths: list[int]) -> int def calcul` |
| docs/longread/quality.md | 79 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def filter_by_length(     reads: list[ReadRecord],     min_l` |
| docs/longread/quality.md | 95 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def detect_adapters(     reads: list[ReadRecord],     adapte` |
| docs/longread/quality.md | 117 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/math/README.md | 46 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.math import ...` |
| docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md | 373 | SyntaxError | Invalid Python syntax: invalid syntax at line 12 | `# Add to __all__: __all__ = [     # ... existing exports ...` |
| docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md | 449 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# tests/test_math_utilities.py    def test_correlation_coeff` |
| docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md | 474 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# tests/test_math_edge_cases.py    def test_large_sample_siz` |
| docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md | 487 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# tests/test_math_integration.py    def test_price_equation_` |
| docs/math/REVIEW_AND_IMPROVEMENT_PLAN.md | 509 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/coalescent.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics` |
| docs/math/coalescent.md | 43 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics.coalescent` |
| docs/math/ddm.md | 13 | ImportError | Cannot import 'ddm' from module 'metainformant.math' | `from metainformant.math import ddm` |
| docs/math/ddm.md | 13 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/demography.md | 27 | ImportError | Cannot import 'exponential_growth_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import exponential_growth` |
| docs/math/demography.md | 27 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 52 | ImportError | Cannot import 'bottleneck_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import bottleneck_effecti` |
| docs/math/demography.md | 52 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 85 | ImportError | Cannot import 'two_epoch_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import two_epoch_effectiv` |
| docs/math/demography.md | 85 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 103 | ImportError | Cannot import 'exponential_growth_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import exponential_growth` |
| docs/math/demography.md | 103 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 103 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/math/demography.md | 126 | ImportError | Cannot import 'bottleneck_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import bottleneck_effecti` |
| docs/math/demography.md | 126 | ImportError | Cannot import 'exponential_growth_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import exponential_growth` |
| docs/math/demography.md | 126 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 155 | ImportError | Cannot import 'exponential_growth_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import exponential_growth` |
| docs/math/demography.md | 155 | ImportError | Cannot import 'two_epoch_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import two_epoch_effectiv` |
| docs/math/demography.md | 155 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 209 | ImportError | Cannot import 'exponential_growth_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import exponential_growth` |
| docs/math/demography.md | 209 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 209 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics.coalescent` |
| docs/math/demography.md | 225 | ImportError | Cannot import 'bottleneck_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import bottleneck_effecti` |
| docs/math/demography.md | 225 | ImportError | Cannot import 'effective_population_size_from_heterozygosity' from module 'metainformant.math.popgen' | `from metainformant.math.popgen import effective_population_s` |
| docs/math/demography.md | 225 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/math/demography.md | 225 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.popgen` |
| docs/math/dynamics.md | 39 | ImportError | Cannot import 'dynamics' from module 'metainformant.math' | `from metainformant.math import dynamics` |
| docs/math/dynamics.md | 39 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 34 | ImportError | Cannot import 'harmonic_mean_effective_size' from module 'metainformant.math' | `from metainformant.math import harmonic_mean_effective_size` |
| docs/math/effective_size.md | 34 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 56 | ImportError | Cannot import 'effective_size_sex_ratio' from module 'metainformant.math' | `from metainformant.math import effective_size_sex_ratio` |
| docs/math/effective_size.md | 56 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 81 | ImportError | Cannot import 'effective_size_from_family_size_variance' from module 'metainformant.math' | `from metainformant.math import effective_size_from_family_si` |
| docs/math/effective_size.md | 81 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 93 | ImportError | Cannot import 'harmonic_mean_effective_size' from module 'metainformant.math' | `from metainformant.math import harmonic_mean_effective_size` |
| docs/math/effective_size.md | 93 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 104 | ImportError | Cannot import 'effective_size_sex_ratio' from module 'metainformant.math' | `from metainformant.math import effective_size_sex_ratio` |
| docs/math/effective_size.md | 104 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 117 | ImportError | Cannot import 'effective_size_from_family_size_variance' from module 'metainformant.math' | `from metainformant.math import effective_size_from_family_si` |
| docs/math/effective_size.md | 117 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/effective_size.md | 133 | ImportError | Cannot import 'effective_size_sex_ratio' from module 'metainformant.math' | `from metainformant.math import effective_size_sex_ratio` |
| docs/math/effective_size.md | 133 | ImportError | Cannot import 'effective_population_size_from_heterozygosity' from module 'metainformant.math' | `from metainformant.math import effective_population_size_fro` |
| docs/math/effective_size.md | 133 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'sir_step' from module 'metainformant.math' | `from metainformant.math import sir_step` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'seir_step' from module 'metainformant.math' | `from metainformant.math import seir_step` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'sis_step' from module 'metainformant.math' | `from metainformant.math import sis_step` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'basic_reproduction_number' from module 'metainformant.math' | `from metainformant.math import basic_reproduction_number` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'effective_reproduction_number' from module 'metainformant.math' | `from metainformant.math import effective_reproduction_number` |
| docs/math/epidemiology.md | 7 | ImportError | Cannot import 'herd_immunity_threshold' from module 'metainformant.math' | `from metainformant.math import herd_immunity_threshold` |
| docs/math/epidemiology.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/epidemiology.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/fst.md | 28 | ImportError | Cannot import 'fst_from_heterozygosity' from module 'metainformant.math' | `from metainformant.math import fst_from_heterozygosity` |
| docs/math/fst.md | 28 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/fst.md | 56 | ImportError | Cannot import 'fst_from_allele_freqs' from module 'metainformant.math' | `from metainformant.math import fst_from_allele_freqs` |
| docs/math/fst.md | 56 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/fst.md | 77 | ImportError | Cannot import 'fst_from_heterozygosity' from module 'metainformant.math' | `from metainformant.math import fst_from_heterozygosity` |
| docs/math/fst.md | 77 | ImportError | Cannot import 'fst_from_allele_freqs' from module 'metainformant.math' | `from metainformant.math import fst_from_allele_freqs` |
| docs/math/fst.md | 77 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/fst.md | 90 | ImportError | Cannot import 'fst_from_heterozygosity' from module 'metainformant.math' | `from metainformant.math import fst_from_heterozygosity` |
| docs/math/fst.md | 90 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/math/fst.md | 90 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/ld.md | 7 | ImportError | Cannot import 'ld_coefficients' from module 'metainformant.math' | `from metainformant.math import ld_coefficients` |
| docs/math/ld.md | 7 | ImportError | Cannot import 'r_squared' from module 'metainformant.math' | `from metainformant.math import r_squared` |
| docs/math/ld.md | 7 | ImportError | Cannot import 'haldane_d_to_c' from module 'metainformant.math' | `from metainformant.math import haldane_d_to_c` |
| docs/math/ld.md | 7 | ImportError | Cannot import 'kosambi_c_to_d' from module 'metainformant.math' | `from metainformant.math import kosambi_c_to_d` |
| docs/math/ld.md | 7 | ImportError | Cannot import 'expected_r2_from_Ne_c' from module 'metainformant.math' | `from metainformant.math import expected_r2_from_Ne_c` |
| docs/math/ld.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/ld.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/popgen.md | 7 | ImportError | Cannot import 'popgen' from module 'metainformant.math' | `from metainformant.math import popgen` |
| docs/math/popgen.md | 7 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/popgen_stats.md | 29 | ImportError | Cannot import 'bootstrap_confidence_interval' from module 'metainformant.math.popgen_stats' | `from metainformant.math.popgen_stats import bootstrap_confid` |
| docs/math/popgen_stats.md | 29 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.popgen_stats` |
| docs/math/popgen_stats.md | 57 | ImportError | Cannot import 'permutation_test' from module 'metainformant.math.popgen_stats' | `from metainformant.math.popgen_stats import permutation_test` |
| docs/math/popgen_stats.md | 57 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.popgen_stats` |
| docs/math/popgen_stats.md | 82 | ImportError | Cannot import 'detect_outliers' from module 'metainformant.math.popgen_stats' | `from metainformant.math.popgen_stats import detect_outliers` |
| docs/math/popgen_stats.md | 82 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.popgen_stats` |
| docs/math/price.md | 14 | ImportError | Cannot import 'price' from module 'metainformant.math' | `from metainformant.math import price` |
| docs/math/price.md | 14 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math` |
| docs/math/selection.md | 5 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.population_genetics` |
| docs/mcp/SPEC.md | 25 | ImportError | Cannot import 'register_tool' from module 'metainformant.mcp' | `from metainformant.mcp import register_tool` |
| docs/mcp/SPEC.md | 25 | AttributeError | 'metainformant' has no attribute 'mcp' (module not found) | `metainformant.mcp` |
| docs/mcp/index.md | 246 | ImportError | Cannot import 'register_tool' from module 'metainformant.mcp' | `from metainformant.mcp import register_tool` |
| docs/mcp/index.md | 246 | ImportError | Cannot import 'run_association' from module 'metainformant.gwas' | `from metainformant.gwas import run_association` |
| docs/mcp/index.md | 246 | AttributeError | 'metainformant' has no attribute 'mcp' (module not found) | `metainformant.mcp` |
| docs/mcp/index.md | 246 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/mcp/index.md | 277 | ImportError | Cannot import 'amalgkit_monitor' from module 'tools' | `from tools import amalgkit_monitor` |
| docs/mcp/index.md | 277 | ImportError | Cannot import 'my_tool' from module 'tools' | `from tools import my_tool` |
| docs/menu/README.md | 31 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core.discovery` |
| docs/menu/README.md | 31 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.navigation` |
| docs/menu/README.md | 31 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.display` |
| docs/menu/discovery.md | 73 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core.discovery` |
| docs/menu/discovery.md | 119 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core.discovery` |
| docs/menu/discovery.md | 119 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core.discovery` |
| docs/menu/index.md | 17 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.core` |
| docs/menu/index.md | 17 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui` |
| docs/menu/navigation.md | 105 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.navigation` |
| docs/menu/navigation.md | 105 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.display` |
| docs/menu/navigation.md | 139 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.navigation` |
| docs/menu/navigation.md | 139 | AttributeError | 'metainformant' has no attribute 'menu' (module not found) | `metainformant.menu.ui.display` |
| docs/metabolomics/ARCHITECTURE.md | 160 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/metabolomics/ARCHITECTURE.md | 160 | ImportError | Cannot import 'correlate_omics_layers' from module 'metainformant.multiomics' | `from metainformant.multiomics import correlate_omics_layers` |
| docs/metabolomics/ARCHITECTURE.md | 160 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/ARCHITECTURE.md | 160 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/metabolomics/ARCHITECTURE.md | 160 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/metabolomics/ARCHITECTURE.md | 160 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/metabolomics/CONFIGURATION.md | 29 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/metabolomics/CONFIGURATION.md | 29 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/metabolomics/CONFIGURATION.md | 144 | ImportError | Cannot import 'load_database' from module 'metainformant.metabolomics.io' | `from metainformant.metabolomics.io import load_database` |
| docs/metabolomics/CONFIGURATION.md | 144 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.io` |
| docs/metabolomics/CONFIGURATION.md | 172 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/metabolomics/CONFIGURATION.md | 217 | ImportError | Cannot import 'Combat' from module 'pycombat' | `from pycombat import Combat` |
| docs/metabolomics/CONFIGURATION.md | 308 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `identify_metabolites(..., ppm_tolerance=0) ValueError: ppm_t` |
| docs/metabolomics/CONFIGURATION.md | 371 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `import metainformant, sys    print(f"METAINFORMANT {metainfo` |
| docs/metabolomics/EXAMPLES.md | 24 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 24 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/metabolomics/EXAMPLES.md | 185 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 261 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 261 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/metabolomics/EXAMPLES.md | 331 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 428 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 428 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/metabolomics/EXAMPLES.md | 428 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.io` |
| docs/metabolomics/EXAMPLES.md | 531 | ImportError | Cannot import 'Combat' from module 'pycombat' | `from pycombat import Combat` |
| docs/metabolomics/EXAMPLES.md | 531 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/EXAMPLES.md | 531 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/metabolomics/GETTING_STARTED.md | 52 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/GETTING_STARTED.md | 260 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/metabolomics/INTEGRATION.md | 21 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/metabolomics/INTEGRATION.md | 21 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/metabolomics/INTEGRATION.md | 44 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 44 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/metabolomics/INTEGRATION.md | 44 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.visualization` |
| docs/metabolomics/INTEGRATION.md | 69 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/metabolomics/INTEGRATION.md | 69 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 69 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/metabolomics/INTEGRATION.md | 94 | SyntaxError | Invalid Python syntax: unmatched ')' at line 5 | `from metainformant.protein import quant as prot from metainf` |
| docs/metabolomics/INTEGRATION.md | 114 | ImportError | Cannot import 'qc_metrics' from module 'metainformant.quality' | `from metainformant.quality import qc_metrics` |
| docs/metabolomics/INTEGRATION.md | 114 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality` |
| docs/metabolomics/INTEGRATION.md | 114 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 135 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 135 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/metabolomics/INTEGRATION.md | 177 | ImportError | Cannot import 'factor_analysis' from module 'metainformant.multiomics' | `from metainformant.multiomics import factor_analysis` |
| docs/metabolomics/INTEGRATION.md | 177 | ImportError | Cannot import 'quant' from module 'metainformant.protein' | `from metainformant.protein import quant` |
| docs/metabolomics/INTEGRATION.md | 177 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/metabolomics/INTEGRATION.md | 177 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 177 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/metabolomics/INTEGRATION.md | 177 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/metabolomics/INTEGRATION.md | 224 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 308 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.io` |
| docs/metabolomics/INTEGRATION.md | 362 | ImportError | Cannot import 'integrate' from module 'metainformant.multiomics' | `from metainformant.multiomics import integrate` |
| docs/metabolomics/INTEGRATION.md | 362 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/INTEGRATION.md | 362 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/metabolomics/INTEGRATION.md | 362 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/metabolomics/INTEGRATION.md | 411 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `rule metabolomics_analysis:     input:         csv="data/met` |
| docs/metabolomics/INTEGRATION.md | 437 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/metabolomics/INTEGRATION.md | 437 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/PERFORMANCE.md | 246 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics` |
| docs/metabolomics/PERFORMANCE.md | 327 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `batch_size = 1000    for i in range(0, n_samples, batch_size` |
| docs/metabolomics/PERFORMANCE.md | 406 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `chunk_rows = 1000    for start in range(0, n_metabolites, ch` |
| docs/metabolomics/PERFORMANCE.md | 414 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `import dask.array as da    X = da.from_array(intensities, ch` |
| docs/metabolomics/README.md | 38 | ImportError | Cannot import 'quantification' from module 'metainformant.metabolomics.analysis' | `from metainformant.metabolomics.analysis import quantificati` |
| docs/metabolomics/README.md | 38 | ImportError | Cannot import 'normalization' from module 'metainformant.metabolomics.analysis' | `from metainformant.metabolomics.analysis import normalizatio` |
| docs/metabolomics/README.md | 38 | ImportError | Cannot import 'mapping' from module 'metainformant.metabolomics.pathways' | `from metainformant.metabolomics.pathways import mapping` |
| docs/metabolomics/README.md | 38 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.io` |
| docs/metabolomics/README.md | 38 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.analysis` |
| docs/metabolomics/README.md | 38 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways` |
| docs/metabolomics/README.md | 38 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.visualization` |
| docs/metabolomics/TROUBLESHOOTING.md | 10 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `>>> from metainformant.metabolomics import analysis ModuleNo` |
| docs/metabolomics/TROUBLESHOOTING.md | 141 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Check actual ppm errors    for mz in observed_mz[:5]:     ` |
| docs/metabolomics/TROUBLESHOOTING.md | 152 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# If positive mode, observed m/z = neutral + 1.007 ([M+H]+) ` |
| docs/metabolomics/TROUBLESHOOTING.md | 222 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `nonzero = dataset.intensities.sum(axis=0) > 0    intensities` |
| docs/metabolomics/TROUBLESHOOTING.md | 258 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Check variances    var_a = X[:, group_a].var(axis=1)    va` |
| docs/metabolomics/TROUBLESHOOTING.md | 334 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `query_set = set(query_metabolites)    for name, members in p` |
| docs/metabolomics/TROUBLESHOOTING.md | 404 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Only keep DB entries near observed range (with buffer)    ` |
| docs/metabolomics/TROUBLESHOOTING.md | 413 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from collections import defaultdict    bins = defaultdict(li` |
| docs/metabolomics/TROUBLESHOOTING.md | 437 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `query_set = set(query_metabolites)    pathway_db_filtered = ` |
| docs/metabolomics/TROUBLESHOOTING.md | 449 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from concurrent.futures import ThreadPoolExecutor    with Th` |
| docs/metabolomics/TROUBLESHOOTING.md | 506 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/metabolomics/TROUBLESHOOTING.md | 554 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `import logging    logging.basicConfig(level=logging.DEBUG)` |
| docs/metagenomics/ARCHITECTURE.md | 400 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.diversity` |
| docs/metagenomics/ARCHITECTURE.md | 436 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.functional` |
| docs/metagenomics/ARCHITECTURE.md | 455 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.comparative` |
| docs/metagenomics/ARCHITECTURE.md | 580 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/metagenomics/ARCHITECTURE.md | 580 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.visualization` |
| docs/metagenomics/ARCHITECTURE.md | 597 | ImportError | Cannot import 'metrics' from module 'metainformant.ecology' | `from metainformant.ecology import metrics` |
| docs/metagenomics/ARCHITECTURE.md | 597 | AttributeError | 'metainformant' has no attribute 'ecology' (module not found) | `metainformant.ecology` |
| docs/metagenomics/ARCHITECTURE.md | 616 | ImportError | Cannot import 'construction' from module 'metainformant.networks' | `from metainformant.networks import construction` |
| docs/metagenomics/ARCHITECTURE.md | 616 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/metagenomics/CAPABILITIES.md | 214 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 1 | `reads = {"r1": "ACGTACGT", "r2": "ACGTACGA", "r3": "ACGTACGA` |
| docs/metagenomics/CAPABILITIES.md | 470 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 1 | `contigs = {"c1": "ACGT...", "c2": "GGCC...", ...} coverage =` |
| docs/metagenomics/CAPABILITIES.md | 543 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 2 | `refs = {"genusA_sp1": "ACGT...", "genusA_sp2": "ACGT..."} ta` |
| docs/metagenomics/GETTING_STARTED.md | 38 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/GETTING_STARTED.md | 38 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.amplicon` |
| docs/metagenomics/GETTING_STARTED.md | 176 | ImportError | Cannot import 'assembly' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import assembly` |
| docs/metagenomics/GETTING_STARTED.md | 176 | ImportError | Cannot import 'binning' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import binning` |
| docs/metagenomics/GETTING_STARTED.md | 176 | ImportError | Cannot import 'profiling' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import profiling` |
| docs/metagenomics/GETTING_STARTED.md | 176 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/GETTING_STARTED.md | 300 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.amplicon` |
| docs/metagenomics/GETTING_STARTED.md | 379 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.diversity` |
| docs/metagenomics/GETTING_STARTED.md | 405 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.comparative` |
| docs/metagenomics/GETTING_STARTED.md | 420 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.visualization` |
| docs/metagenomics/GETTING_STARTED.md | 448 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.shotgun` |
| docs/metagenomics/GETTING_STARTED.md | 448 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.shotgun.binning` |
| docs/metagenomics/GETTING_STARTED.md | 587 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.diversity` |
| docs/metagenomics/README.md | 47 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.amplicon` |
| docs/metagenomics/README.md | 47 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.shotgun` |
| docs/metagenomics/README.md | 47 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.diversity` |
| docs/metagenomics/README.md | 47 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.functional` |
| docs/metagenomics/README.md | 47 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.comparative` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'cluster_otus' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import cluster_otus` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'filter_chimeras' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import filter_chimeras` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'denoise_sequences' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import denoise_sequences` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'merge_paired_reads' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import merge_paired_reads` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'classify_taxonomy' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import classify_taxonomy` |
| docs/metagenomics/amplicon.md | 69 | ImportError | Cannot import 'build_taxonomy_tree' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import build_taxonomy_tree` |
| docs/metagenomics/amplicon.md | 69 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/comparative.md | 60 | ImportError | Cannot import 'differential_abundance' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import differential_abundanc` |
| docs/metagenomics/comparative.md | 60 | ImportError | Cannot import 'clr_transform' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import clr_transform` |
| docs/metagenomics/comparative.md | 60 | ImportError | Cannot import 'indicator_species' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import indicator_species` |
| docs/metagenomics/comparative.md | 60 | ImportError | Cannot import 'effect_size_analysis' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import effect_size_analysis` |
| docs/metagenomics/comparative.md | 60 | ImportError | Cannot import 'biomarker_discovery' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import biomarker_discovery` |
| docs/metagenomics/comparative.md | 60 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'alpha_diversity' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import alpha_diversity` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'beta_diversity' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import beta_diversity` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'rarefaction_curve' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import rarefaction_curve` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'rarefy' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import rarefy` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'permanova' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import permanova` |
| docs/metagenomics/diversity.md | 78 | ImportError | Cannot import 'ordination' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import ordination` |
| docs/metagenomics/diversity.md | 78 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/functional.md | 64 | ImportError | Cannot import 'predict_orfs' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import predict_orfs` |
| docs/metagenomics/functional.md | 64 | ImportError | Cannot import 'annotate_genes' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import annotate_genes` |
| docs/metagenomics/functional.md | 64 | ImportError | Cannot import 'classify_gene_families' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import classify_gene_familie` |
| docs/metagenomics/functional.md | 64 | ImportError | Cannot import 'reconstruct_pathways' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import reconstruct_pathways` |
| docs/metagenomics/functional.md | 64 | ImportError | Cannot import 'compare_pathway_profiles' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import compare_pathway_profi` |
| docs/metagenomics/functional.md | 64 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/metagenomics/index.md | 26 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.amplicon` |
| docs/metagenomics/index.md | 26 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.shotgun` |
| docs/metagenomics/index.md | 26 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics.functional` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'assemble_contigs' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import assemble_contigs` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'calculate_assembly_stats' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import calculate_assembly_st` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'scaffold_contigs' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import scaffold_contigs` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'bin_contigs' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import bin_contigs` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'assess_bin_quality' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import assess_bin_quality` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'build_kmer_index' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import build_kmer_index` |
| docs/metagenomics/shotgun.md | 84 | ImportError | Cannot import 'profile_community' from module 'metainformant.metagenomics' | `from metainformant.metagenomics import profile_community` |
| docs/metagenomics/shotgun.md | 84 | AttributeError | 'metainformant' has no attribute 'metagenomics' (module not found) | `metainformant.metagenomics` |
| docs/ml/BROKEN_EXAMPLES.md | 8 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/ml/BROKEN_EXAMPLES.md | 8 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/BROKEN_EXAMPLES.md | 8 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/ml/BROKEN_EXAMPLES.md | 8 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 34 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.deep_learning.sequences` |
| docs/ml/BROKEN_EXAMPLES.md | 34 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/BROKEN_EXAMPLES.md | 47 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.sequence` |
| docs/ml/BROKEN_EXAMPLES.md | 47 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/BROKEN_EXAMPLES.md | 64 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/ml/BROKEN_EXAMPLES.md | 64 | ImportError | Cannot import 'regression' from module 'metainformant.ml' | `from metainformant.ml import regression` |
| docs/ml/BROKEN_EXAMPLES.md | 64 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/ml/BROKEN_EXAMPLES.md | 64 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 79 | ImportError | Cannot import 'explain' from module 'metainformant.ml' | `from metainformant.ml import explain` |
| docs/ml/BROKEN_EXAMPLES.md | 79 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 92 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.regression` |
| docs/ml/BROKEN_EXAMPLES.md | 114 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/ml/BROKEN_EXAMPLES.md | 114 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/BROKEN_EXAMPLES.md | 114 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/ml/BROKEN_EXAMPLES.md | 114 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 132 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/BROKEN_EXAMPLES.md | 132 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/BROKEN_EXAMPLES.md | 132 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 139 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/BROKEN_EXAMPLES.md | 139 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/BROKEN_EXAMPLES.md | 159 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 180 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.interpretability.explainers` |
| docs/ml/BROKEN_EXAMPLES.md | 180 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/ml/BROKEN_EXAMPLES.md | 207 | ImportError | Cannot import 'explain' from module 'metainformant.ml' | `from metainformant.ml import explain` |
| docs/ml/BROKEN_EXAMPLES.md | 207 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 228 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.interpretability.explainers` |
| docs/ml/BROKEN_EXAMPLES.md | 247 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.interpretability.explainers` |
| docs/ml/BROKEN_EXAMPLES.md | 262 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/ml/BROKEN_EXAMPLES.md | 262 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/BROKEN_EXAMPLES.md | 262 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/ml/BROKEN_EXAMPLES.md | 262 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/BROKEN_EXAMPLES.md | 283 | ImportError | Cannot import 'regression' from module 'metainformant.ml' | `from metainformant.ml import regression` |
| docs/ml/BROKEN_EXAMPLES.md | 283 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/README.md | 40 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.ml import ...` |
| docs/ml/SIGNATURE_MISMATCHES.md | 10 | SyntaxError | Invalid Python syntax: expected ':' at line 10 | `def random_search(     model: Any,     X: Any,     y: Any,  ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 17 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `def random_search(model, X, y, param_distributions, ...)` |
| docs/ml/SIGNATURE_MISMATCHES.md | 24 | SyntaxError | Invalid Python syntax: expected ':' at line 10 | `def random_search(     model_fn: Any,                       ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 44 | SyntaxError | Invalid Python syntax: expected ':' at line 11 | `def bayesian_optimization(     model: Any,     X: Any,     y` |
| docs/ml/SIGNATURE_MISMATCHES.md | 59 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def bayesian_optimization(     objective_fn: Any,           ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 78 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def grid_search(     model: Any,     X: Any,     y: Any,    ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 90 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def grid_search(     model_fn: Any,                         ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 104 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def auto_preprocess(     X: Any,     y: Any | None = None,  ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 116 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 5 | `def auto_preprocess(     X: Any,     y: Any | None = None, )` |
| docs/ml/SIGNATURE_MISMATCHES.md | 133 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def cross_validate_biological(     model: Any,     X: np.nda` |
| docs/ml/SIGNATURE_MISMATCHES.md | 148 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 9 | `def cross_validation_scores(     model: Any,     X: np.ndarr` |
| docs/ml/SIGNATURE_MISMATCHES.md | 162 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 10 | `def cross_validate(     model: Any = None,     X: np.ndarray` |
| docs/ml/SIGNATURE_MISMATCHES.md | 207 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def compute_permutation_importance(     model: Any,     X: A` |
| docs/ml/SIGNATURE_MISMATCHES.md | 220 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def compute_permutation_importance(     model: Any,     X: A` |
| docs/ml/SIGNATURE_MISMATCHES.md | 236 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def compute_lime_explanation(     model: Any,     X: Any,   ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 252 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def feature_interaction(     model: Any,     X: Any,     fea` |
| docs/ml/SIGNATURE_MISMATCHES.md | 265 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def compute_shap_values_kernel(     predict_fn: Any,        ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 283 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def boruta_selection(model, X, y, max_iter=100, random_state` |
| docs/ml/SIGNATURE_MISMATCHES.md | 288 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 7 | `def boruta_selection(     X: Any,     y: Any,     max_iter: ` |
| docs/ml/SIGNATURE_MISMATCHES.md | 302 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def recursive_elimination(model, X, y, n_features=None, cv=5` |
| docs/ml/SIGNATURE_MISMATCHES.md | 307 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def recursive_elimination(     model: Any,     X: Any,     y` |
| docs/ml/SIGNATURE_MISMATCHES.md | 392 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/SIGNATURE_MISMATCHES.md | 392 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/SIGNATURE_MISMATCHES.md | 399 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models` |
| docs/ml/SIGNATURE_MISMATCHES.md | 399 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/SIGNATURE_MISMATCHES.md | 405 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/SIGNATURE_MISMATCHES.md | 405 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'random_search' from module 'metainformant.ml' | `from metainformant.ml import random_search` |
| docs/ml/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'bayesian_optimization' from module 'metainformant.ml' | `from metainformant.ml import bayesian_optimization` |
| docs/ml/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'grid_search' from module 'metainformant.ml' | `from metainformant.ml import grid_search` |
| docs/ml/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'model_selection' from module 'metainformant.ml' | `from metainformant.ml import model_selection` |
| docs/ml/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'auto_preprocess' from module 'metainformant.ml' | `from metainformant.ml import auto_preprocess` |
| docs/ml/VALIDATION_REPORT.md | 44 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/VALIDATION_REPORT.md | 91 | ImportError | Cannot import 'sequences' from module 'metainformant.ml.deep_learning' | `from metainformant.ml.deep_learning import sequences` |
| docs/ml/VALIDATION_REPORT.md | 131 | ImportError | Cannot import 'cross_validate_biological' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import cross_val` |
| docs/ml/VALIDATION_REPORT.md | 131 | ImportError | Cannot import 'bootstrap_validation' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import bootstrap` |
| docs/ml/VALIDATION_REPORT.md | 131 | ImportError | Cannot import 'learning_curve_analysis' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import learning_` |
| docs/ml/VALIDATION_REPORT.md | 131 | ImportError | Cannot import 'permutation_test' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import permutati` |
| docs/ml/VALIDATION_REPORT.md | 131 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.evaluation.validation` |
| docs/ml/VALIDATION_REPORT.md | 175 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.features` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'compute_permutation_importance' from module 'metainformant.ml' | `from metainformant.ml import compute_permutation_importance` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'compute_shap_values_kernel' from module 'metainformant.ml' | `from metainformant.ml import compute_shap_values_kernel` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'compute_lime_explanation' from module 'metainformant.ml' | `from metainformant.ml import compute_lime_explanation` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'partial_dependence' from module 'metainformant.ml' | `from metainformant.ml import partial_dependence` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'feature_interaction' from module 'metainformant.ml' | `from metainformant.ml import feature_interaction` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'boruta_selection' from module 'metainformant.ml' | `from metainformant.ml import boruta_selection` |
| docs/ml/VALIDATION_REPORT.md | 220 | ImportError | Cannot import 'stability_selection' from module 'metainformant.ml' | `from metainformant.ml import stability_selection` |
| docs/ml/VALIDATION_REPORT.md | 220 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/VALIDATION_REPORT.md | 294 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/automl.md | 31 | SyntaxError | Invalid Python syntax: expected ':' at line 10 | `def random_search(     model: Any,     X: Any,     y: Any,  ` |
| docs/ml/automl.md | 48 | SyntaxError | Invalid Python syntax: expected ':' at line 11 | `def bayesian_optimization(     model: Any,     X: Any,     y` |
| docs/ml/automl.md | 66 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def grid_search(     model: Any,     X: Any,     y: Any,    ` |
| docs/ml/automl.md | 81 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def model_selection(     X: Any,     y: Any,     task: str =` |
| docs/ml/automl.md | 96 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def auto_preprocess(     X: Any,     y: Any | None = None,  ` |
| docs/ml/automl.md | 111 | ImportError | Cannot import 'random_search' from module 'metainformant.ml' | `from metainformant.ml import random_search` |
| docs/ml/automl.md | 111 | ImportError | Cannot import 'bayesian_optimization' from module 'metainformant.ml' | `from metainformant.ml import bayesian_optimization` |
| docs/ml/automl.md | 111 | ImportError | Cannot import 'grid_search' from module 'metainformant.ml' | `from metainformant.ml import grid_search` |
| docs/ml/automl.md | 111 | ImportError | Cannot import 'model_selection' from module 'metainformant.ml' | `from metainformant.ml import model_selection` |
| docs/ml/automl.md | 111 | ImportError | Cannot import 'auto_preprocess' from module 'metainformant.ml' | `from metainformant.ml import auto_preprocess` |
| docs/ml/automl.md | 111 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/deep_learning.md | 15 | ImportError | Cannot import 'sequences' from module 'metainformant.ml.deep_learning' | `from metainformant.ml.deep_learning import sequences` |
| docs/ml/deep_learning.md | 15 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/deep_learning.md | 15 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.deep_learning` |
| docs/ml/deep_learning.md | 15 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/deep_learning.md | 49 | SyntaxError | Invalid Python syntax: invalid syntax at line 5 | `from metainformant.ml.deep_learning import sequences  # Gene` |
| docs/ml/deep_learning.md | 100 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/ml/deep_learning.md | 100 | ImportError | Cannot import 'sequences' from module 'metainformant.ml.deep_learning' | `from metainformant.ml.deep_learning import sequences` |
| docs/ml/deep_learning.md | 100 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/ml/deep_learning.md | 100 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.deep_learning` |
| docs/ml/deep_learning.md | 116 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/deep_learning.md | 116 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/evaluation.md | 29 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def train_test_split_biological(     X: np.ndarray,     y: n` |
| docs/ml/evaluation.md | 44 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def cross_validate_biological(     model: Any,     X: np.nda` |
| docs/ml/evaluation.md | 60 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def bootstrap_validation(     model: Any,     X: np.ndarray,` |
| docs/ml/evaluation.md | 75 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def learning_curve_analysis(     model: Any,     X: np.ndarr` |
| docs/ml/evaluation.md | 90 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def permutation_test(     model: Any,     X: np.ndarray,    ` |
| docs/ml/evaluation.md | 104 | ImportError | Cannot import 'cross_validate_biological' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import cross_val` |
| docs/ml/evaluation.md | 104 | ImportError | Cannot import 'bootstrap_validation' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import bootstrap` |
| docs/ml/evaluation.md | 104 | ImportError | Cannot import 'learning_curve_analysis' from module 'metainformant.ml.evaluation.validation' | `from metainformant.ml.evaluation.validation import learning_` |
| docs/ml/evaluation.md | 104 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.evaluation.validation` |
| docs/ml/evaluation.md | 104 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/features.md | 28 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def select_features_univariate(     X: np.ndarray,     y: np` |
| docs/ml/features.md | 42 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def pca_reduction(     X: np.ndarray,     n_components: int ` |
| docs/ml/features.md | 54 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def tsne_reduction(     X: np.ndarray,     n_components: int` |
| docs/ml/features.md | 66 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def umap_reduction(     X: np.ndarray,     n_components: int` |
| docs/ml/features.md | 79 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def ica_reduction(     X: np.ndarray,     n_components: int ` |
| docs/ml/features.md | 92 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.features` |
| docs/ml/index.md | 104 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/index.md | 104 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 132 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 160 | ImportError | Cannot import 'validation' from module 'metainformant.ml' | `from metainformant.ml import validation` |
| docs/ml/index.md | 160 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 193 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/ml/index.md | 193 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/index.md | 193 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/ml/index.md | 193 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 212 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/ml/index.md | 212 | ImportError | Cannot import 'regression' from module 'metainformant.ml' | `from metainformant.ml import regression` |
| docs/ml/index.md | 212 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/ml/index.md | 212 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 235 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/ml/index.md | 235 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/index.md | 235 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/ml/index.md | 235 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 264 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 280 | ImportError | Cannot import 'explain' from module 'metainformant.ml' | `from metainformant.ml import explain` |
| docs/ml/index.md | 280 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 300 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/ml/index.md | 300 | ImportError | Cannot import 'classification' from module 'metainformant.ml' | `from metainformant.ml import classification` |
| docs/ml/index.md | 300 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/ml/index.md | 300 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/index.md | 326 | ImportError | Cannot import 'regression' from module 'metainformant.ml' | `from metainformant.ml import regression` |
| docs/ml/index.md | 326 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/interpretability.md | 31 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def compute_permutation_importance(     model: Any,     X: A` |
| docs/ml/interpretability.md | 46 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def compute_shap_values_kernel(     model: Any,     X: Any, ` |
| docs/ml/interpretability.md | 59 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def compute_lime_explanation(     model: Any,     X: Any,   ` |
| docs/ml/interpretability.md | 73 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def partial_dependence(     model: Any,     X: Any,     feat` |
| docs/ml/interpretability.md | 86 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def feature_interaction(     model: Any,     X: Any,     fea` |
| docs/ml/interpretability.md | 98 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def compute_attention_weights(     model: Any,     X: Any, )` |
| docs/ml/interpretability.md | 109 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def boruta_selection(model, X, y, max_iter=100, random_state` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'compute_permutation_importance' from module 'metainformant.ml' | `from metainformant.ml import compute_permutation_importance` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'compute_shap_values_kernel' from module 'metainformant.ml' | `from metainformant.ml import compute_shap_values_kernel` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'compute_lime_explanation' from module 'metainformant.ml' | `from metainformant.ml import compute_lime_explanation` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'partial_dependence' from module 'metainformant.ml' | `from metainformant.ml import partial_dependence` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'feature_interaction' from module 'metainformant.ml' | `from metainformant.ml import feature_interaction` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'boruta_selection' from module 'metainformant.ml' | `from metainformant.ml import boruta_selection` |
| docs/ml/interpretability.md | 120 | ImportError | Cannot import 'stability_selection' from module 'metainformant.ml' | `from metainformant.ml import stability_selection` |
| docs/ml/interpretability.md | 120 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml` |
| docs/ml/llm_integration.md | 15 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/llm_integration.md | 15 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/llm_integration.md | 34 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/llm_integration.md | 34 | ImportError | Cannot import 'LLMConfig' from module 'metainformant.ml.llm.config' | `from metainformant.ml.llm.config import LLMConfig` |
| docs/ml/llm_integration.md | 34 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/llm_integration.md | 34 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm.config` |
| docs/ml/llm_integration.md | 125 | ImportError | Cannot import 'prompts' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import prompts` |
| docs/ml/llm_integration.md | 125 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/llm_integration.md | 140 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/llm_integration.md | 140 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/ml/llm_integration.md | 140 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/llm_integration.md | 154 | ImportError | Cannot import 'differential_expression' from module 'metainformant.rna' | `from metainformant.rna import differential_expression` |
| docs/ml/llm_integration.md | 154 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/llm_integration.md | 154 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/ml/llm_integration.md | 154 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/llm_integration.md | 168 | ImportError | Cannot import 'OllamaClient' from module 'metainformant.ml.llm' | `from metainformant.ml.llm import OllamaClient` |
| docs/ml/llm_integration.md | 168 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.llm` |
| docs/ml/models.md | 39 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `class BiologicalClassifier:     def __init__(         self, ` |
| docs/ml/models.md | 60 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `class BiologicalRegressor:     def __init__(         self,  ` |
| docs/ml/models.md | 80 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.classification` |
| docs/ml/models.md | 80 | AttributeError | 'metainformant' has no attribute 'ml' (module not found) | `metainformant.ml.models.regression` |
| docs/multiomics/README.md | 37 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.colocalization` |
| docs/multiomics/README.md | 37 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.analysis` |
| docs/multiomics/index.md | 97 | ImportError | Cannot import 'MultiOmicsData' from module 'metainformant.multiomics' | `from metainformant.multiomics import MultiOmicsData` |
| docs/multiomics/index.md | 97 | ImportError | Cannot import 'integrate_omics_data' from module 'metainformant.multiomics' | `from metainformant.multiomics import integrate_omics_data` |
| docs/multiomics/index.md | 97 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 125 | ImportError | Cannot import 'joint_pca' from module 'metainformant.multiomics' | `from metainformant.multiomics import joint_pca` |
| docs/multiomics/index.md | 125 | ImportError | Cannot import 'joint_nmf' from module 'metainformant.multiomics' | `from metainformant.multiomics import joint_nmf` |
| docs/multiomics/index.md | 125 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 159 | ImportError | Cannot import 'canonical_correlation' from module 'metainformant.multiomics' | `from metainformant.multiomics import canonical_correlation` |
| docs/multiomics/index.md | 159 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 181 | ImportError | Cannot import 'MultiOmicsData' from module 'metainformant.multiomics' | `from metainformant.multiomics import MultiOmicsData` |
| docs/multiomics/index.md | 181 | ImportError | Cannot import 'from_rna_expression' from module 'metainformant.multiomics' | `from metainformant.multiomics import from_rna_expression` |
| docs/multiomics/index.md | 181 | ImportError | Cannot import 'joint_pca' from module 'metainformant.multiomics' | `from metainformant.multiomics import joint_pca` |
| docs/multiomics/index.md | 181 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 203 | ImportError | Cannot import 'MultiOmicsData' from module 'metainformant.multiomics' | `from metainformant.multiomics import MultiOmicsData` |
| docs/multiomics/index.md | 203 | ImportError | Cannot import 'from_dna_variants' from module 'metainformant.multiomics' | `from metainformant.multiomics import from_dna_variants` |
| docs/multiomics/index.md | 203 | ImportError | Cannot import 'from_rna_expression' from module 'metainformant.multiomics' | `from metainformant.multiomics import from_rna_expression` |
| docs/multiomics/index.md | 203 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 221 | ImportError | Cannot import 'MultiOmicsData' from module 'metainformant.multiomics' | `from metainformant.multiomics import MultiOmicsData` |
| docs/multiomics/index.md | 221 | ImportError | Cannot import 'from_protein_abundance' from module 'metainformant.multiomics' | `from metainformant.multiomics import from_protein_abundance` |
| docs/multiomics/index.md | 221 | ImportError | Cannot import 'from_rna_expression' from module 'metainformant.multiomics' | `from metainformant.multiomics import from_rna_expression` |
| docs/multiomics/index.md | 221 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/index.md | 241 | ImportError | Cannot import 'MultiOmicsData' from module 'metainformant.multiomics' | `from metainformant.multiomics import MultiOmicsData` |
| docs/multiomics/index.md | 241 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/multiomics/integration.md | 123 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 6 | `def integrate_omics_data(     data_dict: Dict[str, Union[pd.` |
| docs/multiomics/integration.md | 151 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 6 | `def joint_pca(     omics_data: MultiOmicsData,     n_compone` |
| docs/multiomics/integration.md | 181 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 8 | `def joint_nmf(     omics_data: MultiOmicsData,     n_compone` |
| docs/multiomics/integration.md | 213 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 6 | `def canonical_correlation(     omics_data: MultiOmicsData,  ` |
| docs/multiomics/integration.md | 244 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 5 | `def from_dna_variants(     vcf_path: Union[str, Path],     s` |
| docs/multiomics/integration.md | 270 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 6 | `def from_rna_expression(     expression_path: Union[str, Pat` |
| docs/multiomics/integration.md | 297 | SyntaxError | Invalid Python syntax: expected an indented block after function definition on line 1 at line 6 | `def from_protein_abundance(     protein_path: Union[str, Pat` |
| docs/multiomics/methods.md | 33 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def joint_nmf(     data_matrices: dict[str, Any],     k: int` |
| docs/multiomics/methods.md | 47 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def mofa_simple(     data_matrices: dict[str, Any],     k: i` |
| docs/multiomics/methods.md | 62 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def tensor_decomposition(     tensor: list[list[list[float]]` |
| docs/multiomics/methods.md | 77 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def similarity_network_fusion(     networks: list[list[list[` |
| docs/multiomics/methods.md | 92 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def canonical_correlation(     X: Any, Y: Any,     n_compone` |
| docs/multiomics/methods.md | 107 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def multi_omic_clustering(     data_matrices: dict[str, Any]` |
| docs/multiomics/methods.md | 121 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def consensus_clustering(     data: Any,     k_range: range ` |
| docs/multiomics/methods.md | 136 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def multi_view_spectral(     similarity_matrices: list[Any],` |
| docs/multiomics/methods.md | 149 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def evaluate_integration(     labels: list[int],     omic_da` |
| docs/multiomics/methods.md | 161 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.methods.factorization` |
| docs/multiomics/methods.md | 161 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.methods.clustering` |
| docs/multiomics/pathways.md | 35 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def multi_omic_enrichment(     gene_sets: dict[str, list[str` |
| docs/multiomics/pathways.md | 56 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def active_module_detection(     network: dict[str, list[str` |
| docs/multiomics/pathways.md | 71 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def pathway_topology_analysis(     pathway_graph: dict[str, ` |
| docs/multiomics/pathways.md | 88 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def cross_omic_pathway_concordance(     pathway_results: dic` |
| docs/multiomics/pathways.md | 103 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.pathways.enrichment` |
| docs/multiomics/survival.md | 34 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def cox_regression(     time: list[float],     event: list[i` |
| docs/multiomics/survival.md | 49 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def kaplan_meier(     time: list[float],     event: list[int` |
| docs/multiomics/survival.md | 63 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def log_rank_test(     time: list[float],     event: list[in` |
| docs/multiomics/survival.md | 76 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def multi_omic_survival_model(     omic_features: dict[str, ` |
| docs/multiomics/survival.md | 92 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def risk_stratification(     risk_scores: list[float],     t` |
| docs/multiomics/survival.md | 108 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def compute_concordance_index(     risk_scores: list[float],` |
| docs/multiomics/survival.md | 121 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.survival.analysis` |
| docs/networks/README.md | 36 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.networks import ...` |
| docs/networks/VALIDATION_REPORT.md | 25 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def create_network(edges: List[Tuple[str, str]], directed: b` |
| docs/networks/VALIDATION_REPORT.md | 43 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def add_edges_from_interactions(     graph: Any,     interac` |
| docs/networks/VALIDATION_REPORT.md | 65 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def add_edges_from_correlation(     graph: Any,     correlat` |
| docs/networks/VALIDATION_REPORT.md | 226 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def pathway_enrichment(     gene_list: List[str],     pathwa` |
| docs/networks/VALIDATION_REPORT.md | 291 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def load_string_interactions(     interactions_df: Any,     ` |
| docs/networks/VALIDATION_REPORT.md | 367 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `infer_grn(expression_data, gene_names, method="correlation",` |
| docs/networks/community.md | 11 | ImportError | Cannot import 'detect_communities' from module 'metainformant.networks' | `from metainformant.networks import detect_communities` |
| docs/networks/community.md | 11 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/community.md | 49 | ImportError | Cannot import 'modularity' from module 'metainformant.networks' | `from metainformant.networks import modularity` |
| docs/networks/community.md | 49 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/community.md | 70 | ImportError | Cannot import 'community_metrics' from module 'metainformant.networks' | `from metainformant.networks import community_metrics` |
| docs/networks/community.md | 70 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/community.md | 158 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/networks/community.md | 158 | ImportError | Cannot import 'community' from module 'metainformant.networks' | `from metainformant.networks import community` |
| docs/networks/community.md | 158 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/community.md | 158 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/networks/community.md | 176 | ImportError | Cannot import 'community' from module 'metainformant.networks' | `from metainformant.networks import community` |
| docs/networks/community.md | 176 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 11 | ImportError | Cannot import 'create_network' from module 'metainformant.networks' | `from metainformant.networks import create_network` |
| docs/networks/graph.md | 11 | ImportError | Cannot import 'add_edges_from_interactions' from module 'metainformant.networks' | `from metainformant.networks import add_edges_from_interactio` |
| docs/networks/graph.md | 11 | ImportError | Cannot import 'add_edges_from_correlation' from module 'metainformant.networks' | `from metainformant.networks import add_edges_from_correlatio` |
| docs/networks/graph.md | 11 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 44 | ImportError | Cannot import 'add_edges_from_correlation' from module 'metainformant.networks' | `from metainformant.networks import add_edges_from_correlatio` |
| docs/networks/graph.md | 44 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 61 | ImportError | Cannot import 'add_edges_from_interactions' from module 'metainformant.networks' | `from metainformant.networks import add_edges_from_interactio` |
| docs/networks/graph.md | 61 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 79 | ImportError | Cannot import 'network_metrics' from module 'metainformant.networks' | `from metainformant.networks import network_metrics` |
| docs/networks/graph.md | 79 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 105 | ImportError | Cannot import 'centrality_measures' from module 'metainformant.networks' | `from metainformant.networks import centrality_measures` |
| docs/networks/graph.md | 105 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 133 | ImportError | Cannot import 'shortest_paths' from module 'metainformant.networks' | `from metainformant.networks import shortest_paths` |
| docs/networks/graph.md | 133 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 188 | ImportError | Cannot import 'create_network' from module 'metainformant.networks' | `from metainformant.networks import create_network` |
| docs/networks/graph.md | 188 | ImportError | Cannot import 'network_metrics' from module 'metainformant.networks' | `from metainformant.networks import network_metrics` |
| docs/networks/graph.md | 188 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/graph.md | 213 | ImportError | Cannot import 'create_network' from module 'metainformant.networks' | `from metainformant.networks import create_network` |
| docs/networks/graph.md | 213 | ImportError | Cannot import 'centrality_measures' from module 'metainformant.networks' | `from metainformant.networks import centrality_measures` |
| docs/networks/graph.md | 213 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 92 | ImportError | Cannot import 'create_network' from module 'metainformant.networks' | `from metainformant.networks import create_network` |
| docs/networks/index.md | 92 | ImportError | Cannot import 'add_edges_from_interactions' from module 'metainformant.networks' | `from metainformant.networks import add_edges_from_interactio` |
| docs/networks/index.md | 92 | ImportError | Cannot import 'network_metrics' from module 'metainformant.networks' | `from metainformant.networks import network_metrics` |
| docs/networks/index.md | 92 | ImportError | Cannot import 'centrality_measures' from module 'metainformant.networks' | `from metainformant.networks import centrality_measures` |
| docs/networks/index.md | 92 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 120 | ImportError | Cannot import 'detect_communities' from module 'metainformant.networks' | `from metainformant.networks import detect_communities` |
| docs/networks/index.md | 120 | ImportError | Cannot import 'modularity' from module 'metainformant.networks' | `from metainformant.networks import modularity` |
| docs/networks/index.md | 120 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 138 | ImportError | Cannot import 'pathway_enrichment' from module 'metainformant.networks' | `from metainformant.networks import pathway_enrichment` |
| docs/networks/index.md | 138 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 162 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/networks/index.md | 162 | ImportError | Cannot import 'infer_grn' from module 'metainformant.networks' | `from metainformant.networks import infer_grn` |
| docs/networks/index.md | 162 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/networks/index.md | 162 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 182 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/networks/index.md | 182 | ImportError | Cannot import 'community' from module 'metainformant.networks' | `from metainformant.networks import community` |
| docs/networks/index.md | 182 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/index.md | 182 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/networks/pathway.md | 11 | ImportError | Cannot import 'PathwayNetwork' from module 'metainformant.networks' | `from metainformant.networks import PathwayNetwork` |
| docs/networks/pathway.md | 11 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 39 | ImportError | Cannot import 'load_pathway_database' from module 'metainformant.networks' | `from metainformant.networks import load_pathway_database` |
| docs/networks/pathway.md | 39 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 78 | ImportError | Cannot import 'pathway_enrichment' from module 'metainformant.networks' | `from metainformant.networks import pathway_enrichment` |
| docs/networks/pathway.md | 78 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 112 | ImportError | Cannot import 'network_enrichment_analysis' from module 'metainformant.networks' | `from metainformant.networks import network_enrichment_analys` |
| docs/networks/pathway.md | 112 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 139 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/networks/pathway.md | 139 | ImportError | Cannot import 'multi_omics_pathway_analysis' from module 'metainformant.networks' | `from metainformant.networks import multi_omics_pathway_analy` |
| docs/networks/pathway.md | 139 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/networks/pathway.md | 139 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 169 | ImportError | Cannot import 'pathway_activity_inference' from module 'metainformant.networks' | `from metainformant.networks import pathway_activity_inferenc` |
| docs/networks/pathway.md | 169 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 260 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/networks/pathway.md | 260 | ImportError | Cannot import 'pathway_enrichment' from module 'metainformant.networks' | `from metainformant.networks import pathway_enrichment` |
| docs/networks/pathway.md | 260 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/networks/pathway.md | 260 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/pathway.md | 298 | ImportError | Cannot import 'proteomes' from module 'metainformant.protein' | `from metainformant.protein import proteomes` |
| docs/networks/pathway.md | 298 | ImportError | Cannot import 'pathway_enrichment' from module 'metainformant.networks' | `from metainformant.networks import pathway_enrichment` |
| docs/networks/pathway.md | 298 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/networks/pathway.md | 298 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 11 | ImportError | Cannot import 'ProteinNetwork' from module 'metainformant.networks' | `from metainformant.networks import ProteinNetwork` |
| docs/networks/ppi.md | 11 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 31 | ImportError | Cannot import 'load_string_interactions' from module 'metainformant.networks' | `from metainformant.networks import load_string_interactions` |
| docs/networks/ppi.md | 31 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 51 | ImportError | Cannot import 'predict_interactions' from module 'metainformant.networks' | `from metainformant.networks import predict_interactions` |
| docs/networks/ppi.md | 51 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 81 | ImportError | Cannot import 'get_protein_partners' from module 'metainformant.networks' | `from metainformant.networks import get_protein_partners` |
| docs/networks/ppi.md | 81 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 103 | ImportError | Cannot import 'filter_by_confidence' from module 'metainformant.networks' | `from metainformant.networks import filter_by_confidence` |
| docs/networks/ppi.md | 103 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 126 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/networks/ppi.md | 174 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/networks/ppi.md | 174 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/networks/ppi.md | 174 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/networks/ppi.md | 174 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 197 | ImportError | Cannot import 'proteomes' from module 'metainformant.protein' | `from metainformant.protein import proteomes` |
| docs/networks/ppi.md | 197 | ImportError | Cannot import 'ppi' from module 'metainformant.networks' | `from metainformant.networks import ppi` |
| docs/networks/ppi.md | 197 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/networks/ppi.md | 197 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/ppi.md | 197 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/networks/regulatory.md | 11 | ImportError | Cannot import 'GeneRegulatoryNetwork' from module 'metainformant.networks' | `from metainformant.networks import GeneRegulatoryNetwork` |
| docs/networks/regulatory.md | 11 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 32 | ImportError | Cannot import 'infer_grn' from module 'metainformant.networks' | `from metainformant.networks import infer_grn` |
| docs/networks/regulatory.md | 32 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 61 | ImportError | Cannot import 'add_transcription_factor' from module 'metainformant.networks' | `from metainformant.networks import add_transcription_factor` |
| docs/networks/regulatory.md | 61 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 94 | ImportError | Cannot import 'get_targets' from module 'metainformant.networks' | `from metainformant.networks import get_targets` |
| docs/networks/regulatory.md | 94 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 112 | ImportError | Cannot import 'get_regulators' from module 'metainformant.networks' | `from metainformant.networks import get_regulators` |
| docs/networks/regulatory.md | 112 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 129 | ImportError | Cannot import 'filter_by_confidence' from module 'metainformant.networks' | `from metainformant.networks import filter_by_confidence` |
| docs/networks/regulatory.md | 129 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 152 | ImportError | Cannot import 'regulatory_motifs' from module 'metainformant.networks' | `from metainformant.networks import regulatory_motifs` |
| docs/networks/regulatory.md | 152 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 181 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/networks/regulatory.md | 181 | ImportError | Cannot import 'infer_grn' from module 'metainformant.networks' | `from metainformant.networks import infer_grn` |
| docs/networks/regulatory.md | 181 | ImportError | Cannot import 'regulatory_motifs' from module 'metainformant.networks' | `from metainformant.networks import regulatory_motifs` |
| docs/networks/regulatory.md | 181 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/networks/regulatory.md | 181 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/networks/regulatory.md | 207 | ImportError | Cannot import 'motifs' from module 'metainformant.dna' | `from metainformant.dna import motifs` |
| docs/networks/regulatory.md | 207 | ImportError | Cannot import 'add_transcription_factor' from module 'metainformant.networks' | `from metainformant.networks import add_transcription_factor` |
| docs/networks/regulatory.md | 207 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/networks/regulatory.md | 207 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/ontology/README.md | 36 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.ontology import ...` |
| docs/ontology/go.md | 7 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core.go` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'common_ancestors' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import common_ancestors` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'path_to_root' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import path_to_root` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'distance' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import distance` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'find_term_by_name' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import find_term_by_name` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'filter_by_namespace' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import filter_by_namespa` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'get_roots' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import get_roots` |
| docs/ontology/go.md | 24 | ImportError | Cannot import 'get_leaves' from module 'metainformant.ontology.core.go' | `from metainformant.ontology.core.go import get_leaves` |
| docs/ontology/go.md | 24 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core.go` |
| docs/ontology/go.md | 50 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core.go` |
| docs/ontology/obo.md | 33 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def parse_obo(path: str | Path) -> Ontology` |
| docs/ontology/obo.md | 45 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `>>> ontology = parse_obo("data/go.obo") >>> len(ontology) 50` |
| docs/ontology/obo.md | 53 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_obo_format(path: str | Path) -> tuple[bool, Lis` |
| docs/ontology/obo.md | 66 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_obo_statistics(path: str | Path) -> Dict[str, Any]` |
| docs/ontology/obo.md | 93 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core.obo` |
| docs/ontology/pathway_enrichment.md | 25 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def over_representation_analysis(     gene_list: list[str], ` |
| docs/ontology/pathway_enrichment.md | 44 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def compute_enrichment_score(     ranked_list: list[str],   ` |
| docs/ontology/pathway_enrichment.md | 59 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def gsea(     ranked_genes: list[tuple],     gene_sets: dict` |
| docs/ontology/pathway_enrichment.md | 76 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def pathway_network(     enrichment_results: list[dict],    ` |
| docs/ontology/pathway_enrichment.md | 90 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def compare_enrichments(     results_a: list[dict],     resu` |
| docs/ontology/pathway_enrichment.md | 115 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.pathway_enrichment.enrichment` |
| docs/ontology/query.md | 26 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def ancestors(     onto: Ontology, term_id: str, relation_ty` |
| docs/ontology/query.md | 37 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def descendants(     onto: Ontology, term_id: str, relation_` |
| docs/ontology/query.md | 47 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def common_ancestors(     onto: Ontology, term1: str, term2:` |
| docs/ontology/query.md | 57 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def most_informative_common_ancestor(     onto: Ontology, te` |
| docs/ontology/query.md | 69 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def path_to_root(     onto: Ontology, term_id: str, relation` |
| docs/ontology/query.md | 80 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def shortest_path(     onto: Ontology, term1: str, term2: st` |
| docs/ontology/query.md | 91 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def distance(     onto: Ontology, term1: str, term2: str, re` |
| docs/ontology/query.md | 101 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def get_subontology(     onto: Ontology, root_terms: Iterabl` |
| docs/ontology/query.md | 112 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def subgraph(     onto: Ontology, term_ids: List[str], relat` |
| docs/ontology/query.md | 123 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def find_terms_by_name(onto: Ontology, name_pattern: str, ca` |
| docs/ontology/query.md | 135 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_roots(onto: Ontology, relation_type: str = "is_a") -` |
| docs/ontology/query.md | 144 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def information_content(onto: Ontology, term_id: str, corpus` |
| docs/ontology/query.md | 153 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_subontology_stats(onto: Ontology) -> Dict[str, Any]` |
| docs/ontology/query.md | 163 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_ontology_integrity(onto: Ontology) -> Tuple[boo` |
| docs/ontology/query.md | 172 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def clear_cache() -> None def set_cache_enabled(enabled: boo` |
| docs/ontology/query.md | 179 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core.obo` |
| docs/ontology/query.md | 179 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.query.query` |
| docs/pharmacogenomics/CAPABILITIES.md | 1277 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `,         the built-in database from :func:` |
| docs/pharmacogenomics/CAPABILITIES.md | 1306 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Assess polypharmacy risk from multiple concurrent medicati` |
| docs/pharmacogenomics/CAPABILITIES.md | 1321 | SyntaxError | Invalid Python syntax: unmatched ')' at line 1 | `).         - competing_pathways: Dict mapping CYP enzyme to ` |
| docs/pharmacogenomics/CAPABILITIES.md | 1376 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Compute CPIC-style activity score from a diplotype string.` |
| docs/pharmacogenomics/CAPABILITIES.md | 1380 | SyntaxError | Invalid Python syntax: unexpected indent at line 1 | ` and sums the function values of each allele from the provid` |
| docs/pharmacogenomics/CAPABILITIES.md | 1396 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Classify metabolizer phenotype from an activity score.  Us` |
| docs/pharmacogenomics/CAPABILITIES.md | 1413 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Predict metabolizer phenotype from genotype information.  ` |
| docs/pharmacogenomics/CONFIGURATION.md | 17 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/pharmacogenomics/CONFIGURATION.md | 58 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/pharmacogenomics/CONFIGURATION.md | 97 | ImportError | Cannot import 'load_all_allele_definitions' from module 'metainformant.pharmacogenomics.alleles.star_allele' | `from metainformant.pharmacogenomics.alleles.star_allele impo` |
| docs/pharmacogenomics/CONFIGURATION.md | 97 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/CONFIGURATION.md | 108 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.cpic` |
| docs/pharmacogenomics/CONFIGURATION.md | 129 | ImportError | Cannot import 'refresh_pharmgkb_cache' from module 'metainformant.pharmacogenomics.annotations.pharmgkb' | `from metainformant.pharmacogenomics.annotations.pharmgkb imp` |
| docs/pharmacogenomics/CONFIGURATION.md | 129 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.pharmgkb` |
| docs/pharmacogenomics/CONFIGURATION.md | 138 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/pharmacogenomics/EXAMPLES.md | 5 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/EXAMPLES.md | 5 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.diplotype` |
| docs/pharmacogenomics/EXAMPLES.md | 5 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.phenotype` |
| docs/pharmacogenomics/EXAMPLES.md | 5 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.cpic` |
| docs/pharmacogenomics/EXAMPLES.md | 23 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/pharmacogenomics/EXAMPLES.md | 47 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.cpic` |
| docs/pharmacogenomics/EXAMPLES.md | 63 | ImportError | Cannot import 'analyze_drug_gene_interactions' from module 'metainformant.pharmacogenomics.interaction.drug_interactions' | `from metainformant.pharmacogenomics.interaction.drug_interac` |
| docs/pharmacogenomics/EXAMPLES.md | 63 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.interaction.drug_interactions` |
| docs/pharmacogenomics/EXAMPLES.md | 90 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.pathogenicity` |
| docs/pharmacogenomics/EXAMPLES.md | 121 | ImportError | Cannot import 'ClinicalReportBuilder' from module 'metainformant.pharmacogenomics.clinical.reporting' | `from metainformant.pharmacogenomics.clinical.reporting impor` |
| docs/pharmacogenomics/EXAMPLES.md | 121 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.reporting` |
| docs/pharmacogenomics/EXAMPLES.md | 150 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/EXAMPLES.md | 182 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.phenotype` |
| docs/pharmacogenomics/EXAMPLES.md | 203 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.visualization.plots` |
| docs/pharmacogenomics/EXAMPLES.md | 218 | ImportError | Cannot import 'get_drug_phenotypes' from module 'metainformant.pharmacogenomics.annotations.pharmgkb' | `from metainformant.pharmacogenomics.annotations.pharmgkb imp` |
| docs/pharmacogenomics/EXAMPLES.md | 218 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.pharmgkb` |
| docs/pharmacogenomics/EXAMPLES.md | 274 | ImportError | Cannot import 'Client' from module 'dask.distributed' | `from dask.distributed import Client` |
| docs/pharmacogenomics/EXAMPLES.md | 312 | ImportError | Cannot import 'get_cached_pharmgkb_path' from module 'metainformant.pharmacogenomics.annotations.pharmgkb' | `from metainformant.pharmacogenomics.annotations.pharmgkb imp` |
| docs/pharmacogenomics/EXAMPLES.md | 312 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.pharmgkb` |
| docs/pharmacogenomics/EXAMPLES.md | 324 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/EXAMPLES.md | 336 | ImportError | Cannot import '_ACTIVITY_SCORE_TABLES' from module 'metainformant.pharmacogenomics.alleles.activity' | `from metainformant.pharmacogenomics.alleles.activity import ` |
| docs/pharmacogenomics/EXAMPLES.md | 336 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.activity` |
| docs/pharmacogenomics/EXAMPLES.md | 373 | ImportError | Cannot import 'FastAPI' from module 'fastapi' | `from fastapi import FastAPI` |
| docs/pharmacogenomics/EXAMPLES.md | 373 | ImportError | Cannot import 'HTTPException' from module 'fastapi' | `from fastapi import HTTPException` |
| docs/pharmacogenomics/EXAMPLES.md | 373 | ImportError | Cannot import 'BaseModel' from module 'pydantic' | `from pydantic import BaseModel` |
| docs/pharmacogenomics/EXAMPLES.md | 373 | ModuleNotFoundError | Module 'uvicorn' not found in project or standard library | `import uvicorn` |
| docs/pharmacogenomics/EXAMPLES.md | 423 | ModuleNotFoundError | Module 'vcfpy' not found in project or standard library | `import vcfpy` |
| docs/pharmacogenomics/GETTING_STARTED.md | 9 | ImportError | Cannot import 'call_star_alleles' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import call_star_alleles` |
| docs/pharmacogenomics/GETTING_STARTED.md | 9 | ImportError | Cannot import 'determine_diplotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import determine_diploty` |
| docs/pharmacogenomics/GETTING_STARTED.md | 9 | ImportError | Cannot import 'classify_phenotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import classify_phenotyp` |
| docs/pharmacogenomics/GETTING_STARTED.md | 9 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics` |
| docs/pharmacogenomics/GETTING_STARTED.md | 9 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.cpic` |
| docs/pharmacogenomics/GETTING_STARTED.md | 63 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/pharmacogenomics/GETTING_STARTED.md | 80 | ImportError | Cannot import 'Client' from module 'dask.distributed' | `from dask.distributed import Client` |
| docs/pharmacogenomics/GETTING_STARTED.md | 94 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.reporting` |
| docs/pharmacogenomics/GETTING_STARTED.md | 131 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/pharmacogenomics/GETTING_STARTED.md | 159 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/INTEGRATION.md | 8 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.reporting` |
| docs/pharmacogenomics/INTEGRATION.md | 34 | ImportError | Cannot import 'FastAPI' from module 'fastapi' | `from fastapi import FastAPI` |
| docs/pharmacogenomics/INTEGRATION.md | 34 | ImportError | Cannot import 'HTTPException' from module 'fastapi' | `from fastapi import HTTPException` |
| docs/pharmacogenomics/INTEGRATION.md | 34 | ImportError | Cannot import 'BaseModel' from module 'pydantic' | `from pydantic import BaseModel` |
| docs/pharmacogenomics/INTEGRATION.md | 96 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `rule pgx_phenotype:     input:         vcf='data/{sample}.vc` |
| docs/pharmacogenomics/PERFORMANCE.md | 33 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/pharmacogenomics/PERFORMANCE.md | 83 | ImportError | Cannot import 'disk_cache' from module 'metainformant.core.cache' | `from metainformant.core.cache import disk_cache` |
| docs/pharmacogenomics/PERFORMANCE.md | 83 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.cache` |
| docs/pharmacogenomics/PERFORMANCE.md | 152 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics` |
| docs/pharmacogenomics/PERFORMANCE.md | 176 | ModuleNotFoundError | Module 'polars' not found in project or standard library | `import polars` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'call_star_alleles' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import call_star_alleles` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'determine_diplotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import determine_diploty` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'classify_phenotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import classify_phenotyp` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'lookup_drug_gene' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import lookup_drug_gene` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'get_dosing_recommendation' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import get_dosing_recomm` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'generate_clinical_report' from module 'metainformant.pharmacogenomics.clinical' | `from metainformant.pharmacogenomics.clinical import generate` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'export_report' from module 'metainformant.pharmacogenomics.clinical' | `from metainformant.pharmacogenomics.clinical import export_r` |
| docs/pharmacogenomics/README.md | 32 | ImportError | Cannot import 'classify_variant_acmg' from module 'metainformant.pharmacogenomics.clinical' | `from metainformant.pharmacogenomics.clinical import classify` |
| docs/pharmacogenomics/README.md | 32 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics` |
| docs/pharmacogenomics/README.md | 32 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical` |
| docs/pharmacogenomics/README.md | 32 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations` |
| docs/pharmacogenomics/TROUBLESHOOTING.md | 29 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.pharmacogenomics.alleles.star_allele impo` |
| docs/pharmacogenomics/TROUBLESHOOTING.md | 35 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.pharmacogenomics.annotations.cpic import ` |
| docs/pharmacogenomics/TROUBLESHOOTING.md | 42 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.pharmacogenomics.alleles.activity import ` |
| docs/pharmacogenomics/acmg.md | 23 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.pathogenicity` |
| docs/pharmacogenomics/cpic.md | 16 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations.cpic` |
| docs/pharmacogenomics/drug_interactions.md | 8 | ImportError | Cannot import 'analyze_drug_gene_interactions' from module 'metainformant.pharmacogenomics.interaction.drug_interactions' | `from metainformant.pharmacogenomics.interaction.drug_interac` |
| docs/pharmacogenomics/drug_interactions.md | 8 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.interaction.drug_interactions` |
| docs/pharmacogenomics/drug_interactions.md | 95 | ImportError | Cannot import 'MY_DRUG' from module 'metainformant.pharmacogenomics.interaction.data.my_drug' | `from metainformant.pharmacogenomics.interaction.data.my_drug` |
| docs/pharmacogenomics/drug_interactions.md | 95 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.interaction.data.my_drug` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'call_star_alleles' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import call_star_alleles` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'determine_diplotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import determine_diploty` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'classify_phenotype' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import classify_phenotyp` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'lookup_drug_gene' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import lookup_drug_gene` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'get_dosing_recommendation' from module 'metainformant.pharmacogenomics' | `from metainformant.pharmacogenomics import get_dosing_recomm` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'generate_clinical_report' from module 'metainformant.pharmacogenomics.clinical' | `from metainformant.pharmacogenomics.clinical import generate` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'classify_variant_acmg' from module 'metainformant.pharmacogenomics.clinical' | `from metainformant.pharmacogenomics.clinical import classify` |
| docs/pharmacogenomics/index.md | 48 | ImportError | Cannot import 'analyze_drug_gene_interactions' from module 'metainformant.pharmacogenomics.interaction' | `from metainformant.pharmacogenomics.interaction import analy` |
| docs/pharmacogenomics/index.md | 48 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics` |
| docs/pharmacogenomics/index.md | 48 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.annotations` |
| docs/pharmacogenomics/index.md | 48 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical` |
| docs/pharmacogenomics/index.md | 48 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.interaction` |
| docs/pharmacogenomics/reporting.md | 22 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.clinical.reporting` |
| docs/pharmacogenomics/reporting.md | 46 | SyntaxError | Invalid Python syntax: invalid character '‑' (U+2011) at line 3 | `{   'patient': patient_dict,   'genotype_table': list of per` |
| docs/pharmacogenomics/star_alleles.md | 45 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/pharmacogenomics/star_alleles.md | 66 | AttributeError | 'metainformant' has no attribute 'pharmacogenomics' (module not found) | `metainformant.pharmacogenomics.alleles.star_allele` |
| docs/phenotype/README.md | 48 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.phenotype import ...` |
| docs/phenotype/VALIDATION_REPORT.md | 35 | ImportError | Cannot import 'Measurement' from module 'measurement' | `from measurement import Measurement` |
| docs/phenotype/VALIDATION_REPORT.md | 35 | ImportError | Cannot import 'MorphometricProfile' from module 'profile' | `from profile import MorphometricProfile` |
| docs/phenotype/VALIDATION_REPORT.md | 57 | ImportError | Cannot import 'Ethogram' from module 'ethogram' | `from ethogram import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 57 | ImportError | Cannot import 'BehaviorSequence' from module 'sequence' | `from sequence import BehaviorSequence` |
| docs/phenotype/VALIDATION_REPORT.md | 63 | ImportError | Cannot import 'Ethogram' from module 'metainformant.phenotype.behavior' | `from metainformant.phenotype.behavior import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 63 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.behavior.ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 63 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.behavior.sequence` |
| docs/phenotype/VALIDATION_REPORT.md | 63 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.behavior` |
| docs/phenotype/VALIDATION_REPORT.md | 81 | ImportError | Cannot import 'Compound' from module 'compound' | `from compound import Compound` |
| docs/phenotype/VALIDATION_REPORT.md | 81 | ImportError | Cannot import 'ChemicalProfile' from module 'profile' | `from profile import ChemicalProfile` |
| docs/phenotype/VALIDATION_REPORT.md | 228 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.gwas_integration.phewas` |
| docs/phenotype/VALIDATION_REPORT.md | 232 | ImportError | Cannot import 'run_phewas' from module 'metainformant.phenotype.gwas_integration' | `from metainformant.phenotype.gwas_integration import run_phe` |
| docs/phenotype/VALIDATION_REPORT.md | 232 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.gwas_integration` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | ImportError | Cannot import 'Measurement' from module 'metainformant.phenotype' | `from metainformant.phenotype import Measurement` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | ImportError | Cannot import 'Ethogram' from module 'metainformant.phenotype' | `from metainformant.phenotype import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | ImportError | Cannot import 'Compound' from module 'metainformant.phenotype' | `from metainformant.phenotype import Compound` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | ImportError | Cannot import 'TrackingPoint' from module 'metainformant.phenotype' | `from metainformant.phenotype import TrackingPoint` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | ImportError | Cannot import 'AcousticSignal' from module 'metainformant.phenotype' | `from metainformant.phenotype import AcousticSignal` |
| docs/phenotype/VALIDATION_REPORT.md | 319 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/VALIDATION_REPORT.md | 323 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.phenotype.morphological.measurement impor` |
| docs/phenotype/VALIDATION_REPORT.md | 359 | ImportError | Cannot import 'Ethogram' from module 'metainformant.phenotype.behavior' | `from metainformant.phenotype.behavior import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 359 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.behavior` |
| docs/phenotype/VALIDATION_REPORT.md | 363 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.behavior.ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 371 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# behavior/__init__.py   from .ethogram import Ethogram   fr` |
| docs/phenotype/VALIDATION_REPORT.md | 402 | ImportError | Cannot import 'run_phewas' from module 'phewas' | `from phewas import run_phewas` |
| docs/phenotype/VALIDATION_REPORT.md | 402 | ImportError | Cannot import 'genetic_risk_score' from module 'phewas' | `from phewas import genetic_risk_score` |
| docs/phenotype/VALIDATION_REPORT.md | 402 | ImportError | Cannot import 'phenotype_correlation_matrix' from module 'phewas' | `from phewas import phenotype_correlation_matrix` |
| docs/phenotype/VALIDATION_REPORT.md | 402 | ImportError | Cannot import 'phenotype_heritability_screen' from module 'phewas' | `from phewas import phenotype_heritability_screen` |
| docs/phenotype/VALIDATION_REPORT.md | 402 | ImportError | Cannot import 'categorize_phenotypes' from module 'phewas' | `from phewas import categorize_phenotypes` |
| docs/phenotype/VALIDATION_REPORT.md | 412 | SyntaxError | Invalid Python syntax: invalid syntax. Perhaps you forgot a comma? at line 2 | `from . import (analysis, behavior, chemical, data, electroni` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'Measurement' from module 'morphological.measurement' | `from morphological.measurement import Measurement` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'MorphometricProfile' from module 'morphological.profile' | `from morphological.profile import MorphometricProfile` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'Ethogram' from module 'behavior.ethogram' | `from behavior.ethogram import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'BehaviorSequence' from module 'behavior.sequence' | `from behavior.sequence import BehaviorSequence` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'Compound' from module 'chemical.compound' | `from chemical.compound import Compound` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'ChemicalProfile' from module 'chemical.profile' | `from chemical.profile import ChemicalProfile` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'TrackingPoint' from module 'electronic.tracking' | `from electronic.tracking import TrackingPoint` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'Trajectory' from module 'electronic.tracking' | `from electronic.tracking import Trajectory` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'AcousticSignal' from module 'sonic.signal' | `from sonic.signal import AcousticSignal` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'PhenotypePipeline' from module 'workflow.pipeline' | `from workflow.pipeline import PhenotypePipeline` |
| docs/phenotype/VALIDATION_REPORT.md | 420 | ImportError | Cannot import 'PipelineConfig' from module 'workflow.pipeline' | `from workflow.pipeline import PipelineConfig` |
| docs/phenotype/VALIDATION_REPORT.md | 456 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.morphological.measurement` |
| docs/phenotype/VALIDATION_REPORT.md | 527 | ImportError | Cannot import 'Measurement' from module 'metainformant.phenotype' | `from metainformant.phenotype import Measurement` |
| docs/phenotype/VALIDATION_REPORT.md | 527 | ImportError | Cannot import 'Ethogram' from module 'metainformant.phenotype' | `from metainformant.phenotype import Ethogram` |
| docs/phenotype/VALIDATION_REPORT.md | 527 | ImportError | Cannot import 'Compound' from module 'metainformant.phenotype' | `from metainformant.phenotype import Compound` |
| docs/phenotype/VALIDATION_REPORT.md | 527 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/antwiki.md | 7 | ImportError | Cannot import 'antwiki' from module 'metainformant.phenotype' | `from metainformant.phenotype import antwiki` |
| docs/phenotype/antwiki.md | 7 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/antwiki.md | 19 | ImportError | Cannot import 'AntWikiScraper' from module 'metainformant.phenotype.scraper' | `from metainformant.phenotype.scraper import AntWikiScraper` |
| docs/phenotype/antwiki.md | 19 | ImportError | Cannot import 'load_scraper_config' from module 'metainformant.phenotype.scraper' | `from metainformant.phenotype.scraper import load_scraper_con` |
| docs/phenotype/antwiki.md | 19 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.scraper` |
| docs/phenotype/index.md | 32 | ImportError | Cannot import 'Measurement' from module 'metainformant.phenotype' | `from metainformant.phenotype import Measurement` |
| docs/phenotype/index.md | 32 | ImportError | Cannot import 'MorphometricProfile' from module 'metainformant.phenotype' | `from metainformant.phenotype import MorphometricProfile` |
| docs/phenotype/index.md | 32 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 56 | ImportError | Cannot import 'Ethogram' from module 'metainformant.phenotype' | `from metainformant.phenotype import Ethogram` |
| docs/phenotype/index.md | 56 | ImportError | Cannot import 'BehaviorSequence' from module 'metainformant.phenotype' | `from metainformant.phenotype import BehaviorSequence` |
| docs/phenotype/index.md | 56 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 77 | ImportError | Cannot import 'Compound' from module 'metainformant.phenotype' | `from metainformant.phenotype import Compound` |
| docs/phenotype/index.md | 77 | ImportError | Cannot import 'ChemicalProfile' from module 'metainformant.phenotype' | `from metainformant.phenotype import ChemicalProfile` |
| docs/phenotype/index.md | 77 | ImportError | Cannot import 'distance_matrix' from module 'metainformant.phenotype.chemical' | `from metainformant.phenotype.chemical import distance_matrix` |
| docs/phenotype/index.md | 77 | ImportError | Cannot import 'identify_marker_compounds' from module 'metainformant.phenotype.chemical' | `from metainformant.phenotype.chemical import identify_marker` |
| docs/phenotype/index.md | 77 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 77 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.chemical` |
| docs/phenotype/index.md | 88 | ImportError | Cannot import 'TrackingPoint' from module 'metainformant.phenotype' | `from metainformant.phenotype import TrackingPoint` |
| docs/phenotype/index.md | 88 | ImportError | Cannot import 'Trajectory' from module 'metainformant.phenotype' | `from metainformant.phenotype import Trajectory` |
| docs/phenotype/index.md | 88 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 109 | ImportError | Cannot import 'AcousticSignal' from module 'metainformant.phenotype' | `from metainformant.phenotype import AcousticSignal` |
| docs/phenotype/index.md | 109 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 127 | ImportError | Cannot import 'PhenotypePipeline' from module 'metainformant.phenotype' | `from metainformant.phenotype import PhenotypePipeline` |
| docs/phenotype/index.md | 127 | ImportError | Cannot import 'PipelineConfig' from module 'metainformant.phenotype' | `from metainformant.phenotype import PipelineConfig` |
| docs/phenotype/index.md | 127 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/index.md | 145 | ImportError | Cannot import 'run_phewas' from module 'metainformant.phenotype' | `from metainformant.phenotype import run_phewas` |
| docs/phenotype/index.md | 145 | ImportError | Cannot import 'phenotype_correlation_matrix' from module 'metainformant.phenotype' | `from metainformant.phenotype import phenotype_correlation_ma` |
| docs/phenotype/index.md | 145 | ImportError | Cannot import 'genetic_risk_score' from module 'metainformant.phenotype' | `from metainformant.phenotype import genetic_risk_score` |
| docs/phenotype/index.md | 145 | ImportError | Cannot import 'phenotype_heritability_screen' from module 'metainformant.phenotype' | `from metainformant.phenotype import phenotype_heritability_s` |
| docs/phenotype/index.md | 145 | ImportError | Cannot import 'categorize_phenotypes' from module 'metainformant.phenotype' | `from metainformant.phenotype import categorize_phenotypes` |
| docs/phenotype/index.md | 145 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype` |
| docs/phenotype/life_course.md | 23 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.analysis.life_course` |
| docs/phenotype/life_course.md | 45 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.analysis.life_course` |
| docs/phenotype/life_course.md | 60 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def extract_phenotypes_from_events(     event_sequence: Even` |
| docs/phenotype/life_course.md | 73 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def aggregate_temporal_phenotypes(     sequences: List[Event` |
| docs/phenotype/life_course.md | 87 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def analyze_life_course_trajectories(     sequences: List[Ev` |
| docs/phenotype/life_course.md | 100 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def identify_critical_periods(     sequences: List[EventSequ` |
| docs/phenotype/life_course.md | 112 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def predict_life_course_outcomes(     sequences: List[EventS` |
| docs/phenotype/life_course.md | 125 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def identify_trajectory_patterns(     sequences: List[EventS` |
| docs/phenotype/life_course.md | 137 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def analyze_life_course(     sequences: List[EventSequence],` |
| docs/phenotype/life_course.md | 149 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def create_life_course_report(     sequences: List[EventSequ` |
| docs/phenotype/life_course.md | 162 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.analysis.life_course` |
| docs/phenotype/visualization.md | 22 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def plot_trait_distribution(     trait_values: np.ndarray,  ` |
| docs/phenotype/visualization.md | 39 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def plot_trait_correlation_matrix(     trait_data: pd.DataFr` |
| docs/phenotype/visualization.md | 55 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def plot_life_course_trajectory(     life_events: List[Dict[` |
| docs/phenotype/visualization.md | 72 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def plot_morphological_measurements(     measurements: Dict[` |
| docs/phenotype/visualization.md | 88 | SyntaxError | Invalid Python syntax: expected ':' at line 10 | `def plot_behavioral_patterns(     behavioral_data: pd.DataFr` |
| docs/phenotype/visualization.md | 106 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def plot_phenotype_pca(     phenotype_data: np.ndarray,     ` |
| docs/phenotype/visualization.md | 123 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def plot_trait_heritability(     heritability_estimates: Dic` |
| docs/phenotype/visualization.md | 140 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def plot_life_history_comparison(     species_data: Dict[str` |
| docs/phenotype/visualization.md | 155 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def plot_phenotype_network(     phenotype_correlations: np.n` |
| docs/phenotype/visualization.md | 173 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def create_interactive_phenotype_browser(     phenotype_data` |
| docs/phenotype/visualization.md | 187 | AttributeError | 'metainformant' has no attribute 'phenotype' (module not found) | `metainformant.phenotype.visualization.visualization` |
| docs/protein/README.md | 38 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.protein import ...` |
| docs/protein/alphafold.md | 26 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def build_alphafold_url(     uniprot_acc: str,     *,     ve` |
| docs/protein/alphafold.md | 40 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def fetch_alphafold_model(     uniprot_acc: str,     out_dir` |
| docs/protein/alphafold.md | 55 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def batch_download_alphafold_models(     uniprot_accessions:` |
| docs/protein/alphafold.md | 68 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_alphafold_metadata(uniprot_acc: str) -> Dict[str, An` |
| docs/protein/alphafold.md | 77 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def parse_alphafold_confidence(pdb_path: Path) -> List[float` |
| docs/protein/alphafold.md | 86 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_alphafold_structure(pdb_path: Path) -> Dict[str` |
| docs/protein/alphafold.md | 99 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_alphafold_structure_quality(pdb_path: Path) -> Dict[` |
| docs/protein/alphafold.md | 109 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def find_alphafold_models_by_sequence(     sequence: str,   ` |
| docs/protein/alphafold.md | 121 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def search_alphafold_by_keyword(     keyword: str,     max_r` |
| docs/protein/alphafold.md | 132 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_alphafold_coverage() -> Dict[str, Any]` |
| docs/protein/alphafold.md | 141 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein.structure.alphafold` |
| docs/protein/contacts.md | 31 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def calculate_residue_contacts(     coords: np.ndarray,     ` |
| docs/protein/contacts.md | 45 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def identify_hydrogen_bonds(     atoms: List[Dict[str, Any]]` |
| docs/protein/contacts.md | 60 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def identify_salt_bridges(     atoms: List[Dict[str, Any]], ` |
| docs/protein/contacts.md | 74 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def identify_hydrophobic_contacts(     atoms: List[Dict[str,` |
| docs/protein/contacts.md | 87 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def identify_disulfide_bonds(     atoms: List[Dict[str, Any]` |
| docs/protein/contacts.md | 99 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def classify_contact_types(     atoms: List[Dict[str, Any]],` |
| docs/protein/contacts.md | 112 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def analyze_contact_network(     contact_map: np.ndarray, ) ` |
| docs/protein/contacts.md | 126 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def calculate_contact_persistence(     contact_maps: List[np` |
| docs/protein/contacts.md | 137 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def compute_ca_contact_pairs(     coords: list,     threshol` |
| docs/protein/contacts.md | 148 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `>>> coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (10.0, 0.0, ` |
| docs/protein/contacts.md | 156 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein.structure.contacts` |
| docs/protein/index.md | 8 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from pathlib import Path     from metainformant.protein.sequ` |
| docs/protein/index.md | 24 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.protein.uniprot import map_ids_uniprot, f` |
| docs/protein/index.md | 35 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from pathlib import Path     from metainformant.protein.pdb ` |
| docs/protein/index.md | 44 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from pathlib import Path     from metainformant.protein.alph` |
| docs/protein/index.md | 53 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.protein.interpro import fetch_interpro_do` |
| docs/protein/index.md | 61 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `import numpy as np     from pathlib import Path     from met` |
| docs/protein/proteomes.md | 5 | ImportError | Cannot import 'proteomes' from module 'metainformant.protein' | `from metainformant.protein import proteomes` |
| docs/protein/proteomes.md | 5 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/protein/uniprot.md | 25 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def fetch_uniprot_record(uniprot_id: str) -> Dict[str, Any]` |
| docs/protein/uniprot.md | 36 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def fetch_uniprot_fasta(uniprot_id: str) -> Optional[str]` |
| docs/protein/uniprot.md | 44 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def parse_uniprot_fasta_header(header: str) -> Dict[str, str` |
| docs/protein/uniprot.md | 52 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `>>> parse_uniprot_fasta_header("sp|P12345|PROT_HUMAN Protein` |
| docs/protein/uniprot.md | 59 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_uniprot_annotations(uniprot_id: str) -> List[Dict[st` |
| docs/protein/uniprot.md | 67 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def search_uniprot_proteins(     query: str,     max_results` |
| docs/protein/uniprot.md | 80 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_uniprot_taxonomy_info(taxon_id: int) -> Optional[Dic` |
| docs/protein/uniprot.md | 89 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def batch_fetch_uniprot_records(     uniprot_ids: List[str],` |
| docs/protein/uniprot.md | 100 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_uniprot_accession(accession: str) -> bool` |
| docs/protein/uniprot.md | 106 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `>>> validate_uniprot_accession("P12345") True >>> validate_u` |
| docs/protein/uniprot.md | 115 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def map_ids_uniprot(     protein_ids: List[str],     source_` |
| docs/protein/uniprot.md | 130 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein.database.uniprot` |
| docs/quality/README.md | 32 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.quality import ...` |
| docs/quality/contamination.md | 32 | SyntaxError | Invalid Python syntax: expected ':' at line 2 | `class ContaminationDetector:     def __init__(self, referenc` |
| docs/quality/contamination.md | 42 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def detect_microbial_contamination(     self, sequences: Lis` |
| docs/quality/contamination.md | 52 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def detect_cross_species_contamination(     self, sequences:` |
| docs/quality/contamination.md | 63 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def detect_adapter_contamination(     self, sequences: List[` |
| docs/quality/contamination.md | 73 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def detect_duplication_contamination(     self, sequences: L` |
| docs/quality/contamination.md | 82 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def comprehensive_contamination_analysis(     self, sequence` |
| docs/quality/contamination.md | 96 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_rna_contamination(dna_sequences: List[str]) -> Di` |
| docs/quality/contamination.md | 104 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def detect_vector_contamination(     sequences: List[str],  ` |
| docs/quality/contamination.md | 116 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_mycoplasma_contamination(sequences: List[str]) ->` |
| docs/quality/contamination.md | 125 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_rrna_contamination(sequences: List[str]) -> Dict[` |
| docs/quality/contamination.md | 134 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def generate_contamination_report(     contamination_results` |
| docs/quality/contamination.md | 146 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.analysis.contamination` |
| docs/quality/fastq.md | 11 | ImportError | Cannot import 'analyze_fastq_quality' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import analyze_fastq_qualit` |
| docs/quality/fastq.md | 11 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 45 | ImportError | Cannot import 'basic_statistics' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import basic_statistics` |
| docs/quality/fastq.md | 45 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 70 | ImportError | Cannot import 'per_base_quality' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import per_base_quality` |
| docs/quality/fastq.md | 70 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 82 | ImportError | Cannot import 'per_sequence_quality' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import per_sequence_quality` |
| docs/quality/fastq.md | 82 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 94 | ImportError | Cannot import 'sequence_length_distribution' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import sequence_length_dist` |
| docs/quality/fastq.md | 94 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 106 | ImportError | Cannot import 'gc_content_distribution' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import gc_content_distribut` |
| docs/quality/fastq.md | 106 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 118 | ImportError | Cannot import 'adapter_content' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import adapter_content` |
| docs/quality/fastq.md | 118 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 142 | ImportError | Cannot import 'overrepresented_sequences' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import overrepresented_sequ` |
| docs/quality/fastq.md | 142 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 162 | ImportError | Cannot import 'duplication_levels' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import duplication_levels` |
| docs/quality/fastq.md | 162 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 179 | ImportError | Cannot import 'n_content_per_position' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import n_content_per_positi` |
| docs/quality/fastq.md | 179 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 196 | ImportError | Cannot import 'quality_score_distribution' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import quality_score_distri` |
| docs/quality/fastq.md | 196 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/fastq.md | 314 | SyntaxError | Invalid Python syntax: unindent does not match any outer indentation level at line 19 | `def assess_gc_content(results, expected_gc=None):     """Ass` |
| docs/quality/fastq.md | 582 | ImportError | Cannot import 'filter_fastq_quality' from module 'metainformant.dna.fastq' | `from metainformant.dna.fastq import filter_fastq_quality` |
| docs/quality/fastq.md | 582 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.fastq` |
| docs/quality/index.md | 86 | ImportError | Cannot import 'analyze_fastq_quality' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import analyze_fastq_qualit` |
| docs/quality/index.md | 86 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/index.md | 166 | ImportError | Cannot import 'analyze_fastq_quality' from module 'metainformant.quality.fastq' | `from metainformant.quality.fastq import analyze_fastq_qualit` |
| docs/quality/index.md | 166 | ImportError | Cannot import 'process_fastq_file' from module 'metainformant.dna.fastq' | `from metainformant.dna.fastq import process_fastq_file` |
| docs/quality/index.md | 166 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.fastq` |
| docs/quality/index.md | 166 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.fastq` |
| docs/quality/metrics.md | 37 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def calculate_quality_score(     data: Dict[str, Any],     d` |
| docs/quality/metrics.md | 50 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def detect_outliers(     data: List[float],     method: str ` |
| docs/quality/metrics.md | 64 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def calculate_data_integrity_score(     data: Dict[str, Any]` |
| docs/quality/metrics.md | 77 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def compare_quality_metrics(     dataset1: Dict[str, Any],  ` |
| docs/quality/metrics.md | 90 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def generate_quality_report(     quality_data: Dict[str, Any` |
| docs/quality/metrics.md | 104 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def batch_quality_analysis(     file_paths: List[str | Path]` |
| docs/quality/metrics.md | 117 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def calculate_coverage_metrics(     coverage_values: List[fl` |
| docs/quality/metrics.md | 129 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def calculate_duplication_metrics(     duplication_levels: D` |
| docs/quality/metrics.md | 140 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def calculate_gc_metrics(gc_content: List[float]) -> Dict[st` |
| docs/quality/metrics.md | 149 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def calculate_length_metrics(lengths: List[int]) -> Dict[str` |
| docs/quality/metrics.md | 158 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def calculate_quality_metrics(     quality_scores: List[floa` |
| docs/quality/metrics.md | 169 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def calculate_complexity_metrics(     sequences: List[str], ` |
| docs/quality/metrics.md | 181 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality.analysis.metrics` |
| docs/rna/API.md | 29 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def metadata(     params: AmalgkitParams | None = None,     ` |
| docs/rna/API.md | 50 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def integrate(     params: AmalgkitParams | None = None,    ` |
| docs/rna/API.md | 65 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def config(     params: AmalgkitParams | None = None,     **` |
| docs/rna/API.md | 80 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def select(     params: AmalgkitParams | None = None,     **` |
| docs/rna/API.md | 95 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def getfastq(     params: AmalgkitParams | None = None,     ` |
| docs/rna/API.md | 110 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def quant(     params: AmalgkitParams | None = None,     **k` |
| docs/rna/API.md | 125 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def merge(     params: AmalgkitParams | None = None,     **k` |
| docs/rna/API.md | 140 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def cstmm(     params: AmalgkitParams | None = None,     **k` |
| docs/rna/API.md | 155 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def curate(     params: AmalgkitParams | None = None,     **` |
| docs/rna/API.md | 170 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def csca(     params: AmalgkitParams | None = None,     **kw` |
| docs/rna/API.md | 185 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def sanity(     params: AmalgkitParams | None = None,     **` |
| docs/rna/API.md | 204 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_metadata(     params: Mapping[str, Any] | None = Non` |
| docs/rna/API.md | 222 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_integrate(     params: Mapping[str, Any] | None = No` |
| docs/rna/API.md | 238 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_config(     params: Mapping[str, Any] | None = None,` |
| docs/rna/API.md | 254 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_select(     params: Mapping[str, Any] | None = None,` |
| docs/rna/API.md | 270 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_getfastq(     params: Mapping[str, Any] | None = Non` |
| docs/rna/API.md | 288 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_quant(     params: Mapping[str, Any] | None = None, ` |
| docs/rna/API.md | 304 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_merge(     params: Mapping[str, Any] | None = None, ` |
| docs/rna/API.md | 320 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_cstmm(     params: Mapping[str, Any] | None = None, ` |
| docs/rna/API.md | 336 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_curate(     params: Mapping[str, Any] | None = None,` |
| docs/rna/API.md | 352 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_csca(     params: Mapping[str, Any] | None = None,  ` |
| docs/rna/API.md | 368 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def run_sanity(     params: Mapping[str, Any] | None = None,` |
| docs/rna/API.md | 388 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def load_workflow_config(config_file: str | Path) -> Amalgki` |
| docs/rna/API.md | 407 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def plan_workflow(config: AmalgkitWorkflowConfig) -> list[tu` |
| docs/rna/API.md | 426 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def plan_workflow_with_params(     config: AmalgkitWorkflowC` |
| docs/rna/API.md | 441 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def execute_workflow(     config: AmalgkitWorkflowConfig,   ` |
| docs/rna/API.md | 499 | SyntaxError | Invalid Python syntax: expected ':' at line 9 | `def prepare_genome_for_quantification(     genome_dir: Path,` |
| docs/rna/API.md | 523 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def prepare_transcriptome_for_kallisto(     genome_dir: Path` |
| docs/rna/API.md | 543 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def build_kallisto_index(     fasta_path: Path,     index_pa` |
| docs/rna/API.md | 569 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def find_rna_fasta_in_genome_dir(     genome_dir: Path,     ` |
| docs/rna/API.md | 584 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def download_rna_fasta_from_ftp(     ftp_url: str,     genom` |
| docs/rna/API.md | 602 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def download_cds_fasta_from_ftp(     ftp_url: str,     genom` |
| docs/rna/API.md | 619 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def extract_transcripts_from_gff(     gff_path: Path,     ge` |
| docs/rna/API.md | 635 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_expected_index_path(work_dir: Path, species_name: st` |
| docs/rna/API.md | 647 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def verify_genome_status(     genome_dir: Path,     work_dir` |
| docs/rna/API.md | 666 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def orchestrate_genome_setup(     config_dir: Path = Path("c` |
| docs/rna/API.md | 691 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def discover_species_configs(     config_dir: Path = Path("c` |
| docs/rna/API.md | 707 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def run_workflow_for_species(     config_path: Path,     ste` |
| docs/rna/API.md | 731 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def check_workflow_status(     config_path: Path,     *,    ` |
| docs/rna/API.md | 757 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def cleanup_unquantified_samples(     config_path: Path,    ` |
| docs/rna/API.md | 779 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def monitor_workflows(     species_configs: dict[str, Path],` |
| docs/rna/API.md | 804 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_cli_available() -> tuple[bool, str]` |
| docs/rna/API.md | 818 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def ensure_cli_available(     *,     auto_install: bool = Fa` |
| docs/rna/API.md | 835 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def build_cli_args(     params: AmalgkitParams | None,     *` |
| docs/rna/API.md | 857 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def build_amalgkit_command(     subcommand: str,     params:` |
| docs/rna/API.md | 874 | SyntaxError | Invalid Python syntax: expected ':' at line 11 | `def run_amalgkit(     subcommand: str,     params: AmalgkitP` |
| docs/rna/API.md | 912 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def quantify_sample(     sample_id: str,     metadata_rows: ` |
| docs/rna/API.md | 933 | SyntaxError | Invalid Python syntax: expected ':' at line 8 | `def convert_sra_to_fastq(     sample_id: str,     sra_file: ` |
| docs/rna/API.md | 967 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def delete_sample_fastqs(     sample_id: str,     fastq_dir:` |
| docs/rna/API.md | 982 | SyntaxError | Invalid Python syntax: expected ':' at line 12 | `def run_download_quant_workflow(     metadata_path: str | Pa` |
| docs/rna/API.md | 1049 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def count_quantified_samples(config_path: Path) -> tuple[int` |
| docs/rna/API.md | 1063 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_sample_status(config_path: Path, sample_id: str) -> ` |
| docs/rna/API.md | 1082 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def analyze_species_status(config_path: Path) -> dict[str, A` |
| docs/rna/API.md | 1104 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def find_unquantified_samples(config_path: Path) -> list[str` |
| docs/rna/API.md | 1118 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_active_downloads() -> set[str]` |
| docs/rna/API.md | 1132 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_workflow_progress(config_path: Path) -> dict[str, ` |
| docs/rna/API.md | 1154 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def assess_all_species_progress(     config_dir: Path,     *` |
| docs/rna/API.md | 1172 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def initialize_progress_tracking(     config_path: Path,    ` |
| docs/rna/API.md | 1194 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_amalgkit() -> tuple[bool, str]` |
| docs/rna/API.md | 1208 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_sra_toolkit() -> tuple[bool, str]` |
| docs/rna/API.md | 1222 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_kallisto() -> tuple[bool, str]` |
| docs/rna/API.md | 1236 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_metainformant() -> tuple[bool, str]` |
| docs/rna/API.md | 1250 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_virtual_env() -> tuple[bool, str]` |
| docs/rna/API.md | 1264 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_rscript() -> tuple[bool, str]` |
| docs/rna/API.md | 1278 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def check_dependencies() -> dict[str, tuple[bool, str]]` |
| docs/rna/API.md | 1292 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_environment() -> dict[str, Any]` |
| docs/rna/API.md | 1313 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def cleanup_partial_downloads(     config_path: Path,     *,` |
| docs/rna/API.md | 1335 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def fix_abundance_naming(quant_dir: Path, sample_id: str) ->` |
| docs/rna/API.md | 1353 | SyntaxError | Invalid Python syntax: expected ':' at line 3 | `def fix_abundance_naming_for_species(     config_path: Path,` |
| docs/rna/API.md | 1373 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def calculate_translation_efficiency(     rna_df: pd.DataFra` |
| docs/rna/API.md | 1396 | SyntaxError | Invalid Python syntax: expected ':' at line 6 | `def predict_protein_abundance_from_rna(     rna_df: pd.DataF` |
| docs/rna/API.md | 1421 | SyntaxError | Invalid Python syntax: expected ':' at line 4 | `def ribosome_profiling_integration(     rna_df: pd.DataFrame` |
| docs/rna/API.md | 1451 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def validate_all_samples(config: AmalgkitWorkflowConfig) -> ` |
| docs/rna/API.md | 1472 | SyntaxError | Invalid Python syntax: expected ':' at line 5 | `def search_species_with_rnaseq(     search_query: str,     *` |
| docs/rna/API.md | 1496 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def get_genome_info(taxonomy_id: str, species_name: str) -> ` |
| docs/rna/API.md | 1514 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `def generate_config_yaml(     species_name: str,     species` |
| docs/rna/ARCHITECTURE.md | 133 | ImportError | Cannot import 'run_download_quant_workflow' from module 'metainformant.rna.engine.pipeline' | `from metainformant.rna.engine.pipeline import run_download_q` |
| docs/rna/ARCHITECTURE.md | 133 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.pipeline` |
| docs/rna/ARCHITECTURE.md | 166 | ImportError | Cannot import 'step_name' from module 'amalgkit' | `from amalgkit import step_name` |
| docs/rna/ARCHITECTURE.md | 383 | ImportError | Cannot import 'new_step' from module 'amalgkit' | `from amalgkit import new_step` |
| docs/rna/ARCHITECTURE.md | 392 | ImportError | Cannot import 'run' from module 'new_step' | `from new_step import run` |
| docs/rna/CONFIGURATION.md | 41 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.core.configs` |
| docs/rna/CONFIGURATION.md | 52 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/CONFIGURATION.md | 68 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.core.configs` |
| docs/rna/CONFIGURATION.md | 377 | ImportError | Cannot import 'get_recommended_temp_dir' from module 'metainformant.core.disk' | `from metainformant.core.disk import get_recommended_temp_dir` |
| docs/rna/CONFIGURATION.md | 377 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.disk` |
| docs/rna/FILE_PATH_STORAGE.md | 136 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/FILE_PATH_STORAGE.md | 197 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/rna/FILE_PATH_STORAGE.md | 313 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.progress_tracker` |
| docs/rna/FILE_PATH_STORAGE.md | 746 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/ORCHESTRATION.md | 146 | ImportError | Cannot import 'quantify_sample' from module 'metainformant.rna.engine.workflow_steps' | `from metainformant.rna.engine.workflow_steps import quantify` |
| docs/rna/ORCHESTRATION.md | 146 | ImportError | Cannot import 'delete_sample_fastqs' from module 'metainformant.rna.engine.sra_extraction' | `from metainformant.rna.engine.sra_extraction import delete_s` |
| docs/rna/ORCHESTRATION.md | 146 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow_steps` |
| docs/rna/ORCHESTRATION.md | 146 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.sra_extraction` |
| docs/rna/ORCHESTRATION.md | 146 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/ORCHESTRATION.md | 146 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/rna/ORCHESTRATION.md | 175 | ImportError | Cannot import 'cleanup_unquantified_samples' from module 'metainformant.rna.orchestration' | `from metainformant.rna.orchestration import cleanup_unquanti` |
| docs/rna/ORCHESTRATION.md | 175 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.orchestration` |
| docs/rna/README.md | 79 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/README.md | 136 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics.analysis` |
| docs/rna/README.md | 136 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping.colocalization` |
| docs/rna/VALIDATION.md | 42 | ImportError | Cannot import 'validate_all_samples' from module 'metainformant.rna.validation' | `from metainformant.rna.validation import validate_all_sample` |
| docs/rna/VALIDATION.md | 42 | ImportError | Cannot import 'get_sample_pipeline_status' from module 'metainformant.rna.validation' | `from metainformant.rna.validation import get_sample_pipeline` |
| docs/rna/VALIDATION.md | 42 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/VALIDATION.md | 42 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.validation` |
| docs/rna/VALIDATION.md | 304 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/rna/VALIDATION.md | 331 | ImportError | Cannot import 'get_sample_pipeline_status' from module 'metainformant.rna.validation' | `from metainformant.rna.validation import get_sample_pipeline` |
| docs/rna/VALIDATION.md | 331 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.validation` |
| docs/rna/VALIDATION.md | 351 | ImportError | Cannot import 'validate_all_samples' from module 'metainformant.rna.validation' | `from metainformant.rna.validation import validate_all_sample` |
| docs/rna/VALIDATION.md | 351 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/VALIDATION.md | 351 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.validation` |
| docs/rna/amalgkit/R_INSTALLATION.md | 182 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `# Use scikit-learn or pandas for basic QC    import pandas a` |
| docs/rna/amalgkit/amalgkit.md | 51 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.amalgkit` |
| docs/rna/amalgkit/genome_preparation.md | 75 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/genome_preparation.md | 121 | ImportError | Cannot import 'prepare_transcriptome_for_kallisto' from module 'metainformant.rna.genome_prep' | `from metainformant.rna.genome_prep import prepare_transcript` |
| docs/rna/amalgkit/genome_preparation.md | 121 | ImportError | Cannot import 'build_kallisto_index' from module 'metainformant.rna.genome_prep' | `from metainformant.rna.genome_prep import build_kallisto_ind` |
| docs/rna/amalgkit/genome_preparation.md | 121 | ImportError | Cannot import 'prepare_genome_for_quantification' from module 'metainformant.rna.genome_prep' | `from metainformant.rna.genome_prep import prepare_genome_for` |
| docs/rna/amalgkit/genome_preparation.md | 121 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.genome_prep` |
| docs/rna/amalgkit/genome_preparation.md | 411 | SyntaxError | Invalid Python syntax: unexpected indent at line 2 | `from metainformant.rna.genome_prep import build_kallisto_ind` |
| docs/rna/amalgkit/genome_preparation.md | 461 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `prepare_genome_for_quantification(     genome_dir: Path,    ` |
| docs/rna/amalgkit/genome_preparation.md | 483 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `prepare_transcriptome_for_kallisto(     genome_dir: Path,   ` |
| docs/rna/amalgkit/genome_preparation.md | 499 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `build_kallisto_index(     fasta_path: Path,     index_path: ` |
| docs/rna/amalgkit/guide.md | 106 | ImportError | Cannot import 'workflow' from module 'metainformant.rna' | `from metainformant.rna import workflow` |
| docs/rna/amalgkit/guide.md | 106 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/01_metadata.md | 19 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `from metainformant.rna import amalgkit  # High-level wrapper` |
| docs/rna/amalgkit/steps/01_metadata.md | 55 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/01_metadata.md | 423 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/01_metadata.md | 447 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/02_config.md | 29 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/02_config.md | 472 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/03_select.md | 31 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/03_select.md | 151 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `# Minimum quality thresholds min_nspots = 5,000,000         ` |
| docs/rna/amalgkit/steps/03_select.md | 511 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/03_select.md | 539 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/rna/amalgkit/steps/05_integrate.md | 32 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/05_integrate.md | 372 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/06_quant.md | 21 | SyntaxError | Invalid Python syntax: expected ':' at line 7 | `from metainformant.rna import amalgkit  # High-level wrapper` |
| docs/rna/amalgkit/steps/06_quant.md | 59 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 29 | ImportError | Cannot import 'quantify_sample' from module 'metainformant.rna.engine.workflow_steps' | `from metainformant.rna.engine.workflow_steps import quantify` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 29 | ImportError | Cannot import 'delete_sample_fastqs' from module 'metainformant.rna.engine.sra_extraction' | `from metainformant.rna.engine.sra_extraction import delete_s` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 29 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow_steps` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 29 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.sra_extraction` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 52 | ImportError | Cannot import 'cleanup_unquantified_samples' from module 'metainformant.rna.orchestration' | `from metainformant.rna.orchestration import cleanup_unquanti` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 52 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.orchestration` |
| docs/rna/amalgkit/steps/06_quant_advanced.md | 148 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/07_merge.md | 31 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/07_merge.md | 567 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/08_cstmm.md | 30 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/08_cstmm.md | 328 | SyntaxError | Invalid Python syntax: unexpected indent at line 4 | `import pandas as pd        # Add species prefix if needed   ` |
| docs/rna/amalgkit/steps/08_cstmm.md | 480 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/09_curate.md | 39 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/09_curate.md | 562 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/10_csca.md | 31 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/10_csca.md | 495 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/11_sanity.md | 48 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/rna/amalgkit/steps/11_sanity.md | 444 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/amalgkit/steps/11_sanity.md | 461 | SyntaxError | Invalid Python syntax: unindent does not match any outer indentation level at line 16 | `from pathlib import Path  work_dir = Path("output/amalgkit/a` |
| docs/rna/retrieval/ena_downloader.md | 27 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.retrieval` |
| docs/rna/workflow.md | 36 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/workflow.md | 86 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/rna/workflow.md | 122 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.monitoring` |
| docs/rna/workflow.md | 133 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.monitoring` |
| docs/rna/workflow.md | 170 | SyntaxError | Invalid Python syntax: expected an indented block after 'if' statement on line 2 at line 6 | `# Network errors (retry) if "connection" in str(error).lower` |
| docs/rna/workflow.md | 209 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.core.cleanup` |
| docs/rna/workflow.md | 282 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna.engine.workflow` |
| docs/simulation/README.md | 38 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.simulation import ...` |
| docs/simulation/index.md | 50 | ImportError | Cannot import 'generate_random_dna' from module 'metainformant.simulation' | `from metainformant.simulation import generate_random_dna` |
| docs/simulation/index.md | 50 | ImportError | Cannot import 'mutate_sequence' from module 'metainformant.simulation' | `from metainformant.simulation import mutate_sequence` |
| docs/simulation/index.md | 50 | ImportError | Cannot import 'evolve_sequence' from module 'metainformant.simulation' | `from metainformant.simulation import evolve_sequence` |
| docs/simulation/index.md | 50 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/index.md | 65 | ImportError | Cannot import 'simulate_rnaseq_counts' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_rnaseq_counts` |
| docs/simulation/index.md | 65 | ImportError | Cannot import 'simulate_differential_expression' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_differential_e` |
| docs/simulation/index.md | 65 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/index.md | 81 | ImportError | Cannot import 'generate_population_sequences' from module 'metainformant.simulation' | `from metainformant.simulation import generate_population_seq` |
| docs/simulation/index.md | 81 | ImportError | Cannot import 'generate_genotype_matrix' from module 'metainformant.simulation' | `from metainformant.simulation import generate_genotype_matri` |
| docs/simulation/index.md | 81 | ImportError | Cannot import 'simulate_bottleneck_population' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_bottleneck_pop` |
| docs/simulation/index.md | 81 | ImportError | Cannot import 'simulate_selection' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_selection` |
| docs/simulation/index.md | 81 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/index.md | 105 | ImportError | Cannot import 'create_ecosystem' from module 'metainformant.simulation' | `from metainformant.simulation import create_ecosystem` |
| docs/simulation/index.md | 105 | ImportError | Cannot import 'run_simulation' from module 'metainformant.simulation' | `from metainformant.simulation import run_simulation` |
| docs/simulation/index.md | 105 | ImportError | Cannot import 'Agent' from module 'metainformant.simulation' | `from metainformant.simulation import Agent` |
| docs/simulation/index.md | 105 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/index.md | 121 | ImportError | Cannot import 'SimulationConfig' from module 'metainformant.simulation' | `from metainformant.simulation import SimulationConfig` |
| docs/simulation/index.md | 121 | ImportError | Cannot import 'create_simulation_config' from module 'metainformant.simulation' | `from metainformant.simulation import create_simulation_confi` |
| docs/simulation/index.md | 121 | ImportError | Cannot import 'run_simulation_workflow' from module 'metainformant.simulation' | `from metainformant.simulation import run_simulation_workflow` |
| docs/simulation/index.md | 121 | ImportError | Cannot import 'validate_simulation_output' from module 'metainformant.simulation' | `from metainformant.simulation import validate_simulation_out` |
| docs/simulation/index.md | 121 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/index.md | 166 | ImportError | Cannot import 'generate_random_dna' from module 'metainformant.simulation' | `from metainformant.simulation import generate_random_dna` |
| docs/simulation/index.md | 166 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/popgen.md | 36 | ImportError | Cannot import 'generate_population_sequences' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_populat` |
| docs/simulation/popgen.md | 36 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 36 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 66 | ImportError | Cannot import 'generate_two_populations' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_two_pop` |
| docs/simulation/popgen.md | 66 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 66 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 100 | ImportError | Cannot import 'generate_genotype_matrix' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_genotyp` |
| docs/simulation/popgen.md | 100 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 100 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 132 | ImportError | Cannot import 'simulate_bottleneck_population' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import simulate_bottlen` |
| docs/simulation/popgen.md | 132 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 132 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 166 | ImportError | Cannot import 'simulate_population_expansion' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import simulate_populat` |
| docs/simulation/popgen.md | 166 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 166 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 196 | ImportError | Cannot import 'generate_site_frequency_spectrum' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_site_fr` |
| docs/simulation/popgen.md | 196 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 225 | ImportError | Cannot import 'generate_linkage_disequilibrium_data' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_linkage` |
| docs/simulation/popgen.md | 225 | ImportError | Cannot import 'r_squared' from module 'metainformant.math.ld' | `from metainformant.math.ld import r_squared` |
| docs/simulation/popgen.md | 225 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 225 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.ld` |
| docs/simulation/popgen.md | 249 | ImportError | Cannot import 'generate_population_sequences' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_populat` |
| docs/simulation/popgen.md | 249 | ImportError | Cannot import 'calculate_summary_statistics' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import calculate_` |
| docs/simulation/popgen.md | 249 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 249 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/simulation/popgen.md | 271 | ImportError | Cannot import 'simulate_bottleneck_population' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import simulate_bottlen` |
| docs/simulation/popgen.md | 271 | ImportError | Cannot import 'simulate_population_expansion' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import simulate_populat` |
| docs/simulation/popgen.md | 271 | ImportError | Cannot import 'neutrality_test_suite' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import neutrality` |
| docs/simulation/popgen.md | 271 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 271 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/simulation/popgen.md | 306 | ImportError | Cannot import 'generate_genotype_matrix' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_genotyp` |
| docs/simulation/popgen.md | 306 | ImportError | Cannot import 'compute_pca' from module 'metainformant.gwas.structure' | `from metainformant.gwas.structure import compute_pca` |
| docs/simulation/popgen.md | 306 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 306 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.structure` |
| docs/simulation/popgen.md | 327 | ImportError | Cannot import 'generate_two_populations' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_two_pop` |
| docs/simulation/popgen.md | 327 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 327 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/popgen.md | 351 | ImportError | Cannot import 'generate_two_populations' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_two_pop` |
| docs/simulation/popgen.md | 351 | ImportError | Cannot import 'compare_populations' from module 'metainformant.dna.population_analysis' | `from metainformant.dna.population_analysis import compare_po` |
| docs/simulation/popgen.md | 351 | ImportError | Cannot import 'bottleneck_effective_size' from module 'metainformant.math.demography' | `from metainformant.math.demography import bottleneck_effecti` |
| docs/simulation/popgen.md | 351 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 351 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population_analysis` |
| docs/simulation/popgen.md | 351 | AttributeError | 'metainformant' has no attribute 'math' (module not found) | `metainformant.math.demography` |
| docs/simulation/popgen.md | 404 | ImportError | Cannot import 'generate_population_sequences' from module 'metainformant.simulation.popgen' | `from metainformant.simulation.popgen import generate_populat` |
| docs/simulation/popgen.md | 404 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation.popgen` |
| docs/simulation/popgen.md | 404 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna.population` |
| docs/simulation/rna_counts.md | 34 | ImportError | Cannot import 'simulate_counts_negative_binomial' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_counts_negativ` |
| docs/simulation/rna_counts.md | 34 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 69 | ImportError | Cannot import 'simulate_rnaseq_counts' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_rnaseq_counts` |
| docs/simulation/rna_counts.md | 69 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 103 | ImportError | Cannot import 'simulate_differential_expression' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_differential_e` |
| docs/simulation/rna_counts.md | 103 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 137 | ImportError | Cannot import 'simulate_bulk_rnaseq' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_bulk_rnaseq` |
| docs/simulation/rna_counts.md | 137 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 165 | ImportError | Cannot import 'simulate_single_cell_rnaseq' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_single_cell_rn` |
| docs/simulation/rna_counts.md | 165 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 197 | ImportError | Cannot import 'simulate_time_series_expression' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_time_series_ex` |
| docs/simulation/rna_counts.md | 197 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 208 | ImportError | Cannot import 'simulate_spatial_expression' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_spatial_expres` |
| docs/simulation/rna_counts.md | 208 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/rna_counts.md | 228 | ImportError | Cannot import 'simulate_rnaseq_counts' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_rnaseq_counts` |
| docs/simulation/rna_counts.md | 228 | ImportError | Cannot import 'add_technical_noise' from module 'metainformant.simulation' | `from metainformant.simulation import add_technical_noise` |
| docs/simulation/rna_counts.md | 228 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 29 | ImportError | Cannot import 'generate_random_dna' from module 'metainformant.simulation' | `from metainformant.simulation import generate_random_dna` |
| docs/simulation/sequences.md | 29 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 50 | ImportError | Cannot import 'generate_random_protein' from module 'metainformant.simulation' | `from metainformant.simulation import generate_random_protein` |
| docs/simulation/sequences.md | 50 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 61 | ImportError | Cannot import 'generate_coding_sequence' from module 'metainformant.simulation' | `from metainformant.simulation import generate_coding_sequenc` |
| docs/simulation/sequences.md | 61 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 77 | ImportError | Cannot import 'mutate_sequence' from module 'metainformant.simulation' | `from metainformant.simulation import mutate_sequence` |
| docs/simulation/sequences.md | 77 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 98 | ImportError | Cannot import 'evolve_sequence' from module 'metainformant.simulation' | `from metainformant.simulation import evolve_sequence` |
| docs/simulation/sequences.md | 98 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 123 | ImportError | Cannot import 'translate_dna_to_protein' from module 'metainformant.simulation' | `from metainformant.simulation import translate_dna_to_protei` |
| docs/simulation/sequences.md | 123 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 136 | ImportError | Cannot import 'reverse_transcribe_protein_to_dna' from module 'metainformant.simulation' | `from metainformant.simulation import reverse_transcribe_prot` |
| docs/simulation/sequences.md | 136 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 152 | ImportError | Cannot import 'generate_sequence_family' from module 'metainformant.simulation' | `from metainformant.simulation import generate_sequence_famil` |
| docs/simulation/sequences.md | 152 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 166 | ImportError | Cannot import 'analyze_sequence_divergence' from module 'metainformant.simulation' | `from metainformant.simulation import analyze_sequence_diverg` |
| docs/simulation/sequences.md | 166 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 186 | ImportError | Cannot import 'simulate_gene_duplication' from module 'metainformant.simulation' | `from metainformant.simulation import simulate_gene_duplicati` |
| docs/simulation/sequences.md | 186 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/simulation/sequences.md | 202 | ImportError | Cannot import 'calculate_sequence_similarity' from module 'metainformant.simulation' | `from metainformant.simulation import calculate_sequence_simi` |
| docs/simulation/sequences.md | 202 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/singlecell/README.md | 34 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.singlecell import ...` |
| docs/singlecell/celltyping.md | 110 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/celltyping.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/celltyping.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/celltyping.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/celltyping.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/celltyping.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.celltyping.annotation` |
| docs/singlecell/clustering.md | 11 | ImportError | Cannot import 'leiden_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import leiden_clust` |
| docs/singlecell/clustering.md | 11 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 32 | ImportError | Cannot import 'louvain_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import louvain_clus` |
| docs/singlecell/clustering.md | 32 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 54 | ImportError | Cannot import 'kmeans_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import kmeans_clust` |
| docs/singlecell/clustering.md | 54 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 79 | ImportError | Cannot import 'hierarchical_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import hierarchical` |
| docs/singlecell/clustering.md | 79 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 103 | ImportError | Cannot import 'find_markers' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import find_markers` |
| docs/singlecell/clustering.md | 103 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 129 | ImportError | Cannot import 'cluster_composition' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import cluster_comp` |
| docs/singlecell/clustering.md | 129 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 144 | ImportError | Cannot import 'silhouette_scores' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import silhouette_s` |
| docs/singlecell/clustering.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/clustering.md | 219 | ImportError | Cannot import 'compute_neighbors' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/clustering.md | 219 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/clustering.md | 243 | ImportError | Cannot import 'cluster_composition' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import cluster_comp` |
| docs/singlecell/clustering.md | 243 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/config/README.md | 19 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.config import ...` |
| docs/singlecell/differential.md | 97 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/differential.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/differential.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/differential.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/differential.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/differential.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.differential.expression` |
| docs/singlecell/dimensionality.md | 11 | ImportError | Cannot import 'select_hvgs' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import select_h` |
| docs/singlecell/dimensionality.md | 11 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 54 | ImportError | Cannot import 'compute_pca' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/dimensionality.md | 54 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 112 | ImportError | Cannot import 'compute_neighbors' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/dimensionality.md | 112 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 142 | ImportError | Cannot import 'compute_umap' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/dimensionality.md | 142 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 177 | ImportError | Cannot import 'compute_tsne' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/dimensionality.md | 177 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 205 | ImportError | Cannot import 'compute_diffusion_map' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/dimensionality.md | 205 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 321 | ImportError | Cannot import '*' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import *` |
| docs/singlecell/dimensionality.md | 321 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/dimensionality.md | 405 | ImportError | Cannot import 'leiden_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import leiden_clust` |
| docs/singlecell/dimensionality.md | 405 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'load_count_matrix' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import load_coun` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'calculate_qc_metrics' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import calculate` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'select_hvgs' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import select_h` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'compute_pca' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'compute_umap' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import compute_` |
| docs/singlecell/index.md | 111 | ImportError | Cannot import 'leiden_clustering' from module 'metainformant.singlecell.clustering' | `from metainformant.singlecell.clustering import leiden_clust` |
| docs/singlecell/index.md | 111 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/index.md | 111 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/index.md | 111 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.clustering` |
| docs/singlecell/integration.md | 11 | ImportError | Cannot import 'concatenate_datasets' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import concatenate` |
| docs/singlecell/integration.md | 11 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 32 | ImportError | Cannot import 'intersect_datasets' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import intersect_d` |
| docs/singlecell/integration.md | 32 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 51 | ImportError | Cannot import 'union_datasets' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import union_datas` |
| docs/singlecell/integration.md | 51 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 71 | ImportError | Cannot import 'batch_correction_scaling' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import batch_corre` |
| docs/singlecell/integration.md | 71 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 86 | ImportError | Cannot import 'batch_correction_combat' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import batch_corre` |
| docs/singlecell/integration.md | 86 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 105 | ImportError | Cannot import 'batch_correction_harmony' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import batch_corre` |
| docs/singlecell/integration.md | 105 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 125 | ImportError | Cannot import '*' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import *` |
| docs/singlecell/integration.md | 125 | ImportError | Cannot import '*' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import *` |
| docs/singlecell/integration.md | 125 | ImportError | Cannot import '*' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import *` |
| docs/singlecell/integration.md | 125 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/integration.md | 125 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/integration.md | 125 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 189 | ImportError | Cannot import 'integration_metrics' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import integration` |
| docs/singlecell/integration.md | 189 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/integration.md | 210 | ImportError | Cannot import 'plot_integration_assessment' from module 'metainformant.singlecell.integration' | `from metainformant.singlecell.integration import plot_integr` |
| docs/singlecell/integration.md | 210 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.integration` |
| docs/singlecell/preprocessing.md | 28 | ImportError | Cannot import 'SingleCellData' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import SingleCel` |
| docs/singlecell/preprocessing.md | 28 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 47 | ImportError | Cannot import 'load_count_matrix' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import load_coun` |
| docs/singlecell/preprocessing.md | 47 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 78 | ImportError | Cannot import 'calculate_qc_metrics' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import calculate` |
| docs/singlecell/preprocessing.md | 78 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 103 | ImportError | Cannot import 'filter_cells' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import filter_ce` |
| docs/singlecell/preprocessing.md | 103 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 130 | ImportError | Cannot import 'filter_genes' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import filter_ge` |
| docs/singlecell/preprocessing.md | 130 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 147 | ImportError | Cannot import 'normalize_counts' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import normalize` |
| docs/singlecell/preprocessing.md | 147 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 174 | ImportError | Cannot import 'log_transform' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import log_trans` |
| docs/singlecell/preprocessing.md | 174 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 197 | ImportError | Cannot import 'scale_data' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import scale_dat` |
| docs/singlecell/preprocessing.md | 197 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'load_count_matrix' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import load_coun` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'calculate_qc_metrics' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import calculate` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'filter_cells' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import filter_ce` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'filter_genes' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import filter_ge` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'normalize_counts' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import normalize` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'log_transform' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import log_trans` |
| docs/singlecell/preprocessing.md | 221 | ImportError | Cannot import 'scale_data' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import scale_dat` |
| docs/singlecell/preprocessing.md | 221 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/preprocessing.md | 317 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/singlecell/trajectory.md | 11 | ImportError | Cannot import 'compute_pseudotime' from module 'metainformant.singlecell.trajectory' | `from metainformant.singlecell.trajectory import compute_pseu` |
| docs/singlecell/trajectory.md | 11 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.trajectory` |
| docs/singlecell/trajectory.md | 38 | ImportError | Cannot import 'trajectory_analysis' from module 'metainformant.singlecell.trajectory' | `from metainformant.singlecell.trajectory import trajectory_a` |
| docs/singlecell/trajectory.md | 38 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.trajectory` |
| docs/singlecell/trajectory.md | 61 | ImportError | Cannot import 'compute_gene_trends' from module 'metainformant.singlecell.trajectory' | `from metainformant.singlecell.trajectory import compute_gene` |
| docs/singlecell/trajectory.md | 61 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.trajectory` |
| docs/singlecell/trajectory.md | 82 | ImportError | Cannot import 'identify_lineages' from module 'metainformant.singlecell.trajectory' | `from metainformant.singlecell.trajectory import identify_lin` |
| docs/singlecell/trajectory.md | 82 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.trajectory` |
| docs/singlecell/trajectory.md | 102 | ImportError | Cannot import 'find_branch_genes' from module 'metainformant.singlecell.trajectory' | `from metainformant.singlecell.trajectory import find_branch_` |
| docs/singlecell/trajectory.md | 102 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.trajectory` |
| docs/singlecell/trajectory.md | 155 | ImportError | Cannot import '*' from module 'metainformant.singlecell.preprocessing' | `from metainformant.singlecell.preprocessing import *` |
| docs/singlecell/trajectory.md | 155 | ImportError | Cannot import '*' from module 'metainformant.singlecell.dimensionality' | `from metainformant.singlecell.dimensionality import *` |
| docs/singlecell/trajectory.md | 155 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.preprocessing` |
| docs/singlecell/trajectory.md | 155 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.dimensionality` |
| docs/singlecell/velocity.md | 110 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/velocity.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/velocity.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/velocity.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/velocity.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/velocity.md | 144 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.velocity.rna_velocity` |
| docs/singlecell/visualization.md | 11 | ImportError | Cannot import 'plot_qc_metrics' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_qc_m` |
| docs/singlecell/visualization.md | 11 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 27 | ImportError | Cannot import 'plot_qc_scatter' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_qc_s` |
| docs/singlecell/visualization.md | 27 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 47 | ImportError | Cannot import 'plot_embedding' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_embe` |
| docs/singlecell/visualization.md | 47 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 73 | ImportError | Cannot import 'plot_pca' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_pca` |
| docs/singlecell/visualization.md | 73 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 93 | ImportError | Cannot import 'plot_gene_expression' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_gene` |
| docs/singlecell/visualization.md | 93 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 118 | ImportError | Cannot import 'plot_heatmap' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_heat` |
| docs/singlecell/visualization.md | 118 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 138 | ImportError | Cannot import 'plot_clusters' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_clus` |
| docs/singlecell/visualization.md | 138 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 156 | ImportError | Cannot import 'plot_cluster_composition' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_clus` |
| docs/singlecell/visualization.md | 156 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 181 | ImportError | Cannot import 'plot_trajectory' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_traj` |
| docs/singlecell/visualization.md | 181 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 200 | ImportError | Cannot import 'plot_gene_trends' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_gene` |
| docs/singlecell/visualization.md | 200 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 220 | ImportError | Cannot import 'plot_comparison' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_comp` |
| docs/singlecell/visualization.md | 220 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 237 | ImportError | Cannot import 'plot_split' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_spli` |
| docs/singlecell/visualization.md | 237 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 257 | ImportError | Cannot import 'plot_marker_genes' from module 'metainformant.singlecell.visualization' | `from metainformant.singlecell.visualization import plot_mark` |
| docs/singlecell/visualization.md | 257 | AttributeError | 'metainformant' has no attribute 'singlecell' (module not found) | `metainformant.singlecell.visualization` |
| docs/singlecell/visualization.md | 388 | ImportError | Cannot import 'setup_matplotlib_style' from module 'metainformant.visualization.plots' | `from metainformant.visualization.plots import setup_matplotl` |
| docs/singlecell/visualization.md | 388 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/singlecell/visualization.md | 388 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/spatial/CAPABILITIES.md | 846 | SyntaxError | Invalid Python syntax: unexpected indent at line 1 | ` (gene name). If None, uses the         built-in database fr` |
| docs/spatial/CAPABILITIES.md | 892 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Build cell-cell communication network from interaction res` |
| docs/spatial/CAPABILITIES.md | 1052 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Build cell-type reference profiles from single-cell data. ` |
| docs/spatial/CAPABILITIES.md | 1148 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 8) at line 8 | `  Identify tissue niches from cell type composition and spat` |
| docs/spatial/CAPABILITIES.md | 1164 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `: Array of niche assignments (length n_spots),           int` |
| docs/spatial/CAPABILITIES.md | 1274 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 5) at line 5 | `  Anchor-based label transfer from scRNA-seq to spatial data` |
| docs/spatial/CAPABILITIES.md | 1306 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Map scRNA-seq cell type annotations to spatial spots.  Uni` |
| docs/spatial/CAPABILITIES.md | 1400 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Load individual transcript spot coordinates from a MERFISH` |
| docs/spatial/CAPABILITIES.md | 1418 | SyntaxError | Invalid Python syntax: unterminated string literal (detected at line 6) at line 6 | `  Aggregate individual transcript spots to cell-level expres` |
| docs/spatial/CAPABILITIES.md | 1440 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Load a complete MERFISH dataset from a directory.  Expects` |
| docs/spatial/CAPABILITIES.md | 1463 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Find the first matching column name from a list of candida` |
| docs/spatial/CAPABILITIES.md | 1623 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Create a unified SpatialDataset from components.  Args:   ` |
| docs/spatial/CAPABILITIES.md | 1719 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Read per-transcript coordinates from a Xenium transcripts ` |
| docs/spatial/CAPABILITIES.md | 1753 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Load cell boundaries from a parquet file using pandas.  ##` |
| docs/spatial/CAPABILITIES.md | 1759 | SyntaxError | Invalid Python syntax: invalid decimal literal at line 3 | `  Load a complete 10x Xenium dataset from a directory.  Expe` |
| docs/spatial/CAPABILITIES.md | 1784 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Find first matching column name from candidates.  ### Clas` |
| docs/spatial/CAPABILITIES.md | 1839 | SyntaxError | Invalid Python syntax: invalid decimal literal at line 11 | `  Identify tissue niches from local cell type composition.  ` |
| docs/spatial/CAPABILITIES.md | 2041 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Plot pie charts per spatial spot showing cell type fractio` |
| docs/spatial/GETTING_STARTED.md | 17 | ImportError | Cannot import 'load_xenium' from module 'metainformant.spatial' | `from metainformant.spatial import load_xenium` |
| docs/spatial/GETTING_STARTED.md | 17 | ModuleNotFoundError | Module 'spatial.analysis' not found in project or standard library | `import spatial.analysis` |
| docs/spatial/GETTING_STARTED.md | 17 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'load_visium' from module 'metainformant.spatial' | `from metainformant.spatial import load_visium` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'load_merfish' from module 'metainformant.spatial' | `from metainformant.spatial import load_merfish` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'load_xenium' from module 'metainformant.spatial' | `from metainformant.spatial import load_xenium` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/spatial/INTEGRATION.md | 17 | ImportError | Cannot import 'quickplot' from module 'metainformant.visualization' | `from metainformant.visualization import quickplot` |
| docs/spatial/INTEGRATION.md | 17 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/spatial/INTEGRATION.md | 17 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/spatial/INTEGRATION.md | 17 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/spatial/INTEGRATION.md | 55 | ImportError | Cannot import 'submit_batch' from module 'metainformant.cloud' | `from metainformant.cloud import submit_batch` |
| docs/spatial/INTEGRATION.md | 55 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/spatial/INTEGRATION.md | 65 | ImportError | Cannot import 'analyze' from module 'metainformant.spatial' | `from metainformant.spatial import analyze` |
| docs/spatial/INTEGRATION.md | 65 | ImportError | Cannot import 'plot_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import plot_heatmap` |
| docs/spatial/INTEGRATION.md | 65 | ImportError | Cannot import 'plot_timeseries' from module 'metainformant.visualization' | `from metainformant.visualization import plot_timeseries` |
| docs/spatial/INTEGRATION.md | 65 | ImportError | Cannot import 'plot_network' from module 'metainformant.visualization' | `from metainformant.visualization import plot_network` |
| docs/spatial/INTEGRATION.md | 65 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/spatial/INTEGRATION.md | 65 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/spatial/INTEGRATION.md | 84 | SyntaxError | Invalid Python syntax: invalid syntax at line 5 | `class Result:     values: dict            # Primary output d` |
| docs/spatial/TROUBLESHOOTING.md | 17 | ImportError | Cannot import 'load_visium' from module 'metainformant.spatial' | `from metainformant.spatial import load_visium` |
| docs/spatial/TROUBLESHOOTING.md | 17 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial` |
| docs/spatial/TROUBLESHOOTING.md | 62 | ImportError | Cannot import 'benchmark' from module 'metainformant.spatial.performance' | `from metainformant.spatial.performance import benchmark` |
| docs/spatial/TROUBLESHOOTING.md | 62 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.performance` |
| docs/spatial/autocorrelation.md | 9 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.autocorrelation` |
| docs/spatial/autocorrelation.md | 31 | ImportError | Cannot import 'local_moran' from module 'metainformant.spatial.analysis.autocorrelation' | `from metainformant.spatial.analysis.autocorrelation import l` |
| docs/spatial/autocorrelation.md | 31 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.autocorrelation` |
| docs/spatial/autocorrelation.md | 46 | ImportError | Cannot import 'plot_spatial_categorical' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import plot_spatial` |
| docs/spatial/autocorrelation.md | 46 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/autocorrelation.md | 55 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.autocorrelation` |
| docs/spatial/autocorrelation.md | 67 | ImportError | Cannot import 'getis_ord' from module 'metainformant.spatial.analysis.autocorrelation' | `from metainformant.spatial.analysis.autocorrelation import g` |
| docs/spatial/autocorrelation.md | 67 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.autocorrelation` |
| docs/spatial/clustering.md | 18 | ImportError | Cannot import 'spatial_leiden' from module 'metainformant.spatial.analysis.clustering' | `from metainformant.spatial.analysis.clustering import spatia` |
| docs/spatial/clustering.md | 18 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.clustering` |
| docs/spatial/clustering.md | 48 | ImportError | Cannot import 'find_markers' from module 'metainformant.spatial.analysis.clustering' | `from metainformant.spatial.analysis.clustering import find_m` |
| docs/spatial/clustering.md | 48 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.clustering` |
| docs/spatial/clustering.md | 60 | ImportError | Cannot import 'spatial_scatter' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import spatial_scat` |
| docs/spatial/clustering.md | 60 | ImportError | Cannot import 'domain_outlines' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import domain_outli` |
| docs/spatial/clustering.md | 60 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/clustering.md | 60 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/clustering.md | 73 | ImportError | Cannot import 'merge_domains' from module 'metainformant.spatial.analysis.clustering' | `from metainformant.spatial.analysis.clustering import merge_` |
| docs/spatial/clustering.md | 73 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.clustering` |
| docs/spatial/clustering.md | 86 | ImportError | Cannot import 'spatial_kmeans' from module 'metainformant.spatial.analysis.clustering' | `from metainformant.spatial.analysis.clustering import spatia` |
| docs/spatial/clustering.md | 86 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.clustering` |
| docs/spatial/communication.md | 8 | ImportError | Cannot import 'ligand_receptor_score' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import lig` |
| docs/spatial/communication.md | 8 | ImportError | Cannot import 'spatial_ligand_receptor_pairs' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import spa` |
| docs/spatial/communication.md | 8 | ImportError | Cannot import 'permutation_significance' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import per` |
| docs/spatial/communication.md | 8 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.communication` |
| docs/spatial/communication.md | 43 | ImportError | Cannot import 'load_lr_database' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import loa` |
| docs/spatial/communication.md | 43 | ImportError | Cannot import 'spatial_ligand_receptor_pairs' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import spa` |
| docs/spatial/communication.md | 43 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.communication` |
| docs/spatial/communication.md | 64 | ImportError | Cannot import 'lr_network_graph' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import lr_network_g` |
| docs/spatial/communication.md | 64 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/communication.md | 74 | ImportError | Cannot import 'register_custom_lr' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import reg` |
| docs/spatial/communication.md | 74 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.communication` |
| docs/spatial/deconvolution.md | 17 | ImportError | Cannot import 'stereoscope' from module 'metainformant.spatial.analysis.deconvolution' | `from metainformant.spatial.analysis.deconvolution import ste` |
| docs/spatial/deconvolution.md | 17 | ImportError | Cannot import 'train_stereoscope' from module 'metainformant.spatial.analysis.deconvolution' | `from metainformant.spatial.analysis.deconvolution import tra` |
| docs/spatial/deconvolution.md | 17 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.deconvolution` |
| docs/spatial/deconvolution.md | 37 | ImportError | Cannot import 'tangram' from module 'metainformant.spatial.analysis.deconvolution' | `from metainformant.spatial.analysis.deconvolution import tan` |
| docs/spatial/deconvolution.md | 37 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.deconvolution` |
| docs/spatial/integration.md | 7 | ImportError | Cannot import 'spatial_batch_correct' from module 'metainformant.spatial.integration' | `from metainformant.spatial.integration import spatial_batch_` |
| docs/spatial/integration.md | 7 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.integration` |
| docs/spatial/integration.md | 32 | ImportError | Cannot import 'transfer_labels' from module 'metainformant.spatial.integration' | `from metainformant.spatial.integration import transfer_label` |
| docs/spatial/integration.md | 32 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.integration` |
| docs/spatial/integration.md | 48 | ImportError | Cannot import 'register_coordinates' from module 'metainformant.spatial.integration' | `from metainformant.spatial.integration import register_coord` |
| docs/spatial/integration.md | 48 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.integration` |
| docs/spatial/integration.md | 64 | ImportError | Cannot import 'downsample_to_visium' from module 'metainformant.spatial.integration' | `from metainformant.spatial.integration import downsample_to_` |
| docs/spatial/integration.md | 64 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.integration` |
| docs/spatial/io.md | 21 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io.visium` |
| docs/spatial/io.md | 59 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io.xenium` |
| docs/spatial/io.md | 87 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io.merfish` |
| docs/spatial/io.md | 102 | ImportError | Cannot import 'load_slide_seq' from module 'metainformant.spatial.io.slideseq' | `from metainformant.spatial.io.slideseq import load_slide_seq` |
| docs/spatial/io.md | 102 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io.slideseq` |
| docs/spatial/io.md | 117 | ImportError | Cannot import 'load_cosmx' from module 'metainformant.spatial.io.cosmx' | `from metainformant.spatial.io.cosmx import load_cosmx` |
| docs/spatial/io.md | 117 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.io.cosmx` |
| docs/spatial/io.md | 136 | ImportError | Cannot import 'load_myplatform' from module 'myplatform' | `from myplatform import load_myplatform` |
| docs/spatial/neighborhood.md | 12 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.neighborhood` |
| docs/spatial/neighborhood.md | 29 | ImportError | Cannot import 'interaction_index' from module 'metainformant.spatial.analysis.neighborhood' | `from metainformant.spatial.analysis.neighborhood import inte` |
| docs/spatial/neighborhood.md | 29 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.neighborhood` |
| docs/spatial/neighborhood.md | 47 | ImportError | Cannot import 'ligand_receptor_score' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import lig` |
| docs/spatial/neighborhood.md | 47 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.communication` |
| docs/spatial/neighborhood.md | 65 | ImportError | Cannot import 'nearest_neighbor_enrichment' from module 'metainformant.spatial.analysis.neighborhood' | `from metainformant.spatial.analysis.neighborhood import near` |
| docs/spatial/neighborhood.md | 65 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.neighborhood` |
| docs/spatial/neighborhood.md | 80 | ImportError | Cannot import 'stereoscope' from module 'metainformant.spatial.analysis.deconvolution' | `from metainformant.spatial.analysis.deconvolution import ste` |
| docs/spatial/neighborhood.md | 80 | ImportError | Cannot import 'interaction_matrix' from module 'metainformant.spatial.analysis.neighborhood' | `from metainformant.spatial.analysis.neighborhood import inte` |
| docs/spatial/neighborhood.md | 80 | ImportError | Cannot import 'spatial_ligand_receptor_pairs' from module 'metainformant.spatial.analysis.communication' | `from metainformant.spatial.analysis.communication import spa` |
| docs/spatial/neighborhood.md | 80 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.deconvolution` |
| docs/spatial/neighborhood.md | 80 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.neighborhood` |
| docs/spatial/niche.md | 21 | ImportError | Cannot import 'discover_niches' from module 'metainformant.spatial.analysis.niche' | `from metainformant.spatial.analysis.niche import discover_ni` |
| docs/spatial/niche.md | 21 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.niche` |
| docs/spatial/niche.md | 35 | ImportError | Cannot import 'plot_niche_summary' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import plot_niche_s` |
| docs/spatial/niche.md | 35 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/niche.md | 44 | ImportError | Cannot import 'annotate_niche' from module 'metainformant.spatial.analysis.niche' | `from metainformant.spatial.analysis.niche import annotate_ni` |
| docs/spatial/niche.md | 44 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.niche` |
| docs/spatial/niche.md | 60 | ImportError | Cannot import 'find_markers' from module 'metainformant.spatial.analysis.clustering' | `from metainformant.spatial.analysis.clustering import find_m` |
| docs/spatial/niche.md | 60 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.analysis.clustering` |
| docs/spatial/visualization.md | 8 | ImportError | Cannot import 'spatial_scatter' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import spatial_scat` |
| docs/spatial/visualization.md | 8 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/visualization.md | 24 | ImportError | Cannot import 'domain_outlines' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import domain_outli` |
| docs/spatial/visualization.md | 24 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/visualization.md | 40 | ImportError | Cannot import 'expression_heatmap' from module 'metainformant.spatial.visualization' | `from metainformant.spatial.visualization import expression_h` |
| docs/spatial/visualization.md | 40 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization` |
| docs/spatial/visualization.md | 60 | ImportError | Cannot import 'spatial_scatter_plotly' from module 'metainformant.spatial.visualization.interactive' | `from metainformant.spatial.visualization.interactive import ` |
| docs/spatial/visualization.md | 60 | AttributeError | 'metainformant' has no attribute 'spatial' (module not found) | `metainformant.spatial.visualization.interactive` |
| docs/structural_variants/ARCHITECTURE.md | 373 | SyntaxError | Invalid Python syntax: unexpected indent at line 1 | ` where type_weight(DEL)=1.0, TRA=0.95, INV=0.9, DUP=0.7, INS` |
| docs/structural_variants/ARCHITECTURE.md | 568 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CAPABILITIES.md | 1650 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Merge SV calls across samples using reciprocal overlap and` |
| docs/structural_variants/CAPABILITIES.md | 1761 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Build a consensus SV call from a cluster of merged calls. ` |
| docs/structural_variants/CAPABILITIES.md | 1779 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Raise an import error if matplotlib is not available.  ###` |
| docs/structural_variants/CAPABILITIES.md | 1811 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `  Plot read depth coverage with structural variant overlay. ` |
| docs/structural_variants/CONFIGURATION.md | 21 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 183 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 219 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 236 | ImportError | Cannot import 'configure' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import configure` |
| docs/structural_variants/CONFIGURATION.md | 236 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/CONFIGURATION.md | 321 | ImportError | Cannot import 'paths' from module 'metainformant.core' | `from metainformant.core import paths` |
| docs/structural_variants/CONFIGURATION.md | 321 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/CONFIGURATION.md | 331 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `config.set("core.data_dir", "/my/data") # or export METAINFO` |
| docs/structural_variants/CONFIGURATION.md | 376 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/CONFIGURATION.md | 400 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 400 | ImportError | Cannot import 'logging' from module 'metainformant' | `from metainformant import logging` |
| docs/structural_variants/CONFIGURATION.md | 400 | ImportError | Cannot import 'pipeline' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import pipeline` |
| docs/structural_variants/CONFIGURATION.md | 400 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/CONFIGURATION.md | 438 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 438 | ImportError | Cannot import 'get_default_config' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import get_default_co` |
| docs/structural_variants/CONFIGURATION.md | 438 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/CONFIGURATION.md | 448 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/CONFIGURATION.md | 463 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/CONFIGURATION.md | 473 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 38 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/structural_variants/EXAMPLES.md | 38 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 38 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/EXAMPLES.md | 106 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/structural_variants/EXAMPLES.md | 106 | ImportError | Cannot import 'write_vcf' from module 'metainformant.structural_variants.io' | `from metainformant.structural_variants.io import write_vcf` |
| docs/structural_variants/EXAMPLES.md | 106 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.sv_calling` |
| docs/structural_variants/EXAMPLES.md | 106 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/EXAMPLES.md | 106 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/EXAMPLES.md | 106 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.io` |
| docs/structural_variants/EXAMPLES.md | 219 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 318 | ImportError | Cannot import 'write_vcf_from_merged' from module 'metainformant.structural_variants.io' | `from metainformant.structural_variants.io import write_vcf_f` |
| docs/structural_variants/EXAMPLES.md | 318 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/EXAMPLES.md | 318 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/EXAMPLES.md | 318 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.io` |
| docs/structural_variants/EXAMPLES.md | 377 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/EXAMPLES.md | 398 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/EXAMPLES.md | 398 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/EXAMPLES.md | 398 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/EXAMPLES.md | 476 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/EXAMPLES.md | 476 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/EXAMPLES.md | 536 | ImportError | Cannot import 'plot_size_distribution' from module 'metainformant.structural_variants.visualization.plots' | `from metainformant.structural_variants.visualization.plots i` |
| docs/structural_variants/EXAMPLES.md | 536 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.visualization.plots` |
| docs/structural_variants/EXAMPLES.md | 536 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 536 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 610 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/EXAMPLES.md | 610 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/structural_variants/EXAMPLES.md | 706 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/EXAMPLES.md | 706 | AttributeError | 'metainformant' has no attribute 'longread' (module not found) | `metainformant.longread` |
| docs/structural_variants/EXAMPLES.md | 751 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/structural_variants/EXAMPLES.md | 751 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/EXAMPLES.md | 751 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/EXAMPLES.md | 751 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 751 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.sv_calling` |
| docs/structural_variants/EXAMPLES.md | 751 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/EXAMPLES.md | 751 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/EXAMPLES.md | 836 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/structural_variants/EXAMPLES.md | 836 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/structural_variants/EXAMPLES.md | 836 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/EXAMPLES.md | 836 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/EXAMPLES.md | 873 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/EXAMPLES.md | 873 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/EXAMPLES.md | 873 | AttributeError | 'metainformant' has no attribute 'metabolomics' (module not found) | `metainformant.metabolomics.pathways.enrichment` |
| docs/structural_variants/EXAMPLES.md | 944 | ImportError | Cannot import 'differential_splicing' from module 'metainformant.rna' | `from metainformant.rna import differential_splicing` |
| docs/structural_variants/EXAMPLES.md | 944 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/structural_variants/EXAMPLES.md | 944 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/EXAMPLES.md | 975 | ImportError | Cannot import 'write_vcf' from module 'metainformant.structural_variants.io' | `from metainformant.structural_variants.io import write_vcf` |
| docs/structural_variants/EXAMPLES.md | 975 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.io` |
| docs/structural_variants/GETTING_STARTED.md | 30 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/GETTING_STARTED.md | 30 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/GETTING_STARTED.md | 86 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.sv_calling` |
| docs/structural_variants/GETTING_STARTED.md | 143 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/GETTING_STARTED.md | 185 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/GETTING_STARTED.md | 185 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/GETTING_STARTED.md | 223 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/GETTING_STARTED.md | 268 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/GETTING_STARTED.md | 322 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/GETTING_STARTED.md | 365 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/GETTING_STARTED.md | 400 | SyntaxError | Invalid Python syntax: ':' expected after dictionary key at line 7 | `from metainformant.structural_variants.visualization.plots i` |
| docs/structural_variants/GETTING_STARTED.md | 438 | ImportError | Cannot import 'run_full_pipeline' from module 'metainformant.structural_variants.pipeline' | `from metainformant.structural_variants.pipeline import run_f` |
| docs/structural_variants/GETTING_STARTED.md | 438 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.pipeline` |
| docs/structural_variants/GETTING_STARTED.md | 506 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'detect_cnvs' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import detect_cnvs` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'call_svs' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import call_svs` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'annotate_svs' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import annotate_svs` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'config' from module 'metainformant.core' | `from metainformant.core import config` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'logging' from module 'metainformant.core' | `from metainformant.core import logging` |
| docs/structural_variants/INTEGRATION.md | 18 | ImportError | Cannot import 'quickplot' from module 'metainformant.visualization' | `from metainformant.visualization import quickplot` |
| docs/structural_variants/INTEGRATION.md | 18 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/INTEGRATION.md | 18 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/INTEGRATION.md | 18 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/structural_variants/INTEGRATION.md | 56 | ImportError | Cannot import 'submit_batch' from module 'metainformant.cloud' | `from metainformant.cloud import submit_batch` |
| docs/structural_variants/INTEGRATION.md | 56 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud` |
| docs/structural_variants/INTEGRATION.md | 66 | ImportError | Cannot import 'analyze' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import analyze` |
| docs/structural_variants/INTEGRATION.md | 66 | ImportError | Cannot import 'plot_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import plot_heatmap` |
| docs/structural_variants/INTEGRATION.md | 66 | ImportError | Cannot import 'plot_timeseries' from module 'metainformant.visualization' | `from metainformant.visualization import plot_timeseries` |
| docs/structural_variants/INTEGRATION.md | 66 | ImportError | Cannot import 'plot_network' from module 'metainformant.visualization' | `from metainformant.visualization import plot_network` |
| docs/structural_variants/INTEGRATION.md | 66 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/INTEGRATION.md | 66 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/structural_variants/INTEGRATION.md | 85 | SyntaxError | Invalid Python syntax: invalid syntax at line 5 | `class Result:     values: dict            # Primary output d` |
| docs/structural_variants/PERFORMANCE.md | 51 | ImportError | Cannot import 'config' from module 'metainformant' | `from metainformant import config` |
| docs/structural_variants/PERFORMANCE.md | 62 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/structural_variants/PERFORMANCE.md | 101 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/PERFORMANCE.md | 114 | ImportError | Cannot import 'cache' from module 'metainformant.core' | `from metainformant.core import cache` |
| docs/structural_variants/PERFORMANCE.md | 114 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core` |
| docs/structural_variants/PERFORMANCE.md | 166 | ImportError | Cannot import 'pipeline' from module 'metainformant.structural_variants' | `from metainformant.structural_variants import pipeline` |
| docs/structural_variants/PERFORMANCE.md | 166 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants` |
| docs/structural_variants/PERFORMANCE.md | 215 | ImportError | Cannot import 'profile' from module 'memory_profiler' | `from memory_profiler import profile` |
| docs/structural_variants/PERFORMANCE.md | 235 | ImportError | Cannot import 'Client' from module 'dask.distributed' | `from dask.distributed import Client` |
| docs/structural_variants/PERFORMANCE.md | 235 | ImportError | Cannot import 'parallel' from module 'metainformant' | `from metainformant import parallel` |
| docs/structural_variants/PERFORMANCE.md | 356 | ImportError | Cannot import 'setup_logging' from module 'metainformant.core.utils.logging' | `from metainformant.core.utils.logging import setup_logging` |
| docs/structural_variants/PERFORMANCE.md | 356 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.utils.logging` |
| docs/structural_variants/README.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection` |
| docs/structural_variants/README.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation` |
| docs/structural_variants/README.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering` |
| docs/structural_variants/README.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population` |
| docs/structural_variants/README.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.visualization` |
| docs/structural_variants/SPEC.md | 5 | SyntaxError | Invalid Python syntax: expected ':' at line 1 | `def detect_svs(bam: str, reference: str, callers: list[str],` |
| docs/structural_variants/TROUBLESHOOTING.md | 8 | SyntaxError | Invalid Python syntax: invalid character '×' (U+00D7) at line 1 | `Zero SVs detected from 30× WGS` |
| docs/structural_variants/annotation.md | 52 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/annotation.md | 73 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.overlap` |
| docs/structural_variants/annotation.md | 133 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/annotation.md | 156 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/annotation.md | 177 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/annotation.md | 195 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation.functional_impa` |
| docs/structural_variants/detection.md | 42 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/detection.md | 83 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.sv_calling` |
| docs/structural_variants/detection.md | 125 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/detection.md | 151 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.cnv` |
| docs/structural_variants/detection.md | 151 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.sv_calling` |
| docs/structural_variants/detection.md | 151 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection.breakpoints` |
| docs/structural_variants/filtering.md | 37 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/filtering.md | 54 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/filtering.md | 68 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/filtering.md | 84 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/filtering.md | 128 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/filtering.md | 151 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/filtering.md | 172 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.quality_filter` |
| docs/structural_variants/filtering.md | 172 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering.merge` |
| docs/structural_variants/index.md | 25 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.detection` |
| docs/structural_variants/index.md | 25 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.annotation` |
| docs/structural_variants/index.md | 25 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.filtering` |
| docs/structural_variants/population.md | 24 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 50 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 69 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 97 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 115 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 135 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/structural_variants/population.md | 149 | AttributeError | 'metainformant' has no attribute 'structural_variants' (module not found) | `metainformant.structural_variants.population.sv_population` |
| docs/tasks/VALIDATION_REPORT.md | 44 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/tasks/VALIDATION_REPORT.md | 44 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 49 | ImportError | Cannot import 'composition' from module 'metainformant.dna' | `from metainformant.dna import composition` |
| docs/tasks/VALIDATION_REPORT.md | 49 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 52 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 57 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 61 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/tasks/VALIDATION_REPORT.md | 61 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 65 | ImportError | Cannot import 'effects' from module 'metainformant.dna' | `from metainformant.dna import effects` |
| docs/tasks/VALIDATION_REPORT.md | 65 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 67 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/VALIDATION_REPORT.md | 94 | ImportError | Cannot import 'convert' from module 'metainformant.core.io' | `from metainformant.core.io import convert` |
| docs/tasks/VALIDATION_REPORT.md | 94 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/VALIDATION_REPORT.md | 100 | ImportError | Cannot import 'bam_metrics' from module 'metainformant.core.io' | `from metainformant.core.io import bam_metrics` |
| docs/tasks/VALIDATION_REPORT.md | 100 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/VALIDATION_REPORT.md | 102 | ImportError | Cannot import 'merge' from module 'metainformant.variants' | `from metainformant.variants import merge` |
| docs/tasks/VALIDATION_REPORT.md | 102 | AttributeError | 'metainformant' has no attribute 'variants' (module not found) | `metainformant.variants` |
| docs/tasks/VALIDATION_REPORT.md | 104 | ImportError | Cannot import 'gff' from module 'metainformant.core.io' | `from metainformant.core.io import gff` |
| docs/tasks/VALIDATION_REPORT.md | 104 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/VALIDATION_REPORT.md | 106 | ImportError | Cannot import 'fastq' from module 'metainformant.core.io' | `from metainformant.core.io import fastq` |
| docs/tasks/VALIDATION_REPORT.md | 106 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/VALIDATION_REPORT.md | 144 | ImportError | Cannot import 'register_preempt_handler' from module 'metainformant.cloud.interrupt' | `from metainformant.cloud.interrupt import register_preempt_h` |
| docs/tasks/VALIDATION_REPORT.md | 144 | AttributeError | 'metainformant' has no attribute 'cloud' (module not found) | `metainformant.cloud.interrupt` |
| docs/tasks/VALIDATION_REPORT.md | 200 | SyntaxError | Invalid Python syntax: invalid syntax at line 3 | `from metainformant.mcp import register_tool, ToolSpec @regis` |
| docs/tasks/VALIDATION_REPORT.md | 230 | SyntaxError | Invalid Python syntax: invalid syntax at line 2 | `from metainformant.core.caching import memoize_disk @memoize` |
| docs/tasks/VALIDATION_REPORT.md | 238 | ImportError | Cannot import 'parallel_map' from module 'metainformant.core.parallel' | `from metainformant.core.parallel import parallel_map` |
| docs/tasks/VALIDATION_REPORT.md | 238 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.parallel` |
| docs/tasks/VALIDATION_REPORT.md | 271 | SyntaxError | Invalid Python syntax: invalid syntax. Perhaps you forgot a comma? at line 3 | `from joblib import Parallel, delayed from metainformant.dna ` |
| docs/tasks/VALIDATION_REPORT.md | 280 | ImportError | Cannot import 'read_csv_chunks' from module 'metainformant.core.io' | `from metainformant.core.io import read_csv_chunks` |
| docs/tasks/VALIDATION_REPORT.md | 280 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/VALIDATION_REPORT.md | 287 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `ImportError: cannot import name 'apply_qc_filters' from 'met` |
| docs/tasks/VALIDATION_REPORT.md | 305 | ImportError | Cannot import 'run_association' from module 'metainformant.gwas' | `from metainformant.gwas import run_association` |
| docs/tasks/VALIDATION_REPORT.md | 305 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/VALIDATION_REPORT.md | 333 | ImportError | Cannot import 'regenie' from module 'metainformant.gwas' | `from metainformant.gwas import regenie` |
| docs/tasks/VALIDATION_REPORT.md | 333 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/VALIDATION_REPORT.md | 342 | ImportError | Cannot import 'susie' from module 'metainformant.gwas.finemapping' | `from metainformant.gwas.finemapping import susie` |
| docs/tasks/VALIDATION_REPORT.md | 342 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping` |
| docs/tasks/VALIDATION_REPORT.md | 356 | ImportError | Cannot import 'coloc' from module 'metainformant.gwas.coloc' | `from metainformant.gwas.coloc import coloc` |
| docs/tasks/VALIDATION_REPORT.md | 356 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.coloc` |
| docs/tasks/VALIDATION_REPORT.md | 450 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/tasks/VALIDATION_REPORT.md | 484 | ImportError | Cannot import 'deseq2' from module 'metainformant.rna' | `from metainformant.rna import deseq2` |
| docs/tasks/VALIDATION_REPORT.md | 484 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/tasks/VALIDATION_REPORT.md | 492 | ImportError | Cannot import 'orthologs' from module 'metainformant.rna' | `from metainformant.rna import orthologs` |
| docs/tasks/VALIDATION_REPORT.md | 492 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/tasks/VALIDATION_REPORT.md | 507 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/VALIDATION_REPORT.md | 521 | ImportError | Cannot import 'pca_plots' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plots` |
| docs/tasks/VALIDATION_REPORT.md | 521 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/VALIDATION_REPORT.md | 531 | SyntaxError | Invalid Python syntax: positional argument follows keyword argument at line 2 | `from metainformant.visualization import heatmaps fig = heatm` |
| docs/tasks/VALIDATION_REPORT.md | 539 | SyntaxError | Invalid Python syntax: positional argument follows keyword argument at line 3 | `from metainformant.dna import phylogeny tree = phylogeny.rea` |
| docs/tasks/VALIDATION_REPORT.md | 548 | SyntaxError | Invalid Python syntax: positional argument follows keyword argument at line 3 | `from metainformant.networks import graph, visualize ppi = gr` |
| docs/tasks/VALIDATION_REPORT.md | 562 | ImportError | Cannot import 'circos' from module 'metainformant.visualization' | `from metainformant.visualization import circos` |
| docs/tasks/VALIDATION_REPORT.md | 562 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/VALIDATION_REPORT.md | 570 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dashboards` |
| docs/tasks/VALIDATION_REPORT.md | 571 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.plots` |
| docs/tasks/VALIDATION_REPORT.md | 572 | ImportError | Cannot import 'plot_tree' from module 'metainformant.visualization.genomics.trees' | `from metainformant.visualization.genomics.trees import plot_` |
| docs/tasks/VALIDATION_REPORT.md | 572 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.genomics.trees` |
| docs/tasks/VALIDATION_REPORT.md | 573 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks.analysis` |
| docs/tasks/VALIDATION_REPORT.md | 579 | SyntaxError | Invalid Python syntax: positional argument follows keyword argument at line 2 | `from metainformant.visualization import animation fig = anim` |
| docs/tasks/VALIDATION_REPORT.md | 587 | ImportError | Cannot import 'apply_qc_filters' from module 'metainformant.gwas' | `from metainformant.gwas import apply_qc_filters` |
| docs/tasks/VALIDATION_REPORT.md | 587 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/VALIDATION_REPORT.md | 587 | ImportError | Cannot import 'apply_qc_filters' from module 'metainformant.gwas' | `from metainformant.gwas import apply_qc_filters` |
| docs/tasks/VALIDATION_REPORT.md | 587 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/VALIDATION_REPORT.md | 587 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.analysis.quality` |
| docs/tasks/VALIDATION_REPORT.md | 588 | ImportError | Cannot import 'visualize' from module 'metainformant.protein' | `from metainformant.protein import visualize` |
| docs/tasks/VALIDATION_REPORT.md | 588 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/tasks/VALIDATION_REPORT.md | 588 | ImportError | Cannot import 'association_test_linear' from module 'metainformant.gwas.association' | `from metainformant.gwas.association import association_test_` |
| docs/tasks/VALIDATION_REPORT.md | 588 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.association` |
| docs/tasks/VALIDATION_REPORT.md | 588 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/analyze_dna.md | 23 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/tasks/analyze_dna.md | 23 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 41 | ImportError | Cannot import 'composition' from module 'metainformant.dna' | `from metainformant.dna import composition` |
| docs/tasks/analyze_dna.md | 41 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 54 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 72 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 88 | ImportError | Cannot import 'variants' from module 'metainformant.dna' | `from metainformant.dna import variants` |
| docs/tasks/analyze_dna.md | 88 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 113 | ImportError | Cannot import 'composition' from module 'metainformant.dna' | `from metainformant.dna import composition` |
| docs/tasks/analyze_dna.md | 113 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 133 | ImportError | Cannot import 'effects' from module 'metainformant.dna' | `from metainformant.dna import effects` |
| docs/tasks/analyze_dna.md | 133 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/analyze_dna.md | 151 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/data_conversion.md | 19 | ImportError | Cannot import 'convert' from module 'metainformant.core.io' | `from metainformant.core.io import convert` |
| docs/tasks/data_conversion.md | 19 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/data_conversion.md | 51 | ImportError | Cannot import 'convert' from module 'metainformant.core.io' | `from metainformant.core.io import convert` |
| docs/tasks/data_conversion.md | 51 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/data_conversion.md | 69 | ImportError | Cannot import 'bam_metrics' from module 'metainformant.core.io' | `from metainformant.core.io import bam_metrics` |
| docs/tasks/data_conversion.md | 69 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/data_conversion.md | 96 | ImportError | Cannot import 'merge' from module 'metainformant.variants' | `from metainformant.variants import merge` |
| docs/tasks/data_conversion.md | 96 | AttributeError | 'metainformant' has no attribute 'variants' (module not found) | `metainformant.variants` |
| docs/tasks/data_conversion.md | 118 | ImportError | Cannot import 'gff' from module 'metainformant.core.io' | `from metainformant.core.io import gff` |
| docs/tasks/data_conversion.md | 118 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/data_conversion.md | 140 | ImportError | Cannot import 'fastq' from module 'metainformant.core.io' | `from metainformant.core.io import fastq` |
| docs/tasks/data_conversion.md | 140 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/deploy_cloud.md | 135 | SyntaxError | Invalid Python syntax: invalid syntax at line 12 | `# In your pipeline script, register checkpoint handlers from` |
| docs/tasks/mcp_integration.md | 151 | ImportError | Cannot import 'register_tool' from module 'metainformant.mcp' | `from metainformant.mcp import register_tool` |
| docs/tasks/mcp_integration.md | 151 | ImportError | Cannot import 'ToolSpec' from module 'metainformant.mcp' | `from metainformant.mcp import ToolSpec` |
| docs/tasks/mcp_integration.md | 151 | ImportError | Cannot import 'qc' from module 'metainformant.quality' | `from metainformant.quality import qc` |
| docs/tasks/mcp_integration.md | 151 | AttributeError | 'metainformant' has no attribute 'mcp' (module not found) | `metainformant.mcp` |
| docs/tasks/mcp_integration.md | 151 | AttributeError | 'metainformant' has no attribute 'quality' (module not found) | `metainformant.quality` |
| docs/tasks/performance_tuning.md | 21 | ImportError | Cannot import 'memoize_disk' from module 'metainformant.core.caching' | `from metainformant.core.caching import memoize_disk` |
| docs/tasks/performance_tuning.md | 21 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.caching` |
| docs/tasks/performance_tuning.md | 32 | ImportError | Cannot import 'parallel_map' from module 'metainformant.core.parallel' | `from metainformant.core.parallel import parallel_map` |
| docs/tasks/performance_tuning.md | 32 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.parallel` |
| docs/tasks/performance_tuning.md | 135 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/performance_tuning.md | 150 | ImportError | Cannot import 'read_csv_chunks' from module 'metainformant.core.io' | `from metainformant.core.io import read_csv_chunks` |
| docs/tasks/performance_tuning.md | 150 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/tasks/run_gwas.md | 22 | ImportError | Cannot import 'run_association' from module 'metainformant.gwas' | `from metainformant.gwas import run_association` |
| docs/tasks/run_gwas.md | 22 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/run_gwas.md | 116 | ImportError | Cannot import 'regenie' from module 'metainformant.gwas' | `from metainformant.gwas import regenie` |
| docs/tasks/run_gwas.md | 116 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/tasks/run_gwas.md | 144 | ImportError | Cannot import 'susie' from module 'metainformant.gwas.finemapping' | `from metainformant.gwas.finemapping import susie` |
| docs/tasks/run_gwas.md | 144 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas.finemapping` |
| docs/tasks/run_gwas.md | 167 | SyntaxError | Invalid Python syntax: unexpected indent at line 11 | `from metainformant.gwas.coloc import coloc  # GWAS summary s` |
| docs/tasks/run_rna_pipeline.md | 122 | ImportError | Cannot import 'deseq2' from module 'metainformant.rna' | `from metainformant.rna import deseq2` |
| docs/tasks/run_rna_pipeline.md | 122 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/tasks/run_rna_pipeline.md | 144 | ImportError | Cannot import 'orthologs' from module 'metainformant.rna' | `from metainformant.rna import orthologs` |
| docs/tasks/run_rna_pipeline.md | 144 | AttributeError | 'metainformant' has no attribute 'rna' (module not found) | `metainformant.rna` |
| docs/tasks/visualize_results.md | 22 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/visualize_results.md | 39 | ImportError | Cannot import 'pca_plots' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plots` |
| docs/tasks/visualize_results.md | 39 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/visualize_results.md | 50 | ImportError | Cannot import 'heatmaps' from module 'metainformant.visualization' | `from metainformant.visualization import heatmaps` |
| docs/tasks/visualize_results.md | 50 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/visualize_results.md | 71 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/tasks/visualize_results.md | 86 | ImportError | Cannot import 'graph' from module 'metainformant.networks' | `from metainformant.networks import graph` |
| docs/tasks/visualize_results.md | 86 | ImportError | Cannot import 'visualize' from module 'metainformant.networks' | `from metainformant.networks import visualize` |
| docs/tasks/visualize_results.md | 86 | AttributeError | 'metainformant' has no attribute 'networks' (module not found) | `metainformant.networks` |
| docs/tasks/visualize_results.md | 102 | ImportError | Cannot import 'circos' from module 'metainformant.visualization' | `from metainformant.visualization import circos` |
| docs/tasks/visualize_results.md | 102 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/visualize_results.md | 118 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.dashboards` |
| docs/tasks/visualize_results.md | 135 | ImportError | Cannot import 'animation' from module 'metainformant.visualization' | `from metainformant.visualization import animation` |
| docs/tasks/visualize_results.md | 135 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/tasks/visualize_results.md | 152 | ImportError | Cannot import 'visualize' from module 'metainformant.protein' | `from metainformant.protein import visualize` |
| docs/tasks/visualize_results.md | 152 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/testing.md | 308 | ImportError | Cannot import 'some_function' from module 'metainformant.some_module' | `from metainformant.some_module import some_function` |
| docs/testing.md | 308 | AttributeError | 'metainformant' has no attribute 'core' (module not found) | `metainformant.core.io` |
| docs/testing.md | 308 | AttributeError | 'metainformant' has no attribute 'some_module' (module not found) | `metainformant.some_module` |
| docs/testing.md | 357 | ImportError | Cannot import 'map_ids_uniprot' from module 'metainformant.protein.uniprot' | `from metainformant.protein.uniprot import map_ids_uniprot` |
| docs/testing.md | 357 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein.uniprot` |
| docs/visualization/README.md | 49 | SyntaxError | Invalid Python syntax: invalid syntax at line 1 | `from metainformant.visualization import ...` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'manhattan' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'qq' from module 'metainformant.visualization' | `from metainformant.visualization import qq` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'regional' from module 'metainformant.visualization' | `from metainformant.visualization import regional` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'ideogram' from module 'metainformant.visualization' | `from metainformant.visualization import ideogram` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'pca' from module 'metainformant.visualization' | `from metainformant.visualization import pca` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'volcano' from module 'metainformant.visualization' | `from metainformant.visualization import volcano` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'ma' from module 'metainformant.visualization' | `from metainformant.visualization import ma` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'spring_layout' from module 'metainformant.visualization' | `from metainformant.visualization import spring_layout` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'circular' from module 'metainformant.visualization' | `from metainformant.visualization import circular` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'hive' from module 'metainformant.visualization' | `from metainformant.visualization import hive` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'plot_tree' from module 'metainformant.visualization' | `from metainformant.visualization import plot_tree` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'plot_circular' from module 'metainformant.visualization' | `from metainformant.visualization import plot_circular` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'tanglegram' from module 'metainformant.visualization' | `from metainformant.visualization import tanglegram` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'trajectory' from module 'metainformant.visualization' | `from metainformant.visualization import trajectory` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'evolution' from module 'metainformant.visualization' | `from metainformant.visualization import evolution` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'Dashboard' from module 'metainformant.visualization' | `from metainformant.visualization import Dashboard` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'Explorer' from module 'metainformant.visualization' | `from metainformant.visualization import Explorer` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'set_style' from module 'metainformant.visualization' | `from metainformant.visualization import set_style` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'set_palette' from module 'metainformant.visualization' | `from metainformant.visualization import set_palette` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'save_figure' from module 'metainformant.visualization' | `from metainformant.visualization import save_figure` |
| docs/visualization/SPEC.md | 103 | ImportError | Cannot import 'to_caption' from module 'metainformant.visualization' | `from metainformant.visualization import to_caption` |
| docs/visualization/SPEC.md | 103 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/SPEC.md | 124 | ImportError | Cannot import 'manhattan' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan` |
| docs/visualization/SPEC.md | 124 | ImportError | Cannot import 'pca' from module 'metainformant.visualization' | `from metainformant.visualization import pca` |
| docs/visualization/SPEC.md | 124 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/SPEC.md | 124 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/animations.md | 12 | ImportError | Cannot import 'animate_time_series' from module 'metainformant.visualization' | `from metainformant.visualization import animate_time_series` |
| docs/visualization/animations.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 20 | ImportError | Cannot import 'lineplot' from module 'metainformant.visualization' | `from metainformant.visualization import lineplot` |
| docs/visualization/basic.md | 20 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 33 | ImportError | Cannot import 'scatter_plot' from module 'metainformant.visualization' | `from metainformant.visualization import scatter_plot` |
| docs/visualization/basic.md | 33 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 44 | ImportError | Cannot import 'bar_plot' from module 'metainformant.visualization' | `from metainformant.visualization import bar_plot` |
| docs/visualization/basic.md | 44 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 55 | ImportError | Cannot import 'pie_chart' from module 'metainformant.visualization' | `from metainformant.visualization import pie_chart` |
| docs/visualization/basic.md | 55 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 66 | ImportError | Cannot import 'area_plot' from module 'metainformant.visualization' | `from metainformant.visualization import area_plot` |
| docs/visualization/basic.md | 66 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 77 | ImportError | Cannot import 'step_plot' from module 'metainformant.visualization' | `from metainformant.visualization import step_plot` |
| docs/visualization/basic.md | 77 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/basic.md | 88 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/basic.md | 88 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/dimred.md | 12 | ImportError | Cannot import 'pca_plot' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plot` |
| docs/visualization/dimred.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/epigenome.md | 14 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/visualization/epigenome.md | 151 | ImportError | Cannot import 'methylation' from module 'metainformant.epigenome' | `from metainformant.epigenome import methylation` |
| docs/visualization/epigenome.md | 151 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/visualization/epigenome.md | 166 | ImportError | Cannot import 'chipseq' from module 'metainformant.epigenome' | `from metainformant.epigenome import chipseq` |
| docs/visualization/epigenome.md | 166 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/visualization/epigenome.md | 177 | ImportError | Cannot import 'atacseq' from module 'metainformant.epigenome' | `from metainformant.epigenome import atacseq` |
| docs/visualization/epigenome.md | 177 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/visualization/epigenome.md | 224 | ImportError | Cannot import 'methylation' from module 'metainformant.epigenome' | `from metainformant.epigenome import methylation` |
| docs/visualization/epigenome.md | 224 | ImportError | Cannot import 'chipseq' from module 'metainformant.epigenome' | `from metainformant.epigenome import chipseq` |
| docs/visualization/epigenome.md | 224 | AttributeError | 'metainformant' has no attribute 'epigenome' (module not found) | `metainformant.epigenome` |
| docs/visualization/examples.md | 7 | ImportError | Cannot import 'lineplot' from module 'metainformant.visualization' | `from metainformant.visualization import lineplot` |
| docs/visualization/examples.md | 7 | ImportError | Cannot import 'scatter_plot' from module 'metainformant.visualization' | `from metainformant.visualization import scatter_plot` |
| docs/visualization/examples.md | 7 | ImportError | Cannot import 'bar_plot' from module 'metainformant.visualization' | `from metainformant.visualization import bar_plot` |
| docs/visualization/examples.md | 7 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/examples.md | 7 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/examples.md | 31 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization' | `from metainformant.visualization import histogram` |
| docs/visualization/examples.md | 31 | ImportError | Cannot import 'box_plot' from module 'metainformant.visualization' | `from metainformant.visualization import box_plot` |
| docs/visualization/examples.md | 31 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization' | `from metainformant.visualization import qq_plot` |
| docs/visualization/examples.md | 31 | ImportError | Cannot import 'correlation_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import correlation_heatmap` |
| docs/visualization/examples.md | 31 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/examples.md | 55 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan_plot` |
| docs/visualization/examples.md | 55 | ImportError | Cannot import 'volcano_plot' from module 'metainformant.visualization' | `from metainformant.visualization import volcano_plot` |
| docs/visualization/examples.md | 55 | ImportError | Cannot import 'regional_plot' from module 'metainformant.visualization' | `from metainformant.visualization import regional_plot` |
| docs/visualization/examples.md | 55 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/examples.md | 80 | ImportError | Cannot import 'pca_plot' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plot` |
| docs/visualization/examples.md | 80 | ImportError | Cannot import 'umap_plot' from module 'metainformant.visualization' | `from metainformant.visualization import umap_plot` |
| docs/visualization/examples.md | 80 | ImportError | Cannot import 'pca_scree_plot' from module 'metainformant.visualization' | `from metainformant.visualization import pca_scree_plot` |
| docs/visualization/examples.md | 80 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/examples.md | 100 | ImportError | Cannot import 'create_multi_panel' from module 'metainformant.visualization.layout' | `from metainformant.visualization.layout import create_multi_` |
| docs/visualization/examples.md | 100 | ImportError | Cannot import 'add_shared_axis_labels' from module 'metainformant.visualization.layout' | `from metainformant.visualization.layout import add_shared_ax` |
| docs/visualization/examples.md | 100 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.layout` |
| docs/visualization/examples.md | 112 | ImportError | Cannot import 'save_figure_multiformat' from module 'metainformant.visualization.export' | `from metainformant.visualization.export import save_figure_m` |
| docs/visualization/examples.md | 112 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.export` |
| docs/visualization/expression.md | 12 | ImportError | Cannot import 'expression_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import expression_heatmap` |
| docs/visualization/expression.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 12 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan_plot` |
| docs/visualization/genomics.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 31 | ImportError | Cannot import 'volcano_plot' from module 'metainformant.visualization' | `from metainformant.visualization import volcano_plot` |
| docs/visualization/genomics.md | 31 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 49 | ImportError | Cannot import 'regional_plot' from module 'metainformant.visualization' | `from metainformant.visualization import regional_plot` |
| docs/visualization/genomics.md | 49 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 66 | ImportError | Cannot import 'circular_manhattan_plot' from module 'metainformant.visualization' | `from metainformant.visualization import circular_manhattan_p` |
| docs/visualization/genomics.md | 66 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 85 | ImportError | Cannot import 'chromosome_ideogram' from module 'metainformant.visualization' | `from metainformant.visualization import chromosome_ideogram` |
| docs/visualization/genomics.md | 85 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 101 | ImportError | Cannot import 'coverage_plot' from module 'metainformant.visualization' | `from metainformant.visualization import coverage_plot` |
| docs/visualization/genomics.md | 101 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/genomics.md | 115 | ImportError | Cannot import 'variant_plot' from module 'metainformant.visualization' | `from metainformant.visualization import variant_plot` |
| docs/visualization/genomics.md | 115 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'lineplot' from module 'metainformant.visualization' | `from metainformant.visualization import lineplot` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'scatter_plot' from module 'metainformant.visualization' | `from metainformant.visualization import scatter_plot` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization' | `from metainformant.visualization import histogram` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan_plot` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'volcano_plot' from module 'metainformant.visualization' | `from metainformant.visualization import volcano_plot` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'pca_plot' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plot` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'plot_phylo_tree' from module 'metainformant.visualization' | `from metainformant.visualization import plot_phylo_tree` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'animate_time_series' from module 'metainformant.visualization' | `from metainformant.visualization import animate_time_series` |
| docs/visualization/index.md | 39 | ImportError | Cannot import 'Phylo' from module 'Bio' | `from Bio import Phylo` |
| docs/visualization/index.md | 39 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/information.md | 12 | ImportError | Cannot import 'entropy_plot' from module 'metainformant.visualization' | `from metainformant.visualization import entropy_plot` |
| docs/visualization/information.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import man` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'circular_manhattan_plot' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import cir` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'qq_plot_stratified' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import qq_` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'regional_plot_detailed' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import reg` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'pca_plot_gwas' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import pca` |
| docs/visualization/integration.md | 9 | ImportError | Cannot import 'GWAS_VISUALIZATION_AVAILABLE' from module 'metainformant.visualization.gwas_integration' | `from metainformant.visualization.gwas_integration import GWA` |
| docs/visualization/integration.md | 9 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.gwas_integration` |
| docs/visualization/integration.md | 24 | ImportError | Cannot import 'plot_qc_metrics' from module 'metainformant.visualization.singlecell_integration' | `from metainformant.visualization.singlecell_integration impo` |
| docs/visualization/integration.md | 24 | ImportError | Cannot import 'plot_embedding' from module 'metainformant.visualization.singlecell_integration' | `from metainformant.visualization.singlecell_integration impo` |
| docs/visualization/integration.md | 24 | ImportError | Cannot import 'plot_gene_expression' from module 'metainformant.visualization.singlecell_integration' | `from metainformant.visualization.singlecell_integration impo` |
| docs/visualization/integration.md | 24 | ImportError | Cannot import 'SINGLECELL_VISUALIZATION_AVAILABLE' from module 'metainformant.visualization.singlecell_integration' | `from metainformant.visualization.singlecell_integration impo` |
| docs/visualization/integration.md | 24 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.singlecell_integration` |
| docs/visualization/integration.md | 37 | ImportError | Cannot import 'plot_entropy_distribution' from module 'metainformant.visualization.information_integration' | `from metainformant.visualization.information_integration imp` |
| docs/visualization/integration.md | 37 | ImportError | Cannot import 'plot_mutual_information_matrix' from module 'metainformant.visualization.information_integration' | `from metainformant.visualization.information_integration imp` |
| docs/visualization/integration.md | 37 | ImportError | Cannot import 'plot_information_profile' from module 'metainformant.visualization.information_integration' | `from metainformant.visualization.information_integration imp` |
| docs/visualization/integration.md | 37 | ImportError | Cannot import 'INFORMATION_VISUALIZATION_AVAILABLE' from module 'metainformant.visualization.information_integration' | `from metainformant.visualization.information_integration imp` |
| docs/visualization/integration.md | 37 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.information_integration` |
| docs/visualization/integration.md | 50 | ImportError | Cannot import 'plot_event_timeline' from module 'metainformant.visualization.life_events_integration' | `from metainformant.visualization.life_events_integration imp` |
| docs/visualization/integration.md | 50 | ImportError | Cannot import 'plot_event_embeddings' from module 'metainformant.visualization.life_events_integration' | `from metainformant.visualization.life_events_integration imp` |
| docs/visualization/integration.md | 50 | ImportError | Cannot import 'plot_attention_heatmap' from module 'metainformant.visualization.life_events_integration' | `from metainformant.visualization.life_events_integration imp` |
| docs/visualization/integration.md | 50 | ImportError | Cannot import 'LIFE_EVENTS_VISUALIZATION_AVAILABLE' from module 'metainformant.visualization.life_events_integration' | `from metainformant.visualization.life_events_integration imp` |
| docs/visualization/integration.md | 50 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.life_events_integration` |
| docs/visualization/multidim.md | 12 | ImportError | Cannot import 'pairplot_dataframe' from module 'metainformant.visualization' | `from metainformant.visualization import pairplot_dataframe` |
| docs/visualization/multidim.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/networks.md | 12 | ImportError | Cannot import 'network_plot' from module 'metainformant.visualization' | `from metainformant.visualization import network_plot` |
| docs/visualization/networks.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/ontology.md | 14 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/ontology.md | 14 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/ontology.md | 14 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/ontology.md | 108 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/ontology.md | 108 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/ontology.md | 108 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/ontology.md | 124 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/ontology.md | 124 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/ontology.md | 124 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.query` |
| docs/visualization/ontology.md | 124 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/ontology.md | 146 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/ontology.md | 146 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/ontology.md | 146 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/ontology.md | 212 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/ontology.md | 212 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/ontology.md | 212 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.query` |
| docs/visualization/ontology.md | 212 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/plots.md | 28 | ImportError | Cannot import 'lineplot' from module 'metainformant.visualization' | `from metainformant.visualization import lineplot` |
| docs/visualization/plots.md | 28 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 46 | ImportError | Cannot import 'scatter_plot' from module 'metainformant.visualization' | `from metainformant.visualization import scatter_plot` |
| docs/visualization/plots.md | 46 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 54 | ImportError | Cannot import 'heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import heatmap` |
| docs/visualization/plots.md | 54 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 65 | ImportError | Cannot import 'bar_plot' from module 'metainformant.visualization' | `from metainformant.visualization import bar_plot` |
| docs/visualization/plots.md | 65 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 75 | ImportError | Cannot import 'pie_chart' from module 'metainformant.visualization' | `from metainformant.visualization import pie_chart` |
| docs/visualization/plots.md | 75 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 86 | ImportError | Cannot import 'area_plot' from module 'metainformant.visualization' | `from metainformant.visualization import area_plot` |
| docs/visualization/plots.md | 86 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 96 | ImportError | Cannot import 'step_plot' from module 'metainformant.visualization' | `from metainformant.visualization import step_plot` |
| docs/visualization/plots.md | 96 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 110 | ImportError | Cannot import 'volcano_plot' from module 'metainformant.visualization' | `from metainformant.visualization import volcano_plot` |
| docs/visualization/plots.md | 110 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 125 | ImportError | Cannot import 'manhattan_plot' from module 'metainformant.visualization' | `from metainformant.visualization import manhattan_plot` |
| docs/visualization/plots.md | 125 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 141 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization' | `from metainformant.visualization import qq_plot` |
| docs/visualization/plots.md | 141 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 153 | ImportError | Cannot import 'expression_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import expression_heatmap` |
| docs/visualization/plots.md | 153 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 168 | ImportError | Cannot import 'pca_plot' from module 'metainformant.visualization' | `from metainformant.visualization import pca_plot` |
| docs/visualization/plots.md | 168 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 185 | ImportError | Cannot import 'correlation_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import correlation_heatmap` |
| docs/visualization/plots.md | 185 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 206 | ImportError | Cannot import 'plot_venn_diagram' from module 'metainformant.visualization' | `from metainformant.visualization import plot_venn_diagram` |
| docs/visualization/plots.md | 206 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 227 | ImportError | Cannot import 'plot_pairwise_relationships' from module 'metainformant.visualization' | `from metainformant.visualization import plot_pairwise_relati` |
| docs/visualization/plots.md | 227 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/plots.md | 249 | ImportError | Cannot import 'animate_time_series' from module 'metainformant.visualization' | `from metainformant.visualization import animate_time_series` |
| docs/visualization/plots.md | 249 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/protein.md | 14 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/visualization/protein.md | 131 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/visualization/protein.md | 143 | ImportError | Cannot import 'sequences' from module 'metainformant.protein' | `from metainformant.protein import sequences` |
| docs/visualization/protein.md | 143 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/visualization/protein.md | 155 | ImportError | Cannot import 'interpro' from module 'metainformant.protein' | `from metainformant.protein import interpro` |
| docs/visualization/protein.md | 155 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/visualization/protein.md | 201 | ImportError | Cannot import 'sequences' from module 'metainformant.protein' | `from metainformant.protein import sequences` |
| docs/visualization/protein.md | 201 | AttributeError | 'metainformant' has no attribute 'protein' (module not found) | `metainformant.protein` |
| docs/visualization/quality.md | 12 | ImportError | Cannot import 'qc_metrics_plot' from module 'metainformant.visualization' | `from metainformant.visualization import qc_metrics_plot` |
| docs/visualization/quality.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/simulation.md | 14 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/visualization/simulation.md | 163 | ImportError | Cannot import 'sequences' from module 'metainformant.simulation' | `from metainformant.simulation import sequences` |
| docs/visualization/simulation.md | 163 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/visualization/simulation.md | 176 | ImportError | Cannot import 'popgen' from module 'metainformant.simulation' | `from metainformant.simulation import popgen` |
| docs/visualization/simulation.md | 176 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/visualization/simulation.md | 188 | ImportError | Cannot import 'agents' from module 'metainformant.simulation' | `from metainformant.simulation import agents` |
| docs/visualization/simulation.md | 188 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/visualization/simulation.md | 239 | ImportError | Cannot import 'sequences' from module 'metainformant.simulation' | `from metainformant.simulation import sequences` |
| docs/visualization/simulation.md | 239 | ImportError | Cannot import 'popgen' from module 'metainformant.simulation' | `from metainformant.simulation import popgen` |
| docs/visualization/simulation.md | 239 | AttributeError | 'metainformant' has no attribute 'simulation' (module not found) | `metainformant.simulation` |
| docs/visualization/specialized.md | 14 | ImportError | Cannot import 'specialized' from module 'metainformant.visualization' | `from metainformant.visualization import specialized` |
| docs/visualization/specialized.md | 14 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/specialized.md | 128 | ImportError | Cannot import 'integration' from module 'metainformant.multiomics' | `from metainformant.multiomics import integration` |
| docs/visualization/specialized.md | 128 | ImportError | Cannot import 'specialized' from module 'metainformant.visualization' | `from metainformant.visualization import specialized` |
| docs/visualization/specialized.md | 128 | AttributeError | 'metainformant' has no attribute 'multiomics' (module not found) | `metainformant.multiomics` |
| docs/visualization/specialized.md | 128 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/specialized.md | 152 | ImportError | Cannot import 'specialized' from module 'metainformant.visualization' | `from metainformant.visualization import specialized` |
| docs/visualization/specialized.md | 152 | AttributeError | 'metainformant' has no attribute 'gwas' (module not found) | `metainformant.gwas` |
| docs/visualization/specialized.md | 152 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/specialized.md | 171 | ImportError | Cannot import 'visualization' from module 'metainformant.ontology.visualization' | `from metainformant.ontology.visualization import visualizati` |
| docs/visualization/specialized.md | 171 | ImportError | Cannot import 'specialized' from module 'metainformant.visualization' | `from metainformant.visualization import specialized` |
| docs/visualization/specialized.md | 171 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.core` |
| docs/visualization/specialized.md | 171 | AttributeError | 'metainformant' has no attribute 'ontology' (module not found) | `metainformant.ontology.visualization` |
| docs/visualization/specialized.md | 171 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/specialized.md | 240 | ImportError | Cannot import 'specialized' from module 'metainformant.visualization' | `from metainformant.visualization import specialized` |
| docs/visualization/specialized.md | 240 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 12 | ImportError | Cannot import 'histogram' from module 'metainformant.visualization' | `from metainformant.visualization import histogram` |
| docs/visualization/statistical.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 25 | ImportError | Cannot import 'box_plot' from module 'metainformant.visualization' | `from metainformant.visualization import box_plot` |
| docs/visualization/statistical.md | 25 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 38 | ImportError | Cannot import 'violin_plot' from module 'metainformant.visualization' | `from metainformant.visualization import violin_plot` |
| docs/visualization/statistical.md | 38 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 51 | ImportError | Cannot import 'qq_plot' from module 'metainformant.visualization' | `from metainformant.visualization import qq_plot` |
| docs/visualization/statistical.md | 51 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 64 | ImportError | Cannot import 'correlation_heatmap' from module 'metainformant.visualization' | `from metainformant.visualization import correlation_heatmap` |
| docs/visualization/statistical.md | 64 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 78 | ImportError | Cannot import 'density_plot' from module 'metainformant.visualization' | `from metainformant.visualization import density_plot` |
| docs/visualization/statistical.md | 78 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 91 | ImportError | Cannot import 'ridge_plot' from module 'metainformant.visualization' | `from metainformant.visualization import ridge_plot` |
| docs/visualization/statistical.md | 91 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 104 | ImportError | Cannot import 'roc_curve' from module 'metainformant.visualization' | `from metainformant.visualization import roc_curve` |
| docs/visualization/statistical.md | 104 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 117 | ImportError | Cannot import 'precision_recall_curve' from module 'metainformant.visualization' | `from metainformant.visualization import precision_recall_cur` |
| docs/visualization/statistical.md | 117 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 130 | ImportError | Cannot import 'residual_plot' from module 'metainformant.visualization' | `from metainformant.visualization import residual_plot` |
| docs/visualization/statistical.md | 130 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/statistical.md | 144 | ImportError | Cannot import 'leverage_plot' from module 'metainformant.visualization' | `from metainformant.visualization import leverage_plot` |
| docs/visualization/statistical.md | 144 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/styling.md | 12 | SyntaxError | Invalid Python syntax: expression cannot contain assignment, perhaps you meant "=="? at line 3 | `from metainformant.visualization.style import apply_publicat` |
| docs/visualization/styling.md | 31 | ImportError | Cannot import 'get_color_palette' from module 'metainformant.visualization.style' | `from metainformant.visualization.style import get_color_pale` |
| docs/visualization/styling.md | 31 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization.style` |
| docs/visualization/timeseries.md | 12 | ImportError | Cannot import 'time_series_plot' from module 'metainformant.visualization' | `from metainformant.visualization import time_series_plot` |
| docs/visualization/timeseries.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
| docs/visualization/trees.md | 12 | ImportError | Cannot import 'sequences' from module 'metainformant.dna' | `from metainformant.dna import sequences` |
| docs/visualization/trees.md | 12 | ImportError | Cannot import 'plot_phylo_tree' from module 'metainformant.visualization' | `from metainformant.visualization import plot_phylo_tree` |
| docs/visualization/trees.md | 12 | AttributeError | 'metainformant' has no attribute 'dna' (module not found) | `metainformant.dna` |
| docs/visualization/trees.md | 12 | AttributeError | 'metainformant' has no attribute 'visualization' (module not found) | `metainformant.visualization` |
