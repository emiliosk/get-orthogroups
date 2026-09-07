# Ensembl Synteny-Adjusted Orthology Pipeline

> [!NOTE]
> For the complete GitHub package documentation, installation, and usage guide, please see [README.md](README.md).

This pipeline establishes a high-fidelity orthology baseline for the Mammalian RNA Atlas. It leverages curated Ensembl Compara homologies and enriches them with local genomic synteny and intelligent ranking.

## Core Features

### 1. Dual-Track Clustering
The pipeline generates two levels of orthologous groups:
*   **Standard (`OG_`)**: Inclusive clusters capturing all sequence-similar homologs.
*   **High-Quality (`HQ_`)**: Strict clusters using only Ensembl's `is_high_confidence == 1` links.

### 2. Intelligent Ortholog Ranking
For every gene-ortholog pair, the pipeline calculates a directional `ortholog_rank`. This allows researchers to immediately identify the "best" ortholog (Rank 1) in multi-copy relationships.
**Ranking Criteria (Descending Priority):**
1.  Ensembl High Confidence Status
2.  Ensembl GOC Score
3.  Atlas Pairwise GOC Score (custom neighborhood Jaccard)
4.  Minimum Sequence Identity
5.  Ensembl WGA Coverage

### 3. Version & Species Agnostic
The pipeline is designed to be portable. It has been verified to run from scratch on:
*   **Ensembl Release 109** (Standard Atlas baseline)
*   **Ensembl Release 115** (Validated via standalone test)
*   **Arbitrary Species Sets** (Fully flexible via `config.yaml` species nicknames)

## Input Files Overview

*   **`config.yaml`**: Master configuration file defining species nicknames, scientific names, inline species tree, tool paths, and thresholds.
*   **`paths.metadata`** (e.g., `input/gene_metadata/target_gene_annotations.csv`): Single unified metadata table containing target genes, canonical transcript IDs, peptide IDs, and gene symbols across all species.
*   **Species Tree**: Embedded directly as an inline Newick string in `config.yaml` (`species_tree: "..."`) matching species keys (or can point to a file path).
*   **Ensembl FTP Data** (Compara homologies, FASTAs, GFF3 files): Auto-downloaded by Snakemake on the first run.

For full schemas and column requirements, see [USAGE_GUIDE.md](USAGE_GUIDE.md).

## Quick Start

```bash
# 1. Run Pre-flight Diagnostic Check
python preprocessing/00_check_pipeline_setup.py

# 2. Run Pipeline via Snakemake
snakemake --cores 4
# (If resuming an interrupted run, add --rerun-incomplete)

# 3. View Results Summary
python scripts/09_generate_summary_stats.py
```

## Output Structure

*   `ensembl_pipeline_output/Consensus_Master.tsv`: The main cluster table.
*   `ensembl_pipeline_output/pairwise_tables/`: Detailed tables for every species pair, containing all 7 metrics and the new `ortholog_rank`.
*   `ensembl_pipeline_output/Phylogenetic_Trees.tar.gz`: Reconciled gene trees and MSAs for every orthogroup.

For detailed instructions on adding species or upgrading versions, see [USAGE_GUIDE.md](USAGE_GUIDE.md).

