# Methodology: Ensembl Synteny-Adjusted Orthology Pipeline

This document describes the orthology framework developed for the Mammalian RNA Atlas. The pipeline leverages curated Ensembl Compara pairwise homologies and local genomic synteny to produce a robust consensus orthology baseline.

## 1. Core Logic: Dual-Track Transitivity

The pipeline generates two distinct sets of orthologous groups based on pure sequence similarity to allow for inclusive versus strict comparative analysis:

### Track A: Standard Orthogroups (`ens_orthogroup_id`)
- **Logic:** Builds a transitivity network using all Ensembl `ortholog_` relationships.
- **Goal:** Cast a wide net to capture all possible evolutionary lineages and sequence-similar homologs into inclusive clusters.
- **Prefix:** `OG_` (e.g., `OG_00001`).

### Track B: High-Confidence Orthogroups (`ens_hqorthogroup_id`)
- **Logic:** Builds a network using only Ensembl relationships explicitly flagged as **`is_high_confidence == 1`**.
- **Goal:** Provide a curated subset of orthologs that meet Ensembl's strict internal criteria for Gene Order Conservation (GOC) and Whole Genome Alignment (WGA) coverage.
- **Prefix:** `HQ_` (e.g., `HQ_00001`).

---

## 2. Structural Classification (Cardinality Status)

Every gene is assigned a status label for both the Standard and HQ tracks based purely on the cardinality and composition of its orthogroup:

- **`1to1`**: Clean orthogroups containing at most one gene per species.
- **`multi_copy`**: Orthogroups containing multiple genes for at least one species (e.g., tandem duplicates).
- **`singleton`**: Orthogroups containing exactly one gene (no identified orthologs in other species).

---

## 3. Enriched Evidence Tables

In addition to the master consensus, the pipeline produces enriched pairwise orthology tables (e.g., `HUMAN_MOUSE_pairwise_orthologs.tsv`) in the `pairwise_tables/` directory. These tables integrate the raw clustering with both Ensembl and Atlas evidence:

### Ensembl Native Metrics
- **`ens_identity / ens_homology_identity`**: Sequence identity scores for both the query and target genes.
- **`ens_goc_score`**: Ensembl's global Gene Order Conservation metric.
- **`ens_wga_coverage`**: Whole Genome Alignment coverage score.
- **`ens_is_high_confidence`**: Ensembl's binary confidence flag.

### Atlas Local Metrics (Informational)
- **`atlas_goc_score_pairwise`**: Local Jaccard similarity of the +/- 5 gene neighborhood using Standard labels.
- **`atlas_hq_goc_score_pairwise`**: Local Jaccard similarity of the neighborhood using High-Confidence labels.

---

## 4. Intelligent Ortholog Ranking

To resolve `one-to-many` and `many-to-many` relationships (e.g., when a Human gene has multiple orthologs in Mouse), the pipeline assigns a directional **`ortholog_rank`** from the perspective of every species.

This allows researchers to immediately identify the "best" or "ancestral" ortholog (Rank 1) in complex clusters.

**Ranking Criteria (Descending Priority):**
1.  **Ensembl High Confidence**: Prioritize links explicitly curated as high-confidence by Ensembl.
2.  **Ensembl GOC Score**: Prioritize genes with higher global Gene Order Conservation.
3.  **Atlas Pairwise GOC**: Prioritize genes with high local neighborhood support in the Atlas.
4.  **Min(Identity)**: Use the best overall sequence similarity match.
5.  **WGA Coverage**: Final tie-breaker using alignment depth.

---

## 5. Bit-wise Reproducibility (Stable IDs)

Orthogroup IDs (`OG_XXXXX` and `HQ_XXXXX`) are assigned deterministically to ensure bit-wise reproducibility across different machines and runs.

**Logic:**
Before ID assignment, all connected components in the graph are sorted by:
1.  **Size** (Descending).
2.  **Gene ID string** (Ascending, using a sorted list of all Ensembl IDs within the component).

This ensures that `OG_00001` will always refer to the same set of genes, even if the graph traversal or data order varies.

---

## 6. Automated Snakemake Workflow

The pipeline is fully automated and handles the entire lifecycle:

1.  **`download_data`**: Fetches raw TSV homologies, FASTAs, and GFF3 files directly from the Ensembl FTP.
2.  **`prepare_inputs`**: Filters raw sequences to **canonical-only** FASTAs and extracts genomic coordinates from GFF3.
3.  **`ensembl_transitivity`**: Constructs the all-vs-all dual networks (Standard & HQ) and extracts Connected Components.
4.  **`build_consensus`**: Harmonizes IDs and prepares base pairwise tables.
5.  **`synteny_refinement`**: Calculates local informational GOC scores and assigns cardinality status tags.
6.  **`finalize_tables`**: Generates consensus symbols and produces the final `Consensus_Master.tsv`.
7.  **`compute_trees`** (Optional): Builds MAFFT alignments and FastTree/Treerecs reconciled gene trees.
8.  **`generate_summary_stats`**: Compiles overall dataset metrics, species coverage, cardinality breakdowns, and outputs `Summary_Stats.md`.

---

## Related Documentation

* [README.md](README.md) - Pipeline overview and quick start.
* [USAGE_GUIDE.md](USAGE_GUIDE.md) - Configuration, upgrading Ensembl, and adding species.
* [METRICS_AND_RANKING.md](METRICS_AND_RANKING.md) - Metric definitions and ortholog ranking algorithm.
