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
- **`atlas_msa_pairwise_identity`**: Direct pairwise protein sequence identity (%) calculated from MAFFT canonical peptide multiple sequence alignment.

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

## 6. Iterative Micro-Synteny & Sequence-Based Singleton Rescue Engine (v0.2)

Ensembl Compara occasionally leaves true orthologs unassigned as isolated "singletons" due to assembly gaps, strict heuristic cutoffs, or missing pairwise calls between specific species pairs. The pipeline provides a 100% species-agnostic rescue engine ([`scripts/09_singleton_rescue.py`](file:///Users/emilioskarwan/Documents/SciLifeDrive/get-orthogroups/scripts/09_singleton_rescue.py)) to recover these missing links:

### 1. Micro-Synteny Cassette Detection
* For every singleton gene, the engine identifies flanking conserved **macro-anchors** ($L_2, L_1$ upstream and $R_1, R_2$ downstream) within a configurable cassette window (`max_cassette_size: 5`).
* These anchors must belong to established orthogroups shared with the candidate target species.
* **Relative Strand Polarity**: The engine verifies that the candidate singleton and its flanking anchors preserve ancestral chromosomal orientation relative to the target species.

### 2. Dynamic Programming Sequence Verification (`Needleman-Wunsch / BLOSUM62`)
* Micro-synteny alone cannot distinguish true orthologous genes from non-homologous gene insertions or retrotranspositions.
* The engine aligns the candidate protein against the canonical peptide of the missing species in the candidate orthogroup using exact dynamic programming (**Needleman-Wunsch global alignment with the BLOSUM62 matrix**).
* Candidates must pass a score threshold (`min_alignment_score: 50.0`) to confirm sequence homology.

### 3. Iterative Anchor Chaining
* The rescue process runs across multiple iterations (`max_iterations: 3`).
* **Chaining mechanism**: When a singleton $A$ is successfully rescued in iteration 1, it becomes a validated flanking anchor in iteration 2, unlocking adjacent singletons in dense syntenic blocks until graph convergence.

### 4. Human-in-the-Loop Curation Loop
To prevent erroneous merges, proposed rescues are exported to transparent audit tables:
* `rescued_genes_all.tsv`: Complete audit log with flanking anchor context and scores.
* `rescue_curation_candidates.csv`: Comprehensive diagnostic review table.
* `manual_curation_sample.csv`: Pre-formatted proposal table for human review.

Users can curate and filter candidates (e.g., removing conjoined read-through chimeric transcripts) into `manual_curation.csv`. When provided, these validated edges are ingested during Step 3 (`ensembl_transitivity`), seamlessly integrating rescued orthologs into both the transitivity network and final consensus tables.

---

## 7. Phylogenetic Reconstruction & Species-Tree Reconciliation

For all multi-gene orthogroups ($\ge 2$ genes), the pipeline builds multiple sequence alignments and phylogenetic trees that integrate sequence evidence with species-tree reconciliation:

### Alignment (`MAFFT`)
* Canonical peptide sequences are extracted and aligned using **MAFFT** (`--auto`).
* The `--anysymbol` flag is enabled to handle non-standard amino acid characters, providing full native support for mammalian **selenoproteins** (which contain selenocysteine, `U`).

### Maximum-Likelihood Tree Inference (`FastTree`)
* For orthogroups with $\ge 3$ genes, initial phylogenetic trees are inferred using **FastTree** with the **Le-Gascuel (LG)** amino acid replacement matrix (`-lg`).
* FastTree estimates local Shimodaira-Hasegawa (SH)-like branch support values ($0.0 \le \text{support} \le 1.0$).

### Gene Tree–Species Tree Reconciliation (`Treerecs`)
* Sequence-only gene trees frequently exhibit topological discordance due to short sequence lengths, mutational saturation, or alignment noise. To distinguish true evolutionary events from reconstruction artifacts, trees are reconciled against the user-specified species tree using **Treerecs**:
  1. **Parsimonious Re-rooting (`-r`)**: FastTree generates unrooted trees. Treerecs determines the optimal evolutionary root that minimizes total gene duplication and loss penalties under the species tree, correcting arbitrary rooting.
  2. **Branch Support Contraction (`-t 0.8`)**: Branches with SH-like support $< 0.8$ are contracted into polytomies and rearranged parsimoniously to match the species tree topology. This eliminates false duplication and loss calls in clean orthologs while strictly preserving strongly supported ($\ge 0.8$) duplication signals.
  3. **Robust Fallbacks**: For pairwise 2-gene orthologs, a 2-taxon Newick tree is directly generated. If a complex gene family encounters an unresolvable reconciliation error, the pipeline safely falls back to FastTree's raw ML tree.

---

## 8. Automated Snakemake Workflow

The pipeline is fully automated and handles the entire lifecycle:

1.  **`download_data`**: Fetches raw TSV homologies, FASTAs, and GFF3 files directly from the Ensembl FTP.
2.  **`prepare_inputs`**: Filters raw sequences to **canonical-only** FASTAs and extracts genomic coordinates from GFF3.
3.  **`ensembl_transitivity`**: Constructs the all-vs-all dual networks (Standard & HQ) and extracts Connected Components (ingesting curated rescue edges if provided).
4.  **`run_singleton_rescue`** (Optional/Staged): Evaluates micro-synteny cassettes and BLOSUM62 sequence alignments to produce rescue proposals.
5.  **`build_consensus`**: Harmonizes IDs and prepares base pairwise tables.
6.  **`synteny_refinement`**: Calculates local informational GOC scores and assigns cardinality status tags.
7.  **`finalize_tables`**: Generates consensus symbols and produces the final `Consensus_Master.tsv`.
8.  **`compute_trees`** (Optional): Builds MAFFT alignments and FastTree/Treerecs reconciled gene trees with rerooting and branch contraction.
9.  **`generate_summary_stats`**: Compiles overall dataset metrics, species coverage, cardinality breakdowns, and outputs `Summary_Stats.md`.

---

## Related Documentation

* [README.md](README.md) - Pipeline overview and quick start.
* [USAGE_GUIDE.md](USAGE_GUIDE.md) - Configuration, upgrading Ensembl, and adding species.
* [METRICS_AND_RANKING.md](METRICS_AND_RANKING.md) - Metric definitions and ortholog ranking algorithm.
