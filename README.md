# Ensembl Synteny-Adjusted Orthology Pipeline (`get-orthogroups`)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![Ensembl v115](https://img.shields.io/badge/ensembl-v115-orange.svg)](https://www.ensembl.org/)

A Snakemake workflow for mammalian orthology inference. The pipeline integrates curated Ensembl Compara homologies, dual-track transitivity clustering, local synteny (Atlas GOC), and deterministic ortholog ranking to resolve complex one-to-many and many-to-many orthology relationships.

The package is generalized to support arbitrary species sets and Ensembl releases.

---

## Key Features

* **Dual-Track Transitivity Clustering**:
  * **Standard Orthogroups (`OG_XXXXX`)**: Inclusive transitivity network capturing all sequence-similar homologs across evolutionary lineages.
  * **High-Confidence Orthogroups (`HQ_XXXXX`)**: Strict network constructed exclusively from Ensembl's curated `is_high_confidence == 1` relationships. Note that these results might exclude true in-paralogs for certain orthogroups, prioritizing the ancestral gene.
* **Directional Ortholog Ranking (`ortholog_rank`)**:
  * Resolves multi-copy gene arrays by assigning a directional ranking from the perspective of each species.
  * Prioritized multi-tier ranking hierarchy:
    1. `ens_is_high_confidence` (Ensembl curation)
    2. `ens_goc_score` (Ensembl global Gene Order Conservation)
    3. `atlas_goc_score_pairwise` (Our own synteny gene order conservation score based on local Jaccard similarity across a $\pm 5$ gene flanking window)
    4. `min(identity, homology_identity)` (Reciprocal sequence identity)
    5. `ens_wga_coverage` (Whole Genome Alignment coverage)
* **Zero External Preprocessing**:
  * Driven by a single unified master metadata table (`paths.metadata`), eliminating the need for separate transcript tables or fragile external BioMart/REST preprocessing scripts.
* **Bit-wise Reproducibility**:
  * Deterministic orthogroup ID assignment sorted by cluster size and lexicographical Ensembl IDs, guaranteeing identical cluster naming across independent runs and machines.
* **End-to-End Automation**:
  * Automatically fetches raw Compara homologies, protein FASTAs, and GFF3 annotations directly from Ensembl FTP.

---

## Repository Structure

```
get-orthogroups/
├── Snakefile                     # Snakemake workflow definition (Rules 1-7)
├── config.yaml                   # Master configuration (species, assemblies, tool paths)
├── envs/
│   └── orthology_env.yaml        # Conda environment definition
├── example/                      # Self-contained runnable test dataset
│   ├── config_example.yaml
│   └── input/
│       └── gene_metadata/
│           └── target_gene_annotations_example.csv
├── input/                        # Production inputs (user metadata & auto-fetched FTP DB)
│   ├── DB/                       # Auto-populated by Rule 1 (FASTAs, GFF3, homologies)
│   └── gene_metadata/            # Unified master annotation table
├── preprocessing/
│   └── 00_check_pipeline_setup.py # Pre-flight diagnostic & validation check
├── scripts/                      # Core pipeline execution scripts
│   ├── 01_download_ensembl_data.py
│   ├── 02_prepare_orthology_inputs.py
│   ├── 03_ensembl_transitivity.py
│   ├── 04_build_consensus.py
│   ├── 05_synteny_refinement.py
│   ├── 06_finalize_tables.py
│   ├── 07_batch_compute_trees.py
│   ├── 08_plot_trees.R
│   └── 09_generate_summary_stats.py
├── README.md                     # Main repository overview (this document)
├── USAGE_GUIDE.md                # Guide for upgrading Ensembl or adding species
├── METHODS_PIPELINE.md           # Methodological and theoretical documentation
└── METRICS_AND_RANKING.md        # Reference of all output metrics and ranking rules
```

---

## Installation & Setup

### 1. Clone the Repository
```bash
git clone https://github.com/emiliosk/get-orthogroups.git
cd get-orthogroups
```

### 2. Create the Conda Environment
Create and activate the environment using the provided YAML file:
```bash
conda env create -f envs/orthology_env.yaml
conda activate ensembl_orthology
```

### 3. Optional Phylogeny Tools (Tree Computation)
If you plan to run Rule 7 (`compute_trees: true` to generate reconciled gene trees and MSAs):
* **FastTree** and **MAFFT** are installed automatically inside the conda environment.
* **Treerecs** is required for gene/species tree reconciliation. Ensure `treerecs` is in your `PATH` or specify its absolute path under `tools.treerecs` in `config.yaml`.

---

## Quick Start

### Step 1: Run the Pre-flight Diagnostic Check
Verify your environment, Ensembl REST API connectivity, input schemas, and inline species tree:
```bash
python preprocessing/00_check_pipeline_setup.py
```

### Step 2: Test Run with the Bundled Example
A small 6-species test dataset (227 genes) is included in `example/`:
```bash
snakemake --configfile example/config_example.yaml --cores 4
```

### Step 3: Run the Production Pipeline
To run the full pipeline on your target gene annotations:
```bash
snakemake --cores 4 --latency-wait 60 --rerun-incomplete
```

### Step 4: Generate Summary Statistics
Inspect key metrics, core orthogroup counts, and cardinality breakdowns:
```bash
python scripts/09_generate_summary_stats.py
```
This produces a Markdown summary at `ensembl_pipeline_output/Summary_Stats.md`.

---

## Configuration (`config.yaml`)

The pipeline is configured via `config.yaml`:

```yaml
ensembl_version: "115"

# Inline Species Tree (Newick string matching species keys below)
species_tree: "((((MACFA:1,HUMAN:1)1:1,(RATNO:1,MOUSE:1)1:1)1:1,(PIGXX:1,CANLF:1)1:1)1:1);"

paths:
  metadata: "input/gene_metadata/target_gene_annotations.csv"
  raw_homologies_dir: "input/DB/ensembl_homologies_raw"
  genespace_coords: "input/DB/genespace_coords"
  fasta_dir: "input/DB/canonical_transcripts"

tools:
  mafft: "mafft"
  fasttree: "FastTree"
  treerecs: "treerecs"

output_dir: "ensembl_pipeline_output"
symbol_reference_species: "HUMAN"
compute_trees: false

species:
  HUMAN:
    scientific_name: "homo_sapiens"
    assembly: "GRCh38"
  MOUSE:
    scientific_name: "mus_musculus"
    assembly: "GRCm39"
  MACFA:
    scientific_name: "macaca_fascicularis"
    assembly: "Macaca_fascicularis_6.0"
  RATNO:
    scientific_name: "rattus_norvegicus"
    assembly: "GRCr8"
  PIGXX:
    scientific_name: "sus_scrofa"
    assembly: "Sscrofa11.1"
  CANLF:
    scientific_name: "canis_lupus_familiaris"
    assembly: "ROS_Cfam_1.0"

synteny:
  goc_window: 5
  strict_threshold: 0.3
```

### Unified Metadata Schema (`paths.metadata`)
The pipeline expects a single CSV file with the following columns:

| Column | Required | Description |
| :--- | :--- | :--- |
| `ensembl_gene_id` (or `ensembl_id`) | **Yes** | Primary Ensembl Gene Identifier. |
| `species` | **Yes** | Species name or code (e.g. `human`, `mouse`, `macaque`, `rat`, `pig`, `dog`). |
| `ensembl_transcript_id` | **Yes** | Canonical transcript ID used for sequence extraction. |
| `transcript_is_canonical` | Optional | Boolean / Flag (`1` or `True`). Rows with `1` are treated as canonical. |
| `peptide_id` / `peptide_id_version` | Optional | Ensembl protein ID for CDS-to-gene mapping in coordinate BED files. |
| `gene_symbol` | Optional | Standard gene symbol (falls back to `ensembl_gene_id` if missing). |

---

## Output Deliverables

All outputs are saved to `ensembl_pipeline_output/`:

1. **`Consensus_Master.tsv`**: Master cluster table containing all genes across all species, mapped to their Standard (`ens_orthogroup_id`) and High-Quality (`ens_hqorthogroup_id`) orthogroups, cardinality labels (`1to1`, `multi_copy`, `singleton`), and consensus gene symbols.
2. **`pairwise_tables/{SP1}_{SP2}_pairwise_orthologs.tsv`**: Comprehensive pairwise matrices for all species pairs (e.g., 15 tables for 6 species). Each table contains:
   * Reciprocal sequence identities (`ens_identity`, `ens_homology_identity`)
   * Synteny scores (`ens_goc_score`, `atlas_goc_score_pairwise`)
   * Whole genome alignment coverage (`ens_wga_coverage`)
   * Confidence flags (`ens_is_high_confidence`)
   * Directional ranks (`ortholog_rank_{SP1}`, `ortholog_rank_{SP2}`)
3. **`Phylogenetic_Trees.tar.gz`** (Optional): Generated when `compute_trees: true`. Contains MAFFT multiple sequence alignments and Treerecs reconciled phylogenetic trees for all multi-gene orthogroups.

---

## Documentation Links

* **[USAGE_GUIDE.md](USAGE_GUIDE.md)**: Step-by-step instructions for adding new species, upgrading Ensembl versions, and troubleshooting.
* **[METHODS_PIPELINE.md](METHODS_PIPELINE.md)**: In-depth methodological explanation of dual-track transitivity, cardinality tagging, and TreeRecs reconciliation.
* **[METRICS_AND_RANKING.md](METRICS_AND_RANKING.md)**: Formal mathematical definitions of GOC scores and the multi-tier ranking algorithm.

---

## License

Distributed under the MIT License. See `LICENSE` for more information.
