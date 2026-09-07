# Ensembl Synteny-Adjusted Orthology Pipeline (`get-orthogroups`)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![Ensembl v115](https://img.shields.io/badge/ensembl-v115-orange.svg)](https://www.ensembl.org/)

A Snakemake workflow for mammalian orthology inference. The pipeline processes Ensembl Compara pairwise homologies via dual-track transitivity clustering, evaluates local gene-order conservation using a chromosome-bounded neighborhood Jaccard metric, and ranks multi-copy orthologs through a sequential hierarchy of sequence and synteny criteria.

The workflow is configurable for arbitrary species sets and Ensembl releases.

---

## Methodological Overview

* **Dual-Track Transitivity Clustering**:
  * **Standard Orthogroups (`OG_XXXXX`)**: Connected components constructed from all pairwise Ensembl orthology calls. Captures full gene families, including tandem duplications and recent lineage-specific expansions.
  * **High-Confidence Orthogroups (`HQ_XXXXX`)**: Connected components restricted to pairwise links meeting Ensembl's `is_high_confidence == 1` threshold. This track prioritizes ancestral syntenic anchors and may separate secondary in-paralogs into singletons.
* **Local Synteny Scoring (`atlas_goc_score_pairwise`) [Custom Metric]**:
  * Custom metric calculated by this pipeline to evaluate synteny independently of Ensembl's linear 4-gene window.
  * Quantifies neighborhood conservation by computing the Jaccard similarity of orthogroup assignments within a $\pm 5$ gene flanking window along each chromosome.
* **Directional Ortholog Ranking (`ortholog_rank`)**:
  * For one-to-many and many-to-many relationships, candidate orthologs are sorted directionally from the query species perspective to identify the primary ancestral or most conserved counterpart (Rank 1).
  * Sequential ranking criteria (integrating Ensembl Compara annotations with our custom synteny score):
    1. **Ensembl High Confidence** (`ens_is_high_confidence`) — *[Ensembl Compara]*: Curated binary flag based on Ensembl's whole-genome alignment and gene-order thresholds.
    2. **Ensembl Gene Order Conservation** (`ens_goc_score`) — *[Ensembl Compara]*: Ensembl's 4-gene linear flanking synteny metric.
    3. **Custom Local Synteny** (`atlas_goc_score_pairwise`) — *[Custom Pipeline Metric]*: Set-based Jaccard similarity across the $\pm 5$ gene chromosome-bounded neighborhood.
    4. **Reciprocal Sequence Identity** (`min(ens_identity, ens_homology_identity)`) — *[Ensembl Compara]*: Minimum reciprocal protein sequence identity.
    5. **Whole Genome Alignment Coverage** (`ens_wga_coverage`) — *[Ensembl Compara]*: Depth of whole-genome sequence alignment coverage.
* **Unified Metadata Input**:
  * Driven by a single user-provided CSV table (`paths.metadata`) specifying gene identifiers, canonical transcript annotations, and peptide mappings across all target species.
* **Automated Data Retrieval**:
  * Snakemake automatically retrieves required Ensembl Compara homology tables, protein FASTAs, and GFF3 gene coordinates from the Ensembl FTP server on first execution.

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
* **Treerecs** is required for gene/species tree reconciliation:
  * **Linux / Intel Mac**: Install directly via Bioconda:
    ```bash
    conda install -c bioconda treerecs
    ```
  * **Apple Silicon (M1/M2/M3/M4)**: Bioconda provides `treerecs` built for `osx-64` (Intel). Run the automated installer to configure it via macOS Rosetta 2:
    ```bash
    bash scripts/install_treerecs.sh
    ```
    *(Alternatively, manually create an emulated environment: `CONDA_SUBDIR=osx-64 conda create -n treerecs_x86 -c bioconda treerecs -y && ln -sf $(conda info --base)/envs/treerecs_x86/bin/treerecs $CONDA_PREFIX/bin/`)*

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
# For production pipeline results:
python scripts/09_generate_summary_stats.py

# Or for the example test run:
python scripts/09_generate_summary_stats.py example/config_example.yaml
```
This produces a Markdown summary at `ensembl_pipeline_output/Summary_Stats.md` (or `example/output_example/Summary_Stats.md`).

---

## System Requirements & Benchmarks

| Metric | Bundled Example (227 genes) | Full Genome Analysis (~130,000 genes, 6 species) |
| :--- | :--- | :--- |
| **Peak RAM** | < 1 GB | ~3.5 GB (8 GB+ recommended) |
| **Core Pipeline (Rules 1–6)** | ~15–30 seconds | ~2–3 minutes (on 4 CPU cores) |
| **Phylogeny (Rule 7: MSAs & Trees)** | ~5 seconds (13 trees) | ~30–45 minutes (~20,000 trees on 4 cores) |
| **Disk Space** | ~50 MB | ~1.5 GB (raw Ensembl FTP cache + final tables) |

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
   * Reciprocal sequence identities (`ens_identity`, `ens_homology_identity`) — *[Ensembl]*
   * Ensembl Gene Order Conservation (`ens_goc_score`) — *[Ensembl]*
   * Local synteny Jaccard score (`atlas_goc_score_pairwise`) — *[Custom Metric]*
   * Whole genome alignment coverage (`ens_wga_coverage`) — *[Ensembl]*
   * Confidence flags (`ens_is_high_confidence`) — *[Ensembl]*
   * Directional ranks (`ortholog_rank_{SP1}`, `ortholog_rank_{SP2}`) — *[Custom Pipeline Metric]*
3. **`Phylogenetic_Trees.tar.gz`** (Optional): Generated when `compute_trees: true`. Contains MAFFT multiple sequence alignments and Treerecs reconciled phylogenetic trees for all multi-gene orthogroups.

---

## Documentation Links

* **[USAGE_GUIDE.md](USAGE_GUIDE.md)**: Step-by-step instructions for adding new species, upgrading Ensembl versions, and troubleshooting.
* **[METHODS_PIPELINE.md](METHODS_PIPELINE.md)**: In-depth methodological explanation of dual-track transitivity, cardinality tagging, and TreeRecs reconciliation.
* **[METRICS_AND_RANKING.md](METRICS_AND_RANKING.md)**: Formal mathematical definitions of GOC scores and the multi-tier ranking algorithm.

## Citation

If you use this pipeline or data in your research, please cite:

```bibtex
@article{skarwan2026orthology,
  author    = {Skarwan, Emilio},
  title     = {Ensembl Synteny-Adjusted Orthology Pipeline: Graph transitivity and local gene order for mammalian comparative genomics},
  journal   = {Manuscript in preparation},
  year      = {2026},
  url       = {https://github.com/emiliosk/get-orthogroups}
}
```

---

## License

Distributed under the MIT License. See `LICENSE` for more information.
