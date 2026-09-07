# Pipeline Usage & Extension Guide

This guide details how to configure, maintain, and extend the orthology pipeline.

---

## Input Files & Schemas Reference

The pipeline organizes input files into **Required Configurations**, **Auto-Generated Metadata**, and **Auto-Downloaded Data**.

### 1. User-Provided Configuration & Input Files

| File / Setting | Format | Required | Description |
| :--- | :--- | :--- | :--- |
| `config.yaml` | YAML | **Yes** | Master configuration file (version, species list, tool paths, parameters). |
| `paths.metadata` (e.g. `input/gene_metadata/...`) | CSV | **Yes** | **Unified Gene & Transcript Metadata** (see schema below). Defines target genes, canonical transcripts, and peptide mappings. |
| `species_tree` (in `config.yaml`) | Newick | **Yes** | Species tree defined directly inline as a Newick string (or path to a `.nwk` file). **Tip labels MUST match species keys** in `config.yaml`. |

* **Unified Metadata Table Schema** (`paths.metadata`):
  | Column | Required | Description |
  | :--- | :--- | :--- |
  | `ensembl_gene_id` (or `ensembl_id`) | **Yes** | Ensembl Gene ID. |
  | `species` | **Yes** | Species identifier (e.g. `human`, `mouse`, `macaque`, or species codes). |
  | `ensembl_transcript_id` | **Yes** | Transcript ID for sequence extraction. |
  | `transcript_is_canonical` | Optional | Flag (`1` for canonical). If omitted, all rows are treated as canonical. |
  | `peptide_id` / `peptide_id_version` | Optional | Ensembl protein ID for CDS-to-gene mapping in synteny analysis. |
  | `gene_symbol` | Optional | Gene symbol for reporting (falls back to gene ID if missing). |

* **Zero External Preprocessing Needed**:
  With this unified metadata table, Snakemake runs directly end-to-end without requiring any intermediate preprocessing scripts or external BioMart queries!

* **`config.yaml` Schema Example**:
  ```yaml
  ensembl_version: "115"
  species_tree: "((((MACFA:1,HUMAN:1)1:1,(RATNO:1,MOUSE:1)1:1)1:1,(PIGXX:1,CANLF:1)1:1)1:1);"
  paths:
    metadata: "input/gene_metadata/target_gene_annotations.csv"
    raw_homologies_dir: "input/DB/ensembl_homologies_raw"
    genespace_coords: "input/DB/genespace_coords"
    fasta_dir: "input/DB/canonical_transcripts"
  tools:
    mafft: "mafft"
    fasttree: "FastTree"
    treerecs: "/path/to/treerecs"
  output_dir: "ensembl_pipeline_output"
  species:
    HUMAN:
      scientific_name: "homo_sapiens"
      assembly: "GRCh38"
    MOUSE:
      scientific_name: "mus_musculus"
      assembly: "GRCm39"
  ```

### 2. Auto-Downloaded Database Inputs (`input/DB/`)

Snakemake automatically downloads these from Ensembl FTP on the first run (no manual action needed):

* `input/DB/ensembl_homologies_raw/{scientific_name}.tsv.gz`: Raw Ensembl Compara homologies.
* `input/DB/raw_fastas/{species_code}.fa.gz`: Full protein FASTAs.
* `input/DB/gff3/{species_code}.gff3.gz`: Genomic feature annotations for CDS coordinate extraction.

---

## 1. Upgrading Ensembl Version (e.g., to Release 116+)

To move the entire pipeline to a new Ensembl release:

1.  **Edit `config.yaml`**:
    *   Update `ensembl_version: "XXX"`
    *   Assembly strings in the `species` section are **optional**—the pipeline will automatically resolve default assemblies via Ensembl REST API if omitted.
2.  **Provide Updated Unified Metadata**:
    *   Place the new version's gene and transcript annotations in `input/gene_metadata/` (or update `paths.metadata` in `config.yaml`).
3.  **Run Snakemake**:
    *   Snakemake will detect the new version and automatically download the fresh homologies, FASTAs, and GFF3 files from Ensembl FTP.
    ```bash
    snakemake --cores 4
    ```

## 2. Adding a New Species

1.  **Config**: Add the species to `config.yaml` using any nickname key you prefer:
    ```yaml
    species:
      macaque:                           # Custom species nickname (used for table columns & outputs)
        scientific_name: "macaca_fascicularis" # Ensembl scientific name
        # assembly: "Macaca_fascicularis_6.0" (Optional - auto-detected via REST API if omitted)
    ```
2.  **Metadata**: Ensure the new species gene and canonical transcript annotations are included in your master metadata table (`paths.metadata`).
3.  **Species Tree**: Add your species nickname to `species_tree` in `config.yaml`.
4.  **Run**: Execute Snakemake. It will build all pairwise networks automatically based on the species list in `config.yaml`:
    ```bash
    snakemake --cores 4
    ```

## 3. Validation & Diagnostics

### Pre-flight Setup Check
Before starting a long run, verify your tools, species tree, and API connectivity:
```bash
python preprocessing/00_check_pipeline_setup.py
```
This script checks tool execution permissions, Ensembl REST API connectivity, and metadata integrity.

### Results Summary
After a run completes, generate a biological summary report:
```bash
python scripts/09_generate_summary_stats.py
```
This provides orthogroup counts by cardinality (1to1, multi-copy) and species coverage statistics.

## 4. Troubleshooting "Incomplete Files"

If a download is interrupted, Snakemake might flag files as incomplete. To fix:
```bash
snakemake --cleanup-metadata <file_path>
snakemake --unlock
snakemake --rerun-incomplete
```

## 5. Manual Preprocessing

If you need to run the FASTA filtering or Coordinate extraction manually:
```bash
python scripts/02_prepare_orthology_inputs.py
```
*Note: This script dynamically reads `config.yaml` to retrieve target canonical transcript IDs directly from the unified metadata table.*

---

## Related Documentation

* [README.md](README.md) - Pipeline overview, architecture, and quick start.
* [METHODS_PIPELINE.md](METHODS_PIPELINE.md) - Detailed transitivity, synteny, and phylogenetic methods.
* [METRICS_AND_RANKING.md](METRICS_AND_RANKING.md) - Metric definitions and ortholog ranking algorithm.

