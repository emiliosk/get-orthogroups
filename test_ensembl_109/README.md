# Ensembl 109 Test Suite

This directory is an isolated test environment for running the pipeline against **Ensembl Release 109**.

All raw data downloads, canonical databases, and output tables are kept isolated inside `test_ensembl_109/` without conflicting with production runs or other Ensembl releases.

---

## 1. Input Data Setup

Place your input tables into `test_ensembl_109/input/gene_metadata/`:

### Option A: Two Separate Files (Default Configured)

1. **`protein_coding_genes.csv`** (configured as `paths.metadata`):
   * **Required Columns**:
     * `ensembl_gene_id` (or `ensembl_id`): Ensembl Gene ID (e.g., `ENSG00000139618`).
     * `species`: Species name or code (e.g. `human`, `mouse`, `macaque`, `rat`, `pig`, `dog`, or codes `HUMAN`, `MOUSE`, etc.).
   * **Optional Columns**:
     * `gene_symbol`: Gene symbol (e.g., `BRCA2`).

2. **`canonical_transcripts.csv`** (configured as `paths.transcripts`):
   * **Required Columns**:
     * `ensembl_transcript_id` (or `transcript_id`): Canonical transcript ID (e.g., `ENST00000380152`).
   * **Optional Columns**:
     * `transcript_is_canonical`: Flag (`1` or `True`). If omitted, all rows in the file are treated as canonical transcripts.

---

### Option B: Single Combined Metadata File

If you have all annotations in a single file (like the production `target_gene_annotations.csv`), update `test_ensembl_109/config_e109.yaml`:
```yaml
paths:
  metadata: "test_ensembl_109/input/gene_metadata/target_gene_annotations.csv"
  transcripts: "test_ensembl_109/input/gene_metadata/target_gene_annotations.csv"
```

---

## 2. Ensembl 109 Assembly Specifications

In Ensembl Release 109, the genome assemblies are:
* **HUMAN**: `GRCh38`
* **MOUSE**: `GRCm39`
* **MACFA**: `Macaca_fascicularis_6.0`
* **RATNO**: `mRatBN7.2` *(Note: In e115 this was GRCr8; in e109 it is mRatBN7.2)*
* **PIGXX**: `Sscrofa11.1`
* **CANLF**: `ROS_Cfam_1.0`

---

## 3. Running the Test Pipeline

Once your CSV files are in `test_ensembl_109/input/gene_metadata/`, run:

### Step 1: Pre-flight Verification & Data Initialization
Verifies the environment and metadata, downloads raw Ensembl 109 homologies, protein FASTAs, and GFF3 files, and filters canonical transcripts:
```bash
snakemake --configfile test_ensembl_109/config_e109.yaml --cores 4 init
```

### Step 2: Generate Singleton Rescue Suggestions (v0.2)
Runs initial clustering and the iterative singleton rescue engine (up to 3 iterations):
```bash
snakemake --configfile test_ensembl_109/config_e109.yaml --cores 4 rescue
```
Outputs are written to `test_ensembl_109/output/rescue/`.

### Step 3: Human-in-the-Loop Curation (Optional)
Copy and edit the candidate decisions:
```bash
cp test_ensembl_109/output/rescue/manual_curation_sample.csv test_ensembl_109/input/manual_curation.csv
```
Set `curator_decision` to `ACCEPT` or `REJECT`, then set `curation.manual_file: "test_ensembl_109/input/manual_curation.csv"` in `config_e109.yaml`.

### Step 4: Finalize Baseline
Computes final orthogroups, consensus tables, Jaccard synteny scores, pairwise matrices, and summary statistics:
```bash
snakemake --configfile test_ensembl_109/config_e109.yaml --cores 4 baseline
```

### Step 5: Phylogenetic Trees (Optional)
```bash
snakemake --configfile test_ensembl_109/config_e109.yaml --cores 8 trees
```
