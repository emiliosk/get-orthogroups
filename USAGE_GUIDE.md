# Pipeline Usage & Extension Guide

This guide details how to maintain and extend the orthology pipeline.

## 1. Upgrading Ensembl Version (e.g., to Release 116+)

To move the entire pipeline to a new Ensembl release:

1.  **Edit `config.yaml`**:
    *   Update `ensembl_version: "XXX"`
    *   Update assembly names in the `species` section (check Ensembl FTP for the correct strings).
2.  **Generate Transcript Metadata**:
    *   The pipeline requires a file named `input/DB/ensembl_vXXX_multispecies_transcripts.csv`.
    *   Run the robust Python fetcher:
        ```bash
        python preprocessing/04_generate_canonical_transcripts_REST.py
        ```
3.  **Run Snakemake**:
    *   Snakemake will detect the new version and automatically download the fresh homologies, FASTAs, and GFF3 files.
    ```bash
    snakemake --cores 1
    ```

## 2. Adding a New Species

1.  **Metadata**: Ensure the new species gene IDs are present in your `input/gene_metadata/protein_coding_genes.csv`.
2.  **Config**: Add the species to `config.yaml`:
    ```yaml
    NEWSP:
      scientific_name: "scientific_name_here"
      assembly: "Assembly_Name"
      atlas_name: "newspecies"
    ```
3.  **Transcripts**: Re-run the transcript fetcher (step 1.2 above) to include the new species.
4.  **Run**: Execute Snakemake. It will build a 15-pair network (for 6 species) or 21-pair (for 7 species) automatically.

## 3. Validation & Diagnostics

### Pre-flight Setup Check
Before starting a long run, verify your tools and environment:
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

## 4. Manual Preprocessing

If you need to run the FASTA filtering or Coordinate extraction manually:
```bash
python scripts/02_prepare_orthology_inputs.py
```
*Note: This script now dynamically reads the config to find the correct versioned transcript file.*
