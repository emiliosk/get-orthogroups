configfile: "config.yaml"
import os

# Target variables from config
VERSION = config.get("ensembl_version", "109")
OUTPUT_DIR = config["output_dir"]
SPECIES_LIST = [info['scientific_name'] for info in config["species"].values()]
SPECIES_CODES = list(config["species"].keys())

# Toggles from config
COMPUTE_TREES = config.get("compute_trees", False)
RESCUE_SUGGESTIONS = config.get("rescue", {}).get("generate_suggestions", False) or config.get("rescue", {}).get("enabled", False)

# Baseline targets
ALL_TARGETS = [
    f"{OUTPUT_DIR}/Consensus_Master.tsv",
    f"{OUTPUT_DIR}/pairwise_tables",
    f"{OUTPUT_DIR}/Summary_Stats.md"
]

# Conditionally append rescue deliverables
if RESCUE_SUGGESTIONS:
    ALL_TARGETS.extend([
        f"{OUTPUT_DIR}/rescue/rescued_genes_all.tsv",
        f"{OUTPUT_DIR}/rescue/rescue_curation_candidates.csv",
        f"{OUTPUT_DIR}/rescue/manual_curation_sample.csv"
    ])

# Conditionally append tree computations
if COMPUTE_TREES:
    ALL_TARGETS.append(f"{OUTPUT_DIR}/Phylogenetic_Trees.tar.gz")

# --- Default Target Rule ---
rule all:
    input:
        ALL_TARGETS

# --- Dedicated CLI Alias Target Rules ---
rule baseline:
    input:
        f"{OUTPUT_DIR}/Consensus_Master.tsv",
        f"{OUTPUT_DIR}/pairwise_tables",
        f"{OUTPUT_DIR}/Summary_Stats.md"

rule rescue:
    input:
        f"{OUTPUT_DIR}/rescue/rescued_genes_all.tsv",
        f"{OUTPUT_DIR}/rescue/rescue_curation_candidates.csv",
        f"{OUTPUT_DIR}/rescue/manual_curation_sample.csv"

rule trees:
    input:
        f"{OUTPUT_DIR}/Phylogenetic_Trees.tar.gz"

# --- Core Pipeline Rules ---
rule download_data:
    output:
        homologies = expand(config["paths"]["raw_homologies_dir"] + "/{species}.tsv.gz", species=SPECIES_LIST),
        fastas = expand("input/DB/raw_fastas/{code}.fa.gz", code=SPECIES_CODES),
        gffs = expand("input/DB/gff3/{code}.gff3.gz", code=SPECIES_CODES)
    script:
        "scripts/01_download_ensembl_data.py"

rule prepare_inputs:
    input:
        fastas = expand("input/DB/raw_fastas/{code}.fa.gz", code=SPECIES_CODES),
        gffs = expand("input/DB/gff3/{code}.gff3.gz", code=SPECIES_CODES),
        transcripts = config["paths"].get("transcripts", config["paths"]["metadata"]),
        metadata = config["paths"]["metadata"]
    output:
        canonical_fastas = expand(config["paths"]["fasta_dir"] + "/{code}.fa", code=SPECIES_CODES),
        coords = expand(config["paths"]["genespace_coords"] + "/{code}_coords.bed", code=SPECIES_CODES)
    script:
        "scripts/02_prepare_orthology_inputs.py"

# Curation input dependency resolution
curation_file = config.get("curation", {}).get("manual_file", "")
if not curation_file and os.path.exists("input/manual_curation.csv"):
    curation_file = "input/manual_curation.csv"
curation_input = [curation_file] if curation_file and os.path.exists(curation_file) else []

rule ensembl_transitivity:
    input:
        metadata = config["paths"]["metadata"],
        homologies = expand(config["paths"]["raw_homologies_dir"] + "/{species}.tsv.gz", species=SPECIES_LIST),
        curation = curation_input
    output:
        baseline = f"{OUTPUT_DIR}/Ensembl_N0_Clusters.tsv"
    script:
        "scripts/03_ensembl_transitivity.py"

rule run_singleton_rescue:
    input:
        baseline_clusters = f"{OUTPUT_DIR}/Ensembl_N0_Clusters.tsv",
        metadata = config["paths"]["metadata"],
        canonical_fastas = expand(config["paths"]["fasta_dir"] + "/{code}.fa", code=SPECIES_CODES),
        coords = expand(config["paths"]["genespace_coords"] + "/{code}_coords.bed", code=SPECIES_CODES)
    output:
        rescued_tsv = f"{OUTPUT_DIR}/rescue/rescued_genes_all.tsv",
        candidates_csv = f"{OUTPUT_DIR}/rescue/rescue_curation_candidates.csv",
        curation_sample_csv = f"{OUTPUT_DIR}/rescue/manual_curation_sample.csv"
    params:
        coords_dir = config["paths"]["genespace_coords"],
        fasta_dir = config["paths"]["fasta_dir"],
        species_codes = SPECIES_CODES,
        min_alignment_score = config.get("rescue", {}).get("min_alignment_score", 50.0),
        min_consensus_species = config.get("rescue", {}).get("min_consensus_species_unnamed", 2),
        max_cassette_size = config.get("rescue", {}).get("max_cassette_size", 5),
        max_iterations = config.get("rescue", {}).get("max_iterations", 3)
    script:
        "scripts/09_singleton_rescue.py"

rule build_consensus:
    input:
        cluster_hogs = f"{OUTPUT_DIR}/Ensembl_N0_Clusters.tsv",
        metadata = config["paths"]["metadata"]
    output:
        master = temp(f"{OUTPUT_DIR}/Consensus_Master_Initial.tsv"),
        pairwise_dir = temp(directory(f"{OUTPUT_DIR}/pairwise_raw"))
    script:
        "scripts/04_build_consensus.py"

rule synteny_refinement:
    input:
        master = f"{OUTPUT_DIR}/Consensus_Master_Initial.tsv",
        pairwise_dir = f"{OUTPUT_DIR}/pairwise_raw",
        coords = expand(config["paths"]["genespace_coords"] + "/{code}_coords.bed", code=SPECIES_CODES),
        metadata = config["paths"]["metadata"],
        homologies = expand(config["paths"]["raw_homologies_dir"] + "/{species}.tsv.gz", species=SPECIES_LIST)
    output:
        refined_master = temp(f"{OUTPUT_DIR}/Consensus_Master_Refined.tsv"),
        pairwise_goc_dir = directory(f"{OUTPUT_DIR}/pairwise_tables")
    threads:
        config.get("threads", 4)
    params:
        coords_dir = config["paths"]["genespace_coords"],
        homology_dir = config["paths"]["raw_homologies_dir"]
    script:
        "scripts/05_synteny_refinement.py"

rule finalize_tables:
    input:
        refined_master = f"{OUTPUT_DIR}/Consensus_Master_Refined.tsv"
    output:
        final_master = f"{OUTPUT_DIR}/Consensus_Master.tsv"
    script:
        "scripts/06_finalize_tables.py"

rule compute_trees:
    input:
        master = f"{OUTPUT_DIR}/Consensus_Master.tsv",
        canonical_fastas = expand(config["paths"]["fasta_dir"] + "/{code}.fa", code=SPECIES_CODES)
    output:
        tar = f"{OUTPUT_DIR}/Phylogenetic_Trees.tar.gz"
    threads:
        config.get("threads", 4)
    script:
        "scripts/07_batch_compute_trees.py"

rule generate_summary_stats:
    input:
        master = f"{OUTPUT_DIR}/Consensus_Master.tsv",
        pairwise = f"{OUTPUT_DIR}/pairwise_tables",
        tree_tar = f"{OUTPUT_DIR}/Phylogenetic_Trees.tar.gz" if COMPUTE_TREES else []
    output:
        stats = f"{OUTPUT_DIR}/Summary_Stats.md"
    script:
        "scripts/08_generate_summary_stats.py"
