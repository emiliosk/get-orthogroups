configfile: "config.yaml"

# Target variables from config
VERSION = config.get("ensembl_version", "109")
OUTPUT_DIR = config["output_dir"]
SPECIES_LIST = [info['scientific_name'] for info in config["species"].values()]
SPECIES_CODES = list(config["species"].keys())

COMPUTE_TREES = config.get("compute_trees", False)

ALL_TARGETS = [
    f"{OUTPUT_DIR}/Consensus_Master.tsv",
    f"{OUTPUT_DIR}/pairwise_tables",
    f"{OUTPUT_DIR}/Summary_Stats.md"
]
if COMPUTE_TREES:
    ALL_TARGETS.append(f"{OUTPUT_DIR}/Phylogenetic_Trees.tar.gz")

rule all:
    input:
        ALL_TARGETS

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

rule ensembl_transitivity:
    input:
        metadata = config["paths"]["metadata"],
        homologies = expand(config["paths"]["raw_homologies_dir"] + "/{species}.tsv.gz", species=SPECIES_LIST)
    output:
        baseline = f"{OUTPUT_DIR}/Ensembl_N0_Clusters.tsv"
    script:
        "scripts/03_ensembl_transitivity.py"

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
        master = f"{OUTPUT_DIR}/Consensus_Master.tsv"
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

