import pandas as pd

def analyze_consensus(file_path):
    print(f"Analyzing {file_path}...")
    df = pd.read_csv(file_path, sep='\t')
    
    # 1. Standard Orthogroups (ens_orthogroup_id)
    # Exclude singletons and dropped groups
    std_groups = df[~df['ens_orthogroup_id'].str.contains('SNG|_DROP|_UNASSIGNED', na=False)]
    std_counts = std_groups.groupby('ens_orthogroup_id')['species_code'].nunique()
    std_1to1 = std_counts[std_counts == 5].index
    
    # Verify they are truly 1:1 (exactly 5 genes total)
    std_total_counts = std_groups[std_groups['ens_orthogroup_id'].isin(std_1to1)].groupby('ens_orthogroup_id').size()
    std_perfect_1to1 = std_total_counts[std_total_counts == 5].index
    
    # 2. HQ Orthogroups (ens_hqorthogroup_id)
    hq_groups = df[~df['ens_hqorthogroup_id'].str.contains('SNG|_DROP|_UNASSIGNED', na=False)]
    hq_counts = hq_groups.groupby('ens_hqorthogroup_id')['species_code'].nunique()
    hq_1to1 = hq_counts[hq_counts == 5].index
    
    hq_total_counts = hq_groups[hq_groups['ens_hqorthogroup_id'].isin(hq_1to1)].groupby('ens_hqorthogroup_id').size()
    hq_perfect_1to1 = hq_total_counts[hq_total_counts == 5].index
    
    print(f"\n--- Consensus Statistics (1:1:1:1:1 across all species) ---")
    print(f"Total Genes: {len(df)}")
    print(f"Standard 1:1:1:1:1 Groups: {len(std_perfect_1to1)}")
    print(f"HQ 1:1:1:1:1 Groups:       {len(hq_perfect_1to1)}")
    print(f"Difference:               {len(std_perfect_1to1) - len(hq_perfect_1to1)}")
    
    # Detailed Size Distribution
    def get_distribution(col):
        counts = df.groupby(col).size().value_counts().sort_index()
        singletons = df[df[col].str.contains('SNG', na=False)].shape[0]
        return counts, singletons

    std_dist, std_sng = get_distribution('ens_orthogroup_id')
    hq_dist, hq_sng = get_distribution('ens_hqorthogroup_id')

    print("\n--- Group Size Distribution (How many groups have X genes) ---")
    print("Size | Standard Groups | HQ Groups")
    print("-----|-----------------|----------")
    all_sizes = sorted(set(std_dist.index) | set(hq_dist.index))
    for s in all_sizes:
        print(f"{s:4} | {std_dist.get(s, 0):15} | {hq_dist.get(s, 0):9}")
    
    print(f"\nTotal Singleton Genes (Standard): {std_sng}")
    print(f"Total Singleton Genes (HQ):       {hq_sng}")

if __name__ == "__main__":
    import os
    path = "ensembl_pipeline_output/Consensus_Master.tsv"
    if not os.path.exists(path):
        path = "../../ensembl_pipeline_output/Consensus_Master.tsv"
    
    if os.path.exists(path):
        analyze_consensus(path)
    else:
        print(f"Error: Could not find master table at {path}")
