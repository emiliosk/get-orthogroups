import pandas as pd
import numpy as np
import os
import itertools

def generate_consensus_master(cluster_path, metadata_path):
    print("Loading metadata...")
    meta = pd.read_csv(metadata_path)[['ensembl_id', 'gene_symbol', 'species']].drop_duplicates('ensembl_id')
    
    print("Loading orthology results...")
    cluster = pd.read_csv(cluster_path, sep='\t')[['ensembl_id', 'ens_orthogroup_id', 'ens_hqorthogroup_id']]

    master = cluster.copy()
    master['ens_orthogroup_id'] = master['ens_orthogroup_id'].fillna('unassigned')
    master['ens_hqorthogroup_id'] = master['ens_hqorthogroup_id'].fillna('unassigned')
    
    master = pd.merge(master, meta[['ensembl_id', 'gene_symbol']], on='ensembl_id', how='left').drop_duplicates()
    
    return master

def generate_base_pairwise(master_df, sp1, sp2, output_dir):
    # Map species directly from metadata if not present
    if 'species_code' not in master_df.columns:
        species_config = snakemake.config['species']
        standard_to_code = {info['atlas_name']: code for code, info in species_config.items()}
        
        meta = pd.read_csv(snakemake.input.metadata)
        gene_to_sp = meta.set_index('ensembl_id')['species'].map(standard_to_code).to_dict()
        master_df['species_code'] = master_df['ensembl_id'].map(gene_to_sp)

    df1 = master_df[master_df['species_code'] == sp1].copy()
    df2 = master_df[master_df['species_code'] == sp2].copy()
    
    # Only calculate pairwise connections for genes assigned to an Orthogroup
    df1 = df1[df1['ens_orthogroup_id'] != 'unassigned']
    df2 = df2[df2['ens_orthogroup_id'] != 'unassigned']
    
    combined = pd.merge(df1, df2, on='ens_orthogroup_id', suffixes=('_'+sp1, '_'+sp2))
    
    join_cols = [f'ensembl_id_{sp1}', f'ensembl_id_{sp2}', f'gene_symbol_{sp1}', f'gene_symbol_{sp2}']
    final = combined[join_cols + ['ens_orthogroup_id']]
    
    out_path = os.path.join(output_dir, f"{sp1}_{sp2}_pairwise_orthologs.tsv")
    final.to_csv(out_path, sep='\t', index=False)

def main():
    cluster_path = snakemake.input.cluster_hogs
    meta_path = snakemake.input.metadata
    
    out_master = snakemake.output.master
    pairwise_dir = snakemake.output.pairwise_dir
    os.makedirs(pairwise_dir, exist_ok=True)
    
    master = generate_consensus_master(cluster_path, meta_path)
    
    # ensure species_code is in master for pairwise gen
    species_config = snakemake.config['species']
    standard_to_code = {info['atlas_name']: code for code, info in species_config.items()}
    
    meta = pd.read_csv(meta_path)
    gene_to_sp = meta.set_index('ensembl_id')['species'].map(standard_to_code).to_dict()
    master['species_code'] = master['ensembl_id'].map(gene_to_sp)

    master.to_csv(out_master, sep='\t', index=False)
    
    species = master['species_code'].dropna().unique()
    pairs = list(itertools.combinations(species, 2))
    
    print(f"Generating {len(pairs)} base pairwise tables from Ensembl Clusters...")
    for sp1, sp2 in pairs:
        generate_base_pairwise(master, sp1, sp2, pairwise_dir)

if __name__ == "__main__":
    main()
