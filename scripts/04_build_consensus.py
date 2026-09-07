import pandas as pd
import numpy as np
import os
import itertools

def generate_consensus_master(cluster_path, metadata_path):
    print("Loading metadata...")
    meta = pd.read_csv(metadata_path)
    gene_col = 'ensembl_id' if 'ensembl_id' in meta.columns else 'ensembl_gene_id'
    meta['ensembl_id'] = meta[gene_col].astype(str).str.strip()
    
    symbol_col = 'gene_symbol' if 'gene_symbol' in meta.columns else ('external_gene_name' if 'external_gene_name' in meta.columns else None)
    if symbol_col and symbol_col in meta.columns:
        meta['gene_symbol'] = meta[symbol_col].fillna(meta['ensembl_id'])
    else:
        meta['gene_symbol'] = meta['ensembl_id']
        
    meta = meta[['ensembl_id', 'gene_symbol', 'species']].drop_duplicates('ensembl_id')
    
    print("Loading orthology results...")
    cluster = pd.read_csv(cluster_path, sep='\t')[['ensembl_id', 'ens_orthogroup_id', 'ens_hqorthogroup_id']]

    master = cluster.copy()
    master['ens_orthogroup_id'] = master['ens_orthogroup_id'].fillna('unassigned')
    master['ens_hqorthogroup_id'] = master['ens_hqorthogroup_id'].fillna('unassigned')
    
    master = pd.merge(master, meta[['ensembl_id', 'gene_symbol']], on='ensembl_id', how='left').drop_duplicates()
    
    return master

def get_standard_to_code(species_config):
    mapping = {}
    for code, info in species_config.items():
        mapping[code] = code
        mapping[code.lower()] = code
        mapping[code.upper()] = code
        sci = info.get('scientific_name', '')
        if sci:
            mapping[sci] = code
            mapping[sci.lower()] = code
        for extra_key in ['atlas_name', 'common_name']:
            if extra_key in info and info[extra_key]:
                val = str(info[extra_key])
                mapping[val] = code
                mapping[val.lower()] = code

    # Common English species names auto-mapping fallback if not explicitly in config
    common_defaults = {
        'human': 'HUMAN',
        'mouse': 'MOUSE',
        'rat': 'RATNO',
        'macaque': 'MACFA',
        'pig': 'PIGXX',
        'dog': 'CANLF',
        'marmoset': 'CALJA',
        'chimp': 'PANTR',
        'bonobo': 'PANPA',
        'gorilla': 'GORGO',
        'rhesus': 'MACMU'
    }
    for common_name, default_code in common_defaults.items():
        if common_name not in mapping:
            if default_code in species_config:
                mapping[common_name] = default_code
            else:
                for code in species_config.keys():
                    if common_name.lower() in code.lower() or code.lower().startswith(common_name[:3].lower()):
                        mapping[common_name] = code
                        break
    return mapping

def generate_base_pairwise(master_df, sp1, sp2, output_dir):
    # Map species directly from metadata if not present
    if 'species_code' not in master_df.columns:
        species_config = snakemake.config['species']
        standard_to_code = get_standard_to_code(species_config)
        
        meta = pd.read_csv(snakemake.input.metadata)
        gene_col = 'ensembl_id' if 'ensembl_id' in meta.columns else 'ensembl_gene_id'
        meta['ensembl_id'] = meta[gene_col].astype(str).str.strip()
        gene_to_sp = meta.drop_duplicates('ensembl_id').set_index('ensembl_id')['species'].astype(str).str.lower().map(standard_to_code).to_dict()
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
    standard_to_code = get_standard_to_code(species_config)
    
    meta = pd.read_csv(meta_path)
    gene_col = 'ensembl_id' if 'ensembl_id' in meta.columns else 'ensembl_gene_id'
    meta['ensembl_id'] = meta[gene_col].astype(str).str.strip()
    gene_to_sp = meta.drop_duplicates('ensembl_id').set_index('ensembl_id')['species'].astype(str).str.lower().map(standard_to_code).to_dict()
    master['species_code'] = master['ensembl_id'].map(gene_to_sp)

    master.to_csv(out_master, sep='\t', index=False)
    
    species = master['species_code'].dropna().unique()
    pairs = list(itertools.combinations(species, 2))
    
    print(f"Generating {len(pairs)} base pairwise tables from Ensembl Clusters...")
    for sp1, sp2 in pairs:
        generate_base_pairwise(master, sp1, sp2, pairwise_dir)

if __name__ == "__main__":
    main()
