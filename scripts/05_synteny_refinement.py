import pandas as pd
import numpy as np
import os
import itertools
import gzip
import traceback
import sys

def load_data(master_path):
    print("Loading consensus master...")
    master = pd.read_csv(master_path, sep='\t')
    return master

def load_coords(master, coords_dir, metadata_path):
    print("Loading protein-to-gene mapping...")
    meta = pd.read_csv(metadata_path)
    p2g = {}
    if 'peptide_id' in meta.columns:
        p2g.update(meta.dropna(subset=['peptide_id']).set_index('peptide_id')['ensembl_id'].to_dict())
    if 'peptide_id_version' in meta.columns:
        p2g.update(meta.dropna(subset=['peptide_id_version']).set_index('peptide_id_version')['ensembl_id'].to_dict())
    
    gene_to_og = master.set_index('ensembl_id')['ens_orthogroup_id'].to_dict()
    gene_to_hq_og = master.set_index('ensembl_id')['ens_hqorthogroup_id'].to_dict()
    
    print("Loading and sorting genomic coordinates...")
    species_coords = {}
    for f in os.listdir(coords_dir):
        if f.endswith("_coords.bed"):
            sp_code = f.replace("_coords.bed", "")
            path = os.path.join(coords_dir, f)
            coords = pd.read_csv(path, sep='\t', header=None, names=['chr', 'start', 'end', 'p_id'])
            
            coords['chr'] = coords['chr'].astype(str).str.replace(r'\.0$', '', regex=True)
            coords['g_id'] = coords['p_id'].map(p2g)
            coords['og'] = coords['g_id'].map(gene_to_og)
            coords['hq_og'] = coords['g_id'].map(gene_to_hq_og)
            
            coords = coords.dropna(subset=['g_id'])
            coords = coords.sort_values(['chr', 'start']).reset_index(drop=True)
            species_coords[sp_code] = coords
            
    return species_coords

def get_neighborhoods(coords_df, og_col, window_size=5):
    neighborhoods = {}
    valid_coords = coords_df.dropna(subset=[og_col])
    
    for i, row in coords_df.iterrows():
        gene = row['g_id']
        start_idx = max(0, i - window_size)
        end_idx = min(len(coords_df), i + window_size + 1)
        neighbors = coords_df.iloc[start_idx:end_idx][og_col].dropna().tolist()
        neighborhoods[gene] = set(neighbors)
    return neighborhoods

def calculate_pairwise_goc(g1, g2, neighborhoods1, neighborhoods2):
    if g1 not in neighborhoods1 or g2 not in neighborhoods2:
        return 0
    n1 = neighborhoods1[g1]
    n2 = neighborhoods2[g2]
    intersection = len(n1.intersection(n2))
    union = len(n1.union(n2))
    return intersection / union if union > 0 else 0

def load_ensembl_scores(homology_dir, species_config):
    print("Loading Ensembl homology scores...")
    scores_dict = {}
    
    for sp_code, info in species_config.items():
        sci_name = info['scientific_name']
        file_path = os.path.join(homology_dir, f"{sci_name}.tsv.gz")
        
        if os.path.exists(file_path):
            print(f"  Parsing {sci_name} homologies for scores...")
            with gzip.open(file_path, 'rt') as f:
                header = f.readline().strip().split('\t')
                col_ref = header.index('gene_stable_id')
                col_target = header.index('homology_gene_stable_id')
                
                score_cols = ['identity', 'homology_identity', 'goc_score', 'wga_coverage', 'is_high_confidence']
                col_indices = {c: header.index(c) for c in score_cols if c in header}
                
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) < max(col_indices.values()) + 1: continue
                    ref_id = parts[col_ref]
                    target_id = parts[col_target]
                    row_scores = {c: parts[i] for c, i in col_indices.items()}
                    scores_dict[(ref_id, target_id)] = row_scores
    return scores_dict

def main():
    try:
        master_path = snakemake.input.master
        coords_dir = snakemake.params.coords_dir
        metadata_path = snakemake.input.metadata
        homology_dir = snakemake.params.homology_dir
        
        output_master = snakemake.output.refined_master
        output_pairwise_dir = snakemake.output.pairwise_goc_dir
        os.makedirs(output_pairwise_dir, exist_ok=True)
        
        species_config = snakemake.config['species']
        window_size = snakemake.config['synteny']['goc_window']
        
        # Load the pure sequence-based clusters
        master = load_data(master_path)
        
        # Recalculate n_species based on standard OGs just to ensure it's correct
        n_sp = master.groupby('ens_orthogroup_id')['species_code'].nunique()
        master['n_species'] = master['ens_orthogroup_id'].map(n_sp)
        
        # Add Cardinality Status
        def assign_status(col_name):
            status_map = {}
            for og, group in master.groupby(col_name):
                if pd.isna(og) or og == 'unassigned':
                    status_map[og] = 'unassigned'
                elif 'SNG' in str(og):
                    status_map[og] = 'singleton'
                else:
                    n_genes = len(group)
                    n_species = group['species_code'].nunique()
                    if n_genes > n_species:
                        status_map[og] = 'multi_copy'
                    else:
                        status_map[og] = '1to1'
            return master[col_name].map(status_map)
            
        master['ens_orthogroup_status'] = assign_status('ens_orthogroup_id')
        master['ens_hqorthogroup_status'] = assign_status('ens_hqorthogroup_id')
        
        # Pass the master through unchanged structurally
        master.to_csv(output_master, sep='\t', index=False)
        
        # Load Coords and build neighborhoods for informational GOC
        species_coords = load_coords(master, coords_dir, metadata_path)
        og_neighborhoods = {sp: get_neighborhoods(df, 'og', window_size) for sp, df in species_coords.items()}
        hq_neighborhoods = {sp: get_neighborhoods(df, 'hq_og', window_size) for sp, df in species_coords.items()}
        
        # Load Ensembl scores for pairwise tables
        ensembl_scores = load_ensembl_scores(homology_dir, species_config)

        # Generate Refined Pairwise Tables
        print("Generating pairwise tables with Ensembl and Atlas scores...")
        species_codes = list(species_config.keys())
        for sp1, sp2 in itertools.combinations(species_codes, 2):
            df1 = master[master['species_code'] == sp1]
            df2 = master[master['species_code'] == sp2]
            
            merged = pd.merge(df1, df2, on='ens_orthogroup_id', suffixes=('_'+sp1, '_'+sp2))
            if merged.empty: continue

            # Append Ensembl metrics
            def get_ens_scores(row):
                g1 = row[f'ensembl_id_{sp1}']
                g2 = row[f'ensembl_id_{sp2}']
                s = ensembl_scores.get((g1, g2)) or ensembl_scores.get((g2, g1))
                if s:
                    return pd.Series([s.get('identity'), s.get('homology_identity'), s.get('goc_score'), s.get('wga_coverage'), s.get('is_high_confidence')])
                return pd.Series([None, None, None, None, None])

            ens_score_cols = ['ens_identity', 'ens_homology_identity', 'ens_goc_score', 'ens_wga_coverage', 'ens_is_high_confidence']
            merged[ens_score_cols] = merged.apply(get_ens_scores, axis=1)
            
            # Append specific Pairwise Atlas GOC scores
            merged['atlas_goc_score_pairwise'] = merged.apply(
                lambda row: calculate_pairwise_goc(row[f'ensembl_id_{sp1}'], row[f'ensembl_id_{sp2}'], og_neighborhoods[sp1], og_neighborhoods[sp2]), axis=1
            )
            merged['atlas_hq_goc_score_pairwise'] = merged.apply(
                lambda row: calculate_pairwise_goc(row[f'ensembl_id_{sp1}'], row[f'ensembl_id_{sp2}'], hq_neighborhoods[sp1], hq_neighborhoods[sp2]), axis=1
            )

            # --- ADD ORTHOLOG RANKING ---
            # Ensure numeric for ranking
            rank_metrics = ['ens_is_high_confidence', 'ens_goc_score', 'atlas_goc_score_pairwise', 'ens_identity', 'ens_homology_identity', 'ens_wga_coverage']
            for c in rank_metrics:
                merged[c] = pd.to_numeric(merged[c], errors='coerce').fillna(0)

            # Combined identity metric for sorting
            merged['_min_ident'] = merged[['ens_identity', 'ens_homology_identity']].min(axis=1)

            sort_cols = ['ens_is_high_confidence', 'ens_goc_score', 'atlas_goc_score_pairwise', '_min_ident', 'ens_wga_coverage']

            # Rank from sp1 perspective (ranking its sp2 orthologs)
            merged = merged.sort_values([f'ensembl_id_{sp1}'] + sort_cols, ascending=[True] + [False]*len(sort_cols))
            merged[f'ortholog_rank_{sp1}'] = merged.groupby(f'ensembl_id_{sp1}').cumcount() + 1

            # Rank from sp2 perspective (ranking its sp1 orthologs)
            merged = merged.sort_values([f'ensembl_id_{sp2}'] + sort_cols, ascending=[True] + [False]*len(sort_cols))
            merged[f'ortholog_rank_{sp2}'] = merged.groupby(f'ensembl_id_{sp2}').cumcount() + 1

            merged = merged.drop(columns=['_min_ident'])
            # ----------------------------
            
            out_path = os.path.join(output_pairwise_dir, f"{sp1}_{sp2}_pairwise_orthologs.tsv")
            merged.to_csv(out_path, sep='\t', index=False)
            
    except Exception as e:
        print(f"Error in main: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
