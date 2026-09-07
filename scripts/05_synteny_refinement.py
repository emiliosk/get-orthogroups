import pandas as pd
import numpy as np
import os
import itertools
import gzip
import traceback
import sys
from concurrent.futures import ThreadPoolExecutor

def load_data(master_path):
    print("Loading consensus master...")
    master = pd.read_csv(master_path, sep='\t')
    return master

def load_coords(master, coords_dir, metadata_path):
    print("Loading protein-to-gene mapping...")
    meta = pd.read_csv(metadata_path)
    gene_col = 'ensembl_id' if 'ensembl_id' in meta.columns else 'ensembl_gene_id'
    meta['ensembl_id'] = meta[gene_col].astype(str).str.strip()
    p2g = {}
    if 'peptide_id' in meta.columns:
        p2g.update(meta.dropna(subset=['peptide_id']).drop_duplicates('peptide_id').set_index('peptide_id')['ensembl_id'].to_dict())
    if 'peptide_id_version' in meta.columns:
        p2g.update(meta.dropna(subset=['peptide_id_version']).drop_duplicates('peptide_id_version').set_index('peptide_id_version')['ensembl_id'].to_dict())
    
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
            # Sort by chromosome, then physical start coordinate
            coords = coords.sort_values(['chr', 'start']).reset_index(drop=True)
            species_coords[sp_code] = coords
            
    return species_coords

def get_neighborhoods_chromosome_aware(coords_df, og_col, window_size=5):
    """
    Constructs syntenic neighborhoods strictly within chromosome boundaries.
    Uses native Python lists for orders-of-magnitude faster execution.
    """
    neighborhoods = {}
    # Group by chromosome so sliding window NEVER bleeds into adjacent chromosomes
    for chr_name, group in coords_df.groupby('chr', sort=False):
        genes = group['g_id'].tolist()
        ogs = group[og_col].tolist()
        n_genes = len(genes)
        
        for i in range(n_genes):
            start_idx = max(0, i - window_size)
            end_idx = min(n_genes, i + window_size + 1)
            # Collect non-null, assigned orthogroups in the window
            win_ogs = [og for og in ogs[start_idx:end_idx] if pd.notna(og) and og != 'unassigned']
            neighborhoods[genes[i]] = set(win_ogs)
            
    return neighborhoods

def calculate_pairwise_goc(n1, n2):
    """
    Computes Jaccard index |n1 ∩ n2| / |n1 ∪ n2| without allocating union set.
    """
    if not n1 or not n2:
        return 0.0
    inter_len = len(n1.intersection(n2))
    if inter_len == 0:
        return 0.0
    union_len = len(n1) + len(n2) - inter_len
    return inter_len / union_len

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

def process_species_pair(sp1, sp2, master, ensembl_scores, og_neighborhoods, hq_neighborhoods, output_pairwise_dir):
    df1 = master[master['species_code'] == sp1]
    df2 = master[master['species_code'] == sp2]
    
    merged = pd.merge(df1, df2, on='ens_orthogroup_id', suffixes=('_'+sp1, '_'+sp2))
    if merged.empty:
        return sp1, sp2, 0

    g1_list = merged[f'ensembl_id_{sp1}'].tolist()
    g2_list = merged[f'ensembl_id_{sp2}'].tolist()
    
    # 1. Fast Ensembl score lookup (replaces slow apply(axis=1))
    identities = []
    homology_identities = []
    goc_scores = []
    wga_coverages = []
    is_high_confidences = []
    
    for g1, g2 in zip(g1_list, g2_list):
        s = ensembl_scores.get((g1, g2)) or ensembl_scores.get((g2, g1))
        if s:
            identities.append(s.get('identity'))
            homology_identities.append(s.get('homology_identity'))
            goc_scores.append(s.get('goc_score'))
            wga_coverages.append(s.get('wga_coverage'))
            is_high_confidences.append(s.get('is_high_confidence'))
        else:
            identities.append(None)
            homology_identities.append(None)
            goc_scores.append(None)
            wga_coverages.append(None)
            is_high_confidences.append(None)
            
    merged['ens_identity'] = identities
    merged['ens_homology_identity'] = homology_identities
    merged['ens_goc_score'] = goc_scores
    merged['ens_wga_coverage'] = wga_coverages
    merged['ens_is_high_confidence'] = is_high_confidences
    
    # 2. Fast Pairwise Atlas GOC scores (zero set allocation, chromosome-aware)
    sp1_og_neigh = og_neighborhoods.get(sp1, {})
    sp2_og_neigh = og_neighborhoods.get(sp2, {})
    sp1_hq_neigh = hq_neighborhoods.get(sp1, {})
    sp2_hq_neigh = hq_neighborhoods.get(sp2, {})
    
    merged['atlas_goc_score_pairwise'] = [
        calculate_pairwise_goc(sp1_og_neigh.get(g1), sp2_og_neigh.get(g2))
        for g1, g2 in zip(g1_list, g2_list)
    ]
    merged['atlas_hq_goc_score_pairwise'] = [
        calculate_pairwise_goc(sp1_hq_neigh.get(g1), sp2_hq_neigh.get(g2))
        for g1, g2 in zip(g1_list, g2_list)
    ]

    # 3. Ortholog Ranking
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
    
    out_path = os.path.join(output_pairwise_dir, f"{sp1}_{sp2}_pairwise_orthologs.tsv")
    merged.to_csv(out_path, sep='\t', index=False)
    return sp1, sp2, len(merged)

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
        
        # Load Coords and build chromosome-aware neighborhoods for informational GOC
        species_coords = load_coords(master, coords_dir, metadata_path)
        print("Constructing chromosome-aware syntenic neighborhoods...")
        og_neighborhoods = {sp: get_neighborhoods_chromosome_aware(df, 'og', window_size) for sp, df in species_coords.items()}
        hq_neighborhoods = {sp: get_neighborhoods_chromosome_aware(df, 'hq_og', window_size) for sp, df in species_coords.items()}
        
        # Load Ensembl scores for pairwise tables
        ensembl_scores = load_ensembl_scores(homology_dir, species_config)

        # Generate Refined Pairwise Tables (Multi-threaded across pairs)
        print("Generating pairwise tables with Ensembl and Atlas scores (vectorized & parallel)...")
        species_codes = list(species_config.keys())
        pairs = list(itertools.combinations(species_codes, 2))
        
        max_workers = min(4, len(pairs))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(process_species_pair, sp1, sp2, master, ensembl_scores, og_neighborhoods, hq_neighborhoods, output_pairwise_dir)
                for sp1, sp2 in pairs
            ]
            for future in futures:
                sp1, sp2, count = future.result()
                print(f"  [Pairwise] {sp1}_{sp2}: {count} ortholog pairs processed.")
                
        print("Synteny refinement complete.")
            
    except Exception as e:
        print(f"Error in main: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
