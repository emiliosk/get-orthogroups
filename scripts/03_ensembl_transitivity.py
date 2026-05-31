import pandas as pd
import networkx as nx
import gzip
import os

def main():
    meta_path = snakemake.input.metadata
    homology_files = snakemake.input.homologies
    output_path = snakemake.output.baseline
    species_config = snakemake.config['species']
    
    # Create mapping: scientific_name -> CODE (e.g., homo_sapiens -> HUMAN)
    sci_to_code = {info['scientific_name']: code for code, info in species_config.items()}
    # Create mapping: standard_name (atlas_name) -> CODE (e.g., human -> HUMAN)
    standard_to_code = {info['atlas_name']: code for code, info in species_config.items()}

    # 1. Load metadata
    print("Loading project metadata...")
    meta = pd.read_csv(meta_path)
    pc_gene_ids = set(meta['ensembl_id'].tolist())
    print(f"  Loaded {len(pc_gene_ids)} target protein-coding genes.")

    # 2. Build the all-vs-all graphs
    print("Building all-vs-all orthology networks...")
    G_all = nx.Graph()
    G_hq = nx.Graph()
    G_all.add_nodes_from(pc_gene_ids)
    G_hq.add_nodes_from(pc_gene_ids)
    
    total_edges_all = 0
    total_edges_hq = 0
    for file_path in homology_files:
        print(f"  Processing {os.path.basename(file_path)}...")
        with gzip.open(file_path, 'rt') as f:
            header = f.readline().strip().split('\t')
            col_source = header.index('gene_stable_id')
            col_type = header.index('homology_type')
            col_target = header.index('homology_gene_stable_id')
            col_target_sp = header.index('homology_species')
            col_conf = header.index('is_high_confidence') if 'is_high_confidence' in header else -1
            
            line_count = 0
            for line in f:
                line_count += 1
                if line_count % 1000000 == 0:
                    print(f"    ... processed {line_count // 1000000}M lines")
                
                parts = line.strip().split('\t')
                if len(parts) <= col_target_sp: continue
                
                s_id = parts[col_source]
                type_str = parts[col_type]
                t_id = parts[col_target]
                t_sp_str = parts[col_target_sp]
                is_hq = (col_conf != -1 and str(parts[col_conf]) == '1')
                
                # Check if target species is in our list
                if t_sp_str in sci_to_code:
                    # STRICT: Only take orthology relationships
                    if type_str.startswith('ortholog_'):
                        if s_id in pc_gene_ids and t_id in pc_gene_ids:
                            G_all.add_edge(s_id, t_id)
                            total_edges_all += 1
                            if is_hq:
                                G_hq.add_edge(s_id, t_id)
                                total_edges_hq += 1
                            
    print(f"  Total edges (All): {total_edges_all}")
    print(f"  Total edges (HQ):  {total_edges_hq}")
    
    # 3. Clustering (Connected Components)
    print("Clustering into inclusive HOGs...")
    
    def get_og_df(G, prefix):
        components = list(nx.connected_components(G))
        
        # STABILITY FIX: Sort components by (size descending, then by the sorted string of gene IDs)
        # This ensures IDs (OG_00001, etc.) are deterministic even if graph traversal varies.
        def component_key(c):
            return (-len(c), sorted(list(c)))
        
        components.sort(key=component_key)
        
        print(f"  Found {len(components)} unbiased Ensembl groups for {prefix}.")
        results = []
        
        # Use separate counters for regular groups and singletons
        group_count = 1
        singleton_count = 1
        
        for component in components:
            if len(component) > 1:
                og_id = f"{prefix}_{group_count:05d}"
                group_count += 1
            else:
                og_id = f"{prefix}_SNG_{singleton_count:06d}"
                singleton_count += 1
                
            for gene_id in component:
                results.append({'ensembl_id': gene_id, f'{prefix}_id': og_id})
        return pd.DataFrame(results)

    og_all_df = get_og_df(G_all, "OG")
    og_hq_df = get_og_df(G_hq, "HQ")
    
    # 4. Map back to table
    og_merged = pd.merge(og_all_df, og_hq_df, on='ensembl_id', how='outer')
    
    # 5. Merge with metadata and handle singletons
    meta_subset = meta[['ensembl_id', 'species', 'gene_symbol']].drop_duplicates()
    final_df = meta_subset.merge(og_merged, on='ensembl_id', how='left')
    
    def handle_singletons(df, col_name, prefix):
        mask = df[col_name].isna()
        if mask.any():
            singleton_count = 1
            new_ogs = df[col_name].tolist()
            for idx in range(len(new_ogs)):
                if pd.isna(new_ogs[idx]):
                    new_ogs[idx] = f"{prefix}_SNG_{singleton_count:06d}"
                    singleton_count += 1
            df[col_name] = new_ogs
        return df

    final_df = handle_singletons(final_df, 'OG_id', 'OG')
    final_df = handle_singletons(final_df, 'HQ_id', 'HQ')
    
    # Rename columns to match requested names
    final_df = final_df.rename(columns={'OG_id': 'ens_orthogroup_id', 'HQ_id': 'ens_hqorthogroup_id'})
        
    final_df['species_code'] = final_df['species'].map(standard_to_code).fillna(final_df['species'])

    # 6. Save
    print(f"Saving baseline to {output_path}...")
    final_df.to_csv(output_path, sep='\t', index=False)
    print("Done!")

if __name__ == "__main__":
    main()
