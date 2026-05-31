import pandas as pd
from collections import Counter

def generate_symbols(df, group_col, symbol_col_name, ref_sp=None):
    assigned = df[df[group_col] != 'unassigned'].copy()
    assigned['symbol_upper'] = assigned['gene_symbol'].fillna('').astype(str).str.upper()
    group_names = {}
    
    for group_id, group in assigned.groupby(group_col):
        symbols = [s for s in group['symbol_upper'].unique() if s and s != 'NAN']
        if not symbols:
            consensus = f"FAM_{group_id}"
        else:
            ref_symbols = []
            if ref_sp:
                ref_symbols = group[group['species_code'] == ref_sp]['symbol_upper'].dropna().unique()
                ref_symbols = [s for s in ref_symbols if s and s != 'NAN']
            
            if ref_symbols:
                consensus = ref_symbols[0]
            else:
                counts = Counter(group['symbol_upper'].dropna().tolist())
                if '' in counts: del counts['']
                if 'NAN' in counts: del counts['NAN']
                if counts:
                    consensus = counts.most_common(1)[0][0]
                else:
                    consensus = sorted(symbols)[0]
        group_names[group_id] = consensus

    df[symbol_col_name] = df[group_col].map(group_names).fillna('unassigned')
    return df

def main():
    master_path = snakemake.input.refined_master
    output_path = snakemake.output.final_master
    ref_sp = snakemake.config.get('symbol_reference_species', None)
    
    df = pd.read_csv(master_path, sep='\t')
    
    print(f"Generating symbols for Ensembl Orthogroups (Ref Species: {ref_sp})...")
    df = generate_symbols(df, 'ens_orthogroup_id', 'ens_orthogroup_symbol', ref_sp=ref_sp)
    df = generate_symbols(df, 'ens_hqorthogroup_id', 'ens_hqorthogroup_symbol', ref_sp=ref_sp)
    
    # Optional post-processing refinements (Trims)
    refinement_cols = []
    if 'strict_trim' in df.columns:
        print("Generating symbols for Strict Trims...")
        df = generate_symbols(df, 'strict_trim', 'strict_trim_symbol', ref_sp=ref_sp)
        refinement_cols += ['strict_trim', 'strict_trim_symbol', 'strict_trim_status']

    if 'lifeline_trim' in df.columns:
        print("Generating symbols for Lifeline Trims...")
        df = generate_symbols(df, 'lifeline_trim', 'lifeline_trim_symbol', ref_sp=ref_sp)
        refinement_cols += ['lifeline_trim', 'lifeline_trim_symbol', 'lifeline_trim_status']

    if 'lifeline_tree_trim' in df.columns:
        print("Generating symbols for Lifeline Tree Trims...")
        df = generate_symbols(df, 'lifeline_tree_trim', 'lifeline_tree_trim_symbol', ref_sp=ref_sp)
        refinement_cols += ['lifeline_tree_trim', 'lifeline_tree_trim_symbol', 'lifeline_tree_trim_status']

    # Final logic ordering
    base_cols = [
        'ensembl_id', 'gene_symbol', 'species_code', 
        'ens_orthogroup_id', 'ens_orthogroup_symbol', 'ens_orthogroup_status',
        'ens_hqorthogroup_id', 'ens_hqorthogroup_symbol', 'ens_hqorthogroup_status',
        'n_species'
    ]
    
    cols = base_cols + refinement_cols
    df = df[cols]

    df.to_csv(output_path, sep='\t', index=False)
    print(f"Final tables generated at {output_path}")

if __name__ == "__main__":
    main()
