import pandas as pd
import os
import yaml
import sys

def main():
    lines = []
    lines.append("=== Orthology Pipeline Summary Statistics ===\n")
    
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print("Error: config.yaml not found.")
        sys.exit(1)
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    out_dir = config.get('output_dir', 'ensembl_pipeline_output')
    master_path = os.path.join(out_dir, "Consensus_Master.tsv")
    
    if not os.path.exists(master_path):
        print(f"Error: Master table not found at {master_path}")
        sys.exit(1)
        
    df = pd.read_csv(master_path, sep='\t')
    
    lines.append(f"Total Genes: {len(df)}")
    lines.append(f"Total Species: {df['species_code'].nunique()} ({', '.join(df['species_code'].unique())})")
    
    lines.append("\n--- Orthogroup Stats (Standard Track) ---")
    og_stats = df.drop_duplicates('ens_orthogroup_id')['ens_orthogroup_status'].value_counts()
    lines.append(f"Total Orthogroups: {df['ens_orthogroup_id'].nunique()}")
    for status, count in og_stats.items():
        lines.append(f"  {status:12}: {count}")
        
    lines.append("\n--- Orthogroup Stats (High-Confidence Track) ---")
    hq_stats = df.drop_duplicates('ens_hqorthogroup_id')['ens_hqorthogroup_status'].value_counts()
    lines.append(f"Total HQ Groups: {df['ens_hqorthogroup_id'].nunique()}")
    for status, count in hq_stats.items():
        lines.append(f"  {status:12}: {count}")

    lines.append("\n--- Species Coverage (in Standard OGs) ---")
    coverage = df.groupby('ens_orthogroup_id')['species_code'].nunique().value_counts().sort_index(ascending=False)
    for n_sp, count in coverage.items():
        lines.append(f"  Groups with {n_sp} species: {count}")

    # Pairwise Ranking Check
    pairwise_dir = os.path.join(out_dir, "pairwise_tables")
    if os.path.exists(pairwise_dir):
        n_tables = len([f for f in os.listdir(pairwise_dir) if f.endswith('.tsv')])
        lines.append(f"\n--- Exported Pairwise Tables: {n_tables} ---")

    lines.append("\n" + "="*40)
    
    # Print to console
    output_text = "\n".join(lines)
    print(output_text)
    
    # Save to file
    out_file = os.path.join(out_dir, "Summary_Stats.md")
    os.makedirs(out_dir, exist_ok=True)
    with open(out_file, 'w') as f:
        f.write("# Orthology Pipeline Summary Statistics\n\n```text\n" + output_text + "\n```\n")
    print(f"Summary stats saved to {out_file}")

if __name__ == "__main__":
    main()
