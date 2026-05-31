import pandas as pd
import requests
import json
import time
import os
import yaml

def fetch_canonical(gene_ids):
    """Fetch canonical transcript IDs and versions for a list of genes using Ensembl REST API."""
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    results = []
    
    # Process in chunks of 1000 (Ensembl API limit)
    for i in range(0, len(gene_ids), 1000):
        chunk = gene_ids[i:i+1000]
        data = {"ids": chunk}
        print(f"  Fetching chunk {i//1000 + 1} ({len(chunk)} genes)...")
        
        for attempt in range(3):
            try:
                response = requests.post(url, headers=headers, json=data, params={"expand": 1})
                if response.status_code == 200:
                    batch_data = response.json()
                    for gid, gene_info in batch_data.items():
                        if gene_info and 'Transcript' in gene_info:
                            for t in gene_info['Transcript']:
                                if t.get('is_canonical') == 1:
                                    results.append({
                                        'ensembl_gene_id': gid,
                                        'ensembl_transcript_id': t['id'],
                                        'ensembl_transcript_id_version': f"{t['id']}.{t['version']}",
                                        'transcript_is_canonical': 1
                                    })
                    break
                else:
                    print(f"    Error {response.status_code}, retrying...")
                    time.sleep(2)
            except Exception as e:
                print(f"    Request failed: {e}, retrying...")
                time.sleep(2)
    return pd.DataFrame(results)

def main():
    # Load config to get version
    pipeline_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_path = os.path.join(pipeline_dir, "config.yaml")
    
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    
    version = config.get("ensembl_version", "109")
    meta_path = os.path.join(pipeline_dir, config["paths"]["metadata"])
    out_path = os.path.join(pipeline_dir, f"input/DB/ensembl_v{version}_multispecies_transcripts.csv")
    
    print(f"Loading metadata from {meta_path}...")
    meta = pd.read_csv(meta_path)
    
    # Ensure mandatory columns exist
    if 'ensembl_id' not in meta.columns or 'species' not in meta.columns:
        print("Error: Metadata must contain 'ensembl_id' and 'species' columns.")
        return

    gene_ids = meta['ensembl_id'].unique().tolist()
    gene_to_species = meta[['ensembl_id', 'species']].drop_duplicates().set_index('ensembl_id')['species'].to_dict()

    print(f"Fetching canonical transcripts for {len(gene_ids)} genes via Ensembl REST API...")
    trans_df = fetch_canonical(gene_ids)

    if trans_df.empty:
        print("Error: No canonical transcripts found.")
        return

    # Add species information from original metadata
    trans_df['species'] = trans_df['ensembl_gene_id'].map(gene_to_species)

    # Export
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    trans_df.to_csv(out_path, index=False)
    print(f"Successfully saved {len(trans_df)} canonical transcripts to {out_path}")

if __name__ == "__main__":
    main()
