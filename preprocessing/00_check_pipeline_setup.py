import os
import sys
import shutil
import subprocess
import yaml
import requests
import pandas as pd

def check_tools(config):
    print("--- 1. Checking Tools ---")
    compute_trees = config.get("compute_trees", False)
    all_ok = True
    for tool, path in config.get('tools', {}).items():
        if os.path.isabs(path):
            if os.path.exists(path) and os.access(path, os.X_OK):
                print(f"[OK] {tool:10}: {path}")
            else:
                if not compute_trees and tool == "treerecs":
                    print(f"[INFO] {tool:10}: Not found at '{path}', but compute_trees is false (only needed for Rule 7).")
                else:
                    print(f"[FAIL] {tool:10}: Not found or not executable at {path}")
                    all_ok = False
        else:
            env_var = f"{tool.upper()}_BIN"
            found = shutil.which(path) or os.environ.get(env_var)
            if not found and tool == "treerecs" and os.path.exists("/Users/emilioskarwan/anaconda3/envs/intel_phylo/bin/treerecs"):
                found = "/Users/emilioskarwan/anaconda3/envs/intel_phylo/bin/treerecs"
            if found:
                print(f"[OK] {tool:10}: Found ({found})")
            else:
                if not compute_trees and tool == "treerecs":
                    print(f"[INFO] {tool:10}: Not found in PATH, but compute_trees is false (only needed for Rule 7).")
                else:
                    print(f"[FAIL] {tool:10}: Not found in PATH. Check env or config.")
                    all_ok = False
    return all_ok

def check_ensembl():
    print("\n--- 2. Checking Ensembl Connectivity ---")
    url = "https://rest.ensembl.org/info/ping"
    try:
        r = requests.get(url, headers={"Content-Type": "application/json"})
        if r.status_code == 200 and r.json().get('ping') == 1:
            print("[OK] Ensembl REST API is alive.")
            return True
        else:
            print(f"[FAIL] Ensembl REST API returned status {r.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] Could not connect to Ensembl: {e}")
        return False

def check_metadata(config):
    print("\n--- 3. Checking Input Metadata ---")
    path = config['paths']['metadata']
    if not os.path.exists(path):
        print(f"[FAIL] Metadata file not found: {path}")
        return False
    
    try:
        df = pd.read_csv(path)
        gene_col = 'ensembl_id' if 'ensembl_id' in df.columns else ('ensembl_gene_id' if 'ensembl_gene_id' in df.columns else None)
        if not gene_col:
            print("[FAIL] Metadata missing gene ID column ('ensembl_id' or 'ensembl_gene_id')")
            return False
        if 'species' not in df.columns:
            print("[FAIL] Metadata missing required column: 'species'")
            return False
        
        n_unique_genes = df[gene_col].nunique()
        species_list = df['species'].unique().tolist()
        
        extras = []
        if 'ensembl_transcript_id' in df.columns or 'transcript_id' in df.columns:
            extras.append("transcripts")
        if 'transcript_is_canonical' in df.columns:
            canon_count = df['transcript_is_canonical'].astype(str).str.strip().isin(['1', '1.0', 'True', 'true']).sum()
            extras.append(f"{canon_count} canonical")
        if 'peptide_id' in df.columns or 'peptide_id_version' in df.columns:
            extras.append("peptides")
            
        extra_info = f" [{', '.join(extras)}]" if extras else ""
        print(f"[OK] Metadata valid ({len(df)} rows, {n_unique_genes} unique genes, species: {species_list}){extra_info}")
        return True
    except Exception as e:
        print(f"[FAIL] Error reading metadata: {e}")
        return False

def check_species_config(config):
    print("\n--- 4. Checking Species & Assembly Setup ---")
    species_map = config.get("species", {})
    if not species_map:
        print("[FAIL] No species defined in config.yaml under 'species'")
        return False
    
    all_ok = True
    for code, info in species_map.items():
        sci_name = info.get("scientific_name")
        if not sci_name:
            print(f"[FAIL] Species '{code}' is missing 'scientific_name'")
            all_ok = False
            continue
        
        url = f"https://rest.ensembl.org/info/assembly/{sci_name}?content-type=application/json"
        success = False
        for attempt in range(2):
            try:
                r = requests.get(url, timeout=15)
                if r.status_code == 200:
                    rest_assembly = r.json().get("default_coord_system_version")
                    cfg_assembly = info.get("assembly")
                    print(f"[OK] {code:8} ({sci_name}): REST assembly '{rest_assembly}'" + (f" (Config: '{cfg_assembly}')" if cfg_assembly else " (Auto-detect mode)"))
                    success = True
                    break
            except Exception:
                pass
        if not success:
            cfg_assembly = info.get("assembly")
            if cfg_assembly:
                print(f"[OK] {code:8} ({sci_name}): Using config assembly '{cfg_assembly}' (REST lookup skipped/timed out)")
            else:
                print(f"[WARN] Could not reach Ensembl REST for '{sci_name}' to auto-detect assembly.")
            
    # Check species tree (either inline Newick string or file path)
    tree_val = config.get("species_tree") or config.get("paths", {}).get("species_tree")
    if tree_val:
        if isinstance(tree_val, str) and tree_val.strip().startswith("("):
            # Inline Newick string
            import re
            tips = re.findall(r'([A-Za-z0-9_]+)(?::[0-9.]+)?[,\)]', tree_val)
            # Filter out numbers
            tips = [t for t in tips if not t.replace('.', '', 1).isdigit()]
            print(f"[OK] Species tree: Valid inline Newick tree (tips: {tips})")
            missing_tips = [sp for sp in species_map.keys() if sp not in tips]
            if missing_tips:
                print(f"[WARN] Configured species not found in species tree tips: {missing_tips}")
        elif os.path.exists(tree_val):
            print(f"[OK] Species tree file found: {tree_val}")
        else:
            print(f"[WARN] Species tree file not found at '{tree_val}'. (Tree reconciliation requires a valid species tree).")
    else:
        print("[WARN] No species_tree defined in config.yaml.")
            
    return all_ok

def main():
    print("=== Pipeline Pre-flight Setup Check ===\n")
    
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    if not os.path.exists(config_path):
        config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), config_path)
        
    if not os.path.exists(config_path):
        print("Error: config.yaml not found.")
        sys.exit(1)
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    t = check_tools(config)
    e = check_ensembl()
    m = check_metadata(config)
    s = check_species_config(config)
    
    print("\n" + "="*40)
    if all([t, e, m, s]):
        print("RESULT: Pipeline configuration is valid and ready to run!")
    else:
        print("RESULT: Issues detected. Please fix the failures above.")
        sys.exit(1)

if __name__ == "__main__":
    main()

