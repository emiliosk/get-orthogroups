import os
import sys
import shutil
import subprocess
import yaml
import requests
import pandas as pd

def check_tools(config):
    print("--- 1. Checking Tools ---")
    all_ok = True
    for tool, path in config['tools'].items():
        if os.path.isabs(path):
            if os.path.exists(path) and os.access(path, os.X_OK):
                print(f"[OK] {tool:10}: {path}")
            else:
                print(f"[FAIL] {tool:10}: Not found or not executable at {path}")
                all_ok = False
        else:
            found = shutil.which(path)
            if found:
                print(f"[OK] {tool:10}: Found in PATH ({found})")
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
        required = ['ensembl_id', 'species']
        missing = [c for c in required if c not in df.columns]
        if missing:
            print(f"[FAIL] Metadata missing required columns: {missing}")
            return False
        print(f"[OK] Metadata valid ({len(df)} rows, species: {df['species'].unique().tolist()})")
        return True
    except Exception as e:
        print(f"[FAIL] Error reading metadata: {e}")
        return False

def main():
    print("=== Pipeline Pre-flight Setup Check ===\n")
    
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print("Error: config.yaml not found in current directory.")
        sys.exit(1)
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    t = check_tools(config)
    e = check_ensembl()
    m = check_metadata(config)
    
    print("\n" + "="*40)
    if all([t, e, m]):
        print("RESULT: Pipeline is ready to run!")
    else:
        print("RESULT: Issues detected. Please fix the failures above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
