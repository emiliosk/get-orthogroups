import os
import yaml
import requests
from pathlib import Path

# Load config
config_path = "config.yaml"
if not os.path.exists(config_path):
    config_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "config.yaml")

with open(config_path, "r") as f:
    config = yaml.safe_load(f)

VERSION = config.get("ensembl_version", "109")
SPECIES_CONFIG = config["species"]

def get_ensembl_assembly(scientific_name):
    """Fetch default assembly version from Ensembl REST API if not provided in config."""
    url = f"https://rest.ensembl.org/info/assembly/{scientific_name}?content-type=application/json"
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            data = r.json()
            assembly = data.get("default_coord_system_version")
            if assembly:
                print(f"  [REST] Auto-detected assembly for {scientific_name}: {assembly}")
                return assembly
    except Exception as e:
        print(f"  [Warning] Could not fetch assembly via REST API for {scientific_name}: {e}")
    return None

import gzip
import time

def download_file(url, dest, max_retries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        # Verify gzip integrity if it's a .gz file
        if dest.endswith('.gz'):
            try:
                with gzip.open(dest, 'rb') as gz:
                    gz.read(1024)
                print(f"File already exists and is valid: {dest}")
                return True
            except Exception:
                print(f"Existing file {dest} is corrupted. Re-downloading...")
                try:
                    os.remove(dest)
                except Exception:
                    pass
        else:
            print(f"File already exists: {dest}")
            return True
    
    tmp_dest = f"{dest}.part"
    print(f"Downloading {url} to {dest}...")
    
    headers = {"User-Agent": "Mozilla/5.0 (EnsemblOrthologyPipeline/1.0)"}
    for attempt in range(1, max_retries + 1):
        try:
            with requests.get(url, stream=True, timeout=60, headers=headers) as r:
                r.raise_for_status()
                with open(tmp_dest, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
            
            # Verify gzip integrity before accepting file
            if dest.endswith('.gz'):
                with gzip.open(tmp_dest, 'rb') as gz:
                    gz.read(1024)
            
            # Atomic rename upon verified download
            os.replace(tmp_dest, dest)
            return True
        except Exception as e:
            print(f"  [Attempt {attempt}/{max_retries}] Failed downloading {url}: {e}")
            if os.path.exists(tmp_dest):
                try:
                    os.remove(tmp_dest)
                except Exception:
                    pass
            if attempt == max_retries:
                raise RuntimeError(f"Failed to download {url} after {max_retries} attempts: {e}")
            time.sleep(2 * attempt)

def format_ensembl_folder(sci_name):
    """Formats species scientific name for Ensembl directory (e.g. homo_sapiens -> Homo_sapiens)."""
    parts = sci_name.split('_')
    return parts[0].capitalize() + ('_' + '_'.join(parts[1:]) if len(parts) > 1 else '')

def main():
    # 1. Download Homologies
    homology_dir = Path(config["paths"]["raw_homologies_dir"])
    homology_dir.mkdir(parents=True, exist_ok=True)
    
    for sp_key, info in SPECIES_CONFIG.items():
        sci_name = info["scientific_name"]
        url = f"https://ftp.ensembl.org/pub/release-{VERSION}/tsv/ensembl-compara/homologies/{sci_name}/Compara.{VERSION}.protein_default.homologies.tsv.gz"
        dest = homology_dir / f"{sci_name}.tsv.gz"
        download_file(url, str(dest))

    # 2. Download Protein FASTAs (raw)
    fasta_raw_dir = Path("input/DB/raw_fastas")
    fasta_raw_dir.mkdir(parents=True, exist_ok=True)
    
    # 3. Download GFF3 (for coordinates)
    gff_dir = Path("input/DB/gff3")
    gff_dir.mkdir(parents=True, exist_ok=True)

    for sp_key, info in SPECIES_CONFIG.items():
        sci_name = info["scientific_name"]
        folder = sci_name.lower()
        formatted_sp = format_ensembl_folder(sci_name)
        
        assembly = info.get("assembly")
        if not assembly:
            assembly = get_ensembl_assembly(sci_name)
            if not assembly:
                raise ValueError(f"Could not determine assembly for {sp_key} ({sci_name}). Please add 'assembly: <NAME>' to config.yaml.")
        
        prefix = f"{formatted_sp}.{assembly}"
        
        # FASTA download
        fasta_url = f"https://ftp.ensembl.org/pub/release-{VERSION}/fasta/{folder}/pep/{prefix}.pep.all.fa.gz"
        dest_fasta = fasta_raw_dir / f"{sp_key}.fa.gz"
        download_file(fasta_url, str(dest_fasta))

        # GFF3 download
        gff_url = f"https://ftp.ensembl.org/pub/release-{VERSION}/gff3/{folder}/{prefix}.{VERSION}.gff3.gz"
        dest_gff = gff_dir / f"{sp_key}.gff3.gz"
        download_file(gff_url, str(dest_gff))

if __name__ == "__main__":
    main()

