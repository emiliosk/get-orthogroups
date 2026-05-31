import os
import subprocess
import yaml
from pathlib import Path

# Load config
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

VERSION = config.get("ensembl_version", "109")
SPECIES_CONFIG = config["species"]

def download_file(url, dest):
    if os.path.exists(dest):
        print(f"File already exists: {dest}")
        return
    
    print(f"Downloading {url} to {dest}...")
    try:
        subprocess.run(["wget", "-q", "-O", dest, url], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to download {url}: {e}")

def main():
    # 1. Download Homologies
    homology_dir = Path(config["paths"]["raw_homologies_dir"])
    homology_dir.mkdir(parents=True, exist_ok=True)
    
    for sp_key, info in SPECIES_CONFIG.items():
        sci_name = info["scientific_name"]
        # Ensembl 109 TSV homologies:
        url = f"https://ftp.ensembl.org/pub/release-{VERSION}/tsv/ensembl-compara/homologies/{sci_name}/Compara.{VERSION}.protein_default.homologies.tsv.gz"
        dest = homology_dir / f"{sci_name}.tsv.gz"
        download_file(url, str(dest))

    # 2. Download Protein FASTAs (raw)
    fasta_raw_dir = Path("input/DB/raw_fastas")
    fasta_raw_dir.mkdir(parents=True, exist_ok=True)
    
    for sp_key, info in SPECIES_CONFIG.items():
        folder = info["scientific_name"]
        assembly = info["assembly"]
        prefix = f"{folder.capitalize()}.{assembly}"
        
        url = f"https://ftp.ensembl.org/pub/release-{VERSION}/fasta/{folder}/pep/{prefix}.pep.all.fa.gz"
        dest = fasta_raw_dir / f"{sp_key}.fa.gz"
        download_file(url, str(dest))

    # 3. Download GFF3 (for coordinates)
    gff_dir = Path("input/DB/gff3")
    gff_dir.mkdir(parents=True, exist_ok=True)
    
    for sp_key, info in SPECIES_CONFIG.items():
        folder = info["scientific_name"]
        assembly = info["assembly"]
        prefix = f"{folder.capitalize()}.{assembly}"
        
        url = f"https://ftp.ensembl.org/pub/release-{VERSION}/gff3/{folder}/{prefix}.{VERSION}.gff3.gz"
        dest = gff_dir / f"{sp_key}.gff3.gz"
        download_file(url, str(dest))

if __name__ == "__main__":
    main()
