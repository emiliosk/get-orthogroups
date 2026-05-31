import os
import gzip
import pandas as pd
import yaml
from Bio import SeqIO
import re
from pathlib import Path

def main(input_transcripts, raw_fasta_dir, raw_gff_dir, fasta_out_dir, coords_out_dir, species_config):
    print("Loading transcript metadata...")
    trans_df = pd.read_csv(input_transcripts)
    # Ensure IDs are strings and stripped
    trans_df['ensembl_transcript_id'] = trans_df['ensembl_transcript_id'].astype(str).str.strip('"')
    canonical_transcripts = set(trans_df[trans_df['transcript_is_canonical'] == 1]['ensembl_transcript_id'])

    # 1. Filter FASTAs
    Path(fasta_out_dir).mkdir(parents=True, exist_ok=True)
    
    for code in species_config.keys():
        raw_path = Path(raw_fasta_dir) / f"{code}.fa.gz"
        out_path = Path(fasta_out_dir) / f"{code}.fa"
        
        if not raw_path.exists():
            print(f"Skipping {code} FASTA: Raw file not found at {raw_path}")
            continue

        print(f"Filtering {code} FASTA to canonical transcripts...")
        records = []
        with gzip.open(raw_path, "rt") as f:
            for record in SeqIO.parse(f, "fasta"):
                t_match = re.search(r'transcript:([\w.]+)', record.description)
                if t_match:
                    t_id = t_match.group(1).split('.')[0]
                    if t_id in canonical_transcripts:
                        records.append(record)
                else:
                    t_id = record.id.split('.')[0]
                    if t_id in canonical_transcripts:
                        records.append(record)
        
        with open(out_path, "w") as out_f:
            SeqIO.write(records, out_f, "fasta")
        print(f"  Saved {len(records)} sequences to {out_path}")

    # 2. Extract Coordinates (BED) from GFF3
    Path(coords_out_dir).mkdir(parents=True, exist_ok=True)
    
    for code in species_config.keys():
        gff_path = Path(raw_gff_dir) / f"{code}.gff3.gz"
        out_path = Path(coords_out_dir) / f"{code}_coords.bed"
        
        if not gff_path.exists():
            print(f"Skipping {code} GFF: Raw file not found at {gff_path}")
            continue

        print(f"Parsing {code} GFF3 to BED...")
        coords = []
        with gzip.open(gff_path, "rt") as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9: continue
                if parts[2] == 'CDS':
                    attr = parts[8]
                    p_match = re.search(r'ID=CDS:([\w.]+)', attr)
                    if p_match:
                        p_id = p_match.group(1)
                        coords.append([parts[0], int(parts[3]), int(parts[4]), p_id])
        
        if not coords:
            print(f"  Warning: No CDS found in {gff_path}")
            continue

        df = pd.DataFrame(coords, columns=['chr', 'start', 'end', 'p_id'])
        bed_df = df.groupby('p_id').agg({'chr': 'first', 'start': 'min', 'end': 'max'}).reset_index()
        bed_df[['chr', 'start', 'end', 'p_id']].to_csv(out_path, sep='\t', header=False, index=False)
        print(f"  Saved {len(bed_df)} proteins to {out_path}")

if __name__ == "__main__":
    # If run via Snakemake
    if 'snakemake' in globals():
        main(
            input_transcripts = snakemake.input.transcripts,
            raw_fasta_dir = "input/DB/raw_fastas",
            raw_gff_dir = "input/DB/gff3",
            fasta_out_dir = snakemake.config["paths"]["fasta_dir"],
            coords_out_dir = snakemake.config["paths"]["genespace_coords"],
            species_config = snakemake.config["species"]
        )
    else:
        # Fallback for manual run
        with open("config.yaml", "r") as f:
            config = yaml.safe_load(f)
        version = config.get("ensembl_version", "109")
        main(
            input_transcripts = f"input/DB/ensembl_v{version}_multispecies_transcripts.csv",
            raw_fasta_dir = "input/DB/raw_fastas",
            raw_gff_dir = "input/DB/gff3",
            fasta_out_dir = config["paths"]["fasta_dir"],
            coords_out_dir = config["paths"]["genespace_coords"],
            species_config = config["species"]
        )
