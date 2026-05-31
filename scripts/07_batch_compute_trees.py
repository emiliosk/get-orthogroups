import pandas as pd
import os
import subprocess
import shutil
import concurrent.futures
import time
import tarfile

def init_worker(fasta_dir, species_list):
    """Initializer for each worker process to load sequence data into memory."""
    global SEQ_DICT
    global GENE_TO_SPECIES
    SEQ_DICT = {}
    GENE_TO_SPECIES = {}
    
    for sp in species_list:
        fa_path = os.path.join(fasta_dir, f"{sp}.fa")
        if not os.path.exists(fa_path): continue
        
        current_id = None
        current_seq = []
        with open(fa_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_id: SEQ_DICT[current_id] = "".join(current_seq)
                    header = line[1:].strip()
                    if 'gene:' in header:
                        gid = header.split('gene:')[1].split()[0].split('.')[0]
                        current_id = gid
                    else:
                        current_id = header.split()[0].split('.')[0]
                    GENE_TO_SPECIES[current_id] = sp
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_id: SEQ_DICT[current_id] = "".join(current_seq)

def process_group(group_id, gene_ids, msa_dir, tree_dir, tmp_dir, mafft_bin, fasttree_bin, species_tree, treerecs_bin):
    """Worker function for a single group."""
    group_tmp = os.path.join(tmp_dir, group_id)
    try:
        os.makedirs(group_tmp, exist_ok=True)
        
        # 1. Fetch sequences
        found_genes = []
        raw_fa = os.path.join(group_tmp, f"{group_id}_raw.fa")
        with open(raw_fa, 'w') as f:
            for gid in gene_ids:
                seq = SEQ_DICT.get(gid)
                sp = GENE_TO_SPECIES.get(gid, "UNKNOWN")
                if seq:
                    f.write(f">{sp}_{gid}\n{seq}\n")
                    found_genes.append(gid)
        
        if len(found_genes) < 2:
            if os.path.exists(group_tmp): shutil.rmtree(group_tmp)
            return False
            
        # 2. Align with MAFFT
        aln_path = os.path.join(msa_dir, f"{group_id}.fa")
        subprocess.run(f"{mafft_bin} --auto --quiet {raw_fa} > {aln_path}", shell=True, check=True)
        
        # 3. Tree with FastTree and Reconcile with Treerecs
        if len(found_genes) >= 3:
            raw_tree = os.path.join(group_tmp, f"{group_id}_raw.nwk")
            subprocess.run(f"{fasttree_bin} -lg -quiet {aln_path} > {raw_tree}", shell=True, check=True)
            
            # Treerecs reconciliation
            rec_dir = os.path.join(group_tmp, "rec")
            os.makedirs(rec_dir, exist_ok=True)
            cmd = [treerecs_bin, "-g", raw_tree, "-s", species_tree, "-o", rec_dir, "--force", "--align-nodes"]
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            
            res_tree = os.path.join(rec_dir, f"{group_id}_raw.nwk_recs.nwk")
            if os.path.exists(res_tree):
                shutil.copy(res_tree, os.path.join(tree_dir, f"{group_id}.nwk"))
        
        if os.path.exists(group_tmp): shutil.rmtree(group_tmp)
        return True
    except Exception:
        if os.path.exists(group_tmp): shutil.rmtree(group_tmp)
        return False

def main():
    master_path = snakemake.input.master
    fasta_dir = snakemake.config['paths']['fasta_dir']
    species_tree = snakemake.config['paths']['species_tree']
    out_tar = snakemake.output.tar
    
    mafft_bin = snakemake.config['tools'].get('mafft', 'mafft')
    fasttree_bin = snakemake.config['tools'].get('fasttree', 'FastTree')
    treerecs_bin = snakemake.config['tools'].get('treerecs', 'treerecs')
    output_dir = snakemake.config.get('output_dir', 'ensembl_pipeline_output')
    
    out_dir = os.path.join(output_dir, "Phylogenetic_Trees")
    msa_dir = os.path.join(out_dir, "MSAs")
    tree_dir = os.path.join(out_dir, "Gene_Trees")
    tmp_dir = os.path.join(out_dir, "tmp_work")
    
    os.makedirs(msa_dir, exist_ok=True)
    os.makedirs(tree_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    master = pd.read_csv(master_path, sep='\t')
    groups = master[master['ens_orthogroup_id'] != 'unassigned'].groupby('ens_orthogroup_id')['ensembl_id'].apply(list).to_dict()
    
    eligible = {k: v for k, v in groups.items() if len(v) >= 2}
    total = len(eligible)
    print(f"Starting batch tree computation for {total} groups...")
    
    max_workers = os.cpu_count() - 1
    completed = 0
    errors = 0
    start_time = time.time()
    
    species_list = list(snakemake.config['species'].keys())
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, initializer=init_worker, initargs=(fasta_dir, species_list)) as executor:
        futures = {executor.submit(process_group, gid, genes, msa_dir, tree_dir, tmp_dir, mafft_bin, fasttree_bin, species_tree, treerecs_bin): gid for gid, genes in eligible.items()}
        for future in concurrent.futures.as_completed(futures):
            try:
                if future.result(): completed += 1
                else: errors += 1
            except:
                errors += 1
            
            if (completed + errors) % 500 == 0:
                elapsed = time.time() - start_time
                print(f"  Progress: {completed + errors}/{total} | Errors: {errors} | Elapsed: {elapsed/60:.1f}m")

    if os.path.exists(tmp_dir): shutil.rmtree(tmp_dir)
    print(f"Computation Complete. Success: {completed}, Failures: {errors}")
    
    print(f"Archiving results into {out_tar}...")
    # Add a sentinel file so Snakemake knows it's the right dir
    with tarfile.open(out_tar, "w:gz") as tar:
        tar.add(out_dir, arcname="Phylogenetic_Trees")
    
    # Clean up uncompressed dir to save inodes (DISABLED to keep files for visualization)
    # shutil.rmtree(out_dir)
    print("Done!")

if __name__ == "__main__":
    main()
