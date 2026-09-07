import pandas as pd
import os
import subprocess
import shutil
import concurrent.futures
import time
import tarfile
import sys
import yaml

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

def process_group(group_id, gene_ids, msa_dir, tree_dir, tmp_dir, mafft_bin, fasttree_bin, species_tree, treerecs_bin, contraction_threshold=0.8, reroot=True):
    """Worker function for a single orthogroup."""
    group_tmp = os.path.join(tmp_dir, group_id)
    aln_path = os.path.join(msa_dir, f"{group_id}.fa")
    tree_path = os.path.join(tree_dir, f"{group_id}.nwk")
    
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
            return (False, f"Fewer than 2 sequences found ({len(found_genes)})")
            
        # 2. Align with MAFFT (--anysymbol enables selenocysteine 'U' & non-standard residues)
        with open(aln_path, 'w') as out_f:
            res_mafft = subprocess.run(
                [mafft_bin, "--auto", "--anysymbol", "--quiet", raw_fa],
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True
            )
            if res_mafft.returncode != 0:
                if os.path.exists(aln_path):
                    os.remove(aln_path)
                return (False, f"MAFFT error: {res_mafft.stderr.strip()}")
                
        if not os.path.exists(aln_path) or os.path.getsize(aln_path) == 0:
            if os.path.exists(aln_path):
                os.remove(aln_path)
            return (False, "MAFFT produced empty alignment")

        # 3. Tree Inference & Species-Tree Reconciliation
        if len(found_genes) == 2:
            # 2-taxon tree is topologically trivial
            g1, g2 = found_genes[0], found_genes[1]
            sp1 = GENE_TO_SPECIES.get(g1, "UNKNOWN")
            sp2 = GENE_TO_SPECIES.get(g2, "UNKNOWN")
            with open(tree_path, 'w') as f:
                f.write(f"({sp1}_{g1}:0.1,{sp2}_{g2}:0.1);\n")
            return (True, None)

        # 3+ genes: Infer ML tree with FastTree
        raw_tree = os.path.join(group_tmp, f"{group_id}_raw.nwk")
        with open(raw_tree, 'w') as out_f:
            res_ft = subprocess.run(
                [fasttree_bin, "-lg", "-quiet", aln_path],
                stdout=out_f,
                stderr=subprocess.PIPE,
                text=True
            )
            if res_ft.returncode != 0:
                return (False, f"FastTree error: {res_ft.stderr.strip()}")

        if not os.path.exists(raw_tree) or os.path.getsize(raw_tree) == 0:
            return (False, "FastTree produced empty tree")

        # Treerecs gene tree-species tree reconciliation
        rec_dir = os.path.join(group_tmp, "rec")
        os.makedirs(rec_dir, exist_ok=True)
        cmd = [
            treerecs_bin,
            "-g", raw_tree,
            "-s", species_tree,
            "-o", rec_dir,
            "-c", "_",
            "-p", "Y",
            "-f",
            "-q",
            "--output-without-description"
        ]
        if reroot:
            cmd.append("-r")
        if contraction_threshold and str(contraction_threshold).lower() != "none":
            cmd.extend(["-t", str(contraction_threshold)])

        res_tr = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        res_tree = os.path.join(rec_dir, f"{group_id}_raw.nwk_recs.nwk")

        if res_tr.returncode == 0 and os.path.exists(res_tree) and os.path.getsize(res_tree) > 0:
            shutil.copy(res_tree, tree_path)
        else:
            # Fallback to FastTree raw ML tree if Treerecs reconciliation fails
            shutil.copy(raw_tree, tree_path)

        return (True, None)

    except Exception as e:
        # Clean up any partial/corrupt files
        if os.path.exists(aln_path) and os.path.getsize(aln_path) == 0:
            try: os.remove(aln_path)
            except: pass
        if os.path.exists(tree_path) and os.path.getsize(tree_path) == 0:
            try: os.remove(tree_path)
            except: pass
        return (False, str(e))
        
    finally:
        if os.path.exists(group_tmp):
            shutil.rmtree(group_tmp, ignore_errors=True)

def main():
    if 'snakemake' in globals():
        master_path = snakemake.input.master
        fasta_dir = snakemake.config['paths']['fasta_dir']
        species_tree_val = snakemake.config.get('species_tree') or snakemake.config['paths'].get('species_tree')
        out_tar = snakemake.output.tar
        
        mafft_bin = snakemake.config['tools'].get('mafft', 'mafft')
        fasttree_bin = snakemake.config['tools'].get('fasttree', 'FastTree')
        treerecs_bin = snakemake.config['tools'].get('treerecs', 'treerecs')
        if not shutil.which(treerecs_bin) and os.environ.get("TREERECS_BIN"):
            treerecs_bin = os.environ["TREERECS_BIN"]
            
        output_dir = snakemake.config.get('output_dir', 'ensembl_pipeline_output')
        trees_cfg = snakemake.config.get('trees', {})
        contraction_threshold = trees_cfg.get('contraction_threshold', 0.8)
        reroot = trees_cfg.get('reroot', True)
        species_list = list(snakemake.config['species'].keys())
        max_workers = snakemake.threads if hasattr(snakemake, 'threads') else max(1, (os.cpu_count() or 2) - 1)
    else:
        config_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
        if not os.path.exists(config_path):
            print(f"Error: Config file '{config_path}' not found.")
            sys.exit(1)
        with open(config_path, 'r') as f:
            cfg = yaml.safe_load(f)
            
        output_dir = cfg.get('output_dir', 'ensembl_pipeline_output')
        master_path = os.path.join(output_dir, "Consensus_Master.tsv")
        fasta_dir = cfg['paths']['fasta_dir']
        species_tree_val = cfg.get('species_tree') or cfg['paths'].get('species_tree')
        out_tar = os.path.join(output_dir, "Phylogenetic_Trees.tar.gz")
        
        tools_cfg = cfg.get('tools', {})
        mafft_bin = tools_cfg.get('mafft', 'mafft')
        fasttree_bin = tools_cfg.get('fasttree', 'FastTree')
        treerecs_bin = tools_cfg.get('treerecs', 'treerecs')
        if not shutil.which(treerecs_bin) and os.environ.get("TREERECS_BIN"):
            treerecs_bin = os.environ["TREERECS_BIN"]
            
        trees_cfg = cfg.get('trees', {})
        contraction_threshold = trees_cfg.get('contraction_threshold', 0.8)
        reroot = trees_cfg.get('reroot', True)
        species_list = list(cfg.get('species', {}).keys())
        max_workers = cfg.get('threads', max(1, (os.cpu_count() or 2) - 1))

    out_dir = os.path.join(output_dir, "Phylogenetic_Trees")
    msa_dir = os.path.join(out_dir, "MSAs")
    tree_dir = os.path.join(out_dir, "Gene_Trees")
    tmp_dir = os.path.join(out_dir, "tmp_work")
    
    os.makedirs(msa_dir, exist_ok=True)
    os.makedirs(tree_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)
    
    # If species_tree is an inline string, write it to out_dir (isolated from tmp_dir)
    if isinstance(species_tree_val, str) and species_tree_val.strip().startswith('('):
        tree_file = os.path.join(out_dir, "species_tree.nwk")
        with open(tree_file, "w") as f:
            f.write(species_tree_val.strip() + "\n")
        species_tree = tree_file
    else:
        species_tree = species_tree_val
    
    master = pd.read_csv(master_path, sep='\t')
    groups = master[master['ens_orthogroup_id'] != 'unassigned'].groupby('ens_orthogroup_id')['ensembl_id'].apply(list).to_dict()
    
    eligible = {k: v for k, v in groups.items() if len(v) >= 2}
    total = len(eligible)
    print(f"Starting batch tree computation for {total} groups...")
    print(f"  Configuration: contraction_threshold={contraction_threshold}, reroot={reroot}, workers={max_workers}")
    
    completed = 0
    errors = 0
    error_log = []
    start_time = time.time()
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers, initializer=init_worker, initargs=(fasta_dir, species_list)) as executor:
        futures = {
            executor.submit(
                process_group, gid, genes, msa_dir, tree_dir, tmp_dir,
                mafft_bin, fasttree_bin, species_tree, treerecs_bin,
                contraction_threshold, reroot
            ): gid for gid, genes in eligible.items()
        }
        for future in concurrent.futures.as_completed(futures):
            gid = futures[future]
            try:
                success, err_msg = future.result()
                if success:
                    completed += 1
                else:
                    errors += 1
                    error_log.append((gid, err_msg))
            except Exception as e:
                errors += 1
                error_log.append((gid, str(e)))
            
            if (completed + errors) % 500 == 0:
                elapsed = time.time() - start_time
                print(f"  Progress: {completed + errors}/{total} | Errors: {errors} | Elapsed: {elapsed/60:.1f}m")

    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir, ignore_errors=True)
        
    print(f"Computation Complete. Success: {completed}, Failures: {errors}")
    if error_log:
        err_file = os.path.join(out_dir, "tree_computation_errors.log")
        with open(err_file, "w") as ef:
            ef.write("orthogroup_id\terror_message\n")
            for gid, msg in error_log:
                ef.write(f"{gid}\t{msg}\n")
        print(f"Logged {len(error_log)} errors to {err_file}")
    
    print(f"Archiving results into {out_tar}...")
    with tarfile.open(out_tar, "w:gz") as tar:
        tar.add(out_dir, arcname="Phylogenetic_Trees")
    
    # Clean up uncompressed dir to save inodes
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir, ignore_errors=True)
    print("Done!")

if __name__ == "__main__":
    main()
