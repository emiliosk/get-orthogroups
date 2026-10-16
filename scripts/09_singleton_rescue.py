#!/usr/bin/env python3
"""
Step 9: Species-Agnostic Staged Orthology Rescue & Curation Engine (v0.2)
========================================================================
Features:
1. 100% Species-Agnostic: dynamically operates on arbitrary species sets.
2. Iterative Anchor Chaining: runs up to max_iterations (default: 3) so that
   newly rescued genes unlock cassettes for adjacent singletons.
3. Micro-Synteny & Polarity: hierarchical +/-2 macro-anchors (L2, L1, R1, R2)
   with relative strand polarity matching.
4. Pairwise BLOSUM62 Dynamic Programming: exact Needleman-Wunsch alignment
   for high sensitivity without heuristic dropouts.
5. Dual-Hypothesis Testing: side-by-side Synteny vs. Symbol alignment.
6. Three Standardized Deliverables:
   - rescued_genes_all.tsv (Full audit log with paired 4-anchor context)
   - rescue_curation_candidates.csv (28-column diagnostic review table)
   - manual_curation_sample.csv (9-column turnkey starter file with proposed rescues)
"""
import os
import sys
import time
import re
import pandas as pd
from collections import defaultdict
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices

def run_singleton_rescue():
    start_time = time.time()
    print("=== Step 9: Staged Orthology Rescue & Curation Engine (v0.2) ===")
    
    # 1. Resolve inputs and parameters from Snakemake or fallback config
    if 'snakemake' in globals():
        baseline_path = snakemake.input.baseline_clusters
        meta_path = snakemake.input.metadata
        coords_dir = snakemake.params.coords_dir
        fasta_dir = snakemake.params.fasta_dir
        species_list = list(snakemake.params.species_codes)
        
        rescued_tsv = snakemake.output.rescued_tsv
        candidates_csv = snakemake.output.candidates_csv
        curation_sample_csv = snakemake.output.curation_sample_csv
        
        min_alignment_score = float(snakemake.params.min_alignment_score)
        min_consensus_species = int(snakemake.params.min_consensus_species)
        max_cassette_size = int(snakemake.params.max_cassette_size)
        max_iterations = int(snakemake.params.max_iterations)
    else:
        # Standalone CLI execution fallback
        output_dir = "ensembl_pipeline_output"
        baseline_path = os.path.join(output_dir, "Ensembl_N0_Clusters.tsv")
        meta_path = "input/gene_metadata/target_gene_annotations.csv"
        coords_dir = "input/DB/genespace_coords"
        fasta_dir = "input/DB/canonical_transcripts"
        species_list = ['HUMAN', 'MOUSE', 'MACFA', 'RATNO', 'PIGXX', 'CANLF']
        
        rescue_dir = os.path.join(output_dir, "rescue")
        os.makedirs(rescue_dir, exist_ok=True)
        rescued_tsv = os.path.join(rescue_dir, "rescued_genes_all.tsv")
        candidates_csv = os.path.join(rescue_dir, "rescue_curation_candidates.csv")
        curation_sample_csv = os.path.join(rescue_dir, "manual_curation_sample.csv")
        
        min_alignment_score = 50.0
        min_consensus_species = 2
        max_cassette_size = 5
        max_iterations = 3

    print(f"Target Species ({len(species_list)}): {', '.join(species_list)}")
    print(f"Parameters: min_score={min_alignment_score}, min_consensus={min_consensus_species}, max_gap={max_cassette_size}, max_iter={max_iterations}")

    # Ensure parent output directory exists
    os.makedirs(os.path.dirname(os.path.abspath(rescued_tsv)), exist_ok=True)

    # 2. Setup Biopython BLOSUM62 aligner
    print("1. Initializing BLOSUM62 pairwise aligner...")
    matrix = substitution_matrices.load('BLOSUM62')
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = 'global'

    def calc_alignment(seq1, seq2):
        if not seq1 or not seq2:
            return 0.0, 0.0, 0.0
        try:
            score = aligner.score(seq1, seq2)
            a = next(aligner.align(seq1, seq2))
            identities = a.counts().identities
            max_len = max(len(seq1), len(seq2))
            min_len = min(len(seq1), len(seq2))
            id_pct = round((identities / max_len) * 100, 2) if max_len > 0 else 0.0
            cov_pct = round((min_len / max_len) * 100, 2) if max_len > 0 else 0.0
            return round(score, 1), id_pct, cov_pct
        except Exception:
            return 0.0, 0.0, 0.0

    # 3. Load canonical peptides
    print(f"2. Loading canonical peptide FASTAs from {fasta_dir}...")
    gene_pattern = re.compile(r'gene:([A-Za-z0-9_]+)')
    strand_pattern = re.compile(r':(-?1)\s+gene:')
    peptides = {}
    fasta_strands = {}
    for sp in species_list:
        fa_path = os.path.join(fasta_dir, f"{sp}.fa")
        sp_peps = {}
        if os.path.exists(fa_path):
            for rec in SeqIO.parse(fa_path, 'fasta'):
                m = gene_pattern.search(rec.description)
                if m:
                    gid = m.group(1).split('.')[0]
                    sp_peps[gid] = str(rec.seq)
                    sm = strand_pattern.search(rec.description)
                    if sm:
                        fasta_strands[gid] = '+' if sm.group(1) == '1' else '-'
        peptides[sp] = sp_peps
        print(f"   - {sp}: {len(sp_peps)} peptides loaded.")

    # 4. Load baseline clusters & metadata
    print("3. Loading baseline clusters and metadata...")
    clusters = pd.read_csv(baseline_path, sep='\t')
    meta = pd.read_csv(meta_path)
    
    gene_col = 'ensembl_gene_id' if 'ensembl_gene_id' in meta.columns else 'ensembl_id'
    meta['clean_gene_id'] = meta[gene_col].astype(str).str.strip().str.split('.').str[0]

    # Try loading precomputed strands if available, else fasta strands, else default to '+'
    gene_to_strand = {}
    strand_cache = "test/cache/gene_strands.tsv"
    if os.path.exists(strand_cache):
        strands_df = pd.read_csv(strand_cache, sep='\t')
        gene_to_strand = dict(zip(strands_df['gene_id'], strands_df['strand']))
    elif 'strand' in meta.columns:
        gene_to_strand = dict(zip(meta['clean_gene_id'], meta['strand']))
    for gid, st in fasta_strands.items():
        if gid not in gene_to_strand:
            gene_to_strand[gid] = st

    pep_to_gene = {}
    if 'peptide_id' in meta.columns:
        df_p = meta.dropna(subset=['peptide_id'])
        for p, g in zip(df_p['peptide_id'], df_p['clean_gene_id']):
            p_str = str(p).strip()
            pep_to_gene[p_str] = g
            pep_to_gene[p_str.split('.')[0]] = g

    if 'peptide_id_version' in meta.columns:
        df_pv = meta.dropna(subset=['peptide_id_version'])
        for p, g in zip(df_pv['peptide_id_version'], df_pv['clean_gene_id']):
            p_str = str(p).strip()
            pep_to_gene[p_str] = g
            pep_to_gene[p_str.split('.')[0]] = g

    for g in meta['clean_gene_id'].dropna():
        g_str = str(g).strip()
        pep_to_gene[g_str] = g
        pep_to_gene[g_str.split('.')[0]] = g

    # 5. Load coordinates for all configured species
    print("4. Loading genomic coordinates across all species...")
    coords = {}
    for sp in species_list:
        bed_path = os.path.join(coords_dir, f"{sp}_coords.bed")
        if not os.path.exists(bed_path):
            print(f"   [WARNING] Missing coordinate file: {bed_path}")
            continue
        df = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end', 'pep_id'])
        df['pep_id_str'] = df['pep_id'].astype(str).str.strip()
        df['gene_id'] = df['pep_id_str'].map(pep_to_gene)
        mask_missing = df['gene_id'].isna()
        if mask_missing.any():
            df.loc[mask_missing, 'gene_id'] = df.loc[mask_missing, 'pep_id_str'].str.split('.').str[0].map(pep_to_gene)
        df = df.dropna(subset=['gene_id']).drop_duplicates(subset=['gene_id'])
        df['chr'] = df['chr'].astype(str).str.replace(r'\.0$', '', regex=True)
        df['mid'] = (df['start'] + df['end']) / 2
        df['strand'] = df['gene_id'].map(gene_to_strand).fillna('+')
        coords[sp] = df.sort_values(['chr', 'mid']).reset_index(drop=True)

    sp_gene_idx = {sp: {row['gene_id']: idx for idx, row in df.iterrows()} for sp, df in coords.items()}

    def get_relative_orientation(strand1, strand2):
        return 1 if strand1 == strand2 else -1

    def fmt_anchor(row):
        if row is None:
            return "NONE"
        sym = row.get('sym', row['gene_id'])
        og = row.get('og', 'NA')
        strand = row.get('strand', '+')
        return f"{sym} [{og}] ({strand})"

    # State tracking across iterations
    current_clusters = clusters.copy()
    all_rescued_records = []
    all_candidates_records = []
    already_rescued_gids = set()

    # 6. Iterative Anchor-Chaining Loop
    print(f"5. Beginning Iterative Singleton Rescue (up to {max_iterations} iterations)...")
    
    for iteration in range(1, max_iterations + 1):
        print(f"\n--- Iteration {iteration} of {max_iterations} ---")
        
        # Build live graph state for this iteration
        gene_to_og = dict(zip(current_clusters['ensembl_id'], current_clusters['ens_orthogroup_id']))
        gene_to_hq = dict(zip(current_clusters['ensembl_id'], current_clusters.get('ens_hqorthogroup_id', current_clusters['ens_orthogroup_id'])))
        gene_to_sp = dict(zip(current_clusters['ensembl_id'], current_clusters['species_code']))
        gene_to_sym = dict(zip(current_clusters['ensembl_id'], current_clusters['gene_symbol'].fillna(current_clusters['ensembl_id'])))

        og_to_species = defaultdict(set)
        og_to_genes = defaultdict(list)
        for gid, og in gene_to_og.items():
            og_to_species[og].add(gene_to_sp[gid])
            og_to_genes[og].append(gid)

        sym_to_ogs = defaultdict(set)
        for gid, sym in gene_to_sym.items():
            og = gene_to_og.get(gid)
            if og and not og.startswith('OG_SNG_'):
                sym_to_ogs[sym.upper()].add(og)

        for sp in species_list:
            if sp in coords:
                coords[sp]['og'] = coords[sp]['gene_id'].map(gene_to_og)
                coords[sp]['sym'] = coords[sp]['gene_id'].map(gene_to_sym)

        def is_conserved_anchor(og):
            return bool(og and not og.startswith('OG_SNG_') and len(og_to_species.get(og, set())) >= 2)

        def get_og_rep(og):
            if not og or og == 'SINGLETON' or og not in og_to_genes:
                return "NONE", "NONE", "NONE"
            for g in og_to_genes[og]:
                sp_g = gene_to_sp.get(g)
                sym_g = gene_to_sym.get(g, g)
                if sym_g != g:
                    return g, sp_g, sym_g
            g = og_to_genes[og][0]
            return g, gene_to_sp.get(g, 'NONE'), gene_to_sym.get(g, g)

        singletons = current_clusters[
            current_clusters['ens_orthogroup_id'].str.startswith('OG_SNG_') & 
            (~current_clusters['ensembl_id'].isin(already_rescued_gids))
        ].copy()
        
        print(f"   Singletons to evaluate in this iteration: {len(singletons)}")
        round_rescues = 0

        for _, row in singletons.iterrows():
            s_gid = row['ensembl_id']
            s_sp = row['species_code']
            s_sym = str(row['gene_symbol']) if pd.notna(row['gene_symbol']) else s_gid
            has_symbol = bool(pd.notna(row['gene_symbol']) and str(row['gene_symbol']).strip() != s_gid)

            if s_sp not in sp_gene_idx or s_gid not in sp_gene_idx[s_sp]:
                continue

            s_df = coords[s_sp]
            s_idx = sp_gene_idx[s_sp][s_gid]
            s_row = s_df.iloc[s_idx]
            s_chr = s_row['chr']
            s_strand = s_row['strand']
            s_start = s_row['start']
            s_end = s_row['end']
            query_locus = f"chr{s_chr}:{s_start}-{s_end} ({s_strand})"

            q_seq = peptides[s_sp].get(s_gid)
            q_len = len(q_seq) if q_seq else 0

            chr_df = s_df[s_df['chr'] == s_chr]
            chr_min_idx = chr_df.index.min()
            chr_max_idx = chr_df.index.max()
            is_terminal = (s_idx - chr_min_idx < 2) or (chr_max_idx - s_idx < 2)
            boundary_status = "TERMINAL_BOUNDARY" if is_terminal else "INTERSTITIAL"

            left_anchors = []
            for step in range(1, 4):
                idx_l = s_idx - step
                if idx_l >= 0:
                    r_l = s_df.iloc[idx_l]
                    if r_l['chr'] == s_chr and is_conserved_anchor(r_l['og']):
                        left_anchors.append(r_l)
                        if len(left_anchors) >= 2: break

            right_anchors = []
            for step in range(1, 4):
                idx_r = s_idx + step
                if idx_r < len(s_df):
                    r_r = s_df.iloc[idx_r]
                    if r_r['chr'] == s_chr and is_conserved_anchor(r_r['og']):
                        right_anchors.append(r_r)
                        if len(right_anchors) >= 2: break

            q_l1 = left_anchors[0] if len(left_anchors) >= 1 else None
            q_l2 = left_anchors[1] if len(left_anchors) >= 2 else None
            q_r1 = right_anchors[0] if len(right_anchors) >= 1 else None
            q_r2 = right_anchors[1] if len(right_anchors) >= 2 else None

            # Symbol Candidate Check
            symbol_candidate_og = None
            sym_ref_gid, sym_ref_sp, sym_ref_sym = "NONE", "NONE", "NONE"
            sym_score, sym_id, sym_cov = None, None, None

            if has_symbol and s_sym.upper() in sym_to_ogs:
                ogs = sym_to_ogs[s_sym.upper()]
                eligible_symbol_ogs = [og for og in ogs if s_sp not in og_to_species[og]]
                if len(eligible_symbol_ogs) == 1:
                    symbol_candidate_og = eligible_symbol_ogs[0]
                    sym_ref_gid, sym_ref_sp, sym_ref_sym = get_og_rep(symbol_candidate_og)
                    sym_seq = peptides[sym_ref_sp].get(sym_ref_gid) if sym_ref_sp in peptides else None
                    sym_score, sym_id, sym_cov = calc_alignment(q_seq, sym_seq)

            def make_candidate_record(
                conflict_type, recommended_action, curator_notes,
                synteny_og="NONE", synteny_sym="NONE", synteny_ref="NONE", synteny_sp="NONE",
                syn_score=None, syn_id=None, syn_cov=None,
                sug_og=None, sug_gene=None, sug_sp=None, sug_sym=None
            ):
                return {
                    'gene_id': s_gid,
                    'species': s_sp,
                    'gene_symbol': s_sym,
                    'query_locus': query_locus,
                    'query_peptide_len': q_len,
                    'boundary_status': boundary_status,
                    'query_l1_anchor': fmt_anchor(q_l1),
                    'query_r1_anchor': fmt_anchor(q_r1),
                    'synteny_og': synteny_og,
                    'synteny_target_symbol': synteny_sym,
                    'synteny_ref_gene': synteny_ref,
                    'synteny_supporting_species': synteny_sp,
                    'synteny_align_score': syn_score,
                    'synteny_identity_pct': syn_id,
                    'synteny_coverage_pct': syn_cov,
                    'symbol_og': symbol_candidate_og if symbol_candidate_og else "SINGLETON",
                    'symbol_target_symbol': sym_ref_sym if symbol_candidate_og else "NONE",
                    'symbol_ref_gene': f"{sym_ref_gid} ({sym_ref_sp})" if symbol_candidate_og else "NONE",
                    'symbol_align_score': sym_score,
                    'symbol_identity_pct': sym_id,
                    'symbol_coverage_pct': sym_cov,
                    'conflict_type': conflict_type,
                    'recommended_action': recommended_action,
                    'suggested_target_og': sug_og if sug_og else ("SINGLETON" if recommended_action == "LEAVE_SINGLETON" else "REVIEW"),
                    'suggested_target_gene': sug_gene if sug_gene else "NONE",
                    'suggested_target_species': sug_sp if sug_sp else "NONE",
                    'suggested_target_symbol': sug_sym if sug_sym else "NONE",
                    'curator_notes': curator_notes
                }

            # Case A: Boundary / Scaffold tip
            if not left_anchors or not right_anchors:
                if symbol_candidate_og and sym_score and sym_score >= min_alignment_score:
                    if iteration == 1:
                        all_candidates_records.append(make_candidate_record(
                            conflict_type='TERMINAL_SCAFFOLD_END',
                            recommended_action='ACCEPT_SYMBOL',
                            curator_notes=f"Terminal scaffold gene; matched {symbol_candidate_og} ({sym_ref_sym}) with score {sym_score}",
                            sug_og=symbol_candidate_og, sug_gene=sym_ref_gid, sug_sp=sym_ref_sp, sug_sym=sym_ref_sym
                        ))
                continue

            candidate_votes = defaultdict(lambda: {'count': 0, 'species': set(), 'details': []})

            for ref_sp in species_list:
                if ref_sp == s_sp or ref_sp not in coords:
                    continue

                ref_df = coords[ref_sp]
                matched_pair = None
                cass_label = None

                for la_i, la in enumerate(left_anchors):
                    la_ref_genes = [g for g in og_to_genes[la['og']] if gene_to_sp[g] == ref_sp and g in sp_gene_idx[ref_sp]]
                    if not la_ref_genes: continue
                    la_ref_idx = sp_gene_idx[ref_sp][la_ref_genes[0]]
                    la_ref_row = ref_df.iloc[la_ref_idx]

                    for ra_i, ra in enumerate(right_anchors):
                        ra_ref_genes = [g for g in og_to_genes[ra['og']] if gene_to_sp[g] == ref_sp and g in sp_gene_idx[ref_sp]]
                        if not ra_ref_genes: continue
                        ra_ref_idx = sp_gene_idx[ref_sp][ra_ref_genes[0]]
                        ra_ref_row = ref_df.iloc[ra_ref_idx]

                        if la_ref_row['chr'] == ra_ref_row['chr']:
                            cass_label = f"L{la_i+1}-R{ra_i+1}"
                            matched_pair = (la, ra, la_ref_row, ra_ref_row, la_ref_idx, ra_ref_idx, cass_label)
                            break
                    if matched_pair: break

                if not matched_pair: continue

                la_q, ra_q, la_ref, ra_ref, la_ref_idx, ra_ref_idx, cass_label = matched_pair
                ref_min_idx = min(la_ref_idx, ra_ref_idx)
                ref_max_idx = max(la_ref_idx, ra_ref_idx)

                ref_cassette_rows = ref_df.iloc[ref_min_idx + 1: ref_max_idx]
                if len(ref_cassette_rows) > max_cassette_size:
                    continue

                for _, ref_gene in ref_cassette_rows.iterrows():
                    cand_og = ref_gene['og']
                    if not cand_og or cand_og.startswith('OG_SNG_'): continue
                    if s_sp in og_to_species[cand_og]: continue

                    q_rel_polarity = get_relative_orientation(s_strand, la_q['strand'])
                    ref_rel_polarity = get_relative_orientation(ref_gene['strand'], la_ref['strand'])
                    strand_match = (q_rel_polarity == ref_rel_polarity)

                    if strand_match:
                        step_l = -1 if la_ref_idx < ra_ref_idx else 1
                        step_r = 1 if la_ref_idx < ra_ref_idx else -1

                        ref_l2 = None
                        for s in range(1, 4):
                            idx_out = la_ref_idx + s * step_l
                            if 0 <= idx_out < len(ref_df):
                                cand_r = ref_df.iloc[idx_out]
                                if cand_r['chr'] == la_ref['chr'] and is_conserved_anchor(cand_r['og']):
                                    ref_l2 = cand_r; break

                        ref_r2 = None
                        for s in range(1, 4):
                            idx_out = ra_ref_idx + s * step_r
                            if 0 <= idx_out < len(ref_df):
                                cand_r = ref_df.iloc[idx_out]
                                if cand_r['chr'] == ra_ref['chr'] and is_conserved_anchor(cand_r['og']):
                                    ref_r2 = cand_r; break

                        candidate_votes[cand_og]['count'] += 1
                        candidate_votes[cand_og]['species'].add(ref_sp)
                        candidate_votes[cand_og]['details'].append({
                            'ref_sp': ref_sp,
                            'ref_gene_id': ref_gene['gene_id'],
                            'ref_symbol': ref_gene['sym'],
                            'strand_match': strand_match,
                            'cass_label': cass_label,
                            'ref_l2': ref_l2,
                            'ref_l1': la_ref,
                            'ref_r1': ra_ref,
                            'ref_r2': ref_r2
                        })

            # Case B: Synteny found no cassette candidate
            if not candidate_votes:
                if symbol_candidate_og and sym_score and sym_score >= min_alignment_score:
                    if iteration == 1:
                        all_candidates_records.append(make_candidate_record(
                            conflict_type='SYNTENY_MISSED_INVERSION_OR_TANDEM',
                            recommended_action='ACCEPT_SYMBOL',
                            curator_notes=f"Synteny missed due to local rearrangement; aligned to {symbol_candidate_og} ({sym_ref_sym}) with score {sym_score}",
                            sug_og=symbol_candidate_og, sug_gene=sym_ref_gid, sug_sp=sym_ref_sp, sug_sym=sym_ref_sym
                        ))
                continue

            sorted_cands = sorted(candidate_votes.items(), key=lambda x: len(x[1]['species']), reverse=True)
            top_og, top_data = sorted_cands[0]
            chosen_og = None
            chosen_data = None

            if len(sorted_cands) == 1:
                chosen_og, chosen_data = top_og, top_data
            elif len(sorted_cands) > 1:
                second_og, second_data = sorted_cands[1]
                if len(top_data['species']) > len(second_data['species']):
                    chosen_og, chosen_data = top_og, top_data
                else:
                    if symbol_candidate_og == top_og:
                        chosen_og, chosen_data = top_og, top_data
                    elif symbol_candidate_og == second_og:
                        chosen_og, chosen_data = second_og, second_data
                    else:
                        if iteration == 1:
                            best_d = top_data['details'][0]
                            score_t, id_t, cov_t = calc_alignment(q_seq, peptides[best_d['ref_sp']].get(best_d['ref_gene_id']))
                            top_sp_str = ", ".join(sorted(top_data['species']))
                            sec_sp_str = ", ".join(sorted(second_data['species']))
                            rec_act = "REVIEW"
                            sug_og, sug_gid, sug_sp, sug_sym = None, None, None, None
                            if symbol_candidate_og and sym_score and sym_score >= 100.0:
                                rec_act = "ACCEPT_SYMBOL"
                                sug_og, sug_gid, sug_sp, sug_sym = symbol_candidate_og, sym_ref_gid, sym_ref_sp, sym_ref_sym
                                notes = f"Tied synteny between {top_og} and {second_og}; Symbol match {symbol_candidate_og} ({sym_ref_sym}) verified with score {sym_score}"
                            else:
                                notes = f"Tied synteny vote between {top_og} ({top_sp_str}) and {second_og} ({sec_sp_str})"

                            all_candidates_records.append(make_candidate_record(
                                conflict_type='AMBIGUOUS_SYNTENY_TIE',
                                recommended_action=rec_act,
                                curator_notes=notes,
                                synteny_og=f"{top_og}_OR_{second_og}",
                                synteny_sym=f"{best_d['ref_symbol']}_vs_{second_data['details'][0]['ref_symbol']}",
                                synteny_ref=f"{best_d['ref_gene_id']} ({best_d['ref_sp']})",
                                synteny_sp=f"{top_sp_str} / {sec_sp_str}",
                                syn_score=score_t, syn_id=id_t, syn_cov=cov_t,
                                sug_og=sug_og, sug_gene=sug_gid, sug_sp=sug_sp, sug_sym=sug_sym
                            ))
                        continue

            if chosen_og:
                best_detail = chosen_data['details'][0]
                ref_gid = best_detail['ref_gene_id']
                ref_sp = best_detail['ref_sp']
                ref_sym = best_detail['ref_symbol']
                sp_count = len(chosen_data['species'])
                sp_str = ", ".join(sorted(chosen_data['species']))

                ref_seq = peptides[ref_sp].get(ref_gid)
                align_score, id_pct, cov_pct = calc_alignment(q_seq, ref_seq)
                ref_len = len(ref_seq) if ref_seq else 0

                rescue_tier = "SYNTENY_PLUS_SYMBOL" if (has_symbol and symbol_candidate_og == chosen_og) else "PURE_SYNTENY_UNANNOTATED"

                # Check if it meets clean rescue criteria
                is_clean_rescue = (align_score >= min_alignment_score) and (has_symbol or sp_count >= min_consensus_species)
                
                if is_clean_rescue:
                    target_hq = None
                    for g in og_to_genes[chosen_og]:
                        if gene_to_hq.get(g):
                            target_hq = gene_to_hq[g]; break

                    all_rescued_records.append({
                        'query_gene_id': s_gid,
                        'query_species': s_sp,
                        'query_symbol': s_sym,
                        'query_locus': query_locus,
                        'query_l2_anchor': fmt_anchor(q_l2),
                        'query_l1_anchor': fmt_anchor(q_l1),
                        'query_r1_anchor': fmt_anchor(q_r1),
                        'query_r2_anchor': fmt_anchor(q_r2),
                        'rescued_og': chosen_og,
                        'rescued_hq': target_hq,
                        'ref_ortholog_id': ref_gid,
                        'ref_species': ref_sp,
                        'ref_symbol': ref_sym,
                        'ref_l2_anchor': fmt_anchor(best_detail['ref_l2']),
                        'ref_l1_anchor': fmt_anchor(best_detail['ref_l1']),
                        'ref_r1_anchor': fmt_anchor(best_detail['ref_r1']),
                        'ref_r2_anchor': fmt_anchor(best_detail['ref_r2']),
                        'active_cassette_used': best_detail['cass_label'],
                        'supporting_species_count': sp_count,
                        'supporting_species': sp_str,
                        'boundary_status': boundary_status,
                        'alignment_score': align_score,
                        'seq_identity_pct': id_pct,
                        'seq_coverage_pct': cov_pct,
                        'rescue_tier': rescue_tier,
                        'iteration_rescued': iteration,
                        'decision': 'ACCEPT'
                    })
                    already_rescued_gids.add(s_gid)
                    round_rescues += 1
                elif iteration == 1:
                    # Log edge cases in candidates table
                    if align_score < 0.0:
                        is_trunc = (cov_pct < 15.0 and q_len < 100 and ref_len > 500)
                        if is_trunc and symbol_candidate_og == chosen_og:
                            notes = f"Truncated peptide fragment ({q_len} aa vs {ref_len} aa ref); global penalty gives score {align_score}; symbol matches {chosen_og}"
                            rec_act = "REVIEW"
                        else:
                            notes = f"Synteny matched {chosen_og} ({ref_sym}) but BLOSUM62 score is negative ({align_score}); likely pseudogene"
                            rec_act = "LEAVE_SINGLETON"
                        all_candidates_records.append(make_candidate_record(
                            conflict_type='SUSPECTED_RETROGENE', recommended_action=rec_act, curator_notes=notes,
                            synteny_og=chosen_og, synteny_sym=ref_sym, synteny_ref=f"{ref_gid} ({ref_sp})", synteny_sp=sp_str,
                            syn_score=align_score, syn_id=id_pct, syn_cov=cov_pct,
                            sug_og=chosen_og if rec_act != "LEAVE_SINGLETON" else "SINGLETON",
                            sug_gene=ref_gid if rec_act != "LEAVE_SINGLETON" else "NONE",
                            sug_sp=ref_sp if rec_act != "LEAVE_SINGLETON" else "NONE",
                            sug_sym=ref_sym if rec_act != "LEAVE_SINGLETON" else "NONE"
                        ))
                    elif 0.0 <= align_score < min_alignment_score:
                        if symbol_candidate_og and sym_score and sym_score >= min_alignment_score:
                            all_candidates_records.append(make_candidate_record(
                                conflict_type='LOW_SEQUENCE_SIMILARITY', recommended_action='ACCEPT_SYMBOL',
                                curator_notes=f"Synteny matched adjacent gene {chosen_og} ({ref_sym}, score {align_score}), but symbol matches {symbol_candidate_og} ({sym_ref_sym}) with score {sym_score}",
                                synteny_og=chosen_og, synteny_sym=ref_sym, synteny_ref=f"{ref_gid} ({ref_sp})", synteny_sp=sp_str,
                                syn_score=align_score, syn_id=id_pct, syn_cov=cov_pct,
                                sug_og=symbol_candidate_og, sug_gene=sym_ref_gid, sug_sp=sym_ref_sp, sug_sym=sym_ref_sym
                            ))
                        else:
                            all_candidates_records.append(make_candidate_record(
                                conflict_type='LOW_SEQUENCE_SIMILARITY', recommended_action='REVIEW',
                                curator_notes=f"Positive score ({align_score}) but under confidence floor ({min_alignment_score}); needs manual review",
                                synteny_og=chosen_og, synteny_sym=ref_sym, synteny_ref=f"{ref_gid} ({ref_sp})", synteny_sp=sp_str,
                                syn_score=align_score, syn_id=id_pct, syn_cov=cov_pct,
                                sug_og=chosen_og, sug_gene=ref_gid, sug_sp=ref_sp, sug_sym=ref_sym
                            ))
                    elif symbol_candidate_og and symbol_candidate_og != chosen_og:
                        all_candidates_records.append(make_candidate_record(
                            conflict_type='DISCORDANT_ASSIGNMENT', recommended_action='ACCEPT_SYMBOL',
                            curator_notes=f"Synteny matched {chosen_og} ({ref_sym}, score {align_score}), but symbol matches {symbol_candidate_og} ({sym_ref_sym}, score {sym_score})",
                            synteny_og=chosen_og, synteny_sym=ref_sym, synteny_ref=f"{ref_gid} ({ref_sp})", synteny_sp=sp_str,
                            syn_score=align_score, syn_id=id_pct, syn_cov=cov_pct,
                            sug_og=symbol_candidate_og, sug_gene=sym_ref_gid, sug_sp=sym_ref_sp, sug_sym=sym_ref_sym
                        ))
                    elif not has_symbol and sp_count < min_consensus_species:
                        all_candidates_records.append(make_candidate_record(
                            conflict_type='SINGLE_SPECIES_UNANNOTATED', recommended_action='REVIEW',
                            curator_notes=f"Unnamed transcript with only 1 supporting species ({ref_sp}); score {align_score}",
                            synteny_og=chosen_og, synteny_sym=ref_sym, synteny_ref=f"{ref_gid} ({ref_sp})", synteny_sp=sp_str,
                            syn_score=align_score, syn_id=id_pct, syn_cov=cov_pct,
                            sug_og=chosen_og, sug_gene=ref_gid, sug_sp=ref_sp, sug_sym=ref_sym
                        ))

        print(f"   -> Rescued in iteration {iteration}: {round_rescues} singletons")
        
        if round_rescues == 0:
            print(f"Iterative rescue converged early at iteration {iteration} (0 new rescues).")
            break
            
        # Update in-memory graph for next iteration
        for r in all_rescued_records[-round_rescues:]:
            current_clusters.loc[current_clusters['ensembl_id'] == r['query_gene_id'], 'ens_orthogroup_id'] = r['rescued_og']

    # 7. Write the 3 Standardized Deliverables
    print("\n6. Exporting standardized deliverables...")
    
    # 1. rescued_genes_all.tsv
    rescued_cols = [
        'query_gene_id', 'query_species', 'query_symbol', 'query_locus',
        'query_l2_anchor', 'query_l1_anchor', 'query_r1_anchor', 'query_r2_anchor',
        'rescued_og', 'rescued_hq', 'ref_ortholog_id', 'ref_species', 'ref_symbol',
        'ref_l2_anchor', 'ref_l1_anchor', 'ref_r1_anchor', 'ref_r2_anchor',
        'active_cassette_used', 'supporting_species_count', 'supporting_species',
        'boundary_status', 'alignment_score', 'seq_identity_pct', 'seq_coverage_pct',
        'rescue_tier', 'iteration_rescued', 'decision'
    ]
    df_rescued = pd.DataFrame(all_rescued_records) if all_rescued_records else pd.DataFrame(columns=rescued_cols)
    df_rescued.to_csv(rescued_tsv, sep='\t', index=False)
    print(f"   - Saved {len(df_rescued)} clean rescues across {iteration} iterations to {rescued_tsv}")

    # 2. rescue_curation_candidates.csv
    candidate_cols = [
        'gene_id', 'species', 'gene_symbol', 'query_locus', 'query_peptide_len',
        'boundary_status', 'query_l1_anchor', 'query_r1_anchor',
        'synteny_og', 'synteny_target_symbol', 'synteny_ref_gene', 'synteny_supporting_species',
        'synteny_align_score', 'synteny_identity_pct', 'synteny_coverage_pct',
        'symbol_og', 'symbol_target_symbol', 'symbol_ref_gene',
        'symbol_align_score', 'symbol_identity_pct', 'symbol_coverage_pct',
        'conflict_type', 'recommended_action',
        'suggested_target_og', 'suggested_target_gene', 'suggested_target_species', 'suggested_target_symbol',
        'curator_notes'
    ]
    df_candidates = pd.DataFrame(all_candidates_records) if all_candidates_records else pd.DataFrame(columns=candidate_cols)
    df_candidates.to_csv(candidates_csv, index=False)
    print(f"   - Saved {len(df_candidates)} flagged candidates with rich context to {candidates_csv}")

    # 3. manual_curation_sample.csv (9 columns)
    curation_cols = [
        'gene_id', 'species', 'gene_symbol', 'target_gene_id',
        'target_species', 'target_gene_symbol', 'target_orthogroup',
        'curator_decision', 'curation_reason'
    ]
    curation_rows = []
    for _, r in df_rescued.iterrows():
        curation_rows.append({
            'gene_id': r['query_gene_id'],
            'species': r['query_species'],
            'gene_symbol': r['query_symbol'],
            'target_gene_id': r['ref_ortholog_id'],
            'target_species': r['ref_species'],
            'target_gene_symbol': r['ref_symbol'],
            'target_orthogroup': r['rescued_og'],
            'curator_decision': 'ACCEPT',
            'curation_reason': f"Synteny Rescue (Iter {r['iteration_rescued']}): Cassette {r['active_cassette_used']} ({r['supporting_species_count']} sp); flanked by {r['query_l1_anchor']} and {r['query_r1_anchor']}; Score: {r['alignment_score']}"
        })

    for _, r in df_candidates.iterrows():
        if r['recommended_action'] == 'ACCEPT_SYMBOL':
            curation_rows.append({
                'gene_id': r['gene_id'],
                'species': r['species'],
                'gene_symbol': r['gene_symbol'],
                'target_gene_id': r['suggested_target_gene'],
                'target_species': r['suggested_target_species'],
                'target_gene_symbol': r['suggested_target_symbol'],
                'target_orthogroup': r['suggested_target_og'],
                'curator_decision': 'ACCEPT',
                'curation_reason': f"Curated Addition: {r['conflict_type']}; Symbol & Sequence verified ({r['curator_notes']})"
            })

    df_curation_sample = pd.DataFrame(curation_rows) if curation_rows else pd.DataFrame(columns=curation_cols)
    df_curation_sample.to_csv(curation_sample_csv, index=False)
    print(f"   - Saved streamlined {len(df_curation_sample)} proposed rescues (9 columns) to {curation_sample_csv}")
    print(f"\nStep 9 completed successfully in {round(time.time() - start_time, 2)} seconds.")

if __name__ == '__main__':
    run_singleton_rescue()
