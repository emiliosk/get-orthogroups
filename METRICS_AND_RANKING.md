# Metrics & Ranking Logic

The pipeline calculates 8 core metrics to evaluate orthology confidence.

## Metric Definitions

| Metric | Source | Description |
| :--- | :--- | :--- |
| `ens_identity` | Ensembl | % Identity of the reference protein vs target. |
| `ens_homology_identity` | Ensembl | % Identity of the target protein vs reference. |
| `ens_goc_score` | Ensembl | Ensembl's Gene Order Conservation score. |
| `ens_wga_coverage` | Ensembl | Whole Genome Alignment coverage. |
| `ens_is_high_confidence` | Ensembl | Binary flag (1/0) for Ensembl's expert-curated confidence. |
| `atlas_goc_score_pairwise` | Atlas | Custom Jaccard score of the orthogroup-based neighborhood (Window: 5). |
| `atlas_hq_goc_score_pairwise` | Atlas | Jaccard score using only High-Quality (HQ) neighborhood links. |
| `atlas_msa_pairwise_identity` | Atlas | Custom % sequence identity computed from MAFFT canonical peptide MSA across aligned residues. |

## Ortholog Ranking Logic

To resolve `one-to-many` and `many-to-many` relationships, the pipeline assigns a rank from the perspective of each species.

**Column Name:** `ortholog_rank_{SPECIES}`

**Sorting Algorithm:**
For a given gene, all its orthologs in the target species are sorted by the following columns in descending order:
1.  **`ens_is_high_confidence`**: Prioritize curated links.
2.  **`ens_goc_score`**: Prioritize synteny conservation.
3.  **`atlas_goc_score_pairwise`**: Local Atlas neighborhood support.
4.  **`min(identity, homology_identity)`**: Best overall sequence match.
5.  **`ens_wga_coverage`**: Final tie-breaker using alignment depth.

The first gene in this sorted list becomes **Rank 1**.

### Example: Mouse Insulin
*   **Human INS** → **Mouse Ins2**: Rank 1 (High Conf, 100 GOC)
*   **Human INS** → **Mouse Ins1**: Rank 2 (Low Conf, 0 GOC)

---

## Multiple Sequence Alignment (MSA) Conservation Metrics

When phylogenetic tree and MSA computation are enabled (`compute_trees: true`), the pipeline generates **Ensembl-independent**, package-native sequence conservation metrics computed directly from the canonical protein multiple sequence alignments across all species in each orthogroup.

**Output File:** `MSA_Summary_Scores.tsv` (located directly in the output directory and archived in `Phylogenetic_Trees.tar.gz`).

| Metric | Source | Range | Description |
| :--- | :--- | :---: | :--- |
| `msa_mean_pid` | Atlas | 0 – 100% | Mean pairwise sequence identity (%) across all species pairs in the orthogroup. |
| `msa_conserved_sites_pct` | Atlas | 0 – 100% | Percentage of alignment columns that are 100% identical (invariant) across all species without gaps. |
| `msa_aln_length` | Atlas | Integer | Total length of the multiple sequence alignment (in amino acids, including gaps). |
| `msa_gap_pct` | Atlas | 0 – 100% | Overall percentage of gap characters across the entire alignment matrix. |
| `n_seqs` | Atlas | Integer | Number of protein sequences aligned in the orthogroup. |

---

## Related Documentation

* [README.md](README.md) - Pipeline overview and quick start.
* [USAGE_GUIDE.md](USAGE_GUIDE.md) - Adding species, upgrading Ensembl releases, and troubleshooting.
* [METHODS_PIPELINE.md](METHODS_PIPELINE.md) - Dual-track transitivity and phylogenetic reconciliation methodology.
