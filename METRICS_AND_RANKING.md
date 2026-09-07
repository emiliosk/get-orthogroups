# Metrics & Ranking Logic

The pipeline calculates 7 core metrics to evaluate orthology confidence.

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

## Related Documentation

* [README.md](README.md) - Pipeline overview and quick start.
* [USAGE_GUIDE.md](USAGE_GUIDE.md) - Adding species, upgrading Ensembl releases, and troubleshooting.
* [METHODS_PIPELINE.md](METHODS_PIPELINE.md) - Dual-track transitivity and phylogenetic reconciliation methodology.
