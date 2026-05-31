library(biomaRt)
library(dplyr)
#### run this only for 109 pig, since 115 works fine! has one set canonical transcript per gene already. only 109 doesnt.
ensembl_host <- "https://feb2023.archive.ensembl.org"
species <- "sscrofa_gene_ensembl"
target_attributes <- c(
  "ensembl_gene_id", 
  "ensembl_transcript_id", 
  "ensembl_transcript_id_version", 
  "transcript_is_canonical" , 
  "transcript_appris"
)
mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = species,
  host = ensembl_host
)

# Retrieve the data
pig_109 <- getBM(
  attributes = target_attributes,
  mart = mart
)


#current ensembl
ensembl_host <-"https://sep2025.archive.ensembl.org"

mart <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = species,
  host = ensembl_host
)

# Retrieve the data
pig_115 <- getBM(
  attributes = target_attributes,
  mart = mart
)


appris_priority <- c(
  "principal1", "principal2", "principal3", 
  "principal4", "principal5", "alternative1", "alternative2"
)


pig_109_updated <- pig_109 %>%
  left_join(
    pig_115 %>% 
      select(
        ensembl_gene_id, 
        ensembl_transcript_id, 
        canon_115 = transcript_is_canonical,
        appris_115 = transcript_appris
      ),
    by = c("ensembl_gene_id", "ensembl_transcript_id")
  ) %>%
  group_by(ensembl_gene_id) %>%
  mutate(
    # --- 1. Evaluate Ensembl 115 Canonical ---
    has_canon_115 = any(!is.na(canon_115)),
    
    # --- 2. Evaluate Ensembl 115 APPRIS ---
    rank_115 = match(appris_115, appris_priority),
    min_rank_115 = if (all(is.na(rank_115))) NA_integer_ else min(rank_115, na.rm = TRUE),
    has_appris_115 = !is.na(min_rank_115),
    
    # --- 3. Evaluate Ensembl 109 Canonical ---
    has_canon_109 = any(!is.na(transcript_is_canonical)),
    
    # --- 4. Evaluate Ensembl 109 APPRIS ---
    rank_109 = match(transcript_appris, appris_priority),
    min_rank_109 = if (all(is.na(rank_109))) NA_integer_ else min(rank_109, na.rm = TRUE),
    has_appris_109 = !is.na(min_rank_109),
    
    # --- 5. Evaluate Singleton Status ---
    is_singleton = n() == 1,
    
    # --- Cascade Logic for Source Tracking (Allowing Ties) ---
    canonical_source = case_when(
      # 1. Explicit Canonical in BOTH 115 and 109
      has_canon_115 & !is.na(canon_115) & !is.na(transcript_is_canonical) ~ "ensembl_109/115",
      
      # 2. Explicit Canonical in 115 ONLY
      has_canon_115 & !is.na(canon_115) ~ "ensembl_115",
      
      # 3. Explicit Canonical in 109 (If no 115 canonical exists)
      !has_canon_115 & has_canon_109 & !is.na(transcript_is_canonical) ~ "ensembl_109",
      
      # 4. APPRIS 115 (Matches any transcript sharing the minimum rank)
      !has_canon_115 & !has_canon_109 & has_appris_115 & !is.na(rank_115) & rank_115 == min_rank_115 ~ "best_appris_115",
      
      # 5. APPRIS 109 (Matches any transcript sharing the minimum rank)
      !has_canon_115 & !has_canon_109 & !has_appris_115 & has_appris_109 & !is.na(rank_109) & rank_109 == min_rank_109 ~ "best_appris_109",
      
      # 6. Singleton Fallback
      !has_canon_115 & !has_canon_109 & !has_appris_115 & !has_appris_109 & is_singleton ~ "unique_transcript",
      
      TRUE ~ NA_character_
    ),
    
    # --- Cascade Logic for Canonical Flag (1 or NA) ---
    final_is_canonical = case_when(
      !is.na(canonical_source) ~ 1,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(transcript_is_canonical = final_is_canonical) %>%
  select(
    -has_canon_115, -rank_115, -min_rank_115, -has_appris_115,
    -has_canon_109, -rank_109, -min_rank_109, -has_appris_109,
    -is_singleton, -canon_115, -appris_115, -final_is_canonical, -transcript_appris
  ) %>%
  ungroup()


library(dplyr)

# 1. Summarize the canonical counts per gene
canonical_summary <- pig_109_updated %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    # Count how many transcripts in this gene are marked as '1'
    canonical_count = sum(transcript_is_canonical == 1, na.rm = TRUE),
    total_transcripts = n(),
    .groups = 'drop' # Drops the grouping for future operations
  )

# 2. Extract the problem cases into separate data frames for manual inspection
genes_missing_canonical <- canonical_summary %>% 
  filter(canonical_count == 0)

genes_multiple_canonical <- canonical_summary %>% 
  filter(canonical_count > 1)

genes_perfect <- canonical_summary %>% 
  filter(canonical_count == 1)

# 3. Print a quick diagnostic report to the console
cat("=== Diagnostic Report ===\n")
cat("Total unique genes evaluated:         ", nrow(canonical_summary), "\n")
cat("Genes with EXACTLY 1 canonical:       ", nrow(genes_perfect), "\n")
cat("Genes with ZERO canonical transcripts:", nrow(genes_missing_canonical), "\n")
cat("Genes with MULTIPLE canonicals:       ", nrow(genes_multiple_canonical), "\n")


pig_109_updated %>% readr::write_csv("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Documents/SciLifeDrive/mammalian_RNA_atlas/wd_for_gemini/DB/pig_109_canonical_assignment.csv")





# rewrite the pig, that was missing
multispecies_canonical_table <- readr::read_csv("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Documents/SciLifeDrive/mammalian_RNA_atlas/wd_for_gemini/DB/ensembl_v109_multispecies_transcripts.csv")



multispecies_canonical_table %>% 
  filter(species != 'pig') %>% 
  rbind(
    pig_109_updated %>% 
      select(
        ensembl_gene_id, 
        ensembl_transcript_id, ensembl_transcript_id_version,
        transcript_is_canonical
      ) %>% 
      mutate(species = "pig")
  ) %>% 
  readr::write_csv("~/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Documents/SciLifeDrive/mammalian_RNA_atlas/wd_for_gemini/DB/ensembl_v109_multispecies_transcripts_updated.csv")












                            
                                     