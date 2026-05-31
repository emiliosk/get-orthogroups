library(biomaRt)
library(dplyr)
library(readr)

# --- Configuration ---
ensembl_109_host <- "https://feb2023.archive.ensembl.org"
ensembl_115_host <- "https://sep2025.archive.ensembl.org"

species_datasets <- c(
  human   = "hsapiens_gene_ensembl",
  macaque = "mfascicularis_gene_ensembl",
  mouse   = "mmusculus_gene_ensembl",
  pig     = "sscrofa_gene_ensembl",
  rat     = "rnorvegicus_gene_ensembl"
)

target_attributes <- c(
  "ensembl_gene_id", 
  "ensembl_transcript_id", 
  "ensembl_transcript_id_version", 
  "transcript_is_canonical" 
)

appris_priority <- c(
  "principal1", "principal2", "principal3", 
  "principal4", "principal5", "alternative1", "alternative2"
)

# Relative to preprocessing/
export_path <- "../input/DB/test_ensembl_v109_multispecies_transcripts.csv"
dir.create(dirname(export_path), showWarnings = FALSE, recursive = TRUE)

# --- 1. Fetch Data for All Species from Ensembl 109 ---
message("Fetching Ensembl 109 data for all species...")
all_species_data <- lapply(names(species_datasets), function(sp_name) {
  message(sprintf("  Processing %s...", sp_name))
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = species_datasets[[sp_name]], host = ensembl_109_host)
  
  # For Pig, we also need APPRIS for the fallback logic
  attrs <- target_attributes
  if (sp_name == "pig") {
    attrs <- c(attrs, "transcript_appris")
  }
  
  df <- getBM(attributes = attrs, mart = mart)
  df$species <- sp_name
  return(df)
})

# --- 2. Special Logic for Pig (sscrofa) ---
# Pig in 109 has many genes without a canonical transcript flag.
message("Applying fallback logic for Pig canonical assignments...")
pig_109 <- all_species_data[[which(names(species_datasets) == "pig")]]

# Fetch Pig data from Ensembl 115 for better canonical/APPRIS references
mart_115 <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "sscrofa_gene_ensembl", host = ensembl_115_host)
pig_115 <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_is_canonical", "transcript_appris"),
  mart = mart_115
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
    has_canon_115 = any(!is.na(canon_115) & canon_115 == 1),
    rank_115 = match(appris_115, appris_priority),
    min_rank_115 = if (all(is.na(rank_115))) NA_integer_ else min(rank_115, na.rm = TRUE),
    has_appris_115 = !is.na(min_rank_115),
    
    has_canon_109 = any(!is.na(transcript_is_canonical) & transcript_is_canonical == 1),
    rank_109 = match(transcript_appris, appris_priority),
    min_rank_109 = if (all(is.na(rank_109))) NA_integer_ else min(rank_109, na.rm = TRUE),
    has_appris_109 = !is.na(min_rank_109),
    
    is_singleton = n() == 1,
    
    canonical_source = case_when(
      has_canon_115 & !is.na(canon_115) & canon_115 == 1 & !is.na(transcript_is_canonical) & transcript_is_canonical == 1 ~ "ensembl_109/115",
      has_canon_115 & !is.na(canon_115) & canon_115 == 1 ~ "ensembl_115",
      !has_canon_115 & has_canon_109 & !is.na(transcript_is_canonical) & transcript_is_canonical == 1 ~ "ensembl_109",
      !has_canon_115 & !has_canon_109 & has_appris_115 & !is.na(rank_115) & rank_115 == min_rank_115 ~ "best_appris_115",
      !has_canon_115 & !has_canon_109 & !has_appris_115 & has_appris_109 & !is.na(rank_109) & rank_109 == min_rank_109 ~ "best_appris_109",
      !has_canon_115 & !has_canon_109 & !has_appris_115 & !has_appris_109 & is_singleton ~ "unique_transcript",
      TRUE ~ NA_character_
    ),
    
    final_is_canonical = ifelse(!is.na(canonical_source), 1, NA)
  ) %>%
  # Handle ties in APPRIS by taking the first one if multiple share the min rank
  mutate(final_is_canonical = ifelse(final_is_canonical == 1 & row_number() == which(final_is_canonical == 1)[1], 1, NA)) %>%
  ungroup() %>%
  select(ensembl_gene_id, ensembl_transcript_id, ensembl_transcript_id_version, transcript_is_canonical = final_is_canonical, species)

# --- 3. Combine and Export ---
message("Combining all species...")
final_table <- do.call(rbind, lapply(all_species_data, function(df) {
  if (unique(df$species) == "pig") {
    return(pig_109_updated)
  } else {
    return(df %>% select(ensembl_gene_id, ensembl_transcript_id, ensembl_transcript_id_version, transcript_is_canonical, species))
  }
}))

message(sprintf("Exporting final table to %s", export_path))
write_csv(final_table, export_path)

cat("Done!\n")
