
devtools::load_all("/Users/emilioskarwan/Documents/SciLifeDrive/AnnDatR")
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(biomaRt)
library(purrr)

select <- dplyr::select

# Relative to the preprocessing/ folder
data_loc <- "../../processed_data"
export_dir <- "../input/gene_metadata"
selected_level <- 'region'

species = c("mfascicularis" = "macaque", "mmusculus" = 'mouse', "rnorvegicus" = 'rat', 
            "sscrofa" = 'pig', 'hsapiens' = 'human')

dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# Generate species consensus dataset and normalise --------------------------
list_vars <- species %>%   imap(\(name, code){
  
  adata <- AnnDatR$new(
    prefix_name = selected_level,
    layer = "nTPM",
    var_names = "ensembl_id",
    file_dir = file.path(data_loc, name)
  )
  if(nrow(adata$var) != nrow(adata$X)) {
    stop(paste0("Number of genes in var does not match number of genes in X for ", x, ": ", nrow(adata$var), " vs ", nrow(adata$X)))
  }
  adata$var$species <- name
  
  # get protein id's for each gene in each species 
  ensembl_host <- "https://feb2023.archive.ensembl.org"
  dataset_name <- paste0(code, "_gene_ensembl")
  target_attributes <- c(
    "ensembl_gene_id", 
    "ensembl_peptide_id", 
    "ensembl_peptide_id_version"
  )
  mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = dataset_name,
    host = ensembl_host
  )
  peptide_ids <- getBM(
    attributes = target_attributes,
    mart = mart
  ) %>% rename(
    ensembl_id = ensembl_gene_id,
    peptide_id = ensembl_peptide_id,
    peptide_id_version = ensembl_peptide_id_version
  )
  adata$var <- adata$var %>% left_join(peptide_ids, by = "ensembl_id")
  
  return(adata$var %>% select(ensembl_id, gene_symbol, peptide_id, peptide_id_version, species))
})

# rbind list_vars
list_vars <- list_vars %>% 
  map_dfr(\(vars) vars)

# if gene symbol NA, take ensembl_id
list_vars <- list_vars %>% 
  mutate(gene_symbol = case_when(
    is.na(gene_symbol) ~ ensembl_id, 
    .default = gene_symbol
  ))

list_vars %>% write_csv(file.path(export_dir, 'protein_coding_genes.csv'))
