# Install biomaRt if not already available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
}

library(biomaRt)

# Ensembl Release 109 is hosted on the Feb 2023 archive
ensembl_host <- "https://feb2023.archive.ensembl.org"

# Define the Ensembl dataset names for the target species
species_datasets <- c(
  human   = "hsapiens_gene_ensembl",
  macaque = "mfascicularis_gene_ensembl",
  mouse   = "mmusculus_gene_ensembl",
  pig     = "sscrofa_gene_ensembl",
  rat     = "rnorvegicus_gene_ensembl"
)

# Define the exact biomaRt attribute names corresponding to the requested data
target_attributes <- c(
  "ensembl_gene_id", 
  "ensembl_transcript_id", 
  "ensembl_transcript_id_version", 
  "transcript_is_canonical" 
)

# Use lapply to iterate over the species names functionally
species_data_list <- lapply(names(species_datasets), function(species_name) {
  dataset_name <- species_datasets[[species_name]]
  
  message(sprintf("Querying %s dataset (%s) from Ensembl 109...", species_name, dataset_name))
  
  # Connect to the Ensembl version 109 biomart
  mart <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = dataset_name,
    host = ensembl_host
  )
  
  # Retrieve the data
  retrieved_data <- getBM(
    attributes = target_attributes,
    mart = mart
  )
  
  # Append a column to identify the species in the final merged dataset
  retrieved_data$species <- species_name
  
  return(retrieved_data)
})

# Merge the resulting list of data frames into a single master data frame
final_dataset <- do.call(rbind, species_data_list)
rownames(final_dataset) <- NULL

# Save the dataset to a CSV file to maintain the reproducible output
write.csv(final_dataset, "../DB/ensembl_v109_multispecies_transcripts.csv", row.names = FALSE)

