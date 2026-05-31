library(httr)
library(jsonlite)
library(dplyr)
library(readr)

# --- Configuration ---
metadata_path <- "test_ensembl_115_pipeline/input/gene_metadata/protein_coding_genes.csv"
export_path <- "test_ensembl_115_pipeline/input/DB/ensembl_v115_multispecies_transcripts.csv"

# Read the subset of genes
genes_df <- read_csv(metadata_path)
gene_ids <- unique(genes_df$ensembl_id)

message(sprintf("Fetching Ensembl 115 transcript data for %d genes via REST API...", length(gene_ids)))

# Ensembl REST API POST lookup/id
# We can send up to 1000 IDs at a time
fetch_canonical <- function(ids) {
  url <- "https://rest.ensembl.org/lookup/id"
  body <- list(ids = ids)
  
  # Note: expand=1 gives us transcripts
  response <- POST(url, 
                   body = toJSON(body), 
                   add_headers("Content-Type" = "application/json", "Accept" = "application/json"),
                   query = list(expand = 1))
  
  if (status_code(response) != 200) {
    stop(sprintf("REST API request failed with status %d", status_code(response)))
  }
  
  data <- fromJSON(content(response, as = "text"))
  
  # Process results
  results <- lapply(names(data), function(gid) {
    gene_info <- data[[gid]]
    if (is.null(gene_info)) return(NULL)
    
    transcripts <- gene_info$Transcript
    if (is.null(transcripts) || length(transcripts) == 0) return(NULL)
    
    # Find canonical
    # In REST API response, the gene object itself often has 'canonical_transcript' field
    # or we can check the transcripts for 'is_canonical'
    canonical_id <- gene_info$canonical_transcript
    
    # Map to our required format
    transcripts_df <- as.data.frame(transcripts) %>%
      select(ensembl_transcript_id = id, ensembl_transcript_id_version = version, transcript_is_canonical = is_canonical) %>%
      mutate(ensembl_gene_id = gid,
             ensembl_transcript_id_version = as.character(ensembl_transcript_id_version))
    
    # Filter for canonical only as per requirement
    return(transcripts_df %>% filter(transcript_is_canonical == 1 | transcript_is_canonical == "1"))
  })
  
  return(do.call(rbind, results))
}

# Split into chunks of 1000
chunks <- split(gene_ids, ceiling(seq_along(gene_ids)/1000))
final_results <- do.call(rbind, lapply(chunks, fetch_canonical))

# Add species information from the original metadata
final_results <- final_results %>%
  left_join(genes_df %>% select(ensembl_id, species) %>% distinct(), by = c("ensembl_gene_id" = "ensembl_id"))

message(sprintf("Exporting %d canonical transcripts to %s", nrow(final_results), export_path))
write_csv(final_results, export_path)

cat("Done!\n")
