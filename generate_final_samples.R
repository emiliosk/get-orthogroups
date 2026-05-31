source("scripts/08_plot_trees.R")

# Load data
df <- load_master_for_plots()
# Deterministic IDs from the latest production run
og_ids_1to1 <- c("OG_04604", "OG_09643", "OG_12690", "OG_12724", "OG_14916")
og_ids_multi <- c("OG_00012", "OG_00052", "OG_01224", "OG_01342", "OG_01459")
all_ids <- c(og_ids_1to1, og_ids_multi)

# Output directory for samples
output_dir <- "final_samples"

for (id in all_ids) {
  # Determine category for subfolder
  cat_dir <- if (id %in% og_ids_1to1) file.path(output_dir, "1to1") else file.path(output_dir, "multi_copy")
  dir.create(cat_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Plot Tree (Using NULL path to trigger Archive extraction)
  message("Plotting Tree for ", id, "...")
  p_tree <- plot_atlas_tree(id, df, tree_path = NULL)
  if(!is.null(p_tree)) {
    ggsave(file.path(cat_dir, paste0(id, "_tree.pdf")), p_tree, width=10, height=6)
  }
  
  # 2. Plot MSA (Using NULL path to trigger Archive extraction)
  message("Plotting MSA for ", id, "...")
  # Export interactive HTML
  html_out <- file.path(cat_dir, paste0(id, "_msa.html"))
  export_atlas_msa_html(id, df, msa_path = NULL, html_out)
  
  # Also save a static PDF MSA for the sample folder
  p_msa <- plot_atlas_msa(id, df, msa_path = NULL, start_pos=1, end_pos=100)
  if(!is.null(p_msa)) {
    ggsave(file.path(cat_dir, paste0(id, "_msa.pdf")), p_msa, width=12, height=6)
  }
}
