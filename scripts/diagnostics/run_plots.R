source("scripts/diagnostics/plot_archive_trees.R")

# Define families to plot (Modernized IDs)
ids <- c("OG_01059", "OG_00244", "OG_06940", "OG_00744")

# Create output directory
dir.create("plots_diagnostic", showWarnings = FALSE)

for (id in ids) {
  message("Plotting ", id, "...")
  plots <- plot_archived_family(id)
  
  if (!is.null(plots)) {
    # Check if we got a list (tree+msa) or just a tree
    if (is.list(plots) && !is.null(plots$tree)) {
        p_tree <- plots$tree
        p_msa <- plots$msa
    } else {
        p_tree <- plots
        p_msa <- NULL
    }
    
    ggsave(filename = file.path("plots_diagnostic", paste0(id, "_tree.pdf")), 
           plot = p_tree, width = 10, height = 8)
    
    if (!is.null(p_msa)) {
      ggsave(filename = file.path("plots_diagnostic", paste0(id, "_msa.pdf")), 
             plot = p_msa, width = 12, height = 6)
    }
  }
}

message("Done! Plots saved in ensembl_orthology_pipeline/plots_diagnostic/")
