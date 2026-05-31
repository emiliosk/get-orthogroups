source("scripts/diagnostics/plot_archive_trees.R")

# Target Case 2 examples from our previous analysis
# These were SF_00078 (POLR2J) and SF_00149 (LGALS)
# We will generate plots for them from the archive

# Define families to plot (Modernized IDs)
# OG_00744 (POLR2J) and OG_00450 (LGALS)
case2_ids <- c("OG_00744", "OG_00450")

# Create output directory for plots
dir.create("plots_case2", showWarnings = FALSE)

for (id in case2_ids) {
  message("Plotting ", id, "...")
  plots <- plot_archived_family(id)
  
  if (!is.null(plots)) {
    # Save the tree
    if (is.list(plots)) {
        p_tree <- plots$tree
        p_msa <- plots$msa
    } else {
        p_tree <- plots
        p_msa <- NULL
    }
    
    ggsave(filename = file.path("plots_case2", paste0(id, "_tree.pdf")), 
           plot = p_tree, width = 10, height = 8)
    
    if (!is.null(p_msa)) {
      ggsave(filename = file.path("plots_case2", paste0(id, "_msa.pdf")), 
             plot = p_msa, width = 12, height = 6)
    }
  }
}

message("Done! Plots saved in ensembl_orthology_pipeline/plots_case2/")
