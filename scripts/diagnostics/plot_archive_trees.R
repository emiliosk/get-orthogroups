#' Plot Trees from Tar Archive (Modernized Wrapper)
#' 
#' This script wraps the main 08_plot_trees.R logic to provide 
#' a simplified interface for diagnostic plotting from the archive.

# Load the main plotting library
# Note: we assume being run from the pipeline root or diagnostics dir
if(file.exists("scripts/08_plot_trees.R")) {
  source("scripts/08_plot_trees.R")
} else {
  source("../08_plot_trees.R")
}

#' Plot a family tree directly from the archive
#' (This is a compatibility wrapper for older scripts)
plot_archived_family <- function(og_id) {
  
  df <- load_master_for_plots()
  
  # 1. Plot Tree
  p_tree <- plot_atlas_tree(og_id, df)
  
  # 2. Plot MSA
  p_msa <- plot_atlas_msa(og_id, df, start_pos=1, end_pos=100)
  
  if (!is.null(p_msa)) {
    return(list(tree = p_tree, msa = p_msa))
  }
  
  return(p_tree)
}
