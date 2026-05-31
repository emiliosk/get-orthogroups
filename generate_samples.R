source("scripts/08_plot_trees.R")

# Load data
df <- load_master_for_plots()
og_ids <- c("OG_00001", "OG_00003", "OG_00005", "OG_00010")

# Set paths
base_tree_dir <- "sample_plots/Phylogenetic_Trees/Gene_Trees"
base_msa_dir <- "sample_plots/Phylogenetic_Trees/MSAs"
output_dir <- "sample_plots"

for (id in og_ids) {
  # 1. Plot Tree
  tree_path <- file.path(base_tree_dir, paste0(id, ".nwk"))
  if (file.exists(tree_path)) {
    message("Plotting Tree for ", id, "...")
    p_tree <- plot_atlas_tree(id, df, tree_path)
    ggsave(file.path(output_dir, paste0(id, "_tree.pdf")), p_tree, width=10, height=6)
  }
  
  # 2. Plot MSA
  msa_path <- file.path(base_msa_dir, paste0(id, ".fa"))
  if (file.exists(msa_path)) {
    message("Plotting MSA for ", id, "...")
    p_msa <- plot_atlas_msa(id, df, msa_path)
    ggsave(file.path(output_dir, paste0(id, "_msa.pdf")), p_msa, width=10, height=6)
  }
}
