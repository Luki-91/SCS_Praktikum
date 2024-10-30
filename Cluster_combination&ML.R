# Function to combine clusters and run ML pipeline
analyze_combined_clusters <- function(seurat_obj, clusters_to_combine, new_cluster_id = "combined") {
  # Create a copy of the Seurat object to avoid modifying the original
  so_modified <- seurat_obj

  # Get current identities
  current_idents <- Idents(so_modified)

  # Create new identity vector
  new_idents <- as.character(current_idents)

  # Replace specified clusters with new cluster ID
  new_idents[new_idents %in% as.character(clusters_to_combine)] <- new_cluster_id

  # Update identities
  Idents(so_modified) <- new_idents

  # Run ML pipeline on combined cluster
  results <- run_ml_pipeline(so_modified, new_cluster_id)

  return(results)
}

# Example usage:
# Combine clusters 1, 3, and 6
clusters_to_combine <- c("2", "3", "5","7")

# Run analysis
combined_results <- analyze_combined_clusters(
  seurat_obj = SO.har,
  clusters_to_combine = clusters_to_combine,
  new_cluster_id = "cluster_2357"
)

# Extract and save results
extracted_results <- extract_and_save_results(
  combined_results,
  output_prefix = "clusters_2357_combined"
)

# Optional: Print some basic statistics
print_cluster_stats <- function(seurat_obj, clusters_to_combine, new_cluster_id) {
  # Original cells per cluster
  original_stats <- table(
    Idents(seurat_obj)[Idents(seurat_obj) %in% clusters_to_combine]
  )

  # Cells per condition in combined cluster
  combined_stats <- table(
    seurat_obj$Condition[Idents(seurat_obj) %in% clusters_to_combine]
  )

  # Print statistics
  cat("\nOriginal clusters cell counts:\n")
  print(original_stats)

  cat("\nCombined cluster cells per condition:\n")
  print(combined_stats)

  cat("\nTotal cells in combined cluster:", sum(combined_stats), "\n")
}

# Print statistics
print_cluster_stats(SO.har, clusters_to_combine, "cluster_136")
