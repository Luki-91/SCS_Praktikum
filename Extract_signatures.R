library(Seurat)
library(mclust)
library(tidyverse)
library(Matrix)
library(ComplexHeatmap)
library(circlize)

# Fixed function to calculate GMM scores
calculate_gmm_scores <- function(seurat_obj, cluster_id) {
  # Get all genes
  all_genes <- rownames(seurat_obj)

  # Initialize results storage
  results <- data.frame(
    gene = all_genes,
    complexity_score = 0,
    distribution_change = 0,
    variance_ratio = 0,
    mean_change = 0,
    stringsAsFactors = FALSE
  )

  # Filter cells for specific cluster
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")[, cells_cluster]
  metadata <- seurat_obj@meta.data[cells_cluster, ]

  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(all_genes), style = 3)

  # Calculate scores for each gene
  for (i in seq_along(all_genes)) {
    gene <- all_genes[i]

    # Update progress bar
    setTxtProgressBar(pb, i)

    tryCatch({
      # Get expression values
      vehicle_expr <- expr_matrix[gene, metadata$Condition == "vehicle"]
      treated_expr <- expr_matrix[gene, metadata$Condition == "treated"]

      # Skip genes with no variation in either condition
      if (var(vehicle_expr) == 0 || var(treated_expr) == 0) {
        next
      }

      # Only fit GMM if there are enough non-zero values
      if (sum(vehicle_expr > 0) > 10 && sum(treated_expr > 0) > 10) {
        # Fit GMMs with error handling
        vehicle_gmm <- tryCatch(
          Mclust(vehicle_expr, G = 1:5, verbose = FALSE),
          error = function(e) NULL
        )

        treated_gmm <- tryCatch(
          Mclust(treated_expr, G = 1:5, verbose = FALSE),
          error = function(e) NULL
        )

        # Calculate complexity score only if both GMMs fitted successfully
        if (!is.null(vehicle_gmm) && !is.null(treated_gmm)) {
          results$complexity_score[i] <- vehicle_gmm$G - treated_gmm$G
          results$distribution_change[i] <- abs(vehicle_gmm$BIC - treated_gmm$BIC)
        }
      }

      # Calculate basic statistics that don't require GMM
      results$variance_ratio[i] <- log2(var(vehicle_expr) / max(var(treated_expr), 1e-10))
      results$mean_change[i] <- log2(mean(vehicle_expr) / max(mean(treated_expr), 1e-10))

    }, error = function(e) {
      # If there's an error, keep default values of 0
      message(sprintf("Error processing gene %s: %s", gene, e$message))
    })
  }

  # Close progress bar
  close(pb)

  # Remove rows with NaN or Inf values
  results <- results[!apply(results[,-1], 1, function(x) any(is.nan(x) | is.infinite(x))), ]

  return(results)
}

# Function to identify signature genes with robust error handling
extract_signatures <- function(seurat_obj, cluster_id, n_top = 100) {
  # Calculate GMM scores
  message("Calculating GMM scores...")
  gmm_scores <- calculate_gmm_scores(seurat_obj, cluster_id)

  # Normalize scores with error handling
  message("Normalizing scores...")
  gmm_scores_normalized <- gmm_scores %>%
    mutate(across(-gene, function(x) {
      scaled <- scale(x)
      # Replace any NaN or Inf values with 0
      scaled[is.nan(scaled) | is.infinite(scaled)] <- 0
      return(scaled)
    })) %>%
    mutate(composite_score = complexity_score + distribution_change +
             variance_ratio + mean_change)

  # Get top genes
  top_genes <- gmm_scores_normalized %>%
    arrange(desc(composite_score)) %>%
    head(n_top)

  # Additional characterization of top genes
  message("Calculating additional statistics...")
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")[top_genes$gene, cells_cluster]
  metadata <- seurat_obj@meta.data[cells_cluster, ]

  # Calculate additional metrics with error handling
  detailed_stats <- data.frame(
    gene = top_genes$gene,
    composite_score = top_genes$composite_score,
    pct_expressing_vehicle = apply(expr_matrix[, metadata$Condition == "vehicle"], 1,
                                   function(x) mean(x > 0) * 100),
    pct_expressing_treated = apply(expr_matrix[, metadata$Condition == "treated"], 1,
                                   function(x) mean(x > 0) * 100)
  )

  return(detailed_stats)
}

# Fixed plot_signature_heatmap function for S4 objects
plot_signature_heatmap <- function(seurat_obj, signature_genes, cluster_id) {
  message("Extracting data for heatmap...")

  # Get cluster cells
  cluster_cells <- WhichCells(seurat_obj, idents = cluster_id)

  # Subset the Seurat object first
  seurat_subset <- subset(seurat_obj, cells = cluster_cells)

  # Extract expression data for signature genes using proper subsetting
  expr_matrix <- GetAssayData(seurat_subset, slot = "data")
  expr_matrix <- as.matrix(expr_matrix[rownames(expr_matrix) %in% signature_genes, ])

  # Get condition information
  cell_conditions <- seurat_subset$Condition

  # Scale the data
  expr_matrix_scaled <- t(scale(t(expr_matrix)))

  # Create annotation for columns (cells)
  ha <- HeatmapAnnotation(
    df = data.frame(Condition = cell_conditions),
    col = list(Condition = c(
      "vehicle" = "#009900",
      "treated" = "#FF0000"
    ))
  )

  # Calculate log2FC for row annotation
  vehicle_cells <- cell_conditions == "vehicle"
  treated_cells <- cell_conditions == "treated"

  vehicle_means <- rowMeans(expr_matrix[, vehicle_cells, drop = FALSE])
  treated_means <- rowMeans(expr_matrix[, treated_cells, drop = FALSE])
  log2fc <- log2((vehicle_means + 0.1) / (treated_means + 0.1))

  # Create row annotation
  row_ha <- rowAnnotation(
    log2FC = log2fc,
    col = list(log2FC = colorRamp2(
      c(min(log2fc), 0, max(log2fc)),
      c("blue", "white", "red")
    ))
  )

  message("Creating heatmap...")

  # Create heatmap
  heatmap <- Heatmap(
    expr_matrix_scaled,
    name = "Scaled Expression",
    top_annotation = ha,
    right_annotation = row_ha,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_title = paste("Signature Genes -", cluster_id),
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    col = colorRamp2(
      c(-2, 0, 2),
      c("#3366CC", "white", "#CC3366")
    )
  )

  return(heatmap)
}

# Updated analyze_cluster_signatures function
analyze_cluster_signatures <- function(seurat_obj, cluster_id, n_top = 100) {
  # Extract signatures
  message("Extracting signatures...")
  signatures <- extract_signatures(seurat_obj, cluster_id, n_top)

  if (is.null(signatures) || nrow(signatures) == 0) {
    message("No significant signatures found")
    return(NULL)
  }

  # Create visualizations
  message("Creating heatmap...")
  heatmap <- plot_signature_heatmap(seurat_obj, signatures$gene, cluster_id)

  # Calculate expression patterns
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  seurat_subset <- subset(seurat_obj, cells = cells_cluster)
  expr_matrix <- GetAssayData(seurat_subset, slot = "data")
  expr_matrix <- as.matrix(expr_matrix[rownames(expr_matrix) %in% signatures$gene, ])
  metadata <- seurat_subset@meta.data

  # Create pattern groups
  expression_patterns <- data.frame(
    gene = signatures$gene,
    pattern = case_when(
      signatures$pct_expressing_vehicle > 75 &
        signatures$pct_expressing_treated < 25 ~ "Lost in treated",
      signatures$pct_expressing_vehicle > signatures$pct_expressing_treated * 2 ~ "Strongly reduced",
      signatures$pct_expressing_vehicle > signatures$pct_expressing_treated * 1.5 ~ "Moderately reduced",
      TRUE ~ "Other"
    )
  )

  # Save results
  write.csv(signatures,
            file = paste0("cluster_", cluster_id, "_signatures.csv"),
            row.names = FALSE)

  return(list(
    signatures = signatures,
    heatmap = heatmap,
    patterns = expression_patterns
  ))
}

# Example usage
cluster_id <- "1"  # Replace with your cluster of interest
Idents(SO.har) <- SO.har@meta.data$RNA_snn_res.0.4
results <- analyze_cluster_signatures(SO.har, cluster_id)

# Print top genes
print(head(results$signatures, 20))

# Show heatmap
print(results$heatmap)

# Print pattern distribution
print(table(results$patterns$pattern))

# Export results to file
write.csv2(results$signatures, file = paste0("CD8_cluster_", cluster_id, "_signatures.csv"))
write.csv2(results$patterns, file = paste0("CD8_cluster_", cluster_id, "_patterns.csv"))
