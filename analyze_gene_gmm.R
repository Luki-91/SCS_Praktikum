library(Seurat)
library(mclust)
library(ggplot2)
library(tidyr)
library(dplyr)

analyze_gene_gmm <- function(seurat_obj, gene, cluster_id) {
  tryCatch({
    # First, check if the gene exists in the dataset
    if (!gene %in% rownames(seurat_obj)) {
      stop(paste("Gene", gene, "not found in the dataset"))
    }

    # Extract expression data for the specified gene
    expression_data <- GetAssayData(seurat_obj, slot = "data")[gene, ]

    # Create initial data frame
    expr_df <- data.frame(
      expression = as.numeric(expression_data),
      condition = seurat_obj$Condition,
      cluster = as.character(Idents(seurat_obj))
    )

    # Print debugging information
    print(paste("Total cells:", nrow(expr_df)))
    print(paste("Unique clusters:", paste(unique(expr_df$cluster), collapse = ", ")))
    print(paste("Target cluster:", cluster_id))

    # Filter for the specific cluster
    expr_df <- expr_df %>%
      filter(cluster == as.character(cluster_id))

    print(paste("Cells in cluster", cluster_id, ":", nrow(expr_df)))

    # Split data by condition
    vehicle_expr <- expr_df$expression[expr_df$condition == "vehicle"]
    treated_expr <- expr_df$expression[expr_df$condition == "treated"]

    print(paste("Vehicle cells:", length(vehicle_expr)))
    print(paste("Treated cells:", length(treated_expr)))

    # Print expression statistics
    print(paste("Non-zero vehicle cells:", sum(vehicle_expr > 0)))
    print(paste("Non-zero treated cells:", sum(treated_expr > 0)))

    # Fit GMMs only on non-zero values
    vehicle_gmm <- Mclust(vehicle_expr[vehicle_expr > 0], G = 1:2)
    treated_gmm <- Mclust(treated_expr[treated_expr > 0], G = 1:2)

    # Create single-row metrics data frame
    metrics <- data.frame(
      gene = gene,
      n_vehicle = length(vehicle_expr),
      n_treated = length(treated_expr),
      n_nonzero_vehicle = sum(vehicle_expr > 0),
      n_nonzero_treated = sum(treated_expr > 0),
      mean_vehicle = mean(vehicle_expr),
      mean_treated = mean(treated_expr),
      mean_nonzero_vehicle = mean(vehicle_expr[vehicle_expr > 0]),
      mean_nonzero_treated = mean(treated_expr[treated_expr > 0]),
      pct_expressing_vehicle = mean(vehicle_expr > 0) * 100,
      pct_expressing_treated = mean(treated_expr > 0) * 100,
      components_vehicle = ifelse(!is.null(vehicle_gmm), vehicle_gmm$G, NA),
      components_treated = ifelse(!is.null(treated_gmm), treated_gmm$G, NA)
    )

    # Create plot
    plot_data <- expr_df %>%
      mutate(
        condition = factor(condition, levels = c("vehicle", "treated")),
        nonzero = expression > 0
      )

    p <- ggplot(plot_data, aes(x = condition, y = expression)) +
      # Add violin plot
      geom_violin(aes(fill = condition), alpha = 0.3) +
      # Add jittered points, colored by zero/non-zero
      geom_jitter(aes(color = nonzero), width = 0.2, size = 0.5, alpha = 0.5) +
      # Customize appearance
      scale_fill_manual(values = c("vehicle" = "#99CC99", "treated" = "#FF9999")) +
      scale_color_manual(values = c("TRUE" = "#FF0000", "FALSE" = "#CCCCCC")) +
      theme_minimal() +
      labs(
        title = paste("Expression Distribution for", gene),
        subtitle = paste0(
          "Cluster ", cluster_id, "\n",
          "Vehicle: ", round(metrics$pct_expressing_vehicle, 1), "% expressing\n",
          "Treated: ", round(metrics$pct_expressing_treated, 1), "% expressing"
        ),
        x = "Condition",
        y = "Expression Level"
      ) +
      theme(legend.position = "none")

    # Return results
    result <- list(
      plot = p,
      metrics = metrics,
      vehicle_model = vehicle_gmm,
      treated_model = treated_gmm,
      raw_data = plot_data
    )

    class(result) <- "gene_gmm_analysis"
    return(result)

  }, error = function(e) {
    message(sprintf("Error processing gene %s: %s", gene, e$message))
    return(NULL)
  })
}

# Print method for gene_gmm_analysis
print.gene_gmm_analysis <- function(x, ...) {
  cat("Gene GMM Analysis Results\n")
  cat("-------------------------\n")
  cat(sprintf("Expression metrics:\n"))
  print(x$metrics)
  cat("\nUse $plot to view the visualization\n")
  cat("Use $raw_data to access the underlying data\n")
}
