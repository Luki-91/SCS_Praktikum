library(Seurat)
library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)
library(xgboost)
library(e1071)
library(ROCR)
library(pROC)

# Modified prepare_ml_data function with proper Seurat object handling
prepare_ml_data <- function(seurat_obj, cluster_id) {
  message("Preparing data for cluster ", cluster_id)

  # Extract data for specified cluster using proper Seurat functions
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)

  # Get expression matrix using proper Seurat accessor
  expr_matrix <- GetAssayData(seurat_obj, layer = "data", assay = "RNA")
  expr_matrix <- expr_matrix[, cells_cluster, drop = FALSE]

  # Get metadata for the cluster
  metadata <- seurat_obj@meta.data[cells_cluster, , drop = FALSE]

  # Create response variable (0 for treated, 1 for vehicle)
  y <- as.factor(ifelse(metadata$Condition == "treated", 1, 0))

  # Check if we have enough samples
  if (length(unique(y)) != 2) {
    stop("Both conditions (vehicle and treated) must be present in the cluster")
  }

  # Split data into training and testing sets
  set.seed(42)
  trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)

  # Convert expression matrix to dense matrix if sparse
  if (inherits(expr_matrix, "dgCMatrix")) {
    expr_matrix <- as.matrix(expr_matrix)
  }

  # Split the data
  train_cells <- cells_cluster[trainIndex]
  test_cells <- cells_cluster[-trainIndex]

  # Create training and testing matrices
  train_x <- t(expr_matrix[, train_cells, drop = FALSE])
  test_x <- t(expr_matrix[, test_cells, drop = FALSE])

  # Convert to data frames and handle missing values
  train_x <- as.data.frame(train_x)
  test_x <- as.data.frame(test_x)

  # Handle NA/infinite values
  train_x[is.na(train_x)] <- 0
  test_x[is.na(test_x)] <- 0
  train_x[!is.finite(as.matrix(train_x))] <- 0
  test_x[!is.finite(as.matrix(test_x))] <- 0

  # Get feature names (genes)
  feature_names <- rownames(expr_matrix)

  # Ensure column names are set correctly
  colnames(train_x) <- feature_names
  colnames(test_x) <- feature_names

  # Create final training and testing sets
  train_y <- y[trainIndex]
  test_y <- y[-trainIndex]

  # Print some diagnostic information
  message("Data preparation summary:")
  message("Number of genes: ", length(feature_names))
  message("Number of training samples: ", nrow(train_x))
  message("Number of testing samples: ", nrow(test_x))
  message("Class distribution in training: ")
  print(table(train_y))

  return(list(
    train_x = train_x,
    test_x = test_x,
    train_y = train_y,
    test_y = test_y,
    feature_names = feature_names,
    n_features = length(feature_names),
    n_samples = length(cells_cluster)
  ))
}

# Modified train_ml_models function with additional checks
train_ml_models <- function(train_x, train_y, feature_names) {
  models <- list()

  # Check if train_x is valid
  if (is.null(train_x) || ncol(train_x) == 0 || nrow(train_x) == 0) {
    stop("Invalid training data")
  }

  # 1. Elastic Net
  message("Training Elastic Net...")
  tryCatch({
    # Ensure matrix format
    x_matrix <- as.matrix(train_x)

    cv_fit <- cv.glmnet(x = x_matrix,
                        y = train_y,
                        family = "binomial",
                        alpha = 0.5,
                        nfolds = 5)
    models$elastic_net <- cv_fit
  }, error = function(e) {
    message("Error in Elastic Net: ", e$message)
    models$elastic_net <- NULL
  })

  # 2. Random Forest
  message("Training Random Forest...")
  tryCatch({
    # Select top 1000 features based on variance
    var_features <- apply(train_x, 2, var)
    top_features <- names(sort(var_features, decreasing = TRUE))[1:min(1000, ncol(train_x))]
    rf_train_x <- train_x[, top_features, drop = FALSE]

    rf <- randomForest(x = rf_train_x,
                       y = train_y,
                       ntree = 500,
                       importance = TRUE)
    models$random_forest <- list(
      model = rf,
      features = top_features
    )
  }, error = function(e) {
    message("Error in Random Forest: ", e$message)
    models$random_forest <- NULL
  })

  # 3. XGBoost
  message("Training XGBoost...")
  tryCatch({
    if (!is.null(models$random_forest)) {
      xgb_train_x <- train_x[, top_features, drop = FALSE]

      # Convert to matrix
      xgb_train_matrix <- as.matrix(xgb_train_x)

      dtrain <- xgb.DMatrix(data = xgb_train_matrix,
                            label = as.numeric(train_y) - 1)

      params <- list(
        objective = "binary:logistic",
        eta = 0.1,
        max_depth = 6,
        subsample = 0.8,
        colsample_bytree = 0.8,
        eval_metric = "auc"
      )

      xgb_model <- xgb.train(params = params,
                             data = dtrain,
                             nrounds = 100,
                             verbose = 0)

      models$xgboost <- xgb_model
      models$xgboost_features <- top_features
    }
  }, error = function(e) {
    message("Error in XGBoost: ", e$message)
    models$xgboost <- NULL
    models$xgboost_features <- NULL
  })

  return(models)
}

# Modified evaluate_models function to handle corrected feature names
evaluate_models <- function(models, test_x, test_y, feature_names) {
  results <- list()

  # 1. Elastic Net evaluation
  if (!is.null(models$elastic_net)) {
    tryCatch({
      # Ensure column names are correct for elastic net
      x_matrix <- as.matrix(test_x)
      colnames(x_matrix) <- feature_names

      en_pred <- predict(models$elastic_net, newx = x_matrix,
                         s = "lambda.min", type = "response")
      en_roc <- roc(test_y, as.vector(en_pred))

      # Get coefficients with correct feature names
      coef_matrix <- as.matrix(coef(models$elastic_net, s = "lambda.min"))
      importance_values <- abs(coef_matrix[-1, , drop = TRUE])  # Remove intercept
      names(importance_values) <- feature_names

      results$elastic_net <- list(
        auc = auc(en_roc),
        importance = importance_values
      )
    }, error = function(e) {
      message("Error in Elastic Net evaluation: ", e$message)
    })
  }


  # 2. Random Forest evaluation (unchanged)
  if (!is.null(models$random_forest)) {
    tryCatch({
      rf_test_x <- test_x[, models$random_forest$features, drop = FALSE]
      rf_pred <- predict(models$random_forest$model, rf_test_x, type = "prob")[,2]
      rf_roc <- roc(test_y, rf_pred)
      importance_values <- importance(models$random_forest$model)[,4]
      names(importance_values) <- models$random_forest$features
      results$random_forest <- list(
        auc = auc(rf_roc),
        importance = importance_values
      )
    }, error = function(e) {
      message("Error in Random Forest evaluation: ", e$message)
    })
  }

  # 3. XGBoost evaluation with fixed importance calculation
  if (!is.null(models$xgboost) && !is.null(models$xgboost_features)) {
    tryCatch({
      # Prepare test data
      xgb_test_x <- test_x[, models$xgboost_features, drop = FALSE]

      # Make predictions
      xgb_pred <- predict(models$xgboost, as.matrix(xgb_test_x))
      xgb_roc <- roc(test_y, xgb_pred)

      # Calculate feature importance using xgb.importance
      importance_matrix <- xgb.importance(
        feature_names = models$xgboost_features,
        model = models$xgboost
      )

      # Convert importance matrix to named vector
      importance_values <- rep(0, length(models$xgboost_features))
      names(importance_values) <- models$xgboost_features

      if (nrow(importance_matrix) > 0) {
        matched_features <- match(importance_matrix$Feature, names(importance_values))
        importance_values[matched_features] <- importance_matrix$Gain
      }

      results$xgboost <- list(
        auc = auc(xgb_roc),
        importance = importance_values
      )
    }, error = function(e) {
      message("Error in XGBoost evaluation: ", e$message)
      print(str(models$xgboost_features))
    })
  }

  return(results)
}

# Modified aggregate_importance function to handle missing model results
aggregate_importance <- function(eval_results, feature_names) {
  if (length(eval_results) == 0) {
    warning("No evaluation results provided")
    return(NULL)
  }

  # Get the set of features that have importance scores
  used_features <- unique(c(
    names(eval_results$elastic_net$importance),
    names(eval_results$random_forest$importance),
    names(eval_results$xgboost$importance)
  ))

  if (length(used_features) == 0) {
    warning("No features found in model results")
    return(NULL)
  }

  # Create importance matrix
  importance_df <- data.frame(
    feature = used_features,
    elastic_net = 0,
    random_forest = 0,
    xgboost = 0,
    stringsAsFactors = FALSE
  )

  # Fill in importance scores
  if (!is.null(eval_results$elastic_net))
    importance_df$elastic_net[match(names(eval_results$elastic_net$importance), used_features)] <-
    eval_results$elastic_net$importance

  if (!is.null(eval_results$random_forest))
    importance_df$random_forest[match(names(eval_results$random_forest$importance), used_features)] <-
    eval_results$random_forest$importance

  if (!is.null(eval_results$xgboost))
    importance_df$xgboost[match(names(eval_results$xgboost$importance), used_features)] <-
    eval_results$xgboost$importance

  # Scale importance scores
  importance_df[,2:4] <- lapply(importance_df[,2:4], function(x) {
    if (all(x == 0)) return(x)
    scaled <- scale(x)
    scaled[is.na(scaled)] <- 0
    return(as.vector(scaled))
  })

  # Calculate aggregate score
  importance_df$aggregate_score <- rowMeans(importance_df[,2:4])

  # Sort by aggregate score
  importance_df <- importance_df[order(-importance_df$aggregate_score),]

  # Remove rows where all importance scores are 0
  importance_df <- importance_df[rowSums(abs(importance_df[,2:4])) > 0,]

  return(importance_df)
}

# Modified run_ml_pipeline function with better error handling
run_ml_pipeline <- function(seurat_obj, cluster_id, top_n = 100) {
  message("\nStarting ML pipeline for cluster ", cluster_id)

  # Step 1: Prepare data
  message("\nStep 1/5: Preparing data...")
  tryCatch({
    data <- prepare_ml_data(seurat_obj, cluster_id)
  }, error = function(e) {
    message("Error in data preparation: ", e$message)
    return(NULL)
  })

  # Step 2: Train models
  message("\nStep 2/5: Training models...")
  tryCatch({
    models <- train_ml_models(data$train_x, data$train_y, data$feature_names)
  }, error = function(e) {
    message("Error in model training: ", e$message)
    return(NULL)
  })

  # Step 3: Evaluate models
  message("\nStep 3/5: Evaluating models...")
  tryCatch({
    eval_results <- evaluate_models(models, data$test_x, data$test_y, data$feature_names)
  }, error = function(e) {
    message("Error in model evaluation: ", e$message)
    return(NULL)
  })

  # Step 4: Aggregate results
  message("\nStep 4/5: Aggregating results...")
  tryCatch({
    importance_df <- aggregate_importance(eval_results, data$feature_names)
    if (is.null(importance_df) || nrow(importance_df) == 0) {
      message("No significant features found")
      return(NULL)
    }
  }, error = function(e) {
    message("Error in aggregating results: ", e$message)
    return(NULL)
  })

  # Step 5: Calculate expression statistics
  message("\nStep 5/5: Calculating expression statistics...")
  tryCatch({
    # Get top genes (handle case where we have fewer than top_n genes)
    n_genes <- min(top_n, nrow(importance_df))
    top_genes <- head(importance_df$feature, n_genes)

    # Safely get cells and data
    cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
    if (length(cells_cluster) == 0) {
      message("No cells found for cluster ", cluster_id)
      return(NULL)
    }

    # Safely subset expression matrix
    expr_matrix <- GetAssayData(seurat_obj, layer = "data")
    valid_genes <- intersect(top_genes, rownames(expr_matrix))
    if (length(valid_genes) == 0) {
      message("No valid genes found in expression matrix")
      return(NULL)
    }

    expr_matrix <- as.matrix(expr_matrix[valid_genes, cells_cluster, drop = FALSE])
    metadata <- seurat_obj@meta.data[cells_cluster, , drop = FALSE]

    # Calculate statistics
    vehicle_cells <- metadata$Condition == "vehicle"
    treated_cells <- metadata$Condition == "treated"

    if (sum(vehicle_cells) == 0 || sum(treated_cells) == 0) {
      message("Missing cells for one or both conditions")
      return(NULL)
    }

    expr_stats <- data.frame(
      gene = valid_genes,
      importance_score = importance_df$aggregate_score[match(valid_genes, importance_df$feature)],
      mean_vehicle = rowMeans(expr_matrix[, vehicle_cells, drop = FALSE]),
      mean_treated = rowMeans(expr_matrix[, treated_cells, drop = FALSE])
    )

    expr_stats <- expr_stats %>%
      mutate(
        pct_expressing_vehicle = rowMeans(expr_matrix[, vehicle_cells, drop = FALSE] > 0) * 100,
        pct_expressing_treated = rowMeans(expr_matrix[, treated_cells, drop = FALSE] > 0) * 100,
        log2FC = log2((mean_vehicle + 0.1) / (mean_treated + 0.1)),
        diff_pct = pct_expressing_vehicle - pct_expressing_treated
      )

    # Save results
    write.csv2(expr_stats,
              file = paste0("cluster_", cluster_id, "_ml_signatures.csv"),
              row.names = FALSE)

    message("\nAnalysis complete!")

    return(list(
      top_genes = expr_stats,
      importance = importance_df,
      performance = data.frame(
        model = names(eval_results),
        AUC = sapply(eval_results, function(x) x$auc)
      ),
      models = models,
      eval_results = eval_results
    ))

  }, error = function(e) {
    message("Error in calculating expression statistics: ", e$message)
    return(NULL)
  })
}

# Example usage
Idents(SO.har) <- SO.har@meta.data$RNA_snn_res.0.4
cluster_id <- "5"
results <- run_ml_pipeline(SO.har, cluster_id)

# View results
if (!is.null(results)) {
  print("Model Performance:")
  print(results$performance)

  print("\nTop 20 genes:")
  print(head(results$top_genes, 20))
}

# Modified extract_and_save_results function with clearer interpretation
extract_and_save_results <- function(results, output_prefix = "cluster_analysis") {
  # 1. Model Performance Summary
  model_performance <- results$performance
  model_performance$accuracy <- round(model_performance$AUC * 100, 2)
  write.csv2(model_performance,
            file = paste0(output_prefix, "_model_performance.csv"),
            row.names = FALSE)

  # 2. Top Features Summary
  top_features <- merge(
    results$importance,
    results$top_genes,
    by.x = "feature",
    by.y = "gene",
    all.x = TRUE
  )

  # Add interpretation column
  top_features$interpretation <- ifelse(
    top_features$log2FC > 0,
    "Up in treated",
    "Down in treated"
  )

  # Reorder columns for clarity
  top_features <- top_features[, c(
    "feature",
    "aggregate_score",
    "elastic_net",
    "random_forest",
    "xgboost",
    "mean_treated",
    "mean_vehicle",
    "log2FC",
    "pct_expressing_treated",
    "pct_expressing_vehicle",
    "diff_pct",
    "interpretation"
  )]

  # Round numeric columns for readability
  numeric_cols <- sapply(top_features, is.numeric)
  top_features[numeric_cols] <- round(top_features[numeric_cols], 3)

  # Sort by aggregate score
  top_features <- top_features[order(-top_features$aggregate_score), ]

  # Save full results
  write.csv2(top_features,
            file = paste0(output_prefix, "_feature_importance.csv"),
            row.names = FALSE)

  # Create a simplified top 20 features summary
  top_20 <- head(top_features, 20)
  write.csv2(top_20,
            file = paste0(output_prefix, "_top20_features.csv"),
            row.names = FALSE)

  # Return list of data frames for immediate use
  return(list(
    model_performance = model_performance,
    top_features = top_features,
    top_20 = top_20
  ))
}

# Extract and save results
extracted_results <- extract_and_save_results(results, output_prefix = "cluster5")

# View model performance
print("Model Performance:")
print(extracted_results$model_performance)

# View top 20 features
print("\nTop 20 Features:")
print(extracted_results$top_20)

# The saved files will contain:
# 1. cluster2_model_performance.csv:
#    - Model names
#    - AUC scores
#    - Accuracy percentages

# 2. cluster2_feature_importance.csv:
#    - All features with their importance scores
#    - Expression statistics
#    - Fold changes
#    - Expression percentages

# 3. cluster2_top20_features.csv:
#    - Same columns as above but only for top 20 features