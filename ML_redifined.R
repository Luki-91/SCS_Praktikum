```R
library(Seurat)
library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)
library(xgboost)
library(e1071)
library(ROCR)
library(pROC)

# Function to prepare data for ML
prepare_ml_data <- function(seurat_obj, cluster_id) {
  message("Preparing data for cluster ", cluster_id)

  # Extract data for specified cluster
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")[, cells_cluster]
  metadata <- seurat_obj@meta.data[cells_cluster, , drop = FALSE]

  # Convert to dense matrix
  expr_matrix <- as.matrix(expr_matrix)

  # Create response variable (0 for treated, 1 for vehicle)
  y <- as.factor(ifelse(metadata$Condition == "vehicle", 1, 0))

  # Check if we have enough samples
  if (length(unique(y)) != 2) {
    stop("Both conditions (vehicle and treated) must be present in the cluster")
  }

  # Split data into training and testing sets
  set.seed(42)
  trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)

  # Transpose matrix and split into train/test
  expr_matrix_t <- t(expr_matrix)
  train_x <- expr_matrix_t[trainIndex, ]
  test_x <- expr_matrix_t[-trainIndex, ]

  # Handle NA/infinite values
  train_x[is.na(train_x)] <- 0
  test_x[is.na(test_x)] <- 0

  # Replace infinite values
  train_x[!is.finite(as.matrix(train_x))] <- 0
  test_x[!is.finite(as.matrix(test_x))] <- 0

  # Create final data frames
  train_x <- as.data.frame(train_x)
  test_x <- as.data.frame(test_x)
  train_y <- y[trainIndex]
  test_y <- y[-trainIndex]

  return(list(
    train_x = train_x,
    test_x = test_x,
    train_y = train_y,
    test_y = test_y,
    feature_names = colnames(expr_matrix),
    n_features = nrow(expr_matrix),
    n_samples = ncol(expr_matrix)
  ))
}

# Modified train_ml_models function with fixed XGBoost implementation
train_ml_models <- function(train_x, train_y, feature_names) {
  models <- list()

  # 1. Elastic Net
  message("Training Elastic Net...")
  tryCatch({
    cv_fit <- cv.glmnet(as.matrix(train_x), train_y,
                        family = "binomial",
                        alpha = 0.5,
                        nfolds = 5)
    models$elastic_net <- cv_fit
  }, error = function(e) {
    message("Error in Elastic Net: ", e$message)
    models$elastic_net <- NULL
  })

  # 2. Random Forest with reduced features
  message("Training Random Forest...")
  tryCatch({
    # Select top 1000 features
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

  # 3. XGBoost with fixed feature handling
  message("Training XGBoost...")
  tryCatch({
    # Use same features as Random Forest
    if (!is.null(models$random_forest)) {
      xgb_train_x <- train_x[, top_features, drop = FALSE]

      # Create DMatrix without explicit feature names
      dtrain <- xgb.DMatrix(data = as.matrix(xgb_train_x),
                            label = as.numeric(as.character(train_y)) - 1)

      params <- list(
        objective = "binary:logistic",
        eta = 0.1,
        max_depth = 6,
        subsample = 0.8,
        colsample_bytree = 0.8
      )

      xgb_model <- xgb.train(params = params,
                             data = dtrain,
                             nrounds = 100,
                             verbose = 0)

      # Store model and feature names separately
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

# Modified evaluate_models function with updated XGBoost importance calculation
evaluate_models <- function(models, test_x, test_y, feature_names) {
  results <- list()

  # 1. Elastic Net evaluation
  if (!is.null(models$elastic_net)) {
    tryCatch({
      en_pred <- predict(models$elastic_net, newx = as.matrix(test_x),
                         s = "lambda.min", type = "response")
      en_roc <- roc(test_y, as.vector(en_pred))
      coef_matrix <- as.matrix(coef(models$elastic_net, s = "lambda.min"))
      importance_values <- abs(coef_matrix[-1, , drop = TRUE])
      names(importance_values) <- feature_names
      results$elastic_net <- list(
        auc = auc(en_roc),
        importance = importance_values
      )
    }, error = function(e) {
      message("Error in Elastic Net evaluation: ", e$message)
    })
  }

  # 2. Random Forest evaluation
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
      # Prepare test data using the same features as training
      xgb_test_x <- test_x[, models$xgboost_features, drop = FALSE]

      # Make predictions
      xgb_pred <- predict(models$xgboost, as.matrix(xgb_test_x))
      xgb_roc <- roc(test_y, xgb_pred)

      # Calculate feature importance using a different method
      importance_matrix <- xgb.importance(
        model = models$xgboost,
        feature_names = models$xgboost_features
      )

      # Create importance vector
      importance_values <- rep(0, length(models$xgboost_features))
      names(importance_values) <- models$xgboost_features

      # Fill in importance values
      if (nrow(importance_matrix) > 0) {
        importance_values[importance_matrix$Feature] <- importance_matrix$Gain
      }

      results$xgboost <- list(
        auc = auc(xgb_roc),
        importance = importance_values
      )
    }, error = function(e) {
      message("Error in XGBoost evaluation: ", e$message)
    })
  }

  return(results)
}

# Helper function to safely get importance scores
get_safe_importance <- function(model_results, feature_set) {
  # Initialize vector of zeros with names
  imp <- numeric(length(feature_set))
  names(imp) <- feature_set

  # If no model results, return zeros
  if (is.null(model_results) || is.null(model_results$importance)) {
    return(imp)
  }

  # Get importance scores
  model_imp <- model_results$importance

  # Convert to named vector if it isn't already
  if (is.null(names(model_imp))) {
    if (length(model_imp) > length(feature_set)) {
      # Truncate if longer than feature_set
      model_imp <- model_imp[1:length(feature_set)]
      names(model_imp) <- feature_set
    } else if (length(model_imp) < length(feature_set)) {
      # Pad with zeros if shorter
      temp_imp <- numeric(length(feature_set))
      temp_imp[1:length(model_imp)] <- model_imp
      model_imp <- temp_imp
      names(model_imp) <- feature_set
    } else {
      # Same length, just add names
      names(model_imp) <- feature_set
    }
  }

  # Match features and fill importance scores
  common_features <- intersect(names(model_imp), feature_set)
  imp[common_features] <- model_imp[common_features]

  return(imp)
}

# Modified aggregate_importance function with better feature handling
aggregate_importance <- function(eval_results, feature_names) {
  if (length(eval_results) == 0) {
    warning("No evaluation results provided")
    return(NULL)
  }

  # Get the set of features actually used in the models
  used_features <- unique(c(
    names(eval_results$elastic_net$importance),
    names(eval_results$random_forest$importance),
    names(eval_results$xgboost$importance)
  ))

  if (length(used_features) == 0) {
    warning("No features found in model results")
    return(NULL)
  }

  # Create importance matrix with only used features
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

  # Sort by aggregate score and remove zero-importance features
  importance_df <- importance_df[order(-importance_df$aggregate_score),]
  importance_df <- importance_df[rowSums(abs(importance_df[,2:4])) > 0,]

  if (nrow(importance_df) == 0) {
    warning("No features with non-zero importance scores found")
    return(NULL)
  }

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
    expr_matrix <- GetAssayData(seurat_obj, slot = "data")
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
    write.csv(expr_stats,
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
cluster_id <- "2"
results <- run_ml_pipeline(SO.har, cluster_id)

# View results
if (!is.null(results)) {
  print("Model Performance:")
  print(results$performance)

  print("\nTop 20 genes:")
  print(head(results$top_genes, 20))
}
```