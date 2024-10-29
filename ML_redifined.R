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
  expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data")[, cells_cluster])
  metadata <- seurat_obj@meta.data[cells_cluster, , drop = FALSE]

  # Transpose matrix for ML (samples as rows)
  expr_matrix_t <- t(expr_matrix)

  # Create response variable (0 for treated, 1 for vehicle)
  y <- as.factor(ifelse(metadata$Condition == "vehicle", 1, 0))

  # Check if we have enough samples
  if (length(unique(y)) != 2) {
    stop("Both conditions (vehicle and treated) must be present in the cluster")
  }

  # Split data into training and testing sets
  set.seed(42)
  trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)

  # Prepare training and testing sets
  train_x <- expr_matrix_t[trainIndex, ]
  test_x <- expr_matrix_t[-trainIndex, ]
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

# Function to train multiple ML models
train_ml_models <- function(train_x, train_y, feature_names) {
  models <- list()

  # 1. Elastic Net
  message("Training Elastic Net...")
  set.seed(42)
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

  # 2. Random Forest
  message("Training Random Forest...")
  set.seed(42)
  tryCatch({
    rf <- randomForest(x = train_x,
                       y = train_y,
                       ntree = 500,
                       importance = TRUE)
    models$random_forest <- rf
  }, error = function(e) {
    message("Error in Random Forest: ", e$message)
    models$random_forest <- NULL
  })

  # 3. XGBoost
  message("Training XGBoost...")
  set.seed(42)
  tryCatch({
    dtrain <- xgb.DMatrix(data = as.matrix(train_x),
                          label = as.numeric(train_y) - 1)
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
    models$xgboost <- xgb_model
  }, error = function(e) {
    message("Error in XGBoost: ", e$message)
    models$xgboost <- NULL
  })

  return(models)
}

# Function to evaluate models
evaluate_models <- function(models, test_x, test_y, feature_names) {
  results <- list()

  # 1. Elastic Net evaluation
  if (!is.null(models$elastic_net)) {
    tryCatch({
      en_pred <- predict(models$elastic_net, newx = as.matrix(test_x),
                         s = "lambda.min", type = "response")
      en_roc <- roc(test_y, as.vector(en_pred))
      results$elastic_net <- list(
        auc = auc(en_roc),
        importance = abs(coef(models$elastic_net, s = "lambda.min"))[-1]
      )
    }, error = function(e) {
      message("Error in Elastic Net evaluation: ", e$message)
    })
  }

  # 2. Random Forest evaluation
  if (!is.null(models$random_forest)) {
    tryCatch({
      rf_pred <- predict(models$random_forest, test_x, type = "prob")[,2]
      rf_roc <- roc(test_y, rf_pred)
      results$random_forest <- list(
        auc = auc(rf_roc),
        importance = importance(models$random_forest)[,4]
      )
    }, error = function(e) {
      message("Error in Random Forest evaluation: ", e$message)
    })
  }

  # 3. XGBoost evaluation
  if (!is.null(models$xgboost)) {
    tryCatch({
      xgb_pred <- predict(models$xgboost, as.matrix(test_x))
      xgb_roc <- roc(test_y, xgb_pred)
      importance_matrix <- xgb.importance(feature_names = feature_names,
                                          model = models$xgboost)
      results$xgboost <- list(
        auc = auc(xgb_roc),
        importance = setNames(importance_matrix$Gain, importance_matrix$Feature)
      )
    }, error = function(e) {
      message("Error in XGBoost evaluation: ", e$message)
    })
  }

  return(results)
}

# Function to aggregate feature importance
aggregate_importance <- function(eval_results, feature_names) {
  # Initialize importance matrix
  importance_df <- data.frame(
    feature = feature_names,
    elastic_net = 0,
    random_forest = 0,
    xgboost = 0
  )

  # Fill in importance scores with error handling
  if (!is.null(eval_results$elastic_net)) {
    importance_df$elastic_net <- eval_results$elastic_net$importance[feature_names]
  }
  if (!is.null(eval_results$random_forest)) {
    importance_df$random_forest <- eval_results$random_forest$importance[feature_names]
  }
  if (!is.null(eval_results$xgboost)) {
    importance_df$xgboost <- eval_results$xgboost$importance[feature_names]
  }

  # Scale importances
  importance_df[,2:4] <- apply(importance_df[,2:4], 2, function(x) {
    scaled <- scale(x)
    scaled[is.na(scaled)] <- 0
    return(scaled)
  })

  # Calculate aggregate score
  importance_df$aggregate_score <- rowMeans(importance_df[,2:4])

  # Sort by aggregate score
  importance_df <- importance_df[order(-importance_df$aggregate_score),]

  return(importance_df)
}

# Main function to run the ML pipeline
run_ml_pipeline <- function(seurat_obj, cluster_id, top_n = 100) {
  # Prepare data
  message("Preparing data...")
  data <- prepare_ml_data(seurat_obj, cluster_id)

  # Train models
  message("Training models...")
  models <- train_ml_models(data$train_x, data$train_y, data$feature_names)

  # Evaluate models
  message("Evaluating models...")
  eval_results <- evaluate_models(models, data$test_x, data$test_y, data$feature_names)

  # Aggregate results
  message("Aggregating results...")
  importance_df <- aggregate_importance(eval_results, data$feature_names)

  # Get top genes
  top_genes <- head(importance_df, top_n)

  # Calculate expression statistics
  message("Calculating expression statistics...")
  cells_cluster <- WhichCells(seurat_obj, idents = cluster_id)
  expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data")[top_genes$feature, cells_cluster])
  metadata <- seurat_obj@meta.data[cells_cluster, , drop = FALSE]

  # Calculate statistics
  expr_stats <- data.frame(
    gene = top_genes$feature,
    importance_score = top_genes$aggregate_score,
    mean_vehicle = rowMeans(expr_matrix[, metadata$Condition == "vehicle", drop = FALSE]),
    mean_treated = rowMeans(expr_matrix[, metadata$Condition == "treated", drop = FALSE])
  ) %>%
    mutate(
      pct_expressing_vehicle = rowMeans(expr_matrix[, metadata$Condition == "vehicle", drop = FALSE] > 0) * 100,
      pct_expressing_treated = rowMeans(expr_matrix[, metadata$Condition == "treated", drop = FALSE] > 0) * 100,
      log2FC = log2((mean_vehicle + 0.1) / (mean_treated + 0.1)),
      diff_pct = pct_expressing_vehicle - pct_expressing_treated
    )

  # Calculate model performance
  performance <- data.frame(
    model = character(),
    AUC = numeric()
  )

  if (!is.null(eval_results$elastic_net)) {
    performance <- rbind(performance,
                         data.frame(model = "Elastic Net",
                                    AUC = eval_results$elastic_net$auc))
  }
  if (!is.null(eval_results$random_forest)) {
    performance <- rbind(performance,
                         data.frame(model = "Random Forest",
                                    AUC = eval_results$random_forest$auc))
  }
  if (!is.null(eval_results$xgboost)) {
    performance <- rbind(performance,
                         data.frame(model = "XGBoost",
                                    AUC = eval_results$xgboost$auc))
  }

  # Save results
  write.csv(expr_stats,
            file = paste0("cluster_", cluster_id, "_ml_signatures.csv"),
            row.names = FALSE)

  return(list(
    top_genes = expr_stats,
    importance = importance_df,
    performance = performance,
    models = models,
    eval_results = eval_results
  ))
}

# Example usage
# Set proper identities if needed
cluster_id <- "3"
results <- run_ml_pipeline(SO.har, cluster_id)

# View results
if (!is.null(results)) {
  # Print model performance
  print("Model Performance:")
  print(results$performance)

  # Print top genes
  print("\nTop 20 genes:")
  print(head(results$top_genes, 20))

  # Create visualizations
  top_20_genes <- head(results$top_genes$gene, 20)
  print(VlnPlot(SO.har,
                features = head(top_20_genes, 6),
                split.by = "Condition",
                idents = cluster_id,
                ncol = 3,
                pt.size = 0))

  # Compare with GMM results if available
  gmm_results <- extract_signatures(SO.har, cluster_id)
  if (!is.null(gmm_results)) {
    common_genes <- intersect(results$top_genes$gene[1:100],
                              gmm_results$gene[1:100])
    print(paste("\nNumber of genes in common between ML and GMM approaches:",
                length(common_genes)))
    print("Common genes:")
    print(common_genes)
  }
}