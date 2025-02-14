#' Tune Hyperparameters for bHIVE (Swarm/Grid Search)
#'
#' Performs hyperparameter tuning for the bHIVE algorithm over a grid of 
#' hyperparameter values or an externally provided data frame of parameter 
#' combinations. Evaluates each combination using different metrics:
#'
#' - **Classification**: "accuracy", "balanced_accuracy", "f1", "kappa"
#' - **Regression**: "rmse", "mae", "r2"
#' - **Clustering**: "silhouette", "davies_bouldin", or "calinski_harabasz" 
#'
#' **Note**: Some metrics require additional packages or assumptions 
#' (e.g., multi-class classification for "f1" is calculated as a macro-average).
#'
#' @param X A numeric matrix or data frame of input features (rows = 
#' observations, columns = features).
#' @param y Optional. A target vector: factor for classification, numeric for 
#' regression. 
#'   If \code{NULL}, clustering is performed.
#' @param task Character. One of \code{"clustering"}, \code{"classification"},
#'  or \code{"regression"}.
#' @param grid A data frame specifying the hyperparameter combinations. 
#' Should have columns: \code{nAntibodies}, \code{beta}, \code{epsilon}. 
#' (Optionally more if you want to pass other arguments to \code{bHIVE()}.)
#' @param metric Character. Name of the evaluation metric. Options:
#'   \itemize{
#'     \item \strong{Classification}: "accuracy", "balanced_accuracy", "f1", "kappa"
#'     \item \strong{Regression}: "rmse", "mae", "r2"
#'     \item \strong{Clustering}: "silhouette", "davies_bouldin", "calinski_harabasz"
#'   }
#' @param maxIter Integer. Maximum iterations for each \code{bHIVE} run 
#' (default 50).
#' @param BPPARAM Character. A BiocParallel::bpparam() object that can be used 
#' for parallelization. The function supports \code{SerialParam}, \code{MulticoreParam}, 
#' \code{BatchtoolsParam}, and \code{SnowParam}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A list:
#'   \itemize{
#'     \item \code{best_params}: A list (row) of the best hyperparameters.
#'     \item \code{results}: A data frame with the full grid search results, 
#'       including the \code{metric_value} for each combination.
#'   }
#'
#' @examples
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' y <- iris$Species  # classification
#'
#' # Define hyperparameter grid
#' grid <- expand.grid(
#'   nAntibodies = c(10, 20),
#'   beta        = c(3, 5),
#'   epsilon     = c(0.01, 0.05)
#' )
#'
#' # Perform hyperparameter tuning for classification
#' tuning_results <- swarmbHIVE(X = X, 
#'                              y = y, 
#'                              task = "classification", 
#'                              grid = grid, 
#'                              metric = "balanced_accuracy",
#'                              maxIter = 10)
#'
#' # For clustering with silhouette
#' set.seed(42)
#' X_clust <- matrix(rnorm(100 * 5), ncol = 5)
#' grid_clust <- expand.grid(nAntibodies = c(5, 10), 
#'                           beta = c(3, 5), 
#'                           epsilon = c(0.01, 0.05))
#' res_clust <- swarmbHIVE(X_clust, 
#'                         task = "clustering", 
#'                         grid = grid_clust, 
#'                         metric = "silhouette")
#' res_clust$best_params
#' 
#'
#' @importFrom stats dist
#' @importFrom cluster silhouette
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam BatchtoolsParam SerialParam
#' @importFrom clusterCrit intCriteria
#' @export
swarmbHIVE <- function(X, 
                       y = NULL, 
                       task = c("clustering","classification","regression"),
                       grid,
                       metric = NULL,
                       maxIter = 50,
                       BPPARAM = SerialParam(),
                       verbose = TRUE) {
  task <- match.arg(task)
  
  #-------------------
  # 1) Validate input
  #-------------------
  .validate_bHIVE_input(X, y)  # existing check from your code
  
  # If metric is not given, choose a default
  if (is.null(metric)) {
    metric <- switch(
      task,
      "clustering"     = "silhouette",
      "classification" = "accuracy",
      "regression"     = "rmse"
    )
  }
  
  # Allowed metrics by task
  valid_class_metrics <- c("accuracy", "balanced_accuracy", "f1", "kappa")
  valid_reg_metrics   <- c("rmse", "mae", "r2")
  valid_clust_metrics <- c("silhouette", "davies_bouldin", "calinski_harabasz")
  
  # Check if the chosen metric is valid for the task
  if (task == "classification" && ! metric %in% valid_class_metrics) {
    stop(sprintf("Invalid metric '%s' for classification. Valid: %s",
                 metric, paste(valid_class_metrics, collapse=", ")))
  }
  if (task == "regression" && ! metric %in% valid_reg_metrics) {
    stop(sprintf("Invalid metric '%s' for regression. Valid: %s",
                 metric, paste(valid_reg_metrics, collapse=", ")))
  }
  if (task == "clustering" && ! metric %in% valid_clust_metrics) {
    stop(sprintf("Invalid metric '%s' for clustering. Valid: %s",
                 metric, paste(valid_clust_metrics, collapse=", ")))
  }
  
  # If using silhouette or other clustering metrics, pre-compute distance matrix
  # to avoid repeated calls to dist(X).
  dist_mat <- NULL
  if (task == "clustering" && metric %in% c("silhouette","davies_bouldin","calinski_harabasz")) {
    if (verbose) message("Precomputing distance matrix for clustering metrics.")
    dist_mat <- dist(X)
  }
  
  #-------------------------
  # 2) Define metric function
  #-------------------------
  .calc_metric <- function(model, X, y, task, metric, dist_mat) {
    # Extract outputs from model
    if (task == "clustering") {
      cluster_labels <- model$assignments
      if (metric == "silhouette") {
        sil <- silhouette(cluster_labels, dist_mat)
        return(mean(sil[, "sil_width"], na.rm=TRUE))
      } else if (metric %in% c("davies_bouldin","calinski_harabasz")) {
        # Convert labels to 1-based
        part <- as.integer(factor(cluster_labels))
        cc_res <- intCriteria(traj = as.matrix(X), part = part, 
                                           crit = metric)
        return(cc_res[[metric]])
      }
    }
    
    if (task == "classification") {
      predicted_labels <- model$assignments
      actual_labels <- y
      
      if (metric == "accuracy") {
        return(mean(predicted_labels == actual_labels))
      } else if (metric == "balanced_accuracy") {
        # Balanced accuracy across classes
        # For multi-class, we can do macro-average recall
        tbl <- table(actual_labels, predicted_labels)
        # row = actual, col = predicted
        recalls <- diag(prop.table(tbl, margin=1))
        return(mean(recalls, na.rm=TRUE))
      } else if (metric == "f1") {
        # For multi-class, compute macro-F1
        # F1_class_i = 2 * precision_i * recall_i / (precision_i + recall_i)
        tb <- table(actual_labels, predicted_labels)
        # row = actual, col = pred
        f1s <- c()
        for (cl in rownames(tb)) {
          cl <- as.character(cl)
          precision <- tb[cl, cl] / sum(tb[, cl])
          recall    <- tb[cl, cl] / sum(tb[cl, ])
          if (is.na(precision) || is.na(recall) || (precision+recall)==0) {
            f1s <- c(f1s, 0)
          } else {
            f1s <- c(f1s, 2*precision*recall/(precision+recall))
          }
        }
        return(mean(f1s, na.rm=TRUE))
      } else if (metric == "kappa") {
        # Cohen's Kappa
        tb <- table(actual_labels, predicted_labels)
        n  <- sum(tb)
        p0 <- sum(diag(tb)) / n
        # Expected agreement under random chance
        pE <- sum(rowSums(tb) * colSums(tb)) / (n*n)
        return((p0 - pE) / (1 - pE))
      }
    }
    
    if (task == "regression") {
      predicted_values <- model$predictions
      actual_values <- y
      
      if (length(predicted_values) != length(actual_values)) {
        warning("Length mismatch in predicted vs. actual for regression.")
        return(NA_real_)
      }
      
      if (metric == "rmse") {
        return(sqrt(mean((predicted_values - actual_values)^2, na.rm=TRUE)))
      } else if (metric == "mae") {
        return(mean(abs(predicted_values - actual_values), na.rm=TRUE))
      } else if (metric == "r2") {
        ss_res <- sum((actual_values - predicted_values)^2, na.rm=TRUE)
        ss_tot <- sum((actual_values - mean(actual_values, na.rm=TRUE))^2, na.rm=TRUE)
        return(1 - ss_res/ss_tot)
      }
    }
    
    # fallback
    return(NA_real_)
  }
  
  #-----------------------------
  # 3) Function to run one combo
  #-----------------------------
  .evaluate_combo <- function(params_row) {
    # bHIVE with given params
    model <- bHIVE(
      X = X,
      y = y,
      task = task,
      nAntibodies = params_row$nAntibodies,
      beta        = params_row$beta,
      epsilon     = params_row$epsilon,
      maxIter     = maxIter,
      verbose     = FALSE  # override local verbose to reduce console clutter
    )
    # compute metric
    mvalue <- .calc_metric(model, X, y, task, metric, dist_mat)
    
    # Return the row with metric_value
    cbind(params_row, metric_value = mvalue)
  }
  
  #-----------------------------
  # 4) Iterate Over Param Grid
  #-----------------------------
  n_combos <- nrow(grid)
  if (verbose) {
    message(sprintf(
      "Starting swarmbHIVE with %d parameter combinations (task=%s, metric=%s).",
      n_combos, task, metric
    ))
  }

  results_list <- NULL
  results_list <- bplapply(seq_len(n_combos),
      function(i) {
        if (verbose) {
          cat(sprintf(
            "Evaluating combo %d/%d: nAntibodies=%d, beta=%d, epsilon=%.3f\n",
            i, n_combos, grid$nAntibodies[i], grid$beta[i], grid$epsilon[i]
          ))
        }
        .evaluate_combo(grid[i, ])
      },
      BPPARAM = BPPARAM
    )
  
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  #-----------------------------
  # 5) Identify Best Parameters
  #-----------------------------
  if (task == "regression") {
    # For regression metrics, smaller is better (rmse, mae) except R2 is bigger=better
    if (metric %in% c("rmse","mae")) {
      best_idx <- which.min(results_df$metric_value)
    } else if (metric == "r2") {
      best_idx <- which.max(results_df$metric_value)
    }
  } else if (task == "classification") {
    # All classification metrics: bigger is better
    best_idx <- which.max(results_df$metric_value)
  } else {
    # Clustering metrics: bigger is typically better for silhouette, calinski_harabasz
    # for davies_bouldin, smaller is better
    if (metric == "davies_bouldin") {
      best_idx <- which.min(results_df$metric_value)
    } else {
      best_idx <- which.max(results_df$metric_value)
    }
  }
  
  best_params <- results_df[best_idx, , drop=FALSE]
  
  if (verbose) {
    message("Best parameters found:")
    print(best_params)
  }
  
  # Return results
  list(
    best_params = best_params,
    results     = results_df
  )
}
