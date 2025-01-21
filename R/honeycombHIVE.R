#' honeycombHIVE: Multilayered Artificial Immune Network
#'
#' Implements a hierarchical Artificial Immune Network (AI-Net) algorithm with
#' multiple layers, refining antibody populations iteratively to capture finer
#' granularity in data patterns. Each layer processes the data based on the
#' clusters (or assignments) from the previous layer, then collapses them into
#' prototypes (antibodies) for the next iteration.
#'
#' For clustering, each cluster is replaced by a single representative
#' (centroid, medoid, median, or mode). For classification, we also assign a
#' label to each prototype by majority vote; for regression, we assign an
#' average target value to each prototype.
#'
#' @param X A numeric matrix or data frame of input features (rows = observations,
#' columns = features).
#' @param y (Optional) A target vector for classification (factor) or regression
#' (numeric). Ignored if \code{task = "clustering"} or if set to NULL.
#' @param task Character, one of \code{"clustering"}, \code{"classification"},
#' \code{"regression"}.
#' @param layers Integer, how many layers (iterations) to run.
#' @param nAntibodies Integer, how many antibodies (prototypes) to generate
#' initially in each layer.
#' @param minAntibodies Integer, minimal number of antibodies to keep in each
#' layer (to avoid going below a certain threshold).
#' @param epsilon Numeric, threshold parameter passed to bHIVE.
#' @param beta Numeric, hyperparameter controlling selection pressure in bHIVE.
#' @param maxIter Integer, maximum iterations for bHIVE in each layer.
#' @param collapseMethod How to collapse each cluster: \code{"centroid"},
#' \code{"medoid"}, \code{"median"}, \code{"mode"}. Default is \code{"centroid"}.
#' @param minClusterSize If not NULL, minimum size of clusters. Smaller ones can be
#' merged or discarded. Default \code{NULL} = no special handling.
#' @param distance Distance metric used in the \code{dist} function for medoid
#' calculation. Default \code{"euclidean"}.
#' @param verbose Logical, whether to print progress at each layer.
#' @param ... Additional parameters passed to \code{bHIVE}.
#'
#' @return A list of length \code{layers}. Each element (layer) includes:
#' \itemize{
#'   \item \code{antibodies}: The prototype positions in that layer.
#'   \item \code{assignments}: The antibody index assigned to each row of 
#'   that layer's data.
#'   \item \code{membership}: For each \strong{original} row of \code{X}, 
#'   which cluster/antibody it belongs to in this layer.
#'   \item \code{predictions}: For classification/regression, predicted label 
#'   or numeric value
#'         for each \strong{original row in X}.
#'   \item \code{task}: The specified task (clustering, classification, 
#'   regression).
#' }
#'
#' @examples
#' # Example 1: Clustering
#' data(iris)
#' X_iris <- iris[, 1:4]
#' resC <- honeycombHIVE(
#'   X = X_iris,
#'   task = "clustering",
#'   layers = 3,
#'   nAntibodies = 15,
#'   beta = 5,
#'   maxIter = 10
#' )
#' # Final cluster membership of original rows:
#' head(resC[[3]]$membership)
#'
#' # Example 2: Regression
#' set.seed(42)
#' X_reg <- matrix(rnorm(100*4), ncol = 4)
#' y_reg <- rowSums(X_reg[, 1:2]) + rnorm(100)
#' resReg <- honeycombHIVE(
#'   X = X_reg,
#'   y = y_reg,
#'   task = "regression",
#'   layers = 3,
#'   nAntibodies = 10
#' )
#' # Compare final predictions to original y:
#' cor(resReg[[3]]$predictions, y_reg)
#'
#' @importFrom stats median
#' @export
honeycombHIVE <- function(X,
                          y = NULL,
                          task = c("clustering", "classification", "regression"),
                          layers = 3,
                          nAntibodies = 20,
                          minAntibodies = 5,
                          epsilon = 0.05,
                          beta = 5,
                          maxIter = 10,
                          collapseMethod = c("centroid", "medoid", "median", "mode"),
                          minClusterSize = NULL,
                          distance = "euclidean",
                          verbose = TRUE,
                          ...) {
  
  task <- match.arg(task)
  collapseMethod <- match.arg(collapseMethod)
  
  .validate_bHIVE_input(X, y)
  
  X <- as.data.frame(X)
  n_original <- nrow(X)
  
  # Assign rownames if missing
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("row_", seq_len(n_original))
  }
  
  # rowIndices: a list mapping each "prototype" row to the original row indices
  # Initially, each row is its own group
  rowIndices <- lapply(seq_len(n_original), function(i) i)
  names(rowIndices) <- rownames(X)
  
  current_X <- X
  current_y <- if (!is.null(y)) y else NULL
  
  results <- vector("list", length = layers)
  dots <- list(...) # Capture additional arguments to pass to bHIVE
  
  for (layer in seq_len(layers)) {
    if (verbose) {
      cat(sprintf("\n=== Running layer %d / %d (task = %s) ===\n", layer, 
                  layers, task))
    }
    
    # Make sure we don't exceed number of rows
    nAntibodies_layer <- min(nAntibodies, max(minAntibodies, nrow(current_X)))
    
    # Prepare arguments for bHIVE
    bHIVE_args <- list(X = current_X,
                       y = current_y,
                       task = task,
                       nAntibodies = nAntibodies_layer,
                       epsilon = epsilon,
                       beta = beta,
                       maxIter = maxIter,
                       verbose = verbose)
    
    # Add any extra arguments from dots, excluding those already used
    extra_args <- dots[!names(dots) %in% names(bHIVE_args)]
    bHIVE_args <- c(bHIVE_args, extra_args)
    
    # Call bHIVE for this layer
    res_layer <- do.call(bHIVE, bHIVE_args)
    
    # For each original row, we store the cluster membership
    membership <- rep(NA, n_original)  # cluster membership for original data
    cluster_ids <- res_layer$assignments  # cluster assignments for current layer's data
    unique_clusters <- unique(cluster_ids)
    
    # Track how each newly formed prototype maps back to original rows
    subsets_rownames <- split(rownames(current_X), f = cluster_ids)
    
    new_rowIndices <- list()
    cluster_counter <- 1
    
    # Classification/regression, store predictions for all original rows
    # at this layer, so we can compare to y directly.
    if (task %in% c("classification", "regression")) {
      predictions_all_original <- rep(NA, n_original)
    } else {
      predictions_all_original <- NULL
    }
    
    proto_matrix <- NULL  # will accumulate prototypes
    
    # Iterate over each cluster from bHIVE
    for (cid in unique_clusters) {
      rn_in_cluster <- subsets_rownames[[as.character(cid)]]
      # 'orig_indices' = which original row indices belong to this cluster
      orig_indices <- unlist(rowIndices[rn_in_cluster], use.names = FALSE)
      
      # Assign membership
      membership[orig_indices] <- cluster_counter
      
      # Subset for prototype creation
      sub_data <- current_X[rn_in_cluster, , drop = FALSE]
      
      # Create a prototype of that cluster
      if (nrow(sub_data) == 0) {
        proto <- rep(NA, ncol(current_X))
      } else {
        proto <- switch(
          collapseMethod,
          "centroid" = colMeans(sub_data, na.rm = TRUE),
          "medoid" = {
            dmat <- as.matrix(dist(sub_data, method = distance))
            idx_medoid <- which.min(rowSums(dmat))
            as.numeric(sub_data[idx_medoid, , drop = FALSE])
          },
          "median" = apply(sub_data, 2, median, na.rm = TRUE),
          "mode" = {
            apply(sub_data, 2, function(col_data) {
              tabV <- table(col_data)
              as.numeric(names(tabV)[which.max(tabV)])
            })
          }
        )
      }
      
      # For classification/regression, store predictions for *original rows* in this cluster
      if (task == "classification" && !is.null(current_y)) {
        idx_in_current <- match(rn_in_cluster, rownames(current_X))
        if (length(idx_in_current) > 0) {
          y_vals <- current_y[idx_in_current]
          tb <- table(y_vals)
          class_pred <- names(tb)[which.max(tb)]
          # fill the predictions for the original indices
          predictions_all_original[orig_indices] <- class_pred
        }
        
      } else if (task == "regression" && !is.null(current_y)) {
        idx_in_current <- match(rn_in_cluster, rownames(current_X))
        if (length(idx_in_current) > 0) {
          val_pred <- mean(current_y[idx_in_current], na.rm = TRUE)
          # fill the predictions for the original indices
          predictions_all_original[orig_indices] <- val_pred
        }
      }
      
      # Construct new row (prototype) name
      new_row_name <- paste0("Layer", layer, "_Cluster", cluster_counter)
      new_rowIndices[[new_row_name]] <- orig_indices
      
      # Accumulate prototype in proto_matrix
      if (is.null(proto_matrix)) {
        proto_matrix <- matrix(proto, nrow = 1)
        rownames(proto_matrix) <- new_row_name
      } else {
        temp <- matrix(proto, nrow = 1)
        rownames(temp) <- new_row_name
        proto_matrix <- rbind(proto_matrix, temp)
      }
      
      cluster_counter <- cluster_counter + 1
    }
    
    # Optionally discard small clusters
    if (!is.null(minClusterSize)) {
      keep_idx <- sapply(new_rowIndices, length) >= minClusterSize
      if (any(!keep_idx)) {
        if (verbose) {
          cat(sprintf("Discarding %d clusters smaller than minClusterSize (%d)\n", 
                      sum(!keep_idx), minClusterSize))
        }
        proto_matrix <- proto_matrix[keep_idx, , drop = FALSE]
        new_rowIndices <- new_rowIndices[keep_idx]
      }
    }
    
    if (is.null(proto_matrix) || nrow(proto_matrix) == 0) {
      stop("No valid prototypes generated. Please check parameters and input data.")
    }
    
    # Prepare for the next layer
    new_data <- as.data.frame(proto_matrix)
    current_X <- new_data
    rowIndices <- new_rowIndices
    
    # Save the full predictions (for original data) in res_layer
    if (task %in% c("classification", "regression")) {
      res_layer$predictions <- predictions_all_original
    } else {
      res_layer$predictions <- NULL
    }
    
    # Also store membership for each original row
    res_layer$membership <- membership
    
    # Store the updated layer result
    results[[layer]] <- res_layer
    
    if (verbose) {
      cat(sprintf("Layer %d completed. Next layer has %d data points.\n", 
                  layer, nrow(current_X)))
    }
    
    # If no data remain, we can either break or reinit
    if (nrow(current_X) == 0 && layer < layers) {
      warning("No valid clusters remain after layer ", layer, ". Reinitializing with original dataset.\n")
      current_X <- X
      rowIndices <- lapply(seq_len(n_original), function(i) i)
      names(rowIndices) <- rownames(X)
      current_y <- y
      # If classification, ensure factor
      if (task == "classification" && !is.null(current_y)) {
        current_y <- as.factor(current_y)
      }
      
    } else {
      # Rebuild current_y for next layer (classification/regression) 
      # each new prototype has a label/value from the collapsed cluster
      if (task == "classification" && !is.null(current_y)) {
        new_y <- rep(NA, nrow(current_X))
        nrn <- rownames(current_X)
        for (k in seq_along(nrn)) {
          orig_ids <- rowIndices[[nrn[k]]]
          y_vals <- y[orig_ids]
          tb <- table(y_vals)
          new_y[k] <- names(tb)[which.max(tb)]
        }
        current_y <- as.factor(new_y)
        
      } else if (task == "regression" && !is.null(current_y)) {
        new_y <- rep(NA, nrow(current_X))
        nrn <- rownames(current_X)
        for (k in seq_along(nrn)) {
          orig_ids <- rowIndices[[nrn[k]]]
          new_y[k] <- mean(y[orig_ids], na.rm = TRUE)
        }
        current_y <- new_y
      }
    }
  }
  
  return(results)
}
