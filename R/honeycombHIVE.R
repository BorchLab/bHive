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
#' @param ... Additional parameters passed to the placeholder \code{bHIVE} function.
#'
#' @return A list of length \code{layers}. Each element (layer) includes:
#' \itemize{
#'   \item \code{antibodies}: The prototype positions in that layer.
#'   \item \code{assignments}: The antibody index assigned to each row of that layer's data.
#'   \item \code{membership}: For each \strong{original} row of \code{X}, which cluster/antibody it belongs to in this layer.
#'   \item \code{predictions}: For classification/regression, predicted label or numeric value
#'         for each \strong{row in this layer's data}.
#'   \item \code{task}: The specified task (clustering, classification, regression).
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
#' # Check membership + final predictions:
#' head(resReg[[3]]$membership)
#' head(resReg[[3]]$predictions)
#'
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
  
  # Convert to data.frame for uniform handling
  X <- as.data.frame(X)
  n_original <- nrow(X)
  
  # Ensure row names
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("row_", seq_len(n_original))
  }
  
  # rowIndices[[row_name]] = integer vector of the original row indices that formed this row
  rowIndices <- lapply(seq_len(n_original), function(i) i)
  names(rowIndices) <- rownames(X)
  
  current_X <- X
  if (!is.null(y)) {
    current_y <- y
  } else {
    current_y <- NULL
  }
  
  results <- vector("list", length = layers)
  
  for (layer in seq_len(layers)) {
    if (verbose) {
      cat(sprintf("\n=== Running layer %d / %d (task = %s) ===\n", layer, layers, task))
    }
    
    # Adjust number of antibodies
    nAntibodies_layer <- min(nAntibodies, max(minAntibodies, nrow(current_X)))
    
    # Run bHIVE on the *current* data
    res_layer <- bHIVE(X           = current_X,
                       y           = current_y,
                       task        = task,
                       nAntibodies = nAntibodies_layer,
                       epsilon     = epsilon,
                       beta        = beta,
                       maxIter     = maxIter,
                       verbose     = verbose,
                       ...)
    
    # We'll create a membership vector for all original rows
    membership <- rep(NA, n_original)
    
    # "assignments" is the cluster index for each row in current_X
    cluster_ids <- res_layer$assignments
    unique_clusters <- unique(cluster_ids)
    
    # Group row names of current_X by cluster assignment
    subsets_rownames <- split(rownames(current_X), f = cluster_ids)
    
    # Prepare to build new prototypes (antibodies) and new rowIndices after collapse
    new_rowIndices <- list()
    cluster_counter <- 1
    
    # For classification/regression, we'll store predictions for the *current layer data*
    # i.e., one label/value per row in current_X
    predictions_this_layer <- rep(NA, nrow(current_X))
    
    for (cid in unique_clusters) {
      # The row names in current_X belonging to cluster cid
      rn_in_cluster <- subsets_rownames[[as.character(cid)]]
      
      # Combine the original row indices from each
      orig_indices <- unlist(rowIndices[rn_in_cluster], use.names = FALSE)
      
      # Mark membership for those original rows
      membership[orig_indices] <- cluster_counter
      
      # Collapsing method
      sub_data <- current_X[rn_in_cluster, , drop = FALSE]
      if (nrow(sub_data) == 0) {
        # fallback
        proto <- colMeans(current_X, na.rm = TRUE)
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
      
      # If classification or regression, define the label/target for this prototype
      if (task == "classification") {
        # majority vote among sub_data's indices in current_y
        if (!is.null(current_y)) {
          idx_in_current <- match(rn_in_cluster, rownames(current_X))
          if (length(idx_in_current) > 0) {
            y_vals <- current_y[idx_in_current]
            tb <- table(y_vals)
            class_pred <- names(tb)[which.max(tb)]
            predictions_this_layer[idx_in_current] <- class_pred
          }
        }
      } else if (task == "regression") {
        # mean target among sub_data's indices
        if (!is.null(current_y)) {
          idx_in_current <- match(rn_in_cluster, rownames(current_X))
          if (length(idx_in_current) > 0) {
            val_pred <- mean(current_y[idx_in_current], na.rm = TRUE)
            predictions_this_layer[idx_in_current] <- val_pred
          }
        }
      }
      
      # Build a new row name for this cluster
      new_row_name <- paste0("Layer", layer, "_Cluster", cluster_counter)
      new_rowIndices[[new_row_name]] <- orig_indices
      
      if (cluster_counter == 1) {
        proto_matrix <- matrix(proto, nrow=1)
        rownames(proto_matrix) <- new_row_name
      } else {
        temp <- matrix(proto, nrow=1)
        rownames(temp) <- new_row_name
        proto_matrix <- rbind(proto_matrix, temp)
      }
      
      cluster_counter <- cluster_counter + 1
    }
    
    # Possibly discard small clusters
    if (!is.null(minClusterSize)) {
      # Check how many original points are in each new cluster
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
    
    # Build a data.frame from proto_matrix for the next layer
    new_data <- as.data.frame(proto_matrix)
    
    # Update current_X to these new prototypes
    current_X <- new_data
    # Update rowIndices
    rowIndices <- new_rowIndices
    
    # For classification/regression, store predictions for the *current layer's data*
    if (task %in% c("classification","regression")) {
      res_layer$predictions <- predictions_this_layer
    }
    
    # Store membership (length = n_original) with the cluster ID each original row belongs to
    res_layer$membership <- membership
    
    # Save the updated layer result
    results[[layer]] <- res_layer
    
    # Reporting
    if (verbose) {
      cat(sprintf("Layer %d completed. Next layer has %d data points.\n",
                  layer, nrow(current_X)))
      cat(sprintf("Layer %d used %d antibodies.\n",
                  layer, nAntibodies_layer))
    }
    
    # If all data was discarded, reinitialize if more layers remain
    if (nrow(current_X) == 0 && layer < layers) {
      warning("No valid clusters remain after layer ", layer,
              ". Reinitializing with original dataset.\n")
      current_X <- X
      rowIndices <- lapply(seq_len(n_original), function(i) i)
      names(rowIndices) <- rownames(X)
      if (!is.null(y)) current_y <- y
    } else {
      # Update the target for the next layer if classification/regression
      # We'll assign the predicted label/value of each cluster's prototype
      # as the new_y for the next iteration:
      if (task == "classification" && !is.null(current_y)) {
        # For each new row in current_X, compute the single label
        # from the cluster that formed it:
        new_y <- rep(NA, nrow(current_X))
        nrn <- rownames(current_X)
        for (k in seq_along(nrn)) {
          orig_ids <- rowIndices[[nrn[k]]]
          y_vals <- y[orig_ids]
          tb <- table(y_vals)
          new_y[k] <- names(tb)[which.max(tb)]
        }
        current_y <- new_y
      }
      else if (task == "regression" && !is.null(current_y)) {
        # For each new row in current_X, define the numeric label as the mean among old points
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
