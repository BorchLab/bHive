#' honeycombHIVE: Multilayer AIS with optional gradient-based fine-tuning
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
#' @param refine Logical. If TRUE, apply gradient-based refinement to each 
#' layer's antibody positions via \code{refineB()}.
#' @param refineLoss Character specifying the loss function if \code{refine=TRUE}.
#'   e.g. "mse","mae","categorical_crossentropy","binary_crossentropy","huber",
#'   "kullback_leibler", etc.
#' @param refineSteps Integer. Number of gradient steps in 
#' \code{refineB()}.
#' @param refineLR Numeric. Learning rate for gradient updates.
#' @param refinePushAway Logical. If TRUE (classification only), push 
#' prototypes away
#'   from differently labeled points.
#' @param refineHuberDelta Numeric. Delta parameter if using the "huber" loss.
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
#'   or numeric value for each \strong{original row in X}.
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
                          collapseMethod = c("centroid","medoid","median","mode"),
                          minClusterSize = NULL,
                          distance = "euclidean",
                          verbose = TRUE,
                          refine = FALSE,
                          refineLoss = "mse",
                          refineSteps = 5,
                          refineLR = 0.01,
                          refinePushAway = TRUE,
                          refineHuberDelta = 1.0,
                          ...) {
  
  task <- match.arg(task)
  collapseMethod <- match.arg(collapseMethod)
  
  .validate_bHIVE_input(X, y)
  
  # Validate refineLoss based on task
  valid_losses <- list(
    regression = c("mse", "mae", "huber", "kullback_leibler"),
    classification = c("categorical_crossentropy", "binary_crossentropy", "kullback_leibler"),
    clustering = c("mse", "mae", "huber")
  )
  
  if (refine && !refineLoss %in% valid_losses[[task]]) {
    stop(sprintf("Invalid refineLoss '%s' for task '%s'. Supported losses: %s", 
                 refineLoss, task, paste(valid_losses[[task]], collapse = ", ")))
  }
  
  X <- as.data.frame(X)
  n_original <- nrow(X)
  
  # Assign rownames if missing
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("row_", seq_len(n_original))
  }
  
  # rowIndices: mapping each prototype row -> original row indices
  # initially, each row is its own group
  rowIndices <- lapply(seq_len(n_original), function(i) i)
  names(rowIndices) <- rownames(X)
  
  current_X <- X
  current_y <- if (!is.null(y)) y else NULL
  
  results <- vector("list", length=layers)
  dots <- list(...)  # additional args for bHIVE
  
  for (layer in seq_len(layers)) {
    if (verbose) {
      cat(sprintf("\n=== honeycombHIVE_plus: Layer %d / %d (task = %s) ===\n",
                  layer, layers, task))
    }
    
    # 1) bHIVE on current data
    nAntibodies_layer <- min(nAntibodies, max(minAntibodies, nrow(current_X)))
    bHIVE_args <- list(
      X = current_X,
      y = current_y,
      task = task,
      nAntibodies = nAntibodies_layer,
      epsilon = epsilon,
      beta = beta,
      maxIter = maxIter,
      verbose = verbose
    )
    extra_args <- dots[!names(dots) %in% names(bHIVE_args)]
    bHIVE_args <- c(bHIVE_args, extra_args)
    
    res_layer <- do.call(bHIVE, bHIVE_args)
    
    # 2) refineB step
    if (refine) {
      
      # 'assignments' here are for the layer's data, not the original data.
      assignments_layer <- res_layer$assignments
      if (task == "clustering") {
        
        fake_y <- rep(1, nrow(current_X))
        new_A <- refineB(
          A = res_layer$antibodies,
          X = current_X,
          y = fake_y,  # dummy
          assignments = assignments_layer,
          loss = refineLoss,
          task = "classification",  # treat them as same label
          steps = refineSteps,
          lr = refineLR,
          push_away = FALSE,
          huber_delta = refineHuberDelta
        )
      } else {
        # classification or regression
        new_A <- refineB(
          A = res_layer$antibodies,
          X = current_X,
          y = current_y,
          assignments = assignments_layer,
          loss = refineLoss,
          task = task,
          steps = refineSteps,
          lr = refineLR,
          push_away = refinePushAway,
          huber_delta = refineHuberDelta
        )
      }
      
      # update bHIVE result with refined prototypes
      res_layer$antibodies <- new_A
    }
    
    # 3) We proceed to "collapse" each cluster into one prototype for the next layer
    membership <- rep(NA, n_original)
    cluster_ids <- res_layer$assignments
    unique_clusters <- unique(cluster_ids)
    subsets_rownames <- split(rownames(current_X), f = cluster_ids)
    
    new_rowIndices <- list()
    cluster_counter <- 1
    
    # If classification/regression, we store predictions for all original rows
    if (task %in% c("classification","regression")) {
      predictions_all_original <- rep(NA, n_original)
    } else {
      predictions_all_original <- NULL
    }
    
    proto_matrix <- NULL
    
    # Iterate each cluster
    for (cid in unique_clusters) {
      rn_in_cluster <- subsets_rownames[[as.character(cid)]]
      orig_indices <- unlist(rowIndices[rn_in_cluster], use.names=FALSE)
      membership[orig_indices] <- cluster_counter
      
      sub_data <- current_X[rn_in_cluster, , drop=FALSE]
      
      if (nrow(sub_data) == 0) {
        proto <- rep(NA, ncol(current_X))
      } else {
        proto <- switch(
          collapseMethod,
          "centroid" = colMeans(sub_data, na.rm=TRUE),
          "medoid" = {
            dmat <- as.matrix(dist(sub_data, method=distance))
            idx_medoid <- which.min(rowSums(dmat))
            as.numeric(sub_data[idx_medoid, , drop=FALSE])
          },
          "median" = apply(sub_data, 2, median, na.rm=TRUE),
          "mode" = {
            apply(sub_data, 2, function(col_data) {
              tb <- table(col_data)
              as.numeric(names(tb)[which.max(tb)])
            })
          }
        )
      }
      
      if (task=="classification" && !is.null(current_y)) {
        idx_in_current <- match(rn_in_cluster, rownames(current_X))
        if (length(idx_in_current) > 0) {
          y_vals <- current_y[idx_in_current]
          tb <- table(y_vals)
          class_pred <- names(tb)[which.max(tb)]
          predictions_all_original[orig_indices] <- class_pred
        }
      } else if (task=="regression" && !is.null(current_y)) {
        idx_in_current <- match(rn_in_cluster, rownames(current_X))
        if (length(idx_in_current) > 0) {
          val_pred <- mean(current_y[idx_in_current], na.rm=TRUE)
          predictions_all_original[orig_indices] <- ifelse(is.na(val_pred), mean(y, na.rm=TRUE), val_pred)
        }
      }
      
      new_row_name <- paste0("Layer", layer, "_Cluster", cluster_counter)
      new_rowIndices[[new_row_name]] <- orig_indices
      
      if (is.null(proto_matrix)) {
        proto_matrix <- matrix(proto, nrow=1)
        rownames(proto_matrix) <- new_row_name
      } else {
        tmp <- matrix(proto, nrow=1)
        rownames(tmp) <- new_row_name
        proto_matrix <- rbind(proto_matrix, tmp)
      }
      
      cluster_counter <- cluster_counter+1
    }
    
    # Optionally discard small clusters
    if (!is.null(minClusterSize)) {
      keep_idx <- sapply(new_rowIndices, length) >= minClusterSize
      if (any(!keep_idx)) {
        if (verbose) {
          cat(sprintf("Discarding %d clusters smaller than minClusterSize (%d)\n", 
                      sum(!keep_idx), minClusterSize))
        }
        proto_matrix <- proto_matrix[keep_idx, , drop=FALSE]
        new_rowIndices <- new_rowIndices[keep_idx]
      }
    }
    
    if (is.null(proto_matrix) || nrow(proto_matrix)==0) {
      stop("No valid prototypes generated. Check parameters & input data.")
    }
    
    # For the next layer
    new_data <- as.data.frame(proto_matrix)
    current_X <- new_data
    rowIndices <- new_rowIndices
    
    if (task %in% c("classification","regression")) {
      res_layer$predictions <- predictions_all_original
    } else {
      res_layer$predictions <- NULL
    }
    
    res_layer$membership <- membership
    results[[layer]] <- res_layer
    
    if (verbose) {
      cat(sprintf("Layer %d completed. Next layer has %d data points.\n", layer, nrow(current_X)))
    }
    
    # If no data remain, either break or reinit
    if (nrow(current_X) == 0 && layer < layers) {
      warning("No valid clusters remain after layer ", layer,
              ". Reinitializing with original dataset.\n")
      current_X <- X
      rowIndices <- lapply(seq_len(n_original), function(i) i)
      names(rowIndices) <- rownames(X)
      current_y <- y
      if (task=="classification" && !is.null(current_y)) {
        current_y <- as.factor(current_y)
      }
    } else {
      # rebuild current_y for next layer
      if (task=="classification" && !is.null(current_y)) {
        new_y <- rep(NA, nrow(current_X))
        nrn <- rownames(current_X)
        for (k in seq_along(nrn)) {
          orig_ids <- rowIndices[[nrn[k]]]
          tb <- table(y[orig_ids])
          new_y[k] <- names(tb)[which.max(tb)]
        }
        current_y <- as.factor(new_y)
        
      } else if (task=="regression" && !is.null(current_y)) {
        new_y <- rep(NA, nrow(current_X))
        nrn <- rownames(current_X)
        for (k in seq_along(nrn)) {
          orig_ids <- rowIndices[[nrn[k]]]
          new_y[k] <- mean(y[orig_ids], na.rm=TRUE)
        }
        current_y <- new_y
      }
    }
  }
  
  return(results)
}
