#' honeycombHIVE: Multilayer AIS with optional gradient-based fine-tuning
#' 
#' The \code{honeycombHIVE} function implements a multilayer artificial immune 
#' system that iteratively refines a set of prototypes - referred to as 
#' antibodies - to model the structure of the input data. In each layer, the 
#' function first uses the \code{\link{bHIVE}} algorithm to generate or update 
#' antibodies based on the current data representation and task (clustering, 
#' classification, or regression). Optionally, it applies gradient-based 
#' fine-tuning (via \code{\link{refineB}}) to these antibodies, allowing for 
#' advanced refinement through various optimizers (e.g., SGD, Adam, RMSProp) and 
#' customizable loss functions. The final output is a hierarchical set of layers 
#' that encapsulate both the refined prototypes and the corresponding cluster 
#' assignments or predictions for the original observations, making 
#' \code{honeycombHIVE} a versatile tool for adaptive learning and pattern 
#' recognition.
#' 
#' @param X A numeric matrix or data frame of input features 
#' (rows = observations, columns = features).
#' @param y Optional target vector (factor for classification, numeric 
#' for regression).
#' @param task Character, one of "clustering", "classification", or "regression".
#' @param layers Integer, how many layers (AIS iterations) to run.
#' @param nAntibodies Integer, how many antibodies (prototypes) to generate 
#' initially in each layer.
#' @param minAntibodies Integer, minimal number of antibodies to keep in each 
#' layer.
#' @param epsilon Numeric, threshold param for \code{bHIVE} suppression.
#' @param beta Numeric, selection pressure param for \code{bHIVE}.
#' @param maxIter Integer, maximum iterations for \code{bHIVE} each layer.
#' @param collapseMethod One of "centroid","medoid","median","mode".
#' @param minClusterSize Minimum cluster size. Smaller clusters can be 
#' merged/discarded if not NULL.
#' @param distance Distance metric for medoid calculation, e.g. "euclidean".
#' @param verbose Logical, if TRUE prints progress at each layer.
#' @param refine Logical, if TRUE apply gradient-based refinement via 
#' \code{refineB()} to each layer's prototypes.
#' @param refineLoss Character specifying the loss for \code{refineB()} 
#' (e.g. "mse", "mae", etc.).
#' @param refineSteps Integer, number of gradient steps in \code{refineB()}.
#' @param refineLR Numeric, learning rate for gradient updates.
#' @param refinePushAway Logical, if TRUE and classification, push prototypes 
#' away from differently labeled points.
#' @param refineHuberDelta Numeric, delta parameter if using the "huber" loss.
#' @param refineOptimizer Character, one of \code{"sgd", "momentum", "adagrad", 
#' "adam", "rmsprop"} to be passed to \code{refineB()}.
#' @param refineMomentumCoef Numeric, momentum coefficient (if using momentum).
#' @param refineBeta1 Numeric, first moment decay rate (if using Adam).
#' @param refineBeta2 Numeric, second moment decay rate (if using Adam).
#' @param refineRmspropDecay Numeric, decay rate for the moving average of 
#' squared gradients (if using RMSProp).
#' @param refineEpsilon Numeric, a small constant for numerical stability 
#' (used in adaptive optimizers).
#' @param ... Additional arguments passed to \code{bHIVE}.
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
#' 
#' @return A list of length \code{layers}. Each element (layer) includes:
#'   \itemize{
#'     \item \code{antibodies}: The prototypes in that layer.
#'     \item \code{assignments}: Antibody index (in that layer) for each row 
#'     of \code{current_X}.
#'     \item \code{membership}: For each \strong{original} row in \code{X}, 
#'     which cluster/antibody it belongs to in this layer.
#'     \item \code{predictions}: If classification/regression, predicted label 
#'     or numeric value for each original row in \code{X}.
#'     \item \code{task}: The specified task.
#'   }
#'
#' @importFrom stats median setNames
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
                          refine = FALSE,
                          refineLoss = "mse",
                          refineSteps = 5,
                          refineLR = 0.01,
                          refinePushAway = TRUE,
                          refineHuberDelta = 1.0,
                          refineOptimizer = "sgd",
                          refineMomentumCoef = 0.9,
                          refineBeta1 = 0.9,
                          refineBeta2 = 0.999,
                          refineRmspropDecay = 0.9,
                          refineEpsilon = 1e-8,
                          ...) {
  
  # =======================
  # 0) Basic Checks & Setup
  # =======================
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
  
  # Force X into a data.frame for subsetting operations
  X <- as.data.frame(X)
  n_original <- nrow(X)
  
  # Assign rownames if missing
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("row_", seq_len(n_original))
  }
  
  # Create a mapping from original row names to indices.
  rowIndices <- setNames(as.list(seq_len(n_original)), rownames(X))
  
  # If y is provided, assign names so that we can index by row names.
  if (!is.null(y)) {
    if (is.null(names(y))) {
      names(y) <- rownames(X)
    }
  }
  
  # Current data to pass to bHIVE
  current_X <- X
  current_y <- if (!is.null(y)) y else NULL
  
  # Prepare a list to store results from each layer
  results <- vector("list", length = layers)
  dots <- list(...)  # additional args for bHIVE
  
  # ========================
  # 1) Iteration Over Layers
  # ========================
  for (layer_i in seq_len(layers)) {
    if (verbose) {
      message(sprintf("\n=== honeycombHIVE: Layer %d / %d (task=%s) ===",
                      layer_i, layers, task))
    }
    
    # 1a) Decide how many antibodies for this layer
    n_current <- nrow(current_X)  # Number of points (or prototypes) in this layer
    nAntibodies_layer <- min(nAntibodies, max(minAntibodies, n_current))
    
    # 1b) Construct arguments for bHIVE
    bHIVE_args <- list(
      X           = current_X,
      y           = current_y,
      task        = task,
      nAntibodies = nAntibodies_layer,
      epsilon     = epsilon,
      beta        = beta,
      maxIter     = maxIter,
      verbose     = verbose
    )
    # Merge additional arguments not already specified
    extra_args <- dots[setdiff(names(dots), names(bHIVE_args))]
    bHIVE_args <- c(bHIVE_args, extra_args)
    
    # 1c) Run bHIVE on the current layer's data
    res_layer <- do.call(bHIVE, bHIVE_args)
    
    # 1d) Optional refineB step
    if (refine) {
      assignments_layer <- res_layer$assignments
      
      if (task == "clustering") {
        dummy_y <- rep(1, n_current)
        refined_A <- refineB(
          A           = res_layer$antibodies,
          X           = current_X,
          y           = dummy_y,
          assignments = assignments_layer,
          loss        = refineLoss,
          task        = "classification",  
          steps       = refineSteps,
          lr          = refineLR,
          push_away   = FALSE,
          huber_delta = refineHuberDelta,
          verbose     = verbose,
          optimizer   = refineOptimizer,
          momentum_coef = refineMomentumCoef,
          beta1       = refineBeta1,
          beta2       = refineBeta2,
          rmsprop_decay = refineRmspropDecay,
          epsilon     = refineEpsilon
        )
      } else {
        refined_A <- refineB(
          A           = res_layer$antibodies,
          X           = current_X,
          y           = current_y,
          assignments = assignments_layer,
          loss        = refineLoss,
          task        = task,
          steps       = refineSteps,
          lr          = refineLR,
          push_away   = (task == "classification" && refinePushAway),
          huber_delta = refineHuberDelta,
          verbose     = verbose,
          optimizer   = refineOptimizer,
          momentum_coef = refineMomentumCoef,
          beta1       = refineBeta1,
          beta2       = refineBeta2,
          rmsprop_decay = refineRmspropDecay,
          epsilon     = refineEpsilon
        )
      }
      res_layer$antibodies <- refined_A
    }
    
    # ================
    # 2) Collapse Step
    # ================
    cluster_ids <- res_layer$assignments  # Cluster indices from bHIVE
    unique_clusters <- unique(cluster_ids)
    
    # Initialize membership for each original row.
    membership <- rep(NA_integer_, n_original)
    
    # For classification/regression, prepare a vector for predictions.
    predictions_all_original <- if (task %in% c("classification", "regression")) rep(NA, n_original) else NULL
    
    # Pre-allocate lists for new prototypes and updated row indices.
    proto_list <- vector("list", length(unique_clusters))
    new_rowIndices <- vector("list", length(unique_clusters))
    
    cluster_counter <- 1
    subsets_rownames <- split(rownames(current_X), f = cluster_ids)
    
    for (cid in unique_clusters) {
      rn_in_cluster <- subsets_rownames[[as.character(cid)]]
      # Map these rownames back to original indices.
      orig_indices <- unlist(rowIndices[rn_in_cluster], use.names = FALSE)
      n_members <- length(orig_indices)
      
      # Record membership for these original rows.
      membership[orig_indices] <- cluster_counter
      
      # Extract the sub-data for this cluster.
      sub_data <- current_X[rn_in_cluster, , drop = FALSE]
      
      # If the cluster is too small, build an NA prototype.
      if (!is.null(minClusterSize) && n_members < minClusterSize) {
        proto <- rep(NA_real_, ncol(current_X))
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
            apply(sub_data, 2, function(vec) {
              tb <- table(vec)
              as.numeric(names(tb)[which.max(tb)])
            })
          }
        )
      }
      
      # For regression, compute the prediction from current_y.
      if (task == "regression" && !is.null(current_y)) {
        # Use match() to obtain numeric indices from current_X's rownames.
        indices_in_current <- match(rn_in_cluster, rownames(current_X))
        val_pred <- mean(current_y[indices_in_current], na.rm = TRUE)
        if (is.na(val_pred)) {
          val_pred <- mean(y, na.rm = TRUE)  # fallback using original y
        }
        predictions_all_original[orig_indices] <- val_pred
      } else if (task == "classification" && !is.null(current_y)) {
        y_vals <- current_y[rn_in_cluster]
        tb <- table(y_vals)
        class_pred <- names(tb)[which.max(tb)]
        predictions_all_original[orig_indices] <- class_pred
      }
      
      new_rowIndices[[cluster_counter]] <- orig_indices
      proto_list[[cluster_counter]] <- proto
      cluster_counter <- cluster_counter + 1
    }
    
    # Build the prototype matrix.
    proto_matrix <- do.call(rbind, proto_list)
    
    # Remove any prototypes that are entirely NA.
    keep_mask <- !apply(proto_matrix, 1, function(x) all(is.na(x)))
    if (any(!keep_mask)) {
      if (verbose) {
        message(sprintf("Discarding %d collapsed prototypes that are empty/NA.", 
                        sum(!keep_mask)))
      }
      proto_matrix   <- proto_matrix[keep_mask, , drop = FALSE]
      new_rowIndices <- new_rowIndices[keep_mask]
    }
    
    # Name the rows of the prototype matrix.
    row_names <- paste0("Layer", layer_i, "_Cluster", seq_len(nrow(proto_matrix)))
    rownames(proto_matrix) <- row_names
    names(new_rowIndices) <- row_names
    
    if (nrow(proto_matrix) == 0) {
      stop("No valid prototypes generated. Check parameters & input data.")
    }
    
    # Build new data for the next layer.
    new_data <- as.data.frame(proto_matrix)
    
    # Save predictions and membership for this layer.
    if (task %in% c("classification", "regression")) {
      res_layer$predictions <- predictions_all_original
    } else {
      res_layer$predictions <- NULL
    }
    res_layer$membership <- membership
    results[[layer_i]] <- res_layer
    
    if (verbose) {
      message(sprintf("Layer %d completed. Next layer will use %d prototypes.",
                      layer_i, nrow(new_data)))
    }
    
    # =========================
    # 3) Prepare for Next Layer
    # =========================
    current_X <- new_data
    rowIndices <- new_rowIndices
    
    # If no prototypes remain and there are remaining layers, reinitialize.
    if (nrow(current_X) == 0 && layer_i < layers) {
      warning("No valid clusters remain after layer ", layer_i,
              ". Reinitializing with original dataset for next layer.\n")
      current_X <- X
      rowIndices <- setNames(as.list(seq_len(n_original)), rownames(X))
      current_y <- y
      if (task == "classification" && !is.null(current_y)) {
        current_y <- as.factor(current_y)
      }
      next
    }
    
    # Rebuild current_y for next layer based on collapsed clusters.
    if (task == "classification" && !is.null(current_y)) {
      n_new <- nrow(current_X)
      new_y <- character(n_new)
      rn_protos <- rownames(current_X)
      for (ii in seq_len(n_new)) {
        orig_ids <- rowIndices[[rn_protos[ii]]]
        tbl_ii <- table(y[orig_ids])
        new_y[ii] <- names(tbl_ii)[which.max(tbl_ii)]
      }
      current_y <- factor(new_y, levels = levels(y))
    } else if (task == "regression" && !is.null(current_y)) {
      n_new <- nrow(current_X)
      new_y <- numeric(n_new)
      rn_protos <- rownames(current_X)
      for (ii in seq_len(n_new)) {
        orig_ids <- rowIndices[[rn_protos[ii]]]
        new_y[ii] <- mean(y[orig_ids], na.rm = TRUE)
      }
      current_y <- new_y
    }
  }
  return(results)
}
