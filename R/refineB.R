#' refineB: Gradient-based fine-tuning for bHIVE antibodies with multiple loss functions
#'
#' After running bHIVE (or honeycombHIVE), you have a set of final antibody positions (A)
#' in feature space. This function refines those prototypes by iterating over the
#' data points assigned to each antibody and applying small gradient-based updates
#' according to a user-chosen loss function.
#' 
#' The user must provide:
#' - A numeric matrix \code{A} of antibody/prototype positions (nAb x nFeatures).
#' - A numeric matrix/data frame \code{X} of data (nSamples x nFeatures).
#' - Optional \code{y} for classification/regression. If \code{task="clustering"},
#'   \code{y} can be NULL or ignored.
#' - An integer vector \code{assignments} (length = nSamples) giving the
#'   antibody index to which each data point is assigned (range 1..nAb).
#'   
#' ## Available Losses
#'
#' ### Classification (factor y)
#' - **"categorical_crossentropy"**: Pull prototypes toward data points if they share
#'   the antibody's dominant label; push away otherwise.
#' - **"binary_crossentropy"**: Similar to categorical CE, but we interpret factor y
#'   as binary (two classes). Pull for same label, push for different label.
#' - **"kullback_leibler"**: Very rough approach that treats “dominant label vs. others”
#'   as p and q distributions, pushing/pulling prototypes.
#' - **"cosine"**: Interpreted as trying to maximize cosine similarity for same-label points
#'   and minimize for different-label points.
#'
#' ### Regression (numeric y)
#' - **"mse"**: Mean squared error approximation in feature space (pull prototypes
#'   toward assigned points).
#' - **"mae"**: Mean absolute error approach (sign-based pull).
#' - **"poisson"**: Poisson deviance (toy approach that scales the gradient by
#'   (pred - y)/pred if we stored a predicted rate; here we do a naive version).
#' - **"huber"**: Combines L1 and L2 regions, uses a delta cutoff. We adapt it
#'   to a naive per-point gradient in feature space.
#'
#' @param A Numeric matrix (nAb x d) of antibody/prototype positions.
#' @param X Matrix or data frame (nSamples x d) of feature data.
#' @param y Optional. Factor (classification), numeric (regression), or NULL (clustering).
#' @param assignments Integer vector (length=nSamples), specifying the antibody index
#'                    each sample belongs to (1..nAb).
#' @param task One of c("clustering","classification","regression").
#' @param loss One of c("categorical_crossentropy","binary_crossentropy",
#'   "kullback_leibler","cosine","mse","mae","poisson","huber").
#' @param steps Integer. How many gradient steps to run.
#' @param lr Numeric. Learning rate for each update.
#' @param push_away Logical (for classification). Whether to push prototypes away
#'   from differently-labeled samples.
#' @param huber_delta Numeric. The delta threshold if using huber loss.
#' @return Updated matrix A of shape (nAb x d).
#'
#' @details
#' This function does a simple mini-SGD loop over each data point in random order.
#' For classification, it identifies whether the data point shares the antibody's
#' "dominant label" and pulls/pushes accordingly. For regression/clustering, it
#' primarily pulls prototypes toward assigned points. 
#'
#' @export
refineB <- function(A,
                    X,
                    y = NULL,
                    assignments,
                    task = c("clustering", "classification", "regression"),
                    loss = c("categorical_crossentropy", "binary_crossentropy",
                             "kullback_leibler", "cosine",
                             "mse", "mae", "poisson", "huber"),
                    steps = 5,
                    lr = 0.01,
                    push_away = TRUE,
                    huber_delta = 1.0) {
  
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  # Ensure numeric matrix for X and A
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(X)) {
    stop("'X' must be numeric.")
  }
  if (!is.matrix(A) || !is.numeric(A)) {
    A <- as.matrix(A)
  }
  
  nAb <- nrow(A)
  d   <- ncol(A)
  n   <- nrow(X)
  
  if (ncol(X) != d) {
    stop("Number of columns in 'X' must match the number of columns in 'A'.")
  }
  if (length(assignments) != n) {
    stop("Length of 'assignments' must match nrow(X).")
  }
  # Convert assignments to numeric indices
  unique_assignments <- unique(assignments)
  assignments <- match(assignments, unique_assignments)
  nLabels <- length(unique_assignments)
  
  if (nLabels > nAb) {
    if(verbose) {
      message(sprintf("Number of unique assignments (%d) exceeds the number of antibodies (nAb = %d).",
                      nLabels, nAb))
    }
  }
  
  # Pad prototypes if nAb > nLabels
  if (nLabels < nAb) {
    if(verbose) {
      message(sprintf("Number of prototypes (nAb = %d) exceeds unique labels (%d). Extra prototypes will not be assigned data points.\n",
                      nAb, nLabels))
    }
  }
  
  # Validate assignments
  if (any(assignments < 1 | assignments > nAb, na.rm = TRUE)) {
    invalid_assignments <- assignments[assignments < 1 | assignments > nAb]
    stop(sprintf("'assignments' contains invalid values: %s. Valid range is 1..%d.",
                 paste(unique(invalid_assignments), collapse = ", "), nAb))
  }
  
  
  if (task == "regression" && !is.numeric(y)) {
    stop("y must be numeric for regression.")
  }
  
  # Build a samples_by_ab list if needed (e.g., for classification to find dominant labels)
  samples_by_ab <- vector("list", nAb)
  for (i in seq_len(n)) {
    j <- assignments[i]
    samples_by_ab[[j]] <- c(samples_by_ab[[j]], i)
  }
  
  # Classification: find dominant label per antibody
  ab_label <- rep(NA_character_, nAb)
  if (task == "classification") {
    if (!is.factor(y)) {
      y <- as.factor(y)
    }
    for (ab_idx in seq_len(nAb)) {
      pts <- samples_by_ab[[ab_idx]]
      if (length(pts) == 0) {
        ab_label[ab_idx] <- NA
      } else {
        tbl <- table(y[pts])
        ab_label[ab_idx] <- names(tbl)[which.max(tbl)]
      }
    }
  }
  
  # Mini-SGD loop
  for (st in seq_len(steps)) {
    # Shuffle data
    idx_shuf <- sample.int(n)
    for (i in idx_shuf) {
      j <- assignments[i]
      
      # Defensive checks
      if (is.na(j) || j < 1 || j > nAb) {
        next
      }
      
      # x_i must be a numeric vector of length d
      x_i <- as.numeric(X[i, ])
      if (length(x_i) != d) {
        stop(sprintf("x_i is not the correct length: expected %d, got %d", d, length(x_i)))
      }
      
      # same_label logic for classification
      same_label <- NA
      if (task == "classification" && !is.na(ab_label[j])) {
        same_label <- (y[i] == ab_label[j])
      }
      
      # get the gradient from .update_prototype
      grad <- .update_prototype(
        ab_vec      = A[j, ],
        x_i         = x_i,
        same_label  = same_label,
        task        = task,
        loss        = loss,
        push_away   = push_away,
        huber_delta = huber_delta
      )
      
      # sanity check: grad must be length d
      if (length(grad) != d) {
        stop(sprintf(".update_prototype returned a gradient of length %d, expected %d",
                     length(grad), d))
      }
      
      # update
      A[j, ] <- A[j, ] + lr * grad
    }
  }
  
  return(A)
}