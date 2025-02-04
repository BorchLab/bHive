#' refineB: Gradient-based fine-tuning for bHIVE antibodies with multiple loss 
#' functions and optimizers
#'
#' After running \code{bHIVE} (or within \code{honeycombHIVE}), you have a set of 
#' final antibody positions (A) in feature space. This function refines those
#' prototypes by iterating over the data points assigned to each antibody and 
#' applying gradient-based updates using a user-chosen loss function and optimizer.
#' 
#' The user must provide:
#' - A numeric matrix \code{A} of antibody/prototype positions (nAb x nFeatures).
#' - A numeric matrix/data frame \code{X} of data (nSamples x nFeatures).
#' - Optional \code{y} for classification/regression. If \code{task="clustering"},
#'   \code{y} can be NULL or ignored.
#' - An integer or character vector \code{assignments} (length=nSamples) giving the
#'   antibody index (or label) to which each data point is assigned.
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
#' ## Available Optimizers
#' - **"sgd"**: Basic stochastic gradient descent.
#' - **"momentum"**: SGD with momentum.
#' - **"adagrad"**: Adaptive gradient descent.
#' - **"adam"**: Adaptive moment estimation.
#' - **"rmsprop"**: Root Mean Square Propagation.
#'
#' @param A Numeric matrix (nAb x d) of antibody/prototype positions.
#' @param X Matrix or data frame (nSamples x d) of feature data.
#' @param y Optional. Factor (classification), numeric (regression), or NULL (clustering).
#' @param assignments Integer or character vector (length = nSamples), specifying the antibody index
#'                    (or label) each sample belongs to.
#' @param task One of \code{c("clustering", "classification", "regression")}.
#' @param loss One of \code{c("categorical_crossentropy", "binary_crossentropy",
#'   "kullback_leibler", "cosine", "mse", "mae", "poisson", "huber")}.
#' @param steps Integer. How many gradient steps to run.
#' @param lr Numeric. Learning rate for each update.
#' @param push_away Logical (for classification). Whether to push prototypes away
#'   from differently-labeled samples.
#' @param huber_delta Numeric. The delta threshold if using huber loss.
#' @param verbose Logical. If \code{TRUE}, prints progress messages each iteration.
#' @param optimizer One of \code{c("sgd", "momentum", "adagrad", "adam", "rmsprop")}. Specifies the
#'   gradient-based optimization approach.
#' @param momentum_coef Numeric. Momentum coefficient (used if \code{optimizer=="momentum"}).
#' @param beta1 Numeric. Exponential decay rate for the first moment estimates (used in Adam).
#' @param beta2 Numeric. Exponential decay rate for the second moment estimates (used in Adam).
#' @param rmsprop_decay Numeric. Decay factor for the squared gradient moving average (used in RMSProp).
#' @param epsilon Numeric. A small constant for numerical stability.
#'
#' @return Updated matrix \code{A} of shape (nAb x d).
#'
#' @details
#' This function performs gradient-based updates on each antibody using the selected loss function.
#' Depending on the chosen \code{optimizer}, it may use plain SGD, momentum-based updates, Adagrad,
#' Adam, or RMSProp.
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
                    huber_delta = 1.0, 
                    verbose = TRUE,
                    optimizer = c("sgd", "momentum", "adagrad", "adam", "rmsprop"),
                    momentum_coef = 0.9,
                    beta1 = 0.9,
                    beta2 = 0.999,
                    rmsprop_decay = 0.9,
                    epsilon = 1e-8) {
  
  task <- match.arg(task)
  loss <- match.arg(loss)
  optimizer <- match.arg(optimizer)
  
  # Ensure X and A are numeric matrices
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
  
  # Convert assignments to numeric indices if needed (e.g., if assignments are characters)
  unique_assignments <- unique(assignments)
  assignments <- match(assignments, unique_assignments)
  nLabels <- length(unique_assignments)
  
  if (nLabels > nAb) {
    if (verbose) {
      message(sprintf("Number of unique assignments (%d) exceeds the number of 
                      ntibodies (nAb = %d).", nLabels, nAb))
    }
  }
  if (nLabels < nAb && verbose) {
    message(sprintf("Number of prototypes (nAb = %d) exceeds unique labels 
                    (%d). Extra prototypes will not be assigned data points.",
                    nAb, nLabels))
  }
  if (any(assignments < 1 | assignments > nAb, na.rm = TRUE)) {
    invalid_assignments <- assignments[assignments < 1 | assignments > nAb]
    stop(sprintf("'assignments' contains invalid values: %s. Valid range is 1..%d.",
                 paste(unique(invalid_assignments), collapse = ", "), nAb))
  }
  
  if (task == "regression" && !is.numeric(y)) {
    stop("y must be numeric for regression.")
  }
  
  # Build a list of indices per antibody for potential classification tasks.
  samples_by_ab <- vector("list", nAb)
  for (i in seq_len(n)) {
    j <- assignments[i]
    samples_by_ab[[j]] <- c(samples_by_ab[[j]], i)
  }
  
  # For classification: determine the dominant label per antibody.
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
  
  # Initialize optimizer state variables as needed.
  if (optimizer == "momentum") {
    velocity <- matrix(0, nrow = nAb, ncol = d)
  } else if (optimizer == "adagrad") {
    acc <- matrix(0, nrow = nAb, ncol = d)
  } else if (optimizer == "adam") {
    m <- matrix(0, nrow = nAb, ncol = d)
    v <- matrix(0, nrow = nAb, ncol = d)
    t_count <- rep(0, nAb)
  } else if (optimizer == "rmsprop") {
    rmsprop_cache <- matrix(0, nrow = nAb, ncol = d)
  }
  
  # Main refinement loop over steps.
  for (st in seq_len(steps)) {
    if (verbose) {
      message(sprintf("Refinement step %d/%d", st, steps))
    }
    # Shuffle sample indices.
    idx_shuf <- sample.int(n)
    for (i in idx_shuf) {
      j <- assignments[i]
      
      # Skip invalid assignments.
      if (is.na(j) || j < 1 || j > nAb) next
      
      # Get the data point (ensure it's a numeric vector of length d).
      x_i <- as.numeric(X[i, ])
      if (length(x_i) != d) {
        stop(sprintf("x_i is not the correct length: expected %d, got %d", d, length(x_i)))
      }
      
      # For classification, determine if the current data point shares the dominant label.
      same_label <- NA
      if (task == "classification" && !is.na(ab_label[j])) {
        same_label <- (y[i] == ab_label[j])
      }
      
      # Compute the gradient using your helper function.
      grad <- .update_prototype(
        ab_vec      = A[j, ],
        x_i         = x_i,
        same_label  = same_label,
        task        = task,
        loss        = loss,
        push_away   = push_away,
        huber_delta = huber_delta
      )
      
      if (length(grad) != d) {
        stop(sprintf(".update_prototype returned a gradient of length %d, expected %d", length(grad), d))
      }
      
      # Optimizer-specific update:
      if (optimizer == "sgd") {
        A[j, ] <- A[j, ] + lr * grad
      } else if (optimizer == "momentum") {
        velocity[j, ] <- momentum_coef * velocity[j, ] + grad
        A[j, ] <- A[j, ] + lr * velocity[j, ]
      } else if (optimizer == "adagrad") {
        acc[j, ] <- acc[j, ] + grad^2
        A[j, ] <- A[j, ] + lr * grad / (sqrt(acc[j, ]) + epsilon)
      } else if (optimizer == "adam") {
        m[j, ] <- beta1 * m[j, ] + (1 - beta1) * grad
        v[j, ] <- beta2 * v[j, ] + (1 - beta2) * (grad^2)
        t_count[j] <- t_count[j] + 1
        m_hat <- m[j, ] / (1 - beta1^(t_count[j]))
        v_hat <- v[j, ] / (1 - beta2^(t_count[j]))
        A[j, ] <- A[j, ] + lr * m_hat / (sqrt(v_hat) + epsilon)
      } else if (optimizer == "rmsprop") {
        rmsprop_cache[j, ] <- rmsprop_decay * rmsprop_cache[j, ] + (1 - rmsprop_decay) * (grad^2)
        A[j, ] <- A[j, ] + lr * grad / (sqrt(rmsprop_cache[j, ]) + epsilon)
      }
    }
  }
  
  return(A)
}
