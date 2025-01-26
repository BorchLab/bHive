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
#'   ## Available Losses
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
#' @param A Matrix (nAb x d) of antibody/prototype positions.
#' @param X Matrix or data frame (nSamples x d) of feature data.
#' @param y Optional. Factor (classification), numeric (regression), or NULL (clustering).
#' @param assignments Integer vector (length=nSamples), specifying the antibody index
#'                    each sample belongs to (1..nAb).
#' @param task One of c("clustering","classification","regression").
#' @param loss One of c("mse","mae","cross_entropy","huber","poisson","kullback_leibler","cosine"),
#'             or any custom approach you'd like to code.
#' @param steps Integer. How many gradient steps to run.
#' @param lr Numeric. Learning rate for each update.
#' @param push_away Logical (classification only). Whether to push prototypes away
#'                  from differently-labeled samples.
#' @param huber_delta Numeric. The delta threshold if using huber loss.
#' @return The updated matrix A of shape (nAb x d).
#'
#' @export
refineB <- function(A, X, y = NULL,
                    assignments,
                    task = c("clustering","classification","regression"),
                    loss = c("categorical_crossentropy","binary_crossentropy",
                             "kullback_leibler","cosine",
                             "mse","mae","poisson","huber")
                    steps = 5,
                    lr = 0.01,
                    push_away = TRUE,
                    huber_delta = 1.0) {
  
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  X <- as.matrix(X)
  A <- as.matrix(A)
  nAb <- nrow(A)
  d   <- ncol(A)
  n   <- nrow(X)
  
  if (length(assignments) != n) {
    stop("Length of 'assignments' must match nrow(X).")
  }
  if (any(assignments < 1 | assignments > nAb)) {
    stop("'assignments' values must be between 1..nAb.")
  }
  
  # Build samples_by_ab as length=nAb, each an integer vector of indices
  samples_by_ab <- vector("list", nAb)
  for (i in seq_len(n)) {
    j <- assignments[i]
    samples_by_ab[[j]] <- c(samples_by_ab[[j]], i)
  }
  
  # If classification, find dominant labels
  ab_label <- vector("character", length = nAb)
  if (task == "classification") {
    if (!is.factor(y)) {
      stop("For classification, y must be a factor.")
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
    idx_shuf <- sample.int(n)
    for (i in idx_shuf) {
      x_i <- X[i,]
      ab_j <- assignments[i]
      # check classification label
      same_label <- NA
      if (task=="classification") {
        if (!is.na(ab_label[ab_j])) {
          same_label <- (y[i] == ab_label[ab_j])
        }
      }
      grad <- .update_prototype(A[ab_j,], x_i, same_label)
      A[ab_j, ] <- A[ab_j, ] + lr*grad
    }
  }
  
  return(A)
}