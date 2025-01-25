#' refineB: Gradient-based fine-tuning for bHIVE antibodies with multiple loss functions
#'
#' After running bHIVE (or honeycombHIVE), you have a set of final antibody positions (A)
#' in feature space. This function refines those prototypes by iterating over the
#' data points assigned to each antibody and applying small gradient-based updates
#' according to a user-chosen loss function.
#'
#' ## Classification vs. Regression
#' - **Classification**: The function pushes each antibody closer to data points of
#'   its "dominant" label, and optionally pushes it away from points of other labels.
#'   Multiple loss functions are conceptually adapted (naively) to prototype movement
#'   in feature space.
#' - **Regression**: The function pulls each antibody closer to the data points it
#'   represents, based on an approximate gradient for the chosen numeric loss (e.g. MSE, MAE).
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
#' In reality, each loss function is adapted in a **toy** manner to prototype
#' movement in feature space. True gradient-based training might require separate
#' parameters for the predicted label or value. This function demonstrates how you can
#' incorporate a variety of loss functions if you treat the **prototype** as a parameter
#' to be updated.
#'
#' @param A A matrix of antibody positions (nAntibodies x nFeatures).
#' @param X A matrix or data frame of data (nSamples x nFeatures).
#' @param y The target vector. For classification, a factor; for regression, numeric.
#' @param assignments Integer vector of length nSamples giving the antibody index
#'        that each sample is assigned to (i.e. from \code{which.max(affinity)}).
#' @param loss Character. One of:
#'   - classification: "categorical_crossentropy", "binary_crossentropy", "kullback_leibler", "cosine"
#'   - regression: "mse", "mae", "poisson", "huber"
#' @param task Either "classification" or "regression".
#' @param steps How many gradient steps to run.
#' @param lr Learning rate for the gradient updates.
#' @param push_away Logical (classification only). If `TRUE`, prototypes are pushed away
#'        from points with different labels. If `FALSE`, we only pull toward same-label points.
#' @param huber_delta Numeric, used for the "huber" loss threshold.
#'
#' @return A refined (updated) version of A with the same dimensions.
#'
#' @examples
#  data(iris)
#' X <- as.matrix(iris[, 1:4])  # Numeric features only
#'
#' # Runnning bHIVE Classification
#' y <- iris$Species
#' res <- bHIVE(X = X, 
#'               y = y, 
#'               task = "classification", 
#'               nAntibodies = 30, 
#'               beta = 5, 
#'               epsilon = 0.01, 
#'               maxIter = 20, 
#'               k = 3, 
#'               verbose = FALSE)
#'               
#'  A <- res_bhive$antibodies
#'  assignments <- res_bhive$assignments 
#' 
#' #Refining Antibodies/Bees              
#' refineB(A,
#'                  X, 
#'                  y, 
#'                  assignments,
#'                  loss="mse", 
#'                  task="regression", 
#'                  steps=5, 
#'                  lr=0.01)
#'
#' @export
refineB <- function(A, 
                             X, 
                             y, 
                             assignments,
                             loss = c("categorical_crossentropy","binary_crossentropy",
                                      "kullback_leibler","cosine",
                                      "mse","mae","poisson","huber"),
                             task = c("classification","regression"),
                             steps = 5,
                             lr = 0.01,
                             push_away = TRUE,
                             huber_delta = 1.0) {
  
  task <- match.arg(task)
  loss <- match.arg(loss)
  
  # Convert A, X to numeric matrices
  A <- as.matrix(A)
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  nAb <- nrow(A)
  
  if (length(assignments) != n) {
    stop("Length of 'assignments' must match the number of rows in X.")
  }
  if (task == "classification" && !is.factor(y)) {
    warning("Task='classification' but y is not a factor. Converting to factor forcibly.")
    y <- as.factor(y)
  }
  if (task == "regression" && !is.numeric(y)) {
    stop("Task='regression' but y is not numeric.")
  }
  
  classification_loss_funs <- list(
    categorical_crossentropy = .grad_categorical_ce,
    binary_crossentropy      = .grad_binary_ce,
    kullback_leibler        = .grad_kullback_leibler,
    cosine                  = .grad_cosine_class
  )
  
  regression_loss_funs <- list(
    mse     = .grad_mse,
    mae     = .grad_mae,
    poisson = .grad_poisson,
    huber   = .grad_huber
  )
  
 
  # Precompute "dominant label" for classification 
  samples_by_ab <- split(seq_len(n), assignments)
  
  if (task == "classification") {
    # For each antibody, find the dominant label
    # We store in a vector or list: ab_label[[j]]
    ab_label <- vector("list", nAb)
    for (j in seq_len(nAb)) {
      pts_j <- samples_by_ab[[j]]
      if (length(pts_j) == 0) {
        ab_label[[j]] <- NA
      } else {
        y_sub <- y[pts_j]
        tbl <- table(y_sub)
        ab_label[[j]] <- names(tbl)[which.max(tbl)]
      }
    }
  }
  
  # The chosen loss function
  if (task == "classification") {
    if (!loss %in% names(classification_loss_funs)) {
      stop("For classification, loss must be one of: ",
           paste(names(classification_loss_funs), collapse=", "))
    }
    loss_fun <- classification_loss_funs[[loss]]
  } else {
    if (!loss %in% names(regression_loss_funs)) {
      stop("For regression, loss must be one of: ",
           paste(names(regression_loss_funs), collapse=", "))
    }
    loss_fun <- regression_loss_funs[[loss]]
  }
  
  # Gradient Steps
  for (step_i in seq_len(steps)) {
    # Shuffle data points for stochastic order
    idx_shuf <- sample.int(n)
    
    for (i in idx_shuf) {
      ab_j <- assignments[i]
      if (ab_j < 1 || ab_j > nAb) next
      
      x_i <- X[i, ]
      if (task == "classification") {
        dom_label_j <- ab_label[[ab_j]]
        if (is.na(dom_label_j)) {
          # no assigned label => skip
          next
        }
        same_class <- (y[i] == dom_label_j)
        
        # If same_class => "pull"; if different_class => "push" (if push_away=TRUE).
        if (same_class) {
          grad <- loss_fun(x_i, TRUE, A[ab_j, ])
        } else if (push_away) {
          grad <- loss_fun(x_i, FALSE, A[ab_j, ])
        } else {
          grad <- rep(0, d)
        }
        
        A[ab_j, ] <- A[ab_j, ] + lr * grad
        
      } else { 
        # regression
        grad <- loss_fun(x_i, A[ab_j, ])
        # if it's huber, we might pass huber_delta in, but we did it internally for demonstration
        A[ab_j, ] <- A[ab_j, ] + lr * grad
      }
    }
  }
  
  return(A)
}

# Internal Loss Functions

# Regression subfunctions (pull in feature space)
.grad_mse <- function(x_i, ab_vec) {
  # Move ab_vec slightly toward x_i
  return(x_i - ab_vec)
}

.grad_mae <- function(x_i, ab_vec) {
  # sign-based movement
  return(sign(x_i - ab_vec))
}

.grad_poisson <- function(x_i, ab_vec) {
  # Typically Poisson deviance = 2 * [ y log(y/pred) - (y-pred) ], 
  # but we are not storing "pred" as a separate param. 
  # pull ab_vec -> x_i or scale by some factor. 
  return(x_i - ab_vec)
}

.grad_huber <- function(x_i, ab_vec) {
  # Huber loss = if |res| <= delta => 0.5 res^2, else delta(|res|-0.5 delta)
  residual <- x_i - ab_vec
  dist_abs <- sqrt(sum(residual^2))
  if (dist_abs <= huber_delta) {
    # MSE region
    return(residual)
  } else {
    # L1 region
    return(huber_delta * sign(residual))
  }
}

# Classification subfunctions (push/pull)
.grad_categorical_ce <- function(x_i, same_class, ab_vec) {
  # If same_class, pull ab_vec -> x_i
  # If different_class, push ab_vec <- x_i if push_away=TRUE
  if (same_class) {
    return(x_i - ab_vec) 
  } else {
    return(ab_vec - x_i)
  }
}

.grad_binary_ce <- function(x_i, same_class, ab_vec) {
  # Conceptually the same as categorical CE for 2-class. 
  # We'll do the same push/pull approach.
  if (same_class) {
    return(x_i - ab_vec)
  } else {
    return(ab_vec - x_i)
  }
}

.grad_kullback_leibler <- function(x_i, same_class, ab_vec) {
  # Very naive: we treat "p vs q" as a distribution difference. 
  if (same_class) {
    return(x_i - ab_vec)
  } else {
    return(ab_vec - x_i)
  }
}

.grad_cosine_class <- function(x_i, same_class, ab_vec) {
  # We want to maximize cos similarity if same_class, minimize if different_class.
  # cos_sim(a,b) = a.b / (||a|| ||b||).
  if (same_class) {
    return(x_i - ab_vec)
  } else {
    return(ab_vec - x_i)
  }
}
