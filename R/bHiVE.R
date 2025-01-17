#' bHIVE: B-cell-based Hybrid Immune Virtual Evolution 
#'
#' Implements an artificial immune network algorithm for clustering, classification, 
#' and regression tasks. The algorithm evolves a population of "antibodies" via clonal 
#' selection and mutation, applies network suppression to maintain diversity, and assigns
#' data points based on affinity or distance metrics.
#'
#' @param X A numeric matrix or data frame of input features, with rows as observations 
#'   and columns as features.
#' @param y Optional. A target vector. Use for classification (factor) or regression (numeric).
#'   If NULL, clustering will be performed.
#' @param task Character. Specifies the task to perform: \code{"clustering"}, \code{"classification"}, 
#'   or \code{"regression"}. If NULL, it is inferred based on \code{y}.
#' @param nAntibodies Integer. The initial population size of antibodies. 
#' @param beta Numeric. Clone multiplier (controls how many clones are generated for top-matching antibodies).
#' @param epsilon Numeric. Similarity threshold used in network suppression; antibodies closer
#'   than \code{epsilon} are considered redundant.
#' @param maxIter Integer. Maximum number of iterations to run the AI-Net algorithm.
#' @param affinityFunc Character. Specifies the affinity (similarity) function to use for
#'   antibody-data matching. One of \code{"gaussian"}, \code{"laplace"}, \code{"polynomial"}, 
#'   \code{"cosine"}, or \code{"hamming"}.
#' @param distFunc Character. Specifies the distance function for clustering and suppression. 
#'   One of \code{"euclidean"}, \code{"manhattan"}, \code{"minkowski"}, 
#'   \code{"cosine"}, \code{"mahalanobis"}, or \code{"hamming"}.
#' @param affinityParams A list of optional parameters for the chosen affinity or distance function.
#'   \itemize{
#'     \item \code{alpha} (for RBF or Laplace kernel),
#'     \item \code{c}, \code{p} (for polynomial kernel or Minkowski distance),
#'     \item \code{Sigma} (for Mahalanobis distance).
#'   }
#' @param mutationDecay Numeric. Factor by which the mutation rate decays each iteration 
#'   (should be \eqn{\le 1.0}). Default is 1.0 (no decay).
#' @param mutationMin Numeric. Minimum mutation rate, preventing the mutation scale from 
#'   shrinking to zero. 
#' @param maxClones Numeric. Maximum number of clones per top-matching antibody; defaults to \code{Inf}.
#' @param stopTolerance Numeric. If the change in the number of antibodies (repertoire size)
#'   is \eqn{\le stopTolerance} for consecutive iterations, this may trigger the 
#'   \code{noImprovementLimit}.
#' @param noImprovementLimit Integer. Stops the algorithm early if there is no further
#'   improvement in antibody count (beyond \code{stopTolerance}) for this many consecutive
#'   iterations. Default is \code{Inf}, meaning no early stop based on improvement.
#' @param initMethod Character. Method for initializing antibodies. Can be:
#'   \itemize{
#'     \item \code{"sample"} - randomly selects rows from \code{X} as initial antibodies.
#'     \item \code{"random"} - samples Gaussian noise using \code{X}'s column means/sds.
#'   }
#' @param k Integer. Number of top-matching antibodies (by affinity) to consider cloning for each data point.
#' @param verbose Logical. If \code{TRUE}, prints progress messages each iteration.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{antibodies}: A matrix of the final antibody vectors in the repertoire.
#'     \item \code{assignments}: Cluster assignments (for clustering), predicted labels (classification), 
#'       or predicted values (regression).
#'     \item \code{task}: The type of task performed (\code{"clustering"}, \code{"classification"}, \code{"regression"}).
#'   }
#'
#' @examples
#' # Example 1: Clustering with the Iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])  # Numeric features only
#' res <- bHIVE(X = X, 
#'              task = "clustering", 
#'              nAntibodies = 30, 
#'              beta = 5, 
#'              epsilon = 0.01, 
#'              maxIter = 20, 
#'              k = 3, 
#'              verbose = FALSE)
#' table(res$assignments)
#'
#' # Example 2: Classification with Iris species
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
#' table(res$assignments, y)
#'
#' # Example 3: Regression
#' y <- as.numeric(iris$Sepal.Length)
#' res <- bHIVE(X = X, 
#'              y = y, 
#'              task = "regression", 
#'              nAntibodies = 30, 
#'              beta = 5, 
#'              epsilon = 0.01, 
#'              maxIter = 20, 
#'              k = 3, 
#'              verbose = FALSE)
#' cor(res$assignments, y)
#' 
#' @export
bHIVE <- function(X, 
                  y = NULL, 
                  task = NULL, 
                  nAntibodies = 20, 
                  beta = 5, 
                  epsilon = 0.01, 
                  maxIter = 50,
                  affinityFunc = "gaussian", 
                  distFunc = "euclidean",
                  affinityParams = list(alpha = 1, 
                                        c = 1, 
                                        p = 2, 
                                        Sigma = NULL),
                  mutationDecay = 1.0, 
                  mutationMin = 0.01, 
                  maxClones = Inf,
                  stopTolerance = 0.0, 
                  noImprovementLimit = Inf,
                  initMethod = c("sample", "random"), 
                  k = 3, 
                  verbose = TRUE) {
  
  # Infer task if not explicitly provided
  if (is.null(task)) {
    task <- if (is.null(y)) {
      "clustering"
    } else if (is.factor(y)) {
      "classification"
    } else {
      "regression"
    }
  }
  
  # Validate task input
  if (!task %in% c("clustering", "classification", "regression")) {
    stop("Invalid task. Choose from 'clustering', 'classification', or 'regression'.")
  }
  
  initMethod <- match.arg(initMethod)
  n <- nrow(X)
  d <- ncol(X)
  
  if (initMethod == "sample") {
    A <- X[sample(1:n, size = nAntibodies, replace = TRUE), ]
  } else {
    # "random": e.g., random normal around the data mean
    xMean <- colMeans(X)
    xSd <- apply(X, 2, sd) + 1e-8
    A <- matrix(rnorm(nAntibodies * d, mean = 0, sd = 1), nrow = nAntibodies, ncol = d)
    A <- sweep(A, 2, xSd, `*`)  # scale
    A <- sweep(A, 2, xMean, `+`) # shift
  }
  
  affFn <- switch(
    affinityFunc,
    "gaussian"   = .affinity_RBF_custom,
    "laplace"    = .affinity_laplace_custom,
    "polynomial" = .affinity_poly_custom,
    "cosine"     = .affinity_cosine_custom,
    "hamming"    = .affinity_hamming_custom,
    stop("Invalid affinityFunc provided.")
  )
  
  distFn <- switch(
    distFunc,
    "euclidean"   = .dist_euclidean_custom,
    "manhattan"   = .dist_manhattan_custom,
    "minkowski"   = .dist_minkowski_custom,
    "cosine"      = .dist_cosine_custom,
    "mahalanobis" = .dist_mahalanobis_custom,
    "hamming"     = .dist_hamming_custom,
    stop("Invalid distFunc provided.")
  )
  
  # Track iteration-based stopping
  noImprovementCount <- 0
  prevAntibodyCount  <- nrow(A)
  
  for (iter in 1:maxIter) {
    
    # For each data point, clone + mutate top matching antibodies
    for (i in 1:n) {
      x <- X[i, ]
      
      # 1) Compute affinity to all antibodies
      aff_values <- apply(A, 1, function(a) affFn(x, a, affinityParams))
      
      # 2) Identify the top k matching antibodies
      top_idx <- order(aff_values, decreasing = TRUE)[1:k]
      max_aff <- max(aff_values)
      
      # 3) Clone and mutate top k
      for (j in top_idx) {
        f_j <- aff_values[j]
        # # clones (bounded by maxClones)
        nClones <- floor(beta * (f_j / max_aff))
        nClones <- min(nClones, maxClones)  # cap the total clones if desired
        
        if (nClones > 0) {
          for (cloneId in 1:nClones) {
            # Mutation scale (inverse to affinity), decayed by iteration
            # Also clamped by mutationMin
            baseMutation <- (1.0 - f_j) * (mutationDecay^(iter - 1))
            mutation_rate <- max(baseMutation, mutationMin)
            
            # Create mutated clone
            mutated <- A[j, ] + rnorm(d, mean = 0, sd = mutation_rate)
            
            # Evaluate new clone's affinity
            f_mutated <- affFn(x, mutated, affinityParams)
            # If better, replace parent
            if (f_mutated > f_j) {
              A[j, ] <- mutated
              f_j <- f_mutated
            }
          }
        }
      }
    }
    
    # Network Suppression: remove antibodies that are too similar 
    keep <- rep(TRUE, nrow(A))
    for (u in 1:(nrow(A) - 1)) {
      if (!keep[u]) next
      for (v in (u+1):nrow(A)) {
        if (!keep[v]) next
        # Use the chosen distance function to measure similarity
        dist_uv <- distFn(A[u, ], A[v, ], affinityParams) 
        if (dist_uv < epsilon) {
          # remove the second for simplicity
          keep[v] <- FALSE
        }
      }
    }
    A <- A[keep, , drop = FALSE]
    
    # Check Convergence / Early Stopping 
    currentCount <- nrow(A)
    changeInCount <- abs(currentCount - prevAntibodyCount)
    
    if (changeInCount <= stopTolerance) {
      noImprovementCount <- noImprovementCount + 1
    } else {
      noImprovementCount <- 0
    }
    prevAntibodyCount <- currentCount
    
    if (noImprovementCount >= noImprovementLimit) {
      if (verbose) {
        cat("Early stopping due to no improvement for", noImprovementCount, "iters\n")
      }
      break
    }
    
    if (verbose) {
      cat(sprintf("Iteration: %d | #Antibodies: %d | noImprovementCount: %d\n",
                  iter, currentCount, noImprovementCount))
    }
  }
  
  # Task-specific assignments
  if (task == "clustering") {
    cluster_assignments <- apply(X, 1, function(x) {
      distances <- apply(A, 1, function(a) distFn(x, a, affinityParams))
      which.min(distances)
    })
    assignments <- cluster_assignments
  } else if (task == "classification") {
    antibody_classes <- sample(levels(y), nAntibodies, replace = TRUE) # Random initial assignment
    predicted_labels <- apply(X, 1, function(x) {
      affinities <- apply(A, 1, function(a) distFn(x, a, affinityParams))
      antibody_classes[which.max(affinities)]
    })
    assignments <- predicted_labels
  } else if (task == "regression") {
    antibody_values <- runif(nAntibodies, min(y), max(y)) # Random initial values
    predicted_values <- apply(X, 1, function(x) {
      affinities <- apply(A, 1, function(a) distFn(x, a, affinityParams))
      sum(affinities * antibody_values) / sum(affinities)
    })
    assignments <- predicted_values
  }
  
  list(antibodies = A, assignments = assignments, task = task)
}

#Afinity Functions
.affinity_RBF_custom <- function(x, y, params) {
  dist2 <- sum((x - y)^2)
  exp(-params$alpha * dist2)
}
.affinity_laplace_custom <- function(x, y, params) {
  dist1 <- sum(abs(x - y))
  exp(-params$alpha * dist1)
}
.affinity_poly_custom <- function(x, y, params) {
  (sum(x * y) + params$c)^params$p
}
.affinity_cosine_custom <- function(x, y, params) {
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (denom == 0) return(0)
  sum(x * y) / denom
}
.affinity_hamming_custom <- function(x, y, params) {
  x_bin <- as.integer(x)
  y_bin <- as.integer(y)
  matches <- sum(x_bin == y_bin)
  matches / length(x)
}
  
# Distance Functions 
.dist_euclidean_custom <- function(x, y, params) {
    sqrt(sum((x - y)^2))
}
  
.dist_manhattan_custom <- function(x, y, params) {
    sum(abs(x - y))
}
  
.dist_minkowski_custom <- function(x, y, params) {
    p <- params$p
    sum(abs(x - y)^p)^(1/p)
}
  
.dist_cosine_custom <- function(x, y, params) {
    # 1 - Cosine Similarity
    cs <- .affinity_cosine_custom(x, y, params)
    1 - cs
}
  
.dist_mahalanobis_custom <- function(x, y, params) {
    # Sigma must be in params$Sigma
    if (is.null(params$Sigma)) {
      stop("Mahalanobis distance requires a covariance matrix Sigma in distParams.")
    }
    diff <- x - y
    invSigma_diff <- solve(params$Sigma, diff)   # or precompute solve(params$Sigma) once
    sqrt(sum(diff * invSigma_diff))
}
  
.dist_hamming_custom <- function(x, y, params) {
    x_bin <- as.integer(x)
    y_bin <- as.integer(y)
    sum(x_bin != y_bin)
}
