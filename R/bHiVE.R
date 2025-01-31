#' bHIVE: B-cell Hybrid Immune Variant Engine
#'
#' Implements an artificial immune network algorithm for clustering, classification, 
#' and regression tasks. The algorithm evolves a population of "antibodies" 
#' via clonal selection and mutation, applies network suppression to maintain 
#' diversity, and assigns data points based on affinity or distance metrics.
#'
#' @param X A numeric matrix or data frame of input features, with rows as 
#' observations 
#'   and columns as features.
#' @param y Optional. A target vector. Use for classification (factor) or 
#' regression (numeric).
#'   If NULL, clustering will be performed.
#' @param task Character. Specifies the task to perform: \code{"clustering"}, 
#' \code{"classification"}, or \code{"regression"}. If NULL, it is inferred 
#' based on \code{y}.
#' @param nAntibodies Integer. The initial population size of antibodies. 
#' @param beta Numeric. Clone multiplier (controls how many clones are 
#' generated for top-matching antibodies).
#' @param epsilon Numeric. Similarity threshold used in network suppression; 
#' antibodies closer than \code{epsilon} are considered redundant.
#' @param maxIter Integer. Maximum number of iterations to run the AI-Net 
#' algorithm.
#' @param affinityFunc Character. Specifies the affinity (similarity) function 
#' to use for antibody-data matching. One of \code{"gaussian"}, \code{"laplace"}, 
#' \code{"polynomial"}, \code{"cosine"}, or \code{"hamming"}.
#' @param distFunc Character. Specifies the distance function for clustering 
#' and suppression. One of \code{"euclidean"}, \code{"manhattan"}, 
#' \code{"minkowski"}, \code{"cosine"}, \code{"mahalanobis"}, 
#' or \code{"hamming"}.
#' @param affinityParams A list of optional parameters for the chosen affinity 
#' or distance function.
#'   \itemize{
#'     \item \code{alpha} (for RBF or Laplace kernel),
#'     \item \code{c}, \code{p} (for polynomial kernel or Minkowski distance),
#'     \item \code{Sigma} (for Mahalanobis distance).
#'   }
#' @param mutationDecay Numeric. Factor by which the mutation rate decays each 
#' iteration (should be \eqn{\le 1.0}). Default is 1.0 (no decay).
#' @param mutationMin Numeric. Minimum mutation rate, preventing the mutation 
#' scale from shrinking to zero. 
#' @param maxClones Numeric. Maximum number of clones per top-matching antibody;
#'  defaults to \code{Inf}.
#' @param stopTolerance Numeric. If the change in the number of antibodies 
#' (repertoire size) is \eqn{\le stopTolerance} for consecutive iterations, 
#' this may trigger the \code{noImprovementLimit}.
#' @param noImprovementLimit Integer. Stops the algorithm early if there is no 
#' further improvement in antibody count (beyond \code{stopTolerance}) for this
#' many consecutive iterations. Default is \code{Inf}, meaning no early stop 
#' based on improvement.
#' @param initMethod Character. Method for initializing antibodies. Can be:
#'   \itemize{
#'     \item \code{"sample"} - randomly selects rows from \code{X} as initial 
#'     antibodies.
#'     \item \code{"random"} - samples Gaussian noise using \code{X}'s column 
#'     means/sds.
#'     \item \code{"random_uniform"} - samples uniformly in [min, max] of each 
#'     column.
#'     \item \code{"kmeans++"} - tries a kmeans++-like initialization for 
#'     coverage.
#'   }
#' @param k Integer. Number of top-matching antibodies (by affinity) to 
#' consider cloning for each data point.
#' @param verbose Logical. If \code{TRUE}, prints progress messages each 
#' iteration.
#'
#' @return A list:
#'   \itemize{
#'     \item \code{antibodies}: Final antibody vectors (nAntibodies x nFeatures).
#'     \item \code{assignments}: 
#'         - For clustering: integer cluster IDs in [1..#Antibodies].
#'         - For classification: predicted labels.
#'         - For regression: integer cluster index (in [1..#Antibodies]) if used in synergy with \code{refineB}.
#'     \item \code{predictions}: Only for regression, the numeric predictions per row.
#'     \item \code{task}: The chosen task.
#'   }
#'
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
#' @importFrom stats rnorm runif sd
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
                  initMethod = c("sample", "random", "random_uniform", "kmeans++"), 
                  k = 3, 
                  verbose = TRUE) {
  #TODO Parallel or batch approach to affinity matrix
  # ====================================
  # 0) Basic Validation & Task Inference
  # ====================================
  .validate_bHIVE_input(X, y)  
  
  if (is.null(task)) {
    if (is.null(y)) {
      task <- "clustering"
    } else if (is.factor(y)) {
      task <- "classification"
    } else {
      task <- "regression"
    }
  }
  task <- match.arg(task, c("clustering","classification","regression"))
  
  initMethod <- match.arg(initMethod, c("sample", "random", "random_uniform","kmeans++"))
  
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  # ===================
  # 1. Antibody Initialization
  # ===================
  A <- switch(
    initMethod,
    "sample" = { X[sample.int(n, 
                              size = nAntibodies, 
                              replace = TRUE), , drop=FALSE] },
    "random" = { xMean <- colMeans(X)
                 xSd   <- apply(X, 2, sd) + 1e-8
                 mat   <- matrix(rnorm(nAntibodies * d), nrow = nAntibodies)
                 mat   <- sweep(mat, 2, xSd, `*`)
                 sweep(mat, 2, xMean, `+`)},
    "random_uniform" = { xMin <- apply(X, 2, min)
                         xMax <- apply(X, 2, max)
                         mat  <- matrix(runif(nAntibodies * d), nrow = nAntibodies)
                         for (col_i in seq_len(d)) {
                           range_i <- xMax[col_i] - xMin[col_i]
                           mat[, col_i] <- xMin[col_i] + range_i * mat[, col_i]
                         }
                        mat},
    "kmeans++" = {
      .init_kmeanspp(X, nAntibodies)  # Assume you have a kmeans++ init function
    }
  )
  if (!is.matrix(A) || nrow(A) == 0) {
    stop("Initialization of antibodies failed. Check nAntibodies and input X.")
  }
  
  m <- nrow(A) 
  
  # ===================
  # 2. Pick Affinity & Distance
  # ===================
  affFn <- switch(affinityFunc,
                  "gaussian"   = .affinity_RBF_custom,
                  "laplace"    = .affinity_laplace_custom,
                  "polynomial" = .affinity_poly_custom,
                  "cosine"     = .affinity_cosine_custom,
                  "hamming"    = .affinity_hamming_custom,
                  stop("Invalid affinityFunc."))
  distFn <- switch(distFunc,
                   "euclidean"   = .dist_euclidean_custom,
                   "manhattan"   = .dist_manhattan_custom,
                   "minkowski"   = .dist_minkowski_custom,
                   "cosine"      = .dist_cosine_custom,
                   "mahalanobis" = .dist_mahalanobis_custom,
                   "hamming"     = .dist_hamming_custom,
                   stop("Invalid distFunc."))
  
  # ===================
  # 3. Task-Specific Setup
  # ===================
  if (task == "classification") {
    classes <- levels(y)
    nClasses <- length(classes)
    class_counts <- matrix(0, nrow = m, ncol = nClasses)
    colnames(class_counts) <- classes
  } else if (task == "regression") {
    sum_y   <- numeric(m)
    sum_aff <- numeric(m)
    overall_mean <- mean(y, na.rm=TRUE)  # only compute once
  }
  
  # For early stopping
  noImproveCount <- 0
  prevCount <- m
  
  # ======================
  # 4) Main Iteration Loop
  # ======================
  
  for (iter in seq_len(maxIter)) {
    # reset counters
    if (task=="classification") {
      class_counts[] <- 0
    } else if (task=="regression") {
      sum_y[] <- 0
      sum_aff[] <- 0
    }
    
    # For each data point
    for (i in seq_len(n)) {
      x_i <- X[i,]
      
      # compute affinity
      aff_values <- numeric(m)
      for (j in seq_len(m)) {
        aff_values[j] <- affFn(x_i, A[j, ], affinityParams)
      }
      
      # If max affinity is zero or NA => skip
      max_aff <- max(aff_values, na.rm = TRUE)
      if (max_aff == 0 || is.na(max_aff)) next
      
      #Identify top k antibodies
      k2 <- min(k, m)
      top_idx <- sort.int(aff_values, decreasing=TRUE, index.return=TRUE)$ix[1:k2]
      
      # classification/regression counters
      if (task == "classification") {
        y_class <- as.character(y[i])
        class_col <- match(y_class, colnames(class_counts))
        # Weighted vote
        for (jj in top_idx) {
          class_counts[jj, class_col] <- class_counts[jj, class_col] + aff_values[jj]
        }
      } else if (task == "regression") {
        # Weighted sum
        y_val <- y[i]
        for (jj in top_idx) {
          sum_y[jj]   <- sum_y[jj] + y_val * aff_values[jj]
          sum_aff[jj] <- sum_aff[jj] + aff_values[jj]
        }
      }
      
      # clone/mutate
      for (jj in top_idx) {
        f_j <- aff_values[jj]
        nClones <- min(maxClones, floor(beta * (f_j / max_aff)))
        if (nClones <= 0) next
        
        for (clone_id in seq_len(nClones)) {
          # decayed mutation rate for iteration
          mutation_rate <- max((1.0 - f_j) * mutationDecay^(iter - 1), mutationMin)
          # propose mutated antibody
          mutated   <- A[jj, ] + rnorm(d, mean=0, sd=mutation_rate)
          f_mutated <- affFn(x_i, mutated, affinityParams)
          if (f_mutated > f_j) {
            A[jj, ] <- mutated  # keep improvement
          }
        }
      }
    }
    
    # update classification/regression
    if (task=="classification") {
      # each antibody's label is the class with largest class_counts row
      antibody_classes <- apply(class_counts, 1, function(row) {
        if (all(row==0)) {
          # fallback
          colnames(class_counts)[sample(ncol(class_counts),1)]
        } else {
          colnames(class_counts)[which.max(row)]
        }
      })
    } else if (task=="regression") {
      antibody_values <- ifelse(sum_aff>0, sum_y/sum_aff, overall_mean)
    }
    
    # ======================
    # 5) Network Suppression
    # =======================
    keep <- rep(TRUE, m)
    for (u in seq_len(m - 1)) {
      if (!keep[u]) next
      for (v in seq.int(u+1, m)) {
        if (!keep[v]) next
        #TODO Possible point of improvement approximate with RANN or Annoy
        dist_uv <- distFn(A[u, ], A[v, ], affinityParams)
        if (dist_uv < epsilon) {
          keep[v] <- FALSE
        }
      }
    }
    # Subset the antibodies
    A <- A[keep, , drop=FALSE]
    m_new <- nrow(A)
    
    if (task=="classification") {
      class_counts <- class_counts[keep,,drop=FALSE]
    } else if (task=="regression") {
      #antibody_values <- antibody_values[keep]
      antibody_values <- ifelse(keep, antibody_values, mean(y, na.rm = TRUE))
    }
    
    #If suppressed everything => abort
    if (m_new == 0) {
      stop("All antibodies were suppressed. Increase nAntibodies or decrease epsilon.")
    }
    
    # For next iteration
    m <- m_new
    
    # ========================
    # 6) Early Stopping Check
    # ========================
    changeCount <- abs(m - prevCount)
    if (changeCount <= stopTolerance) {
      noImproveCount <- noImproveCount + 1
    } else {
      noImproveCount <- 0
    }
    prevCount <- m
    
    if (noImproveCount >= noImprovementLimit) {
      if (verbose) {
        cat("Early stopping: no improvement for", noImproveCount, "iterations.\n")
      }
      break
    }
    
    if (verbose) {
      cat(sprintf("Iteration %d | #Antibodies: %d | noImproveCount: %d\n",
                  iter, m, noImproveCount))
    }
  }
  
  # =====================
  # 7) Final assignments
  # =====================
  if (task == "clustering") {
    # Find nearest antibody by distance => cluster IDs
    assignments <- integer(n)
    for (i in seq_len(n)) {
      x_i <- X[i, ]
      min_dist <- Inf
      best_j <- NA
      for (j in seq_len(m)) {
        d_j <- distFn(x_i, A[j, ], affinityParams)
        if (d_j < min_dist) {
          min_dist <- d_j
          best_j <- j
        }
      }
      assignments[i] <- best_j
    }
    # optionally re-label cluster IDs from 1..m
    assignments <- as.numeric(factor(assignments))
    
    res <- list(
      antibodies  = A,
      assignments = assignments,
      task        = task
    )
    
  } else if (task == "classification") {
    # For each data row, choose the antibody with largest affinity
    assignments <- character(n)
    for (i in seq_len(n)) {
      x_i <- X[i, ]
      best_aff <- -Inf
      best_j <- 1L
      for (j in seq_len(m)) {
        f_j <- affFn(x_i, A[j, ], affinityParams)
        if (f_j > best_aff) {
          best_aff <- f_j
          best_j <- j
        }
      }
      assignments[i] <- antibody_classes[best_j]
    }
    
    res <- list(
      antibodies  = A,
      assignments = assignments,
      task        = task
    )
    
  } else {
    # Regression
    # predictions for each data row => weighted average of antibody_values
    predictions <- numeric(n)
    for (i in seq_len(n)) {
      x_i <- X[i, ]
      aff_vec <- numeric(m)
      for (j in seq_len(m)) {
        aff_vec[j] <- affFn(x_i, A[j, ], affinityParams)
      }
      s_aff <- sum(aff_vec, na.rm=TRUE)
      if (s_aff == 0) {
        predictions[i] <- overall_mean
      } else {
        predictions[i] <- sum(aff_vec * antibody_values, na.rm=TRUE) / s_aff
      }
    }
    
    # produce integer cluster index => nearest by distance
    cluster_assign <- integer(n)
    for (i in seq_len(n)) {
      x_i <- X[i, ]
      min_dist <- Inf
      best_j <- 1L
      for (j in seq_len(m)) {
        d_j <- distFn(x_i, A[j, ], affinityParams)
        if (d_j < min_dist) {
          min_dist <- d_j
          best_j <- j
        }
      }
      cluster_assign[i] <- best_j
    }
    
    res <- list(
      antibodies  = A,
      assignments = cluster_assign,  # integer cluster IDs
      predictions = predictions,     # numeric predicted values
      task        = task
    )
  }
  
  return(res)
}


# ---------------------------
# HELPER FUNCTIONS
# ---------------------------

# Simple kmeans++ style initialization
.init_kmeanspp <- function(X, nCenters) {
  X <- as.matrix(X)
  n <- nrow(X)
  d <- ncol(X)
  
  # 1) choose one center uniformly at random
  centers <- matrix(0, nrow = nCenters, ncol = d)
  idx <- sample(n, 1)
  centers[1, ] <- X[idx, ]
  
  # 2) For each data point x, compute D(x) = min distance to any chosen center
  # 3) Choose a new data point at random weighted by D(x)^2
  if (nCenters > 1) {
    for (cId in 2:nCenters) {
      dists <- sapply(1:n, function(i) {
        min(rowSums((centers[1:(cId-1), , drop = FALSE] - X[i, ])^2))
      })
      probs <- dists / sum(dists)
      idx <- sample(n, 1, prob = probs)
      centers[cId, ] <- X[idx, ]
    }
  }
  centers
}


# AFFINITY FUNCTIONS
.affinity_RBF_custom <- function(x, y, params) {
  # Gaussian (RBF) kernel: exp(-alpha * ||x-y||^2)
  dist2 <- sum((x - y)^2)
  exp(-params$alpha * dist2)
}

.affinity_laplace_custom <- function(x, y, params) {
  # Laplace kernel: exp(-alpha * ||x-y||_1)
  dist1 <- sum(abs(x - y))
  exp(-params$alpha * dist1)
}

.affinity_poly_custom <- function(x, y, params) {
  # Polynomial kernel: (x.y + c)^p
  (sum(x * y) + params$c)^params$p
}

.affinity_cosine_custom <- function(x, y, params) {
  # Cosine similarity: (x.y)/(||x||*||y||)
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (denom == 0) return(0)
  sum(x * y) / denom
}

.affinity_hamming_custom <- function(x, y, params) {
  # Hamming similarity:
  x_bin <- as.integer(x)
  y_bin <- as.integer(y)
  matches <- sum(x_bin == y_bin)
  matches / length(x)
}

# DISTANCE FUNCTIONS
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
  # 1 - Cosine similarity
  cs <- .affinity_cosine_custom(x, y, params)
  1 - cs
}

.dist_mahalanobis_custom <- function(x, y, params) {
  # Sigma must be in params$Sigma
  if (is.null(params$Sigma)) {
    stop("Mahalanobis distance requires a covariance matrix Sigma in distParams.")
  }
  diff <- x - y
  invSigma_diff <- solve(params$Sigma, diff)  
  sqrt(sum(diff * invSigma_diff))
}

.dist_hamming_custom <- function(x, y, params) {
  x_bin <- as.integer(x)
  y_bin <- as.integer(y)
  sum(x_bin != y_bin)
}
