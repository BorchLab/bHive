# -------------------------------------------------------
# Example of a minimal aiNet-like algorithm in R with clustering
# -------------------------------------------------------

#TODO For Mahalanobis distance, you’ll need to provide a covariance matrix relevant to your data.
#TODO For Hamming metrics, ensure your data is in the proper binary or integer format.
aiNet_with_clustering <- function(X, 
                                  nAntibodies = 20,    # Initial population size
                                  beta = 5,            # Clone multiplier
                                  epsilon = 0.01,      # Similarity threshold
                                  maxIter = 50,        # Maximum iterations
                                  alpha = 1,           # Rate parameter in affinity function
                                  affinityFunc = "gaussian",
                                  distFunc = "euclidean",
                                  k = 3,               # number of top-matching antibodies
                                  verbose = TRUE) {
  
  affinityFunc <- switch(affinityFunc,
                        "gaussian" = .affinity_RBF,
                        "cosine"  = .affinity_cosine,
                        "laplace"      = .affinity_laplace,
                        "ploynomial"  = .affinity_poly,
                        "hamming"  =  .affinity_hamming,
                        stop("Invalid affinityFunc provided"))
  
  distFunc <- switch(distFunc,
                     "euclidean" = .dist_euclidean,
                     "manhatten"  = .dist_manhatten,
                     "minkowski"      = .dist_minkowski,
                     "cosine"  = .dist_cosine,
                     "mahalanobis"  =  .dist_mahalanobis,
                     "hamming" = .dist_hamming,
                     stop("Invalid distFunc provided"))
  
  set.seed(123)
  n <- nrow(X)
  d <- ncol(X)
  
  # Initialize antibody set A randomly from the data
  A <- X[sample(1:n, size = nAntibodies, replace = TRUE), ]
  
  for (iter in 1:maxIter) {
    
    # For each data point, clone and mutate top matching antibodies
    for (i in 1:n) {
      x <- X[i, ]
      
      # Compute affinity of x to all antibodies in A
      aff_values <- apply(A, 1, function(a) affinityFunc(x, a))
      
      # Identify the top k matching antibodies
      top_idx <- order(aff_values, decreasing = TRUE)[1:k]
      
      # For each of these top antibodies, clone and mutate
      for (j in top_idx) {
        f_j <- aff_values[j]
        nClones <- floor(beta * (f_j / max(aff_values)))
        
        if (nClones > 0) {
          for (cloneId in 1:nClones) {
            # Mutation scale inversely proportional to affinity
            mutation_rate <- 1.0 - f_j
            
            # Create mutated clone
            mutated <- A[j, ] + rnorm(d, mean = 0, sd = mutation_rate)
            
            # Evaluate new clone's affinity
            f_mutated <- affinity(x, mutated)
            if (f_mutated > f_j) {
              # If clone is better, replace the parent
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
        dist_uv <- sqrt(sum((A[u, ] - A[v, ])^2))
        if (dist_uv < epsilon) {
          # Remove the lower-affinity one; 
          # or just remove the second for simplicity
          keep[v] <- FALSE
        }
      }
    }
    A <- A[keep, , drop = FALSE]
    
    if (verbose) {
      cat("Iteration:", iter, "| #Antibodies:", nrow(A), "\n")
    }
  }
  
  # --------------------------------------------------
  # Clustering step: Assign each data point in X to 
  # the nearest antibody in the final set A
  # --------------------------------------------------
  
  # For convenience, define a distance function 
  # (rather than the affinity).
  
  cluster_assignments <- integer(n)  # store cluster labels for each point
  for (i in 1:n) {
    dists <- apply(A, 1, function(a) distFunc(X[i, ], a))
    cluster_assignments[i] <- which.min(dists)  # 1-based index of the nearest antibody
  }
  
  # Return both the final antibody set and the cluster assignments
  return(list(
    antibodies = A,
    cluster = cluster_assignments
  ))
}

.affinity_RBF <- function(x, y) {
  dist2 <- sum((x - y)^2)
  exp(-alpha * dist2)
}

.affinity_cosine <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

.affinity_laplace <- function(x, y, alpha = 1) {
  dist1 <- sum(abs(x - y))
  exp(-alpha * dist1)
}

.affinity_poly <- function(x, y, c = 1, p = 2) {
  (sum(x * y) + c)^p
}

.affinity_hamming <- function(x, y) {
  # x and y should be integer or logical vectors of the same length
  # Convert them to 0/1 if logical:
  x_bin <- as.integer(x)
  y_bin <- as.integer(y)
  
  # Number of matches
  matches <- sum(x_bin == y_bin)
  # total length
  matches / length(x)
}

.dist_euclidean <- function(x, y) {
  sqrt(sum((x - y)^2))
}

.dist_manhattan <- function(x, y) {
  sum(abs(x - y))
}

.dist_minkowski <- function(x, y, p = 2) {
  sum(abs(x - y)^p)^(1 / p)
}

dist_cosine <- function(x, y) {
  # 1 - Cosine Similarity
  cos_sim <- affinity_cosine(x, y)
  1 - cos_sim
}

.dist_mahalanobis <- function(x, y, Sigma) {
  # x, y: numeric vectors
  # Sigma: covariance matrix (must be invertible)
  
  diff <- x - y
  # Solve for Sigma^-1 * diff
  invSigma_diff <- solve(Sigma, diff)
  sqrt(sum(diff * invSigma_diff))
}

.dist_hamming <- function(x, y) {
  x_bin <- as.integer(x)
  y_bin <- as.integer(y)
  sum(x_bin != y_bin)
}

data(iris)
X <- as.matrix(iris[, 1:4])   # 4D data (Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)

# Run the aiNet
res <- aiNet_with_clustering(
  X,
  nAntibodies = 20,
  beta = 5,
  epsilon = 0.1,
  maxIter = 20,
  alpha = 1,
  k = 3, # top 3 matching antibodies
  verbose = TRUE
)

# Extract results
antibodies_final <- res$antibodies
cluster_labels   <- res$cluster
