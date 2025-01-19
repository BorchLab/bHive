#' Tune Hyperparameters for bHIVE
#'
#' Performs hyperparameter tuning for the bHIVE algorithm using grid search. 
#' The function evaluates multiple combinations of hyperparameters to find the 
#' best parameters based on a user-defined evaluation metric.
#'
#' @param X A numeric matrix or data frame of input features, with rows as 
#' observations 
#'   and columns as features.
#' @param y Optional. A target vector. Use for classification (factor) or 
#' regression (numeric).
#'   If NULL, clustering will be performed.
#' @param task Character. Specifies the task to perform: \code{"clustering"},
#'  \code{"classification"}, 
#'   or \code{"regression"}.
#' @param grid A data frame specifying the hyperparameter combinations to 
#' evaluate. 
#'   Should include columns \code{nAntibodies}, \code{beta}, and \code{epsilon}.
#' @param metric Character. The evaluation metric to use for selecting the best 
#' hyperparameters. 
#'   Options are:
#'   \itemize{
#'     \item \code{"accuracy"}: For classification tasks.
#'     \item \code{"rmse"}: For regression tasks.
#'     \item \code{"silhouette"}: For clustering tasks.
#'   }
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{best_params}: A list of the best hyperparameters.
#'     \item \code{results}: A data frame with the grid search results, 
#'     including the performance metric for each combination.
#'   }
#'
#' @examples
#' # Load Iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])  # Features
#' y <- iris$Species            # Labels for classification
#'
#' # Define hyperparameter grid
#' grid <- expand.grid(
#'   nAntibodies = c(10, 20, 30),
#'   beta = c(3, 5, 7),
#'   epsilon = c(0.01, 0.05, 0.1)
#' )
#'
#' # Perform hyperparameter tuning for classification
#' tuning_results <- swarmbHIVE(X = X, 
#'                             y = y,
#'                              task = "classification", 
#'                              grid = grid, 
#'                              metric = "accuracy")
#'
#' # Best parameters
#' print(tuning_results$best_params)
#' @importFrom magrittr %>%
#' @importFrom cluster silhouette
#' @importFrom stats dist
#' @export
swarmbHIVE <- function(X, 
                       y = NULL, 
                       task = "clustering", 
                       grid, metric = "silhouette", 
                       verbose = TRUE) {
  
  # Validate input data
  .validate_bHIVE_input(X, y)
  
  # Validate task
  if (!task %in% c("clustering", "classification", "regression")) {
    stop("Invalid task. Choose from 'clustering', 'classification', or 
         'regression'.")
  }
  
  # Validate metric
  valid_metrics <- c("accuracy", "rmse", "silhouette")
  if (!metric %in% valid_metrics) {
    stop("Invalid metric. Choose from 'accuracy', 'rmse', or 'silhouette'.")
  }
  
  # Initialize results data frame
  results <- data.frame(
    nAntibodies = numeric(),
    beta = numeric(),
    epsilon = numeric(),
    metric_value = numeric()
  )
  
  # Grid search
  for (i in seq_len(nrow(grid))) {
    params <- grid[i, ]
    
    if (verbose) {
      cat(sprintf(
        "Evaluating combination: nAntibodies = %d, beta = %d, epsilon = %.3f\n",
        params$nAntibodies, params$beta, params$epsilon
      ))
    }
    
    # Run bHIVE with current parameters
    model <- bHIVE(
      X = X,
      y = y,
      task = task,
      nAntibodies = params$nAntibodies,
      beta = params$beta,
      epsilon = params$epsilon,
      maxIter = 50,
      verbose = FALSE
    )
    
    # Evaluate performance
    if (task == "clustering") {
      cluster_labels <- model$assignments
      metric_value <- if (metric == "silhouette") {
        silhouette(cluster_labels, dist(X))[, 3] %>% mean(na.rm = TRUE)
      } else {
        stop("Invalid metric for clustering.")
      }
    } else if (task == "classification") {
      predicted_labels <- model$assignments
      metric_value <- if (metric == "accuracy") {
        mean(predicted_labels == y)
      } else {
        stop("Invalid metric for classification.")
      }
    } else if (task == "regression") {
      predicted_values <- model$assignments
      metric_value <- if (metric == "rmse") {
        sqrt(mean((predicted_values - y)^2))
      } else {
        stop("Invalid metric for regression.")
      }
    }
    
    # Store results
    results <- rbind(results, cbind(params, metric_value))
  }
  
  # Find the best parameters
  best_idx <- if (metric == "rmse") {
    which.min(results$metric_value)
  } else {
    which.max(results$metric_value)
  }
  best_params <- grid[best_idx, ]
  
  # Return results
  list(
    best_params = best_params,
    results = results
  )
}
