#' B-cell-based Hybrid Immune Virtual Evolution (bHIVE) for caret
#'
#' A wrapper for integrating the B-cell Hybrid Immune Variant Engine 
#' (bHIVE) algorithm with the \code{caret} package. Supports both classification 
#' and regression tasks, providing compatibility with \code{caret::train()} 
#' for model training and validation.
#'
#' @name bHIVEmodel
#' @format A list containing the components required for integration with 
#' the \code{caret} package.
#'
#' @details 
#' The \code{bHIVEmodel} wrapper facilitates the use of bHIVE for classification 
#' and regression. It defines the model label, parameter grid, fitting function, 
#' and prediction methods to conform to the \code{caret} model specification.
#'
#' @section Components:
#' \describe{
#'   \item{\code{label}}{Character string. Identifies the model as 
#'   "B-cell-based Hybrid Immune Virtual Evolution".}
#'   \item{\code{library}}{Character string. Specifies the R package containing 
#'   the bHIVE implementation. Default is "customPackage".}
#'   \item{\code{type}}{Character vector. Specifies the supported tasks: 
#'   "Classification" and "Regression".}
#'   \item{\code{parameters}}{A \code{data.frame} describing the tunable 
#'   parameters: 
#'   \itemize{
#'     \item \code{parameter}: Name of the parameter.
#'     \item \code{class}: Data type of the parameter ("numeric").
#'     \item \code{label}: Short description of the parameter.
#'   }}
#'   \item{\code{grid}}{Function. Generates a grid of tuning parameters for 
#'   hyperparameter optimization.}
#'   \item{\code{fit}}{Function. Trains the bHIVE model using specified 
#'   hyperparameters and task type.}
#'   \item{\code{predict}}{Function. Generates predictions for new data 
#'   (classification labels or regression values).}
#'   \item{\code{prob}}{Function. Calculates class probabilities for 
#'   classification tasks.}
#' }
#'
#' @section Parameters:
#' \describe{
#'   \item{\code{nAntibodies}}{Number of initial antibodies in the 
#'   bHIVE algorithm.}
#'   \item{\code{beta}}{Clone multiplier. Controls the number of clones 
#'   generated for top-matching antibodies.}
#'   \item{\code{epsilon}}{Similarity threshold for antibody suppression. 
#'   Smaller values encourage more diversity in the repertoire.}
#' }
#'
#' @section Functions:
#' \itemize{
#'   \item \code{grid(x, y, len)}: Generates a grid of tuning parameters. 
#'   Accepts:
#'     \itemize{
#'       \item \code{x}: Feature matrix or data frame.
#'       \item \code{y}: Target vector (factor for classification, numeric for 
#'       regression).
#'       \item \code{len}: Number of grid points for each parameter.
#'     }
#'   \item \code{fit(x, y, wts, param, lev, last, classProbs, ...)}: Trains 
#'   the bHIVE model. Key arguments:
#'     \itemize{
#'       \item \code{x}: Feature matrix or data frame.
#'       \item \code{y}: Target vector.
#'       \item \code{param}: List of hyperparameters (\code{nAntibodies}, 
#'       \code{beta}, \code{epsilon}).
#'       \item \code{...}: Additional arguments passed to the bHIVE function.
#'     }
#'   \item \code{predict(modelFit, newdata, submodels)}: Generates predictions 
#'   for new data.
#'     \itemize{
#'       \item \code{modelFit}: Trained bHIVE model.
#'       \item \code{newdata}: New feature data for prediction.
#'     }
#'   \item \code{prob(modelFit, newdata, submodels)}: Calculates class 
#'   probabilities (classification only).
#' }
#'
#' @section Example Usage:
#' \preformatted{
#' library(caret)
#' 
#' # Simulated classification dataset
#' set.seed(123)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' y <- factor(sample(c("Class1", "Class2"), 100, replace = TRUE))
#' 
#' # Train bHIVE model using caret
#' trainControl <- trainControl(method = "cv", number = 5, classProbs = TRUE)
#' tunedModel <- train(
#'   x = X,
#'   y = y,
#'   method = bHIVEmodel,
#'   trControl = trainControl,
#'   tuneLength = 3
#' )
#' 
#' # Predictions
#' predictions <- predict(tunedModel, newdata = X)
#' probabilities <- predict(tunedModel, newdata = X, type = "prob")
#' }
#'
#' @seealso \code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
#'
#' @export
bHIVEmodel <- list(
  label = "Artificial Immune Network (bHIVE)",
  library = "bHIVE",
  type = c("Regression", "Classification", "Clustering"),
  parameters = data.frame(
    parameter = c("nAntibodies", "beta", "epsilon"),
    class = c("numeric", "numeric", "numeric"),
    label = c("Number of Antibodies", "Clone Multiplier", "Epsilon Threshold")
  ),
  grid = function(x, y, len = NULL) {
    expand.grid(
      nAntibodies = seq(10, 50, length.out = len),
      beta = seq(1, 10, length.out = len),
      epsilon = seq(0.01, 0.1, length.out = len)
    )
  },
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    bHIVE(
      X = x,
      y = y,
      task = if (is.factor(y)) "classification" else if (is.numeric(y)) 
        "regression" else "clustering",
      nAntibodies = param$nAntibodies,
      beta = param$beta,
      epsilon = param$epsilon,
      ...
    )
  },
  predict = function(modelFit, newdata, submodels = NULL) {
    # Extract the 'assignments' from the modelFit object
    if ("assignments" %in% names(modelFit)) {
      predictions <- apply(newdata, 1, function(row) {
        # Find the closest antibody for clustering/classification/regression
        distances <- apply(modelFit$antibodies, 1, function(a) {
          sum((row - a)^2)  # Euclidean distance
        })
        closest <- which.min(distances)
        if (modelFit$task == "classification") {
          return(modelFit$assignments[closest])
        } else if (modelFit$task == "regression") {
          return(modelFit$assignments[closest])
        } else {
          return(closest)
        }
      })
    } else {
      stop("Invalid model object: missing 'assignments'")
    }
    
    return(predictions)
  },
  prob = function(modelFit, newdata, submodels = NULL) {
    stop("Probabilities are not supported for this model.")
  }
)
