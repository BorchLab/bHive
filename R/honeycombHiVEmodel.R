#' Mulilayered honeycombHIVE for caret
#'
#' A \code{caret} wrapper for the \code{\link{honeycombHIVE}} function,
#' enabling seamless integration with the \code{caret} package for
#' hyperparameter tuning, cross-validation, and performance evaluation.
#'
#' @name honeycombHIVEmodel
#' @section Parameters:
#'   \itemize{
#'     \item \code{nAntibodies}: Number of initial antibodies in the network.
#'     \item \code{beta}: Clone multiplier controlling the number of clones per
#'           antibody.
#'     \item \code{epsilon}: Threshold for network suppression to remove
#'           redundant antibodies.
#'     \item \code{layers}: Number of hierarchical layers for iterative
#'           refinement.
#'   }
#'
#' @section Supported Tasks:
#'   \itemize{
#'     \item \code{"Regression"}: Predicts numeric target values.
#'     \item \code{"Classification"}: Assigns class labels to input observations.
#'     \item \code{"Clustering"}: Groups data points based on similarity
#'           (though typically \code{caret} is used for supervised tasks).
#'   }
#'
#' @return A \code{caret} model definition list. Pass it to
#'   \code{\link[caret]{train}} for model training and evaluation.
#'
#' @examples
#' \dontrun{
#'   library(caret)
#'   # Example: Classification with Iris
#'   data(iris)
#'   X <- as.matrix(iris[, 1:4])
#'   y <- iris$Species
#'
#'   train_control <- trainControl(method = "cv", number = 5)
#'   set.seed(42)
#'   model <- train(
#'     x = X,
#'     y = y,
#'     method = honeycombHIVEmodel,
#'     trControl = train_control,
#'     tuneGrid = expand.grid(
#'       nAntibodies = c(10, 20),
#'       beta = c(3, 5),
#'       epsilon = c(0.01, 0.05),
#'       layers = c(1, 2)
#'     )
#'   )
#'   print(model)
#' }
#'
#' @seealso \code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
#' @export
honeycombHIVEmodel <- list(
  label = "Multilayered Artificial Immune Network (honeycombHIVE)",
  library = c("bHIVE"),  
  type = c("Regression", "Classification", "Clustering"),
  
  # Hyperparameters recognized by 'caret' for tuning
  parameters = data.frame(
    parameter = c("nAntibodies", "beta", "epsilon", "layers"),
    class = c("numeric", "numeric", "numeric", "integer"),
    label = c(
      "Number of Antibodies",
      "Clone Multiplier (beta)",
      "Epsilon Threshold",
      "Number of Layers"
    )
  ),
  
  # Function to generate a default tuning grid if the user doesn't provide one
  grid = function(x, y, len = NULL) {
    expand.grid(
      nAntibodies = seq(10, 50, length.out = len),
      beta = seq(1, 10, length.out = len),
      epsilon = seq(0.01, 0.1, length.out = len),
      layers = seq(1, 3, length.out = len)
    )
  },
  
  # The core "fit" function: trains honeycombHIVE on the given data
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    task <- if (is.factor(y)) {
      "classification"
    } else if (is.numeric(y)) {
      "regression"
    } else {
      "clustering"
    }
    
    # Train the honeycombHIVE model
    results <- honeycombHIVE(
      X           = x,
      y           = y,
      task        = task,
      nAntibodies = param$nAntibodies,
      beta        = param$beta,
      epsilon     = param$epsilon,
      layers      = param$layers,
      ...
    )
    
    # The model object we return to caret is "modelFit".
    # We'll store:
    # 1) The entire results list from honeycombHIVE
    # 2) The final layer's prototypes, plus any needed info for prediction
    final_layer <- results[[length(results)]]
    
    modelFit <- list(
      task        = task,
      results     = results,      # all layers
      finalLayer  = final_layer,  # convenience
      obsLevels   = if (task == "classification") lev else NULL
    )
    
    if (task == "classification") {
      modelFit$levels <- lev
    } else {
      modelFit$levels <- NULL
    }
    
    class(modelFit) <- "honeycombHIVEmodel"
    return(modelFit)
  },
  
  # Prediction function for new data
  predict = function(modelFit, newdata, submodels = NULL) {
    task <- modelFit$task
    finalLayer <- modelFit$finalLayer
    
    # If the honeycombHIVE is purely clustering:
    if (task == "clustering") {
      stop("Prediction not defined for clustering in this caret wrapper.")
    }

    # 1) Get final antibodies from the final layer
    finalAb <- finalLayer$antibodies
    if (is.null(finalAb)) {
      stop("No antibodies found in final layer. Check your honeycombHIVE output.")
    }
    
    # 2) Optionally, get the label/value for each antibody
    if (task == "classification") {
      # We'll assume finalLayer$antibodyLabels is a factor with length = nrow(finalAb)
      antibodyLabels <- finalLayer$antibodyLabels
      if (is.null(antibodyLabels)) {
        stop("No final antibody labels found. Check how your honeycombHIVE stores them.")
      }
    } else if (task == "regression") {
      # We'll assume finalLayer$antibodyValues is a numeric vector, length = nrow(finalAb)
      antibodyValues <- finalLayer$antibodyValues
      if (is.null(antibodyValues)) {
        stop("No final antibody regression values found. Check how your honeycombHIVE stores them.")
      }
    }
    
    # 3) Compute distances from newdata to finalAb
    newdata <- as.matrix(newdata)  # ensure matrix format
    finalAb <- as.matrix(finalAb)
    preds <- vector("list", nrow(newdata))  # store index of nearest antibody
    
    for (i in seq_len(nrow(newdata))) {
      diffs <- sweep(finalAb, 2, newdata[i,], FUN = "-")
      dists <- rowSums(diffs^2)  # Euclidean squared distance
      idx_min <- which.min(dists)
      preds[[i]] <- idx_min
    }
    # Flatten
    nearest_idx <- unlist(preds)
    
    # 4) Return the label or value for each row
    if (task == "classification") {
      # Then we pick antibodyLabels[nearest_idx]
      out <- antibodyLabels[nearest_idx]
      # Ensure factor with correct levels
      # caretaker typically wants predicted factor with modelFit$levels
      out <- factor(out, levels = modelFit$levels)
      return(out)
      
    } else if (task == "regression") {
      # Numeric
      out <- antibodyValues[nearest_idx]
      return(out)
    }
  },
  
  # Probability predictions for classification
  prob = function(modelFit, newdata, submodels = NULL) {
    stop("Probability predictions are not supported for honeycombHIVE 
         in this wrapper.")
  }
)
