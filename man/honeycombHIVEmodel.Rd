% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/honeycombHiVEmodel.R
\docType{data}
\name{honeycombHIVEmodel}
\alias{honeycombHIVEmodel}
\title{Mulilayered honeycombHIVE for caret}
\format{
An object of class \code{list} of length 8.
}
\usage{
honeycombHIVEmodel
}
\value{
A \code{caret} model definition list. Pass it to
  \code{\link[caret]{train}} for model training and evaluation.
}
\description{
A \code{caret} wrapper for the \code{\link{honeycombHIVE}} function,
enabling seamless integration with the \code{caret} package for
hyperparameter tuning, cross-validation, and performance evaluation.
}
\section{Parameters}{

  \itemize{
    \item \code{nAntibodies}: Number of initial antibodies in the network.
    \item \code{beta}: Clone multiplier controlling the number of clones per
          antibody.
    \item \code{epsilon}: Threshold for network suppression to remove
          redundant antibodies.
    \item \code{layers}: Number of hierarchical layers for iterative
          refinement.
    \item \code{refineOptimizer}: Optimizer for gradient-based refinement
          (e.g. "sgd", "momentum", "adagrad", "adam", "rmsprop").
    \item \code{refineSteps}: Number of gradient update steps in refinement.
    \item \code{refineLR}: Learning rate for refinement.
    \item \code{refineHuberDelta}: Delta parameter used if the "huber" loss is chosen.
  }
}

\section{Supported Tasks}{

  \itemize{
    \item \code{"Regression"}: Predicts numeric target values.
    \item \code{"Classification"}: Assigns class labels to input observations.
    \item \code{"Clustering"}: Groups data points based on similarity
          (though typically \code{caret} is used for supervised tasks).
  }
}

\examples{
\dontrun{
  library(caret)
  # Example: Classification with Iris
  data(iris)
  X <- as.matrix(iris[, 1:4])
  y <- iris$Species

  train_control <- trainControl(method = "cv", number = 5)
  set.seed(42)
  model <- train(
    x = X,
    y = y,
    method = honeycombHIVEmodel,
    trControl = train_control,
    tuneGrid = expand.grid(
      nAntibodies = c(10, 20),
      beta = c(3, 5),
      epsilon = c(0.01, 0.05),
      layers = c(1, 2),
      refineOptimizer = "adam",
      refineSteps = 5,
      refineLR = 0.01,
      refineHuberDelta = 1.0
    )
  )
  print(model)
}

}
\seealso{
\code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
}
\keyword{datasets}
