#' @title ImmuneAlgorithm
#' @description Abstract R6 base class for all immune-inspired algorithms.
#' Subclasses must implement the \code{fit} method.
#'
#' @param config Named list of hyperparameters.
#' @param modules Named list of module instances.
#'
#' @examples
#' # ImmuneAlgorithm is abstract; use AINet for concrete instances
#' algo <- ImmuneAlgorithm$new()
#' print(algo)
#'
#' @importFrom R6 R6Class
#' @export
ImmuneAlgorithm <- R6::R6Class(
  "ImmuneAlgorithm",
  public = list(

    #' @field repertoire An \code{\link{ImmuneRepertoire}} object.
    repertoire = NULL,

    #' @field config Named list of algorithm hyperparameters.
    config = NULL,

    #' @field modules Named list of injected module instances
    #'   (SHMEngine, IdiotypicNetwork, GerminalCenter, etc.).
    modules = NULL,

    #' @field history List of per-iteration metrics.
    history = NULL,

    #' @field result The result from the last call to \code{fit()}.
    result = NULL,

    #' @description Create a new ImmuneAlgorithm.
    #' @param config Named list of hyperparameters.
    #' @param modules Named list of module instances.
    initialize = function(config = list(), modules = list()) {
      self$config  <- config
      self$modules <- modules
      self$history <- list()
    },

    #' @description Fit the algorithm to data. Must be overridden by subclasses.
    #' @param X Numeric matrix (n x d).
    #' @param y Optional factor target for classification.
    #' @param task Character: "clustering" or "classification".
    #' @param ... Additional arguments.
    #' @return The algorithm object (invisibly), with \code{result} populated.
    fit = function(X, y = NULL, task = NULL, ...) {
      stop("ImmuneAlgorithm$fit() is abstract and must be overridden by subclasses.")
    },

    #' @description Predict on new data using the trained repertoire.
    #' @param newdata Numeric matrix (n_new x d).
    #' @return Predictions (class labels, numeric values, or cluster IDs).
    predict = function(newdata) {
      if (is.null(self$result)) {
        stop("Model has not been fitted yet. Call $fit() first.")
      }
      newdata <- as.matrix(newdata)

      task <- self$result$task
      # Prefer the reported prototypes (consolidated / forced-K centroids for
      # clustering; pruned antibodies for classification) so test-time labels
      # match the training partition. Fall back to the raw repertoire.
      A    <- self$result$antibodies %||% self$repertoire$as_matrix()
      cfg  <- self$config

      # Apply the scaling learned at fit() so new data lives in the same space
      # as the trained antibodies. No-op when scale = "none".
      if (!is.null(cfg$scaling)) {
        newdata <- .bhive_apply_scaling(newdata, cfg$scaling)
      }

      alpha <- cfg$affinityParams$alpha %||% 1
      c_p   <- cfg$affinityParams$c %||% 1
      p_p   <- cfg$affinityParams$p %||% 2
      Sigma_inv <- if (!is.null(cfg$affinityParams$Sigma)) {
        solve(cfg$affinityParams$Sigma)
      } else {
        matrix(0, 0, 0)
      }

      fa <- final_assignment_cpp(
        newdata, A,
        cfg$affinityFunc %||% "gaussian",
        cfg$distFunc %||% "euclidean",
        switch(task, clustering = 0L, classification = 1L),
        alpha, c_p, p_p, Sigma_inv
      )

      if (task == "clustering") {
        return(fa$assignments)
      } else {
        classes <- self$result$antibody_classes
        return(classes[fa$best_antibody_idx])
      }
    },

    #' @description Print summary of the algorithm.
    #' @param ... Not used.
    print = function(...) {
      cat(sprintf("<%s>\n", class(self)[1]))
      if (!is.null(self$repertoire)) {
        cat(sprintf("  Repertoire: %d antibodies x %d features\n",
                    self$repertoire$size(), self$repertoire$n_features()))
      }
      if (!is.null(self$result)) {
        cat(sprintf("  Task: %s\n", self$result$task))
        cat(sprintf("  Iterations: %d\n", length(self$history)))
      } else {
        cat("  (not yet fitted)\n")
      }
      invisible(self)
    },

    #' @description Get a summary of the fitting history.
    #' @return Data frame of per-iteration metrics.
    summary = function() {
      if (length(self$history) == 0) return(NULL)
      do.call(rbind, lapply(seq_along(self$history), function(i) {
        h <- self$history[[i]]
        data.frame(iteration = i, n_antibodies = h$n_antibodies,
                   stringsAsFactors = FALSE)
      }))
    }
  )
)
