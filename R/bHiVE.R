#' bHIVE: B-cell Hybrid Immune Variant Engine
#'
#' Implements an artificial immune network algorithm for clustering and
#' classification tasks. The algorithm evolves a population of "antibodies"
#' via clonal selection and mutation, applies network suppression to maintain
#' diversity, and assigns data points based on affinity or distance metrics.
#'
#' @param X A numeric matrix or data frame of input features, with rows as
#' observations and columns as features.
#' @param y Optional. A factor target vector for classification. If NULL,
#' clustering will be performed.
#' @param task Character. Specifies the task to perform: \code{"clustering"} or
#' \code{"classification"}. If NULL, it is inferred based on \code{y}.
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
#' @param scale Character. Per-feature input scaling: \code{"none"} (default),
#' \code{"zscore"}, \code{"robust"} (median/IQR), or \code{"arcsinh"} (CyTOF).
#' Passed to \code{\link{AINet}}; makes \code{epsilon} and distances behave
#' consistently across datasets of different magnitude.
#' @param targetK Integer or NULL. If set, force the clustering result to
#' exactly \code{targetK} clusters via a seeded K-means refinement (the immune
#' network supplies the seeds). NULL (default) keeps the emergent cluster count.
#' @param epsilonQuantile Numeric in (0, 1) or NULL. If set, the suppression
#' threshold adapts each iteration to this quantile of pairwise antibody
#' distances instead of the fixed \code{epsilon}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages each
#' iteration.
#' @param ... Additional arguments forwarded to \code{\link{AINet}}, including
#' immunology modules (\code{shm}, \code{idiotypic}, \code{germinalCenter},
#' \code{microenvironment}, \code{activation}, \code{memory},
#' \code{classSwitcher}, \code{init}) and \code{consolidate} /
#' \code{consolidationSteps}.
#'
#' @return A list:
#'   \itemize{
#'     \item \code{antibodies}: Final antibody vectors (nAntibodies x nFeatures).
#'     \item \code{assignments}:
#'         - For clustering: integer cluster IDs in [1..#Antibodies].
#'         - For classification: predicted labels.
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
                  scale = c("none", "zscore", "robust", "arcsinh"),
                  targetK = NULL,
                  epsilonQuantile = NULL,
                  verbose = TRUE,
                  ...) {
  # ====================================
  # 0) Basic Validation & Task Inference
  # ====================================
  # bHIVE is now a thin functional wrapper over the AINet R6 engine. This is
  # the single code path: the C++ clonal-selection/suppression backends,
  # Lloyd consolidation, scaling, target-K, and the composable immunology
  # modules all live in AINet, so swarmbHIVE() and honeycombHIVE() (which call
  # bHIVE) inherit every one of them. The previous pure-R loop duplicated a
  # slower, module-free subset of this and has been retired.
  .validate_bHIVE_input(X, y)

  if (is.null(task)) {
    task <- if (is.null(y)) "clustering" else "classification"
  }
  task <- match.arg(task, c("clustering", "classification"))
  initMethod <- match.arg(initMethod,
                          c("sample", "random", "random_uniform", "kmeans++"))
  scale <- match.arg(scale)

  model <- AINet$new(
    nAntibodies        = nAntibodies,
    beta               = beta,
    epsilon            = epsilon,
    maxIter            = maxIter,
    k                  = k,
    affinityFunc       = affinityFunc,
    distFunc           = distFunc,
    affinityParams     = affinityParams,
    mutationDecay      = mutationDecay,
    mutationMin        = mutationMin,
    maxClones          = maxClones,
    stopTolerance      = stopTolerance,
    noImprovementLimit = noImprovementLimit,
    initMethod         = initMethod,
    scale              = scale,
    targetK            = targetK,
    epsilonQuantile    = epsilonQuantile,
    verbose            = verbose,
    ...
  )
  model$fit(X, y = y, task = task)

  # Preserve the historical return contract: a plain list with antibodies,
  # assignments and task. Carry the back-transformed prototypes and the
  # fitted model along for callers that want them (NULL-safe for old code).
  res <- model$result
  res$model <- model
  res
}
