% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bHiVE.R
\name{bHIVE}
\alias{bHIVE}
\title{bHIVE: B-cell Hybrid Immune Variant Engine}
\usage{
bHIVE(
  X,
  y = NULL,
  task = NULL,
  nAntibodies = 20,
  beta = 5,
  epsilon = 0.01,
  maxIter = 50,
  affinityFunc = "gaussian",
  distFunc = "euclidean",
  affinityParams = list(alpha = 1, c = 1, p = 2, Sigma = NULL),
  mutationDecay = 1,
  mutationMin = 0.01,
  maxClones = Inf,
  stopTolerance = 0,
  noImprovementLimit = Inf,
  initMethod = c("sample", "random", "random_uniform", "kmeans++"),
  k = 3,
  verbose = TRUE
)
}
\arguments{
\item{X}{A numeric matrix or data frame of input features, with rows as 
observations and columns as features.}

\item{y}{Optional. A target vector. Use for classification (factor) or 
regression (numeric). If NULL, clustering will be performed.}

\item{task}{Character. Specifies the task to perform: \code{"clustering"}, 
\code{"classification"}, or \code{"regression"}. If NULL, it is inferred 
based on \code{y}.}

\item{nAntibodies}{Integer. The initial population size of antibodies.}

\item{beta}{Numeric. Clone multiplier (controls how many clones are 
generated for top-matching antibodies).}

\item{epsilon}{Numeric. Similarity threshold used in network suppression; 
antibodies closer than \code{epsilon} are considered redundant.}

\item{maxIter}{Integer. Maximum number of iterations to run the AI-Net 
algorithm.}

\item{affinityFunc}{Character. Specifies the affinity (similarity) function 
to use for antibody-data matching. One of \code{"gaussian"}, \code{"laplace"}, 
\code{"polynomial"}, \code{"cosine"}, or \code{"hamming"}.}

\item{distFunc}{Character. Specifies the distance function for clustering 
and suppression. One of \code{"euclidean"}, \code{"manhattan"}, 
\code{"minkowski"}, \code{"cosine"}, \code{"mahalanobis"}, 
or \code{"hamming"}.}

\item{affinityParams}{A list of optional parameters for the chosen affinity 
or distance function.
  \itemize{
    \item \code{alpha} (for RBF or Laplace kernel),
    \item \code{c}, \code{p} (for polynomial kernel or Minkowski distance),
    \item \code{Sigma} (for Mahalanobis distance).
  }}

\item{mutationDecay}{Numeric. Factor by which the mutation rate decays each 
iteration (should be \eqn{\le 1.0}). Default is 1.0 (no decay).}

\item{mutationMin}{Numeric. Minimum mutation rate, preventing the mutation 
scale from shrinking to zero.}

\item{maxClones}{Numeric. Maximum number of clones per top-matching antibody;
defaults to \code{Inf}.}

\item{stopTolerance}{Numeric. If the change in the number of antibodies 
(repertoire size) is \eqn{\le stopTolerance} for consecutive iterations, 
this may trigger the \code{noImprovementLimit}.}

\item{noImprovementLimit}{Integer. Stops the algorithm early if there is no 
further improvement in antibody count (beyond \code{stopTolerance}) for this
many consecutive iterations. Default is \code{Inf}, meaning no early stop 
based on improvement.}

\item{initMethod}{Character. Method for initializing antibodies. Can be:
\itemize{
  \item \code{"sample"} - randomly selects rows from \code{X} as initial 
  antibodies.
  \item \code{"random"} - samples Gaussian noise using \code{X}'s column 
  means/sds.
  \item \code{"random_uniform"} - samples uniformly in [min, max] of each 
  column.
  \item \code{"kmeans++"} - tries a kmeans++-like initialization for 
  coverage.
}}

\item{k}{Integer. Number of top-matching antibodies (by affinity) to 
consider cloning for each data point.}

\item{verbose}{Logical. If \code{TRUE}, prints progress messages each 
iteration.}
}
\value{
A list:
  \itemize{
    \item \code{antibodies}: Final antibody vectors (nAntibodies x nFeatures).
    \item \code{assignments}: 
        - For clustering: integer cluster IDs in [1..#Antibodies].
        - For classification: predicted labels.
        - For regression: integer cluster index (in [1..#Antibodies]) if used in synergy with \code{refineB}.
    \item \code{predictions}: Only for regression, the numeric predictions per row.
    \item \code{task}: The chosen task.
  }
}
\description{
Implements an artificial immune network algorithm for clustering, classification, 
and regression tasks. The algorithm evolves a population of "antibodies" 
via clonal selection and mutation, applies network suppression to maintain 
diversity, and assigns data points based on affinity or distance metrics.
}
\examples{
# Example 1: Clustering with the Iris dataset
data(iris)
X <- as.matrix(iris[, 1:4])  # Numeric features only
res <- bHIVE(X = X, 
             task = "clustering", 
             nAntibodies = 30, 
             beta = 5, 
             epsilon = 0.01, 
             maxIter = 20, 
             k = 3, 
             verbose = FALSE)
table(res$assignments)

# Example 2: Classification with Iris species
y <- iris$Species
res <- bHIVE(X = X, 
              y = y, 
              task = "classification", 
              nAntibodies = 30, 
              beta = 5, 
              epsilon = 0.01, 
              maxIter = 20, 
              k = 3, 
              verbose = FALSE)
table(res$assignments, y)

# Example 3: Regression
y <- as.numeric(iris$Sepal.Length)
res <- bHIVE(X = X, 
             y = y, 
             task = "regression", 
             nAntibodies = 30, 
             beta = 5, 
             epsilon = 0.01, 
             maxIter = 20, 
             k = 3, 
             verbose = FALSE)
cor(res$assignments, y)

}
